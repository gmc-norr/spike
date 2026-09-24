//! BAM/CRAM read pair extraction and ReadPool building.
//!
//! Extracts paired-end reads from a BAM/CRAM region, stores them in their original
//! sequencing orientation (FASTQ order).

use anyhow::{Context, Result};
use std::collections::HashMap;

use crate::stats::FragmentDist;
use crate::types::{ReadPair, ReadPool};

/// Check if a path refers to a CRAM file (by extension).
pub fn is_cram(path: &str) -> bool {
    path.ends_with(".cram")
}

/// Build a noodles FASTA repository for CRAM decoding.
pub fn build_fasta_repository(ref_path: &str) -> Result<noodles::fasta::Repository> {
    let fasta_reader = noodles::fasta::io::indexed_reader::Builder::default()
        .build_from_path(ref_path)
        .with_context(|| format!("failed to open FASTA index for CRAM decoding: {}", ref_path))?;
    let adapter = noodles::fasta::repository::adapters::IndexedReader::new(fasta_reader);
    Ok(noodles::fasta::Repository::new(adapter))
}

/// Extract all properly-paired read pairs from a genomic region.
///
/// Supports both BAM and CRAM formats (detected by file extension).
/// CRAM files require `ref_path` to be `Some`.
///
/// Two-pass approach:
/// 1. First pass: collect read1 records (is_first_segment) with seq, qual, pos, tlen
/// 2. Second pass: collect read2 records, match by name
///
/// Reads aligned in reverse complement are reverse-complemented back to FASTQ orientation.
pub fn extract_read_pairs(
    alignment_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
    ref_path: Option<&str>,
) -> Result<Vec<ReadPair>> {
    if is_cram(alignment_path) {
        let rp = ref_path.ok_or_else(|| {
            anyhow::anyhow!("CRAM input requires a reference FASTA (--reference)")
        })?;
        extract_read_pairs_cram(alignment_path, chrom, start, end, min_mapq, rp)
    } else {
        extract_read_pairs_bam(alignment_path, chrom, start, end, min_mapq)
    }
}

/// BAM-specific read pair extraction (existing implementation).
fn extract_read_pairs_bam(
    bam_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
) -> Result<Vec<ReadPair>> {
    log::info!(
        "Extracting read pairs from {}:{}-{} (MAPQ >= {})",
        chrom,
        start,
        end,
        min_mapq,
    );

    // Pass 1: collect both read1 and read2 records in the target region.
    let mut read1_map: HashMap<String, PartialRead> = HashMap::new();
    let mut read2_map: HashMap<String, PartialRead> = HashMap::new();
    let mut pass1_read1_count = 0usize;
    let mut pass1_read2_count = 0usize;
    let mut max_abs_tlen = 0u64;
    let mut missing_qual_count = 0usize;

    {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM: {}", bam_path))?;
        let header = reader.read_header()?;

        let start_pos = safe_noodles_position(start + 1);
        let end_pos = safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();

            if !passes_filters_bam(&flags, min_mapq, &record) {
                continue;
            }
            if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
                continue;
            }
            if !flags.is_first_segment() && !flags.is_last_segment() {
                continue;
            }

            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };
            let partial =
                match parse_partial_from_bam_record(&record, &flags, &mut missing_qual_count) {
                    Some(p) => p,
                    None => continue,
                };
            max_abs_tlen = max_abs_tlen.max(partial.tlen.unsigned_abs() as u64);

            if flags.is_first_segment() {
                read1_map.insert(name, partial);
                pass1_read1_count += 1;
            } else {
                read2_map.insert(name, partial);
                pass1_read2_count += 1;
            }
        }
    }

    // Pair records already complete in pass 1.
    let mut pairs: Vec<ReadPair> = Vec::new();
    let pass1_names: Vec<String> = read1_map.keys().cloned().collect();
    let mut pass1_paired = 0usize;
    for name in pass1_names {
        if let (Some(read1), Some(read2)) = (read1_map.remove(&name), read2_map.remove(&name)) {
            pairs.push(build_pair_from_partials(name, read1, read2, chrom));
            pass1_paired += 1;
        }
    }

    log::info!(
        "Pass 1: collected {} read1 + {} read2 records ({} pairs already complete)",
        pass1_read1_count,
        pass1_read2_count,
        pass1_paired,
    );

    // Pass 2: query a wider region to recover missing mates.
    let wider_padding = max_abs_tlen.saturating_add(200).max(1000);

    {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM for pass 2: {}", bam_path))?;
        let header = reader.read_header()?;

        let wider_start = start.saturating_sub(wider_padding);
        let wider_end = end.saturating_add(wider_padding);

        let start_pos = safe_noodles_position(wider_start + 1);
        let end_pos = safe_noodles_position(wider_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();

            if !passes_filters_bam(&flags, min_mapq, &record) {
                continue;
            }
            if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
                continue;
            }
            if !flags.is_first_segment() && !flags.is_last_segment() {
                continue;
            }

            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };

            if flags.is_first_segment() {
                if read1_map.contains_key(&name) {
                    continue;
                }
                if let Some(read2) = read2_map.remove(&name) {
                    let Some(read1) =
                        parse_partial_from_bam_record(&record, &flags, &mut missing_qual_count)
                    else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            } else {
                if read2_map.contains_key(&name) {
                    continue;
                }
                if let Some(read1) = read1_map.remove(&name) {
                    let Some(read2) =
                        parse_partial_from_bam_record(&record, &flags, &mut missing_qual_count)
                    else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            }
        }
    }

    let unmatched = read1_map.len() + read2_map.len();
    log_extraction_result(missing_qual_count, unmatched, pairs.len(), chrom, start, end);
    Ok(pairs)
}

/// CRAM-specific read pair extraction.
///
/// Converts CRAM records to RecordBuf for uniform field access.
fn extract_read_pairs_cram(
    cram_path: &str,
    chrom: &str,
    start: u64,
    end: u64,
    min_mapq: u8,
    ref_path: &str,
) -> Result<Vec<ReadPair>> {
    log::info!(
        "Extracting read pairs (CRAM) from {}:{}-{} (MAPQ >= {})",
        chrom,
        start,
        end,
        min_mapq,
    );

    let repository = build_fasta_repository(ref_path)?;

    // Pass 1: collect both read1 and read2 records in the target region.
    let mut read1_map: HashMap<String, PartialRead> = HashMap::new();
    let mut read2_map: HashMap<String, PartialRead> = HashMap::new();
    let mut pass1_read1_count = 0usize;
    let mut pass1_read2_count = 0usize;
    let mut max_abs_tlen = 0u64;
    let mut missing_qual_count = 0usize;

    {
        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository.clone())
            .build_from_path(cram_path)
            .with_context(|| format!("failed to open CRAM: {}", cram_path))?;
        let header = reader.read_header()?;

        let start_pos = safe_noodles_position(start + 1);
        let end_pos = safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record
                .try_into_alignment_record(&header)
                .with_context(|| "failed to convert CRAM record to alignment record")?;

            let flags = buf.flags();

            if !passes_filters_buf(min_mapq, &buf) {
                continue;
            }
            if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
                continue;
            }
            if !flags.is_first_segment() && !flags.is_last_segment() {
                continue;
            }

            let name = match buf.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };
            let partial =
                match parse_partial_from_record_buf(&buf, &flags, &mut missing_qual_count) {
                    Some(p) => p,
                    None => continue,
                };
            max_abs_tlen = max_abs_tlen.max(partial.tlen.unsigned_abs() as u64);

            if flags.is_first_segment() {
                read1_map.insert(name, partial);
                pass1_read1_count += 1;
            } else {
                read2_map.insert(name, partial);
                pass1_read2_count += 1;
            }
        }
    }

    // Pair records already complete in pass 1.
    let mut pairs: Vec<ReadPair> = Vec::new();
    let pass1_names: Vec<String> = read1_map.keys().cloned().collect();
    let mut pass1_paired = 0usize;
    for name in pass1_names {
        if let (Some(read1), Some(read2)) = (read1_map.remove(&name), read2_map.remove(&name)) {
            pairs.push(build_pair_from_partials(name, read1, read2, chrom));
            pass1_paired += 1;
        }
    }

    log::info!(
        "Pass 1 (CRAM): collected {} read1 + {} read2 records ({} pairs already complete)",
        pass1_read1_count,
        pass1_read2_count,
        pass1_paired,
    );

    // Pass 2: query a wider region to recover missing mates.
    let wider_padding = max_abs_tlen.saturating_add(200).max(1000);

    {
        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(cram_path)
            .with_context(|| format!("failed to open CRAM for pass 2: {}", cram_path))?;
        let header = reader.read_header()?;

        let wider_start = start.saturating_sub(wider_padding);
        let wider_end = end.saturating_add(wider_padding);

        let start_pos = safe_noodles_position(wider_start + 1);
        let end_pos = safe_noodles_position(wider_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let cram_record = rec_result?;
            let buf = cram_record
                .try_into_alignment_record(&header)
                .with_context(|| "failed to convert CRAM record to alignment record")?;

            let flags = buf.flags();

            if !passes_filters_buf(min_mapq, &buf) {
                continue;
            }
            if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
                continue;
            }
            if !flags.is_first_segment() && !flags.is_last_segment() {
                continue;
            }

            let name = match buf.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };

            if flags.is_first_segment() {
                if read1_map.contains_key(&name) {
                    continue;
                }
                if let Some(read2) = read2_map.remove(&name) {
                    let Some(read1) =
                        parse_partial_from_record_buf(&buf, &flags, &mut missing_qual_count)
                    else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            } else {
                if read2_map.contains_key(&name) {
                    continue;
                }
                if let Some(read1) = read1_map.remove(&name) {
                    let Some(read2) =
                        parse_partial_from_record_buf(&buf, &flags, &mut missing_qual_count)
                    else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            }
        }
    }

    let unmatched = read1_map.len() + read2_map.len();
    log_extraction_result(missing_qual_count, unmatched, pairs.len(), chrom, start, end);
    Ok(pairs)
}

/// Build a ReadPool from extracted read pairs.
///
/// Sorts pairs by ref_start, then name. The name tie-break makes the order
/// independent of extraction order (HashMap order), so a fixed seed gives
/// the same output on every run.
pub fn build_read_pool(pairs: Vec<ReadPair>, frag_dist: FragmentDist) -> ReadPool {
    let mut pairs = pairs;
    pairs.sort_by(|a, b| a.ref_start.cmp(&b.ref_start).then_with(|| a.name.cmp(&b.name)));

    log::info!("Built read pool: {} pairs", pairs.len());

    ReadPool { pairs, frag_dist }
}

// --- Internal helpers ---

struct PartialRead {
    seq: Vec<u8>,
    qual: Vec<u8>,
    pos: u64,
    tlen: i32,
}

/// Compute the fragment end position from template length and read positions.
fn compute_ref_end(read1: &PartialRead, read2: &PartialRead, ref_start: u64) -> u64 {
    let tlen = if read1.tlen != 0 {
        read1.tlen
    } else {
        read2.tlen
    };
    if tlen > 0 {
        ref_start + tlen as u64
    } else if tlen < 0 {
        ref_start + (tlen as i64).unsigned_abs()
    } else {
        // tlen == 0: estimate from positions + read length.
        let r1_end = read1.pos + read1.seq.len() as u64;
        let r2_end = read2.pos + read2.seq.len() as u64;
        r1_end.max(r2_end)
    }
}

/// Build a complete ReadPair from read1/read2 partial records.
fn build_pair_from_partials(
    name: String,
    read1: PartialRead,
    read2: PartialRead,
    chrom: &str,
) -> ReadPair {
    let ref_start = read1.pos.min(read2.pos);
    let ref_end = compute_ref_end(&read1, &read2, ref_start);
    let insert_size = if read1.tlen != 0 {
        (read1.tlen as i64).unsigned_abs() as i64
    } else {
        (read2.tlen as i64).unsigned_abs() as i64
    };

    ReadPair {
        name,
        seq1: read1.seq,
        qual1: read1.qual,
        seq2: read2.seq,
        qual2: read2.qual,
        ref_start,
        ref_end,
        insert_size,
        chrom: chrom.to_string(),
    }
}

/// Whether raw (pre-Phred+33) quality-score bytes indicate that no quality
/// data was stored for this record.
///
/// SAM's QUAL field is all-or-nothing per record: it is either `*` (missing)
/// or a full string the same length as SEQ — never a mix of the two within
/// one record. BAM/CRAM encode "missing" as every per-base byte set to 0xFF
/// (255); a real Phred score never reaches 255 (max is 93), so this check is
/// unambiguous.
fn quality_is_missing(raw_qual: &[u8]) -> bool {
    !raw_qual.is_empty() && raw_qual.iter().all(|&b| b == 0xff)
}

/// Parse sequence/quality/position fields from a BAM record into a partial read.
///
/// Records whose quality is entirely missing (SAM `*`, encoded as every
/// per-base byte 0xFF) are skipped, and `missing_qual_count` is incremented.
fn parse_partial_from_bam_record(
    record: &noodles::bam::Record,
    flags: &noodles::sam::alignment::record::Flags,
    missing_qual_count: &mut usize,
) -> Option<PartialRead> {
    let pos = match record.alignment_start() {
        Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
        _ => return None,
    };
    if record.mate_alignment_start().is_none() {
        return None;
    }
    let tlen = record.template_length();

    let raw_qual = record.quality_scores();
    let raw_qual = raw_qual.as_ref();
    if quality_is_missing(raw_qual) {
        *missing_qual_count += 1;
        return None;
    }
    let mut seq: Vec<u8> = record.sequence().iter().collect();
    let mut qual: Vec<u8> = raw_qual.iter().map(|s| s.wrapping_add(33)).collect();

    if flags.is_reverse_complemented() {
        reverse_complement(&mut seq);
        qual.reverse();
    }

    Some(PartialRead {
        seq,
        qual,
        pos,
        tlen,
    })
}

/// Parse sequence/quality/position fields from a CRAM RecordBuf into a partial read.
///
/// Records whose quality is entirely missing (SAM `*`, encoded as every
/// per-base byte 0xFF) are skipped, and `missing_qual_count` is incremented.
fn parse_partial_from_record_buf(
    buf: &noodles::sam::alignment::RecordBuf,
    flags: &noodles::sam::alignment::record::Flags,
    missing_qual_count: &mut usize,
) -> Option<PartialRead> {
    let pos = match buf.alignment_start() {
        Some(p) => usize::from(p).saturating_sub(1) as u64,
        None => return None,
    };
    if buf.mate_alignment_start().is_none() {
        return None;
    }
    let tlen = buf.template_length();

    let raw_qual = buf.quality_scores();
    let raw_qual = raw_qual.as_ref();
    if quality_is_missing(raw_qual) {
        *missing_qual_count += 1;
        return None;
    }
    let mut seq: Vec<u8> = buf.sequence().as_ref().to_vec();
    let mut qual: Vec<u8> = raw_qual.iter().map(|s| s.wrapping_add(33)).collect();

    if flags.is_reverse_complemented() {
        reverse_complement(&mut seq);
        qual.reverse();
    }

    Some(PartialRead {
        seq,
        qual,
        pos,
        tlen,
    })
}

/// Log extraction results (shared between BAM and CRAM paths).
fn log_extraction_result(
    missing_qual_count: usize,
    unmatched: usize,
    pair_count: usize,
    chrom: &str,
    start: u64,
    end: u64,
) {
    if missing_qual_count > 0 {
        // M14: a record whose quality is entirely missing (SAM `*`) is
        // dropped rather than silently turned into space characters, which
        // used to poison the learned quality model. Surface the count so a
        // user can tell this is happening to their donor BAM.
        log::warn!(
            "{} record(s) in {}:{}-{} had no quality scores (SAM '*') and were skipped",
            missing_qual_count,
            chrom,
            start,
            end,
        );
    }
    if unmatched > 0 {
        log::debug!(
            "{} records had no matching mate in the queried windows",
            unmatched,
        );
    }
    log::info!(
        "Extracted {} complete read pairs from {}:{}-{}",
        pair_count,
        chrom,
        start,
        end,
    );
}

/// Filter check for BAM records.
fn passes_filters_bam(
    flags: &noodles::sam::alignment::record::Flags,
    min_mapq: u8,
    record: &noodles::bam::Record,
) -> bool {
    if flags.is_unmapped()
        || flags.is_secondary()
        || flags.is_supplementary()
        || flags.is_duplicate()
        || flags.is_qc_fail()
    {
        return false;
    }
    let mq: u8 = match record.mapping_quality() {
        Some(q) => u8::from(q),
        None => 0,
    };
    mq >= min_mapq
}

/// Filter check for RecordBuf (used by CRAM path).
fn passes_filters_buf(min_mapq: u8, buf: &noodles::sam::alignment::RecordBuf) -> bool {
    let flags = buf.flags();
    if flags.is_unmapped()
        || flags.is_secondary()
        || flags.is_supplementary()
        || flags.is_duplicate()
        || flags.is_qc_fail()
    {
        return false;
    }
    let mq: u8 = match buf.mapping_quality() {
        Some(q) => u8::from(q),
        None => 0,
    };
    mq >= min_mapq
}

/// Convert a coordinate to a noodles 1-based Position.
///
/// Saturates to [`noodles::core::Position::MIN`] (= 1) on underflow and to
/// `usize::MAX` on overflow. Both extremes are valid noodles positions, so
/// this function never panics.
pub fn safe_noodles_position(pos: u64) -> noodles::core::Position {
    let pos_usize = usize::try_from(pos).unwrap_or(usize::MAX).max(1);
    // pos_usize >= 1 is guaranteed by .max(1) above; new() only fails for 0.
    noodles::core::Position::new(pos_usize).unwrap_or(noodles::core::Position::MIN)
}

/// Reverse-complement a DNA sequence in place.
/// Complement a single DNA base: A↔T, C↔G.
pub fn complement_base(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other, // N stays N
    }
}

pub fn reverse_complement(seq: &mut [u8]) {
    seq.reverse();
    for base in seq.iter_mut() {
        *base = match *base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            other => other, // N stays N
        };
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_read_pool_order_does_not_depend_on_input_order() {
        // Extraction yields pairs in HashMap order, which changes between runs.
        // Pairs with the same start must still end up in one fixed order, or
        // the same --seed suppresses different reads.
        let pair = |name: &str, start: u64| ReadPair {
            name: name.to_string(),
            seq1: vec![],
            qual1: vec![],
            seq2: vec![],
            qual2: vec![],
            ref_start: start,
            ref_end: start + 400,
            insert_size: 400,
            chrom: "chr1".to_string(),
        };
        let names = |input: Vec<ReadPair>| -> Vec<String> {
            let pool = build_read_pool(input, FragmentDist::from_stats(400.0, 80.0));
            pool.pairs.into_iter().map(|p| p.name).collect()
        };

        let run1 = names(vec![pair("b", 10), pair("a", 10), pair("c", 5)]);
        let run2 = names(vec![pair("a", 10), pair("c", 5), pair("b", 10)]);

        assert_eq!(run1, vec!["c", "a", "b"]);
        assert_eq!(run2, vec!["c", "a", "b"]);
    }

    fn partial(pos: u64, len: usize, tlen: i32) -> PartialRead {
        PartialRead {
            seq: vec![b'A'; len],
            qual: vec![b'!' + 30; len],
            pos,
            tlen,
        }
    }

    #[test]
    fn test_reverse_complement() {
        let mut seq = b"ACGTNN".to_vec();
        reverse_complement(&mut seq);
        assert_eq!(&seq, b"NNACGT");
    }

    #[test]
    fn test_reverse_complement_empty() {
        let mut seq = Vec::new();
        reverse_complement(&mut seq);
        assert!(seq.is_empty());
    }

    #[test]
    fn test_is_cram() {
        assert!(is_cram("sample.cram"));
        assert!(!is_cram("sample.bam"));
        assert!(!is_cram("sample.cram.bai"));
    }

    #[test]
    fn test_compute_ref_end_uses_positive_tlen() {
        let r1 = partial(100, 150, 400);
        let r2 = partial(320, 150, -400);
        let ref_end = compute_ref_end(&r1, &r2, 100);
        assert_eq!(ref_end, 500);
    }

    #[test]
    fn test_compute_ref_end_uses_negative_tlen_abs() {
        let r1 = partial(300, 150, -420);
        let r2 = partial(100, 150, 420);
        let ref_end = compute_ref_end(&r1, &r2, 100);
        assert_eq!(ref_end, 520);
    }

    #[test]
    fn test_compute_ref_end_falls_back_when_tlen_zero() {
        let r1 = partial(100, 150, 0);
        let r2 = partial(260, 150, 0);
        let ref_end = compute_ref_end(&r1, &r2, 100);
        assert_eq!(ref_end, 410); // max(100+150, 260+150)
    }

    #[test]
    fn test_build_pair_from_partials_sets_expected_fields() {
        let r1 = partial(300, 150, -420);
        let mut r2 = partial(100, 150, 420);
        r2.seq.fill(b'T');
        let pair = build_pair_from_partials("read1".to_string(), r1, r2, "chr1");

        assert_eq!(pair.name, "read1");
        assert_eq!(pair.chrom, "chr1");
        assert_eq!(pair.ref_start, 100);
        assert_eq!(pair.ref_end, 520);
        assert_eq!(pair.insert_size, 420);
        assert_eq!(pair.seq1.len(), 150);
        assert_eq!(pair.seq2.len(), 150);
    }

    // --- M14: missing base qualities must not silently become spaces ---

    #[test]
    fn test_quality_is_missing_detects_all_0xff() {
        // BAM/CRAM encode a SAM `*` (no quality stored) as every per-base
        // raw byte set to 0xFF. That's the exact input that used to produce
        // `s.wrapping_add(33) == 32` (space) for every base.
        assert!(quality_is_missing(&[0xff; 10]));
    }

    #[test]
    fn test_quality_is_missing_false_for_real_quality() {
        // Real Phred scores never reach 255 (max encodable is 93), so any
        // record with real quality data must not be flagged as missing.
        assert!(!quality_is_missing(&[30, 30, 40, 2, 0]));
    }

    #[test]
    fn test_quality_is_missing_false_for_empty() {
        // An empty quality slice (e.g. a zero-length read) isn't the
        // "missing" sentinel; don't misclassify it.
        assert!(!quality_is_missing(&[]));
    }

    fn record_buf_with_quality(qual: Vec<u8>) -> noodles::sam::alignment::RecordBuf {
        let seq_len = qual.len();
        noodles::sam::alignment::RecordBuf::builder()
            .set_alignment_start(safe_noodles_position(100))
            .set_mate_alignment_start(safe_noodles_position(500))
            .set_template_length(400)
            .set_sequence(noodles::sam::alignment::record_buf::Sequence::from(vec![
                b'A';
                seq_len
            ]))
            .set_quality_scores(noodles::sam::alignment::record_buf::QualityScores::from(
                qual,
            ))
            .build()
    }

    #[test]
    fn test_parse_partial_from_record_buf_skips_missing_quality_and_counts_it() {
        let buf = record_buf_with_quality(vec![0xff; 10]);
        let flags = noodles::sam::alignment::record::Flags::empty();
        let mut missing_qual_count = 0usize;

        let result = parse_partial_from_record_buf(&buf, &flags, &mut missing_qual_count);

        assert!(result.is_none(), "record with all-0xFF quality must be skipped");
        assert_eq!(missing_qual_count, 1);
    }

    #[test]
    fn test_parse_partial_from_record_buf_keeps_real_quality() {
        let buf = record_buf_with_quality(vec![30; 10]);
        let flags = noodles::sam::alignment::record::Flags::empty();
        let mut missing_qual_count = 0usize;

        let result = parse_partial_from_record_buf(&buf, &flags, &mut missing_qual_count);

        let partial = result.expect("record with real quality must be kept");
        assert_eq!(partial.qual, vec![30u8 + 33; 10]);
        assert_eq!(missing_qual_count, 0);
    }
}
