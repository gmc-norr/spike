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
            let partial = match parse_partial_from_bam_record(&record, &flags) {
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
                    let Some(read1) = parse_partial_from_bam_record(&record, &flags) else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            } else {
                if read2_map.contains_key(&name) {
                    continue;
                }
                if let Some(read1) = read1_map.remove(&name) {
                    let Some(read2) = parse_partial_from_bam_record(&record, &flags) else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            }
        }
    }

    let unmatched = read1_map.len() + read2_map.len();
    log_extraction_result(unmatched, pairs.len(), chrom, start, end);
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
            let partial = match parse_partial_from_record_buf(&buf, &flags) {
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
                    let Some(read1) = parse_partial_from_record_buf(&buf, &flags) else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            } else {
                if read2_map.contains_key(&name) {
                    continue;
                }
                if let Some(read1) = read1_map.remove(&name) {
                    let Some(read2) = parse_partial_from_record_buf(&buf, &flags) else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            }
        }
    }

    let unmatched = read1_map.len() + read2_map.len();
    log_extraction_result(unmatched, pairs.len(), chrom, start, end);
    Ok(pairs)
}

/// Build a ReadPool from extracted read pairs.
///
/// Sorts pairs by ref_start.
pub fn build_read_pool(pairs: Vec<ReadPair>, frag_dist: FragmentDist) -> ReadPool {
    let mut pairs = pairs;
    pairs.sort_by_key(|p| p.ref_start);

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

/// Parse sequence/quality/position fields from a BAM record into a partial read.
fn parse_partial_from_bam_record(
    record: &noodles::bam::Record,
    flags: &noodles::sam::alignment::record::Flags,
) -> Option<PartialRead> {
    let pos = match record.alignment_start() {
        Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
        _ => return None,
    };
    if record.mate_alignment_start().is_none() {
        return None;
    }
    let tlen = record.template_length();

    let mut seq: Vec<u8> = record.sequence().iter().collect();
    let mut qual: Vec<u8> = record
        .quality_scores()
        .as_ref()
        .iter()
        .map(|s| s.wrapping_add(33))
        .collect();

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
fn parse_partial_from_record_buf(
    buf: &noodles::sam::alignment::RecordBuf,
    flags: &noodles::sam::alignment::record::Flags,
) -> Option<PartialRead> {
    let pos = match buf.alignment_start() {
        Some(p) => usize::from(p).saturating_sub(1) as u64,
        None => return None,
    };
    if buf.mate_alignment_start().is_none() {
        return None;
    }
    let tlen = buf.template_length();

    let mut seq: Vec<u8> = buf.sequence().as_ref().to_vec();
    let mut qual: Vec<u8> = buf
        .quality_scores()
        .as_ref()
        .iter()
        .map(|s| s.wrapping_add(33))
        .collect();

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
fn log_extraction_result(unmatched: usize, pair_count: usize, chrom: &str, start: u64, end: u64) {
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

/// Convert a 0-based half-open end coordinate to a noodles 1-based Position.
/// Noodles positions must be >= 1, so we clamp.
pub fn safe_noodles_position(pos: u64) -> noodles::core::Position {
    let pos_usize = usize::try_from(pos).unwrap_or(usize::MAX).max(1);
    noodles::core::Position::try_from(pos_usize).expect("position must be >= 1")
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
}
