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

    // Pass 1: collect read1 records.
    let mut read1_map: HashMap<String, PartialPair> = HashMap::new();

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
            if !flags.is_first_segment() {
                continue;
            }
            if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
                continue;
            }

            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };

            let pos = match record.alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => continue,
            };

            let mate_pos = match record.mate_alignment_start() {
                Some(Ok(p)) => usize::from(p).saturating_sub(1) as u64,
                _ => continue,
            };

            let tlen = record.template_length();
            let is_reverse = flags.is_reverse_complemented();

            let mut seq: Vec<u8> = record.sequence().iter().collect();
            let mut qual: Vec<u8> = record
                .quality_scores()
                .as_ref()
                .iter()
                .map(|s| s.wrapping_add(33))
                .collect();

            // Reverse-complement back to FASTQ orientation if needed.
            if is_reverse {
                reverse_complement(&mut seq);
                qual.reverse();
            }

            read1_map.insert(
                name,
                PartialPair {
                    seq1: seq,
                    qual1: qual,
                    pos,
                    mate_pos,
                    tlen,
                },
            );
        }
    }

    log::info!("Pass 1: collected {} read1 records", read1_map.len());

    // Pass 2: collect read2 records, match by name.
    let mut pairs: Vec<ReadPair> = Vec::with_capacity(read1_map.len());

    {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(bam_path)
            .with_context(|| format!("failed to open BAM for pass 2: {}", bam_path))?;
        let header = reader.read_header()?;

        // Query a wider region to capture mates that may be slightly outside.
        let wider_start = start.saturating_sub(1000);
        let wider_end = end.saturating_add(1000);

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
            if !flags.is_last_segment() {
                continue;
            }

            let name = match record.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };

            if let Some(partial) = read1_map.remove(&name) {
                let is_reverse = flags.is_reverse_complemented();

                let mut seq: Vec<u8> = record.sequence().iter().collect();
                let mut qual: Vec<u8> = record
                    .quality_scores()
                    .as_ref()
                    .iter()
                    .map(|s| s.wrapping_add(33))
                    .collect();

                if is_reverse {
                    reverse_complement(&mut seq);
                    qual.reverse();
                }

                let ref_start = partial.pos.min(partial.mate_pos);
                let ref_end = compute_ref_end(&partial, ref_start, seq.len());

                pairs.push(ReadPair {
                    name,
                    seq1: partial.seq1,
                    qual1: partial.qual1,
                    seq2: seq,
                    qual2: qual,
                    ref_start,
                    ref_end,
                    insert_size: partial.tlen as i64,
                    chrom: chrom.to_string(),
                });
            }
        }
    }

    log_extraction_result(&read1_map, pairs.len(), chrom, start, end);
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

    // Pass 1: collect read1 records.
    let mut read1_map: HashMap<String, PartialPair> = HashMap::new();

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
            if !flags.is_first_segment() {
                continue;
            }
            if !flags.is_properly_segmented() || flags.is_mate_unmapped() {
                continue;
            }

            let name = match buf.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };

            let pos = match buf.alignment_start() {
                Some(p) => usize::from(p).saturating_sub(1) as u64,
                None => continue,
            };

            let mate_pos = match buf.mate_alignment_start() {
                Some(p) => usize::from(p).saturating_sub(1) as u64,
                None => continue,
            };

            let tlen = buf.template_length();
            let is_reverse = flags.is_reverse_complemented();

            let mut seq: Vec<u8> = buf.sequence().as_ref().to_vec();
            let mut qual: Vec<u8> = buf
                .quality_scores()
                .as_ref()
                .iter()
                .map(|s| s.wrapping_add(33))
                .collect();

            if is_reverse {
                reverse_complement(&mut seq);
                qual.reverse();
            }

            read1_map.insert(
                name,
                PartialPair {
                    seq1: seq,
                    qual1: qual,
                    pos,
                    mate_pos,
                    tlen,
                },
            );
        }
    }

    log::info!("Pass 1 (CRAM): collected {} read1 records", read1_map.len());

    // Pass 2: collect read2 records, match by name.
    let mut pairs: Vec<ReadPair> = Vec::with_capacity(read1_map.len());

    {
        let mut reader = noodles::cram::io::indexed_reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(cram_path)
            .with_context(|| format!("failed to open CRAM for pass 2: {}", cram_path))?;
        let header = reader.read_header()?;

        let wider_start = start.saturating_sub(1000);
        let wider_end = end.saturating_add(1000);

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
            if !flags.is_last_segment() {
                continue;
            }

            let name = match buf.name() {
                Some(n) => String::from_utf8_lossy(n.as_ref()).into_owned(),
                None => continue,
            };

            if let Some(partial) = read1_map.remove(&name) {
                let is_reverse = flags.is_reverse_complemented();

                let mut seq: Vec<u8> = buf.sequence().as_ref().to_vec();
                let mut qual: Vec<u8> = buf
                    .quality_scores()
                    .as_ref()
                    .iter()
                    .map(|s| s.wrapping_add(33))
                    .collect();

                if is_reverse {
                    reverse_complement(&mut seq);
                    qual.reverse();
                }

                let ref_start = partial.pos.min(partial.mate_pos);
                let ref_end = compute_ref_end(&partial, ref_start, seq.len());

                pairs.push(ReadPair {
                    name,
                    seq1: partial.seq1,
                    qual1: partial.qual1,
                    seq2: seq,
                    qual2: qual,
                    ref_start,
                    ref_end,
                    insert_size: partial.tlen as i64,
                    chrom: chrom.to_string(),
                });
            }
        }
    }

    log_extraction_result(&read1_map, pairs.len(), chrom, start, end);
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

struct PartialPair {
    seq1: Vec<u8>,
    qual1: Vec<u8>,
    pos: u64,
    mate_pos: u64,
    tlen: i32,
}

/// Compute the fragment end position from template length and read positions.
fn compute_ref_end(partial: &PartialPair, ref_start: u64, r2_seq_len: usize) -> u64 {
    if partial.tlen > 0 {
        ref_start + partial.tlen as u64
    } else if partial.tlen < 0 {
        ref_start + (partial.tlen as i64).unsigned_abs()
    } else {
        // tlen == 0: estimate from positions + read length.
        let r1_end = partial.pos + partial.seq1.len() as u64;
        let r2_end = partial.mate_pos + r2_seq_len as u64;
        r1_end.max(r2_end)
    }
}

/// Log extraction results (shared between BAM and CRAM paths).
fn log_extraction_result(
    read1_map: &HashMap<String, PartialPair>,
    pair_count: usize,
    chrom: &str,
    start: u64,
    end: u64,
) {
    let unmatched = read1_map.len();
    if unmatched > 0 {
        log::debug!(
            "{} read1 records had no matching read2 in the query region",
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
fn passes_filters_buf(
    min_mapq: u8,
    buf: &noodles::sam::alignment::RecordBuf,
) -> bool {
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
    noodles::core::Position::try_from((pos as usize).max(1))
        .expect("position must be >= 1")
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
}
