//! BAM read pair extraction and ReadPool building.
//!
//! Extracts paired-end reads from a BAM region, stores them in their original
//! sequencing orientation (FASTQ order).

use anyhow::{Context, Result};
use std::collections::HashMap;

use crate::stats::FragmentDist;
use crate::types::{ReadPair, ReadPool};

/// Extract all properly-paired read pairs from a genomic region.
///
/// Two-pass approach:
/// 1. First pass: collect read1 records (is_first_segment) with seq, qual, pos, tlen
/// 2. Second pass: collect read2 records, match by name
///
/// Reads aligned in reverse complement are reverse-complemented back to FASTQ orientation.
pub fn extract_read_pairs(
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

        let start_pos =
            noodles::core::Position::try_from(1_usize.max(start as usize + 1)).unwrap();
        let end_pos = noodles::core::Position::try_from(end as usize).unwrap();
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();

            if !passes_filters(&flags, min_mapq, &record) {
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

        let start_pos =
            noodles::core::Position::try_from(1_usize.max(wider_start as usize + 1)).unwrap();
        let end_pos = noodles::core::Position::try_from(wider_end as usize).unwrap();
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let query = reader.query(&header, &region)?;

        for rec_result in query {
            let record = rec_result?;
            let flags = record.flags();

            if !passes_filters(&flags, min_mapq, &record) {
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
                let ref_end = if partial.tlen > 0 {
                    ref_start + partial.tlen as u64
                } else if partial.tlen < 0 {
                    ref_start + (-partial.tlen) as u64
                } else {
                    // tlen == 0: estimate from positions + read length.
                    let r1_end = partial.pos + partial.seq1.len() as u64;
                    let r2_end = partial.mate_pos + seq.len() as u64;
                    r1_end.max(r2_end)
                };

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

    let unmatched = read1_map.len();
    if unmatched > 0 {
        log::debug!(
            "{} read1 records had no matching read2 in the query region",
            unmatched,
        );
    }

    log::info!(
        "Extracted {} complete read pairs from {}:{}-{}",
        pairs.len(),
        chrom,
        start,
        end,
    );

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

fn passes_filters(
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

/// Reverse-complement a DNA sequence in place.
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
}
