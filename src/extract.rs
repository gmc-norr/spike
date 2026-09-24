//! BAM/CRAM read pair extraction and ReadPool building.
//!
//! Extracts paired-end reads from a BAM/CRAM region, stores them in their original
//! sequencing orientation (FASTQ order).

use anyhow::{Context, Result};
use std::collections::{BTreeSet, HashMap, HashSet};

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

/// What extraction found in one region.
pub struct Extraction {
    /// Complete pairs, usable as donor material.
    pub pairs: Vec<ReadPair>,
    /// Names of the pairs dropped because a mate's stored quality is unusable
    /// (M14). spike cannot replace them — it has no quality string to write —
    /// but merge.sh still removes them from the merged BAM, so they cannot sit
    /// inside a simulated event as un-suppressible reference support.
    pub unusable_qual_names: BTreeSet<String>,
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
) -> Result<Extraction> {
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
) -> Result<Extraction> {
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
    let mut tally = UnusableQualTally::default();

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
            let partial = match parse_partial_from_bam_record(&record, &flags, &name, &mut tally)
            {
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
                        parse_partial_from_bam_record(&record, &flags, &name, &mut tally)
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
                        parse_partial_from_bam_record(&record, &flags, &name, &mut tally)
                    else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            }
        }
    }

    let unmatched = read1_map.len() + read2_map.len();
    log_extraction_result(&tally, unmatched, pairs.len(), chrom, start, end);
    Ok(Extraction {
        pairs,
        unusable_qual_names: tally.pair_names(),
    })
}

/// Pick the CRAM index entries whose slice can hold a record overlapping
/// `interval` on `reference_sequence_id`.
///
/// A `.crai` entry covers one slice, and the slice's records all lie within
/// `alignment_start .. alignment_start + alignment_span`, so a slice outside
/// the interval cannot hold a record inside it. Entries that name no
/// reference sequence are kept: `-1` is noodles' `UNMAPPED`, a slice of
/// unplaced reads, and it carries no coordinates to judge. (A container
/// holding several contigs is *not* this case — htslib writes one `.crai`
/// line per contig, each with that contig's own start and span, all pointing
/// at the same offset.) Keeping them changes nothing either way:
/// `Query::read_next_container` skips every entry whose reference id is not
/// the queried one.
fn select_crai_entries(
    index: &[noodles::cram::crai::Record],
    reference_sequence_id: usize,
    interval: noodles::core::region::Interval,
) -> Vec<noodles::cram::crai::Record> {
    index
        .iter()
        .filter(|entry| {
            let Some(id) = entry.reference_sequence_id() else {
                return true;
            };
            if id != reference_sequence_id {
                return false;
            }
            let Some(start) = entry.alignment_start() else {
                return true;
            };
            let last = usize::from(start).saturating_add(entry.alignment_span().saturating_sub(1));
            let end = noodles::core::Position::new(last).unwrap_or(start);
            interval.intersects((start..=end).into())
        })
        .cloned()
        .collect()
}

/// Open an indexed CRAM reader that only has to visit `region`'s containers.
///
/// noodles-cram 0.74's `Query::read_next_container` compares an index entry's
/// reference id and nothing else, so a query seeks to and fully decodes every
/// container on the chromosome (M15). The index it queries is the one the
/// reader was built with, so pruning that index to the slices that can overlap
/// `region` is enough to skip the rest. Record-level filtering is unchanged:
/// noodles still returns only the records intersecting the region.
fn open_cram_reader_for_region(
    cram_path: &str,
    repository: &noodles::fasta::Repository,
    region: &noodles::core::Region,
) -> Result<(
    noodles::cram::io::IndexedReader<std::fs::File>,
    noodles::sam::Header,
)> {
    // The header is what maps the region's name to a reference id, so it has
    // to be read before the index can be pruned.
    let mut reader = noodles::cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .build_from_path(cram_path)?;
    let header = reader.read_header()?;

    let Some(reference_sequence_id) = header.reference_sequences().get_index_of(region.name())
    else {
        // Unknown contig: leave the full index in place and let query() raise
        // its own "invalid reference sequence name".
        return Ok((reader, header));
    };
    let index = select_crai_entries(reader.index(), reference_sequence_id, region.interval());

    let mut reader = noodles::cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository.clone())
        .set_index(index)
        .build_from_path(cram_path)?;
    let header = reader.read_header()?;

    Ok((reader, header))
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
) -> Result<Extraction> {
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
    let mut tally = UnusableQualTally::default();

    {
        let start_pos = safe_noodles_position(start + 1);
        let end_pos = safe_noodles_position(end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) = open_cram_reader_for_region(cram_path, &repository, &region)
            .with_context(|| format!("failed to open CRAM: {}", cram_path))?;
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
            let partial = match parse_partial_from_record_buf(&buf, &flags, &name, &mut tally) {
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
        let wider_start = start.saturating_sub(wider_padding);
        let wider_end = end.saturating_add(wider_padding);

        let start_pos = safe_noodles_position(wider_start + 1);
        let end_pos = safe_noodles_position(wider_end);
        let region = noodles::core::Region::new(chrom, start_pos..=end_pos);
        let (mut reader, header) = open_cram_reader_for_region(cram_path, &repository, &region)
            .with_context(|| format!("failed to open CRAM for pass 2: {}", cram_path))?;
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
                        parse_partial_from_record_buf(&buf, &flags, &name, &mut tally)
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
                        parse_partial_from_record_buf(&buf, &flags, &name, &mut tally)
                    else {
                        continue;
                    };
                    pairs.push(build_pair_from_partials(name, read1, read2, chrom));
                }
            }
        }
    }

    let unmatched = read1_map.len() + read2_map.len();
    log_extraction_result(&tally, unmatched, pairs.len(), chrom, start, end);
    Ok(Extraction {
        pairs,
        unusable_qual_names: tally.pair_names(),
    })
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

/// Drop read pairs already seen under the same name, keeping the first.
///
/// Extraction windows can share reads -- a `--region` used by both sides of a
/// fusion, or a fragment whose mate lands in a neighbouring window -- and the
/// same fragment must not enter the pool twice.
pub fn dedup_pairs_by_name(pairs: &mut Vec<ReadPair>) {
    let mut seen: HashSet<String> = HashSet::with_capacity(pairs.len());
    pairs.retain(|p| seen.insert(p.name.clone()));
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

/// Highest raw (pre-Phred+33) quality score SAM allows: Q93, which becomes
/// `~` (126) once shifted into Phred+33.
const MAX_RAW_QUAL: u8 = 93;

/// Whether a record stores no quality at all (SAM `*`).
///
/// SAM's QUAL field is all-or-nothing per record — `*`, or a full string the
/// same length as SEQ — but the two containers encode `*` differently and
/// spike sees both:
///   * BAM hands over the raw array with every per-base byte set to 0xFF.
///   * CRAM never does: noodles normalises that array to an *empty* slice
///     before spike sees it (`read_quality_scores_stored_as_array` ends with
///     `if buf.iter().all(|&n| n == MISSING) { buf.clear(); }`), so on the
///     CRAM path an empty quality beside a non-empty sequence *is* the `*`.
///
/// `seq_len` is what separates that from a genuinely empty record.
fn quality_is_missing(raw_qual: &[u8], seq_len: usize) -> bool {
    if raw_qual.is_empty() {
        return seq_len > 0;
    }
    raw_qual.iter().all(|&b| b == 0xff)
}

/// The first raw quality byte this record may not carry, if any.
///
/// SAM caps Phred at Q93, so a raw byte above that is malformed — a partially
/// 0xFF record, or corruption. `wrapping_add(33)` turns it into a byte outside
/// the printable Phred+33 range, which poisons the learned quality profile
/// even when the record is only ever used as a suppressed donor and never
/// written. Enforce the bound here rather than assume the SAM spec holds.
fn quality_out_of_range(raw_qual: &[u8]) -> Option<u8> {
    raw_qual.iter().copied().find(|&b| b > MAX_RAW_QUAL)
}

/// Records dropped because their stored quality is unusable, keyed by
/// (read name, first-segment flag).
///
/// Keyed rather than counted because pass 2 queries a superset of pass 1's
/// window and re-parses any record whose mate pass 1 kept, so the same record
/// is offered to the tally twice and a plain counter double-counts it.
#[derive(Default)]
struct UnusableQualTally {
    /// Records with no quality stored at all (SAM `*`).
    missing: HashSet<(String, bool)>,
    /// Records whose quality carries a raw byte above Q93.
    out_of_range: HashSet<(String, bool)>,
}

impl UnusableQualTally {
    fn add_missing(&mut self, name: &str, is_first_segment: bool) {
        self.missing.insert((name.to_string(), is_first_segment));
    }

    fn add_out_of_range(&mut self, name: &str, is_first_segment: bool) {
        self.out_of_range.insert((name.to_string(), is_first_segment));
    }

    /// Names of the pairs these records belong to, deduplicated across mates.
    fn pair_names(&self) -> BTreeSet<String> {
        self.missing
            .iter()
            .chain(self.out_of_range.iter())
            .map(|(name, _)| name.clone())
            .collect()
    }
}

/// Parse sequence/quality/position fields from a BAM record into a partial read.
///
/// Records whose stored quality is unusable — missing entirely (SAM `*`), or
/// carrying a raw byte above Q93 — are skipped and recorded in `tally`.
fn parse_partial_from_bam_record(
    record: &noodles::bam::Record,
    flags: &noodles::sam::alignment::record::Flags,
    name: &str,
    tally: &mut UnusableQualTally,
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
    let mut seq: Vec<u8> = record.sequence().iter().collect();
    if quality_is_missing(raw_qual, seq.len()) {
        tally.add_missing(name, flags.is_first_segment());
        return None;
    }
    if quality_out_of_range(raw_qual).is_some() {
        tally.add_out_of_range(name, flags.is_first_segment());
        return None;
    }
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
/// Records whose stored quality is unusable — missing entirely (SAM `*`), or
/// carrying a raw byte above Q93 — are skipped and recorded in `tally`.
fn parse_partial_from_record_buf(
    buf: &noodles::sam::alignment::RecordBuf,
    flags: &noodles::sam::alignment::record::Flags,
    name: &str,
    tally: &mut UnusableQualTally,
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
    let mut seq: Vec<u8> = buf.sequence().as_ref().to_vec();
    if quality_is_missing(raw_qual, seq.len()) {
        tally.add_missing(name, flags.is_first_segment());
        return None;
    }
    if quality_out_of_range(raw_qual).is_some() {
        tally.add_out_of_range(name, flags.is_first_segment());
        return None;
    }
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
    tally: &UnusableQualTally,
    unmatched: usize,
    pair_count: usize,
    chrom: &str,
    start: u64,
    end: u64,
) {
    // M14: a record whose stored quality is unusable is dropped rather than
    // silently turned into space characters, which used to poison the learned
    // quality model. Surface the counts so a user can tell this is happening
    // to their donor BAM. The window named is the one being simulated; pass 2
    // can add a few records from just outside it.
    if !tally.missing.is_empty() {
        log::warn!(
            "{} record(s) considered for the {}:{}-{} donor pool had no quality \
             scores (SAM '*'); their pairs were dropped from the pool and merge.sh \
             removes them from the merged BAM",
            tally.missing.len(),
            chrom,
            start,
            end,
        );
    }
    if !tally.out_of_range.is_empty() {
        log::warn!(
            "{} record(s) considered for the {}:{}-{} donor pool carry a raw \
             quality score above Q{}; their pairs were dropped from the pool and \
             merge.sh removes them from the merged BAM",
            tally.out_of_range.len(),
            chrom,
            start,
            end,
            MAX_RAW_QUAL,
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

    #[test]
    fn test_dedup_pairs_by_name_keeps_one_copy_of_a_pair_in_two_windows() {
        // A pair that sits in two extraction windows -- the --region shared by
        // both sides of a fusion (M8), or a fragment straddling the boundary
        // between neighbouring windows -- is extracted once per window and
        // must still enter the pool once.
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

        let mut pairs = vec![pair("a", 10), pair("b", 20), pair("a", 10)];
        dedup_pairs_by_name(&mut pairs);

        let names: Vec<String> = pairs.into_iter().map(|p| p.name).collect();
        assert_eq!(names, vec!["a", "b"]);
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

    fn crai_record(
        reference_sequence_id: Option<usize>,
        start: usize,
        span: usize,
        offset: u64,
    ) -> noodles::cram::crai::Record {
        noodles::cram::crai::Record::new(
            reference_sequence_id,
            noodles::core::Position::new(start),
            span,
            offset,
            0,
            0,
        )
    }

    fn query_interval(start: usize, end: usize) -> noodles::core::region::Interval {
        (noodles::core::Position::new(start).unwrap()..=noodles::core::Position::new(end).unwrap())
            .into()
    }

    #[test]
    fn test_select_crai_entries_skips_slices_outside_the_region() {
        // Five 1000 bp slices on reference 19 plus one on another contig. A
        // query of 2500-3200 can only be served by the 2001-3000 and 3001-4000
        // slices; noodles-cram 0.74 seeks to and decodes all five, because
        // `Query::read_next_container` compares only the reference id (M15).
        let index = vec![
            crai_record(Some(19), 1, 1000, 100),
            crai_record(Some(19), 1001, 1000, 200),
            crai_record(Some(19), 2001, 1000, 300),
            crai_record(Some(19), 3001, 1000, 400),
            crai_record(Some(19), 4001, 1000, 500),
            crai_record(Some(7), 2001, 1000, 600),
        ];
        let selected = select_crai_entries(&index, 19, query_interval(2500, 3200));
        let offsets: Vec<u64> = selected.iter().map(|r| r.offset()).collect();
        assert_eq!(offsets, vec![300, 400]);
    }

    #[test]
    fn test_select_crai_entries_keeps_the_slice_ending_on_the_first_base() {
        // A slice that ends exactly on the first queried base still holds
        // records inside the region; one ending a base earlier cannot.
        let index = vec![
            crai_record(Some(19), 1, 2499, 100),
            crai_record(Some(19), 1, 2500, 200),
        ];
        let selected = select_crai_entries(&index, 19, query_interval(2500, 3200));
        let offsets: Vec<u64> = selected.iter().map(|r| r.offset()).collect();
        assert_eq!(offsets, vec![200]);
    }

    #[test]
    fn test_select_crai_entries_keeps_the_slice_starting_on_the_last_base() {
        // The right edge, pinned at one base: a slice whose first base is the
        // region's last base can hold a record inside the region.
        let index = vec![crai_record(Some(19), 3200, 1000, 100)];
        let selected = select_crai_entries(&index, 19, query_interval(2500, 3200));
        let offsets: Vec<u64> = selected.iter().map(|r| r.offset()).collect();
        assert_eq!(offsets, vec![100]);
    }

    #[test]
    fn test_select_crai_entries_drops_the_slice_starting_one_base_past_the_region() {
        // The same edge from the other side: a slice starting one base after
        // the region ends cannot hold a record inside it, however long it is.
        let index = vec![crai_record(Some(19), 3201, 1000, 200)];
        let selected = select_crai_entries(&index, 19, query_interval(2500, 3200));
        let offsets: Vec<u64> = selected.iter().map(|r| r.offset()).collect();
        assert_eq!(offsets, Vec::<u64>::new());
    }

    #[test]
    fn test_select_crai_entries_keeps_entries_without_a_reference_id() {
        // A `.crai` entry with no reference id is `-1`, a slice of unplaced
        // reads; it has no coordinates on this reference, so there is nothing
        // to judge it by and it stays. The second entry (same contig, well
        // outside the region) is what makes this a test rather than a
        // restatement of "keep everything".
        let index = vec![
            crai_record(None, 1, 0, 100),
            crai_record(Some(19), 4001, 1000, 200),
        ];
        let selected = select_crai_entries(&index, 19, query_interval(2500, 3200));
        let offsets: Vec<u64> = selected.iter().map(|r| r.offset()).collect();
        assert_eq!(offsets, vec![100]);
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

    // --- M14: unusable base qualities must not silently become spaces ---

    #[test]
    fn test_quality_is_missing_detects_all_0xff() {
        // BAM encodes a SAM `*` (no quality stored) as every per-base raw
        // byte set to 0xFF. That's the exact input that used to produce
        // `s.wrapping_add(33) == 32` (space) for every base.
        assert!(quality_is_missing(&[0xff; 10], 10));
    }

    #[test]
    fn test_quality_is_missing_false_for_real_quality() {
        // Real Phred scores never reach 255 (max encodable is 93), so any
        // record with real quality data must not be flagged as missing.
        assert!(!quality_is_missing(&[30, 30, 40, 2, 0], 5));
    }

    #[test]
    fn test_quality_is_missing_true_for_empty_qual_with_sequence() {
        // CRAM's encoding of `*`: noodles clears the all-0xFF buffer before
        // spike sees it, so a CRAM donor never reaches the all-0xFF branch.
        // An empty quality beside a 151 bp sequence is a missing quality --
        // exactly the record that used to reach the FASTQ writer as a 151 bp
        // SEQ line with a zero-length QUAL line.
        assert!(quality_is_missing(&[], 151));
    }

    #[test]
    fn test_quality_is_missing_false_for_empty_record() {
        // A record with neither sequence nor quality is not the `*` case.
        assert!(!quality_is_missing(&[], 0));
    }

    #[test]
    fn test_quality_out_of_range_flags_byte_above_q93() {
        // A partially-0xFF or otherwise malformed record: SAM caps Phred at
        // Q93, so anything above that would wrap past `~` on +33.
        assert_eq!(quality_out_of_range(&[30, 30, 200, 30]), Some(200));
        assert_eq!(quality_out_of_range(&[30, 94, 30]), Some(94));
    }

    #[test]
    fn test_quality_out_of_range_accepts_q0_through_q93() {
        assert_eq!(quality_out_of_range(&[0, 1, 40, 93]), None);
    }

    #[test]
    fn test_unusable_qual_tally_counts_each_record_once() {
        // Pass 2 queries a superset of pass 1's window and re-parses any
        // record whose mate pass 1 kept, so the same record reaches the
        // tally twice and must still be counted once.
        let mut tally = UnusableQualTally::default();
        tally.add_missing("readA", true);
        tally.add_missing("readA", true);
        assert_eq!(tally.missing.len(), 1, "the same record must be tallied once");

        // The two mates of one pair are two distinct records, but one name.
        tally.add_missing("readA", false);
        assert_eq!(tally.missing.len(), 2);
        assert_eq!(
            tally.pair_names().into_iter().collect::<Vec<_>>(),
            vec!["readA".to_string()],
        );
    }

    #[test]
    fn test_unusable_qual_tally_pair_names_covers_both_reasons() {
        let mut tally = UnusableQualTally::default();
        tally.add_missing("noQual", true);
        tally.add_out_of_range("badByte", false);
        assert_eq!(
            tally.pair_names().into_iter().collect::<Vec<_>>(),
            vec!["badByte".to_string(), "noQual".to_string()],
        );
    }

    fn record_buf_with_seq_and_qual(
        seq_len: usize,
        qual: Vec<u8>,
    ) -> noodles::sam::alignment::RecordBuf {
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

    fn record_buf_with_quality(qual: Vec<u8>) -> noodles::sam::alignment::RecordBuf {
        record_buf_with_seq_and_qual(qual.len(), qual)
    }

    #[test]
    fn test_parse_partial_from_record_buf_skips_missing_quality_and_counts_it() {
        let buf = record_buf_with_quality(vec![0xff; 10]);
        let flags = noodles::sam::alignment::record::Flags::empty();
        let mut tally = UnusableQualTally::default();

        let result = parse_partial_from_record_buf(&buf, &flags, "allFF", &mut tally);

        assert!(result.is_none(), "record with all-0xFF quality must be skipped");
        assert_eq!(tally.missing.len(), 1);
    }

    #[test]
    fn test_parse_partial_from_record_buf_skips_empty_cram_quality() {
        // The CRAM shape of a `*` record: 151 bases, zero quality bytes.
        let buf = record_buf_with_seq_and_qual(151, Vec::new());
        let flags = noodles::sam::alignment::record::Flags::empty();
        let mut tally = UnusableQualTally::default();

        let result = parse_partial_from_record_buf(&buf, &flags, "cramStar", &mut tally);

        assert!(
            result.is_none(),
            "CRAM record with an empty quality buffer must be skipped"
        );
        assert!(tally.pair_names().contains("cramStar"));
    }

    #[test]
    fn test_parse_partial_from_record_buf_skips_out_of_range_quality_byte() {
        // Not expressible in SAM text, but nothing stops a malformed BAM/CRAM
        // from storing it, and it poisons the quality profile if it survives.
        let mut qual = vec![30u8; 10];
        qual[4] = 200;
        let buf = record_buf_with_quality(qual);
        let flags = noodles::sam::alignment::record::Flags::empty();
        let mut tally = UnusableQualTally::default();

        let result = parse_partial_from_record_buf(&buf, &flags, "badByte", &mut tally);

        assert!(
            result.is_none(),
            "record with a raw quality byte above Q93 must be skipped"
        );
        assert_eq!(tally.out_of_range.len(), 1);
        assert!(tally.pair_names().contains("badByte"));
    }

    #[test]
    fn test_parse_partial_from_record_buf_keeps_real_quality() {
        let buf = record_buf_with_quality(vec![30; 10]);
        let flags = noodles::sam::alignment::record::Flags::empty();
        let mut tally = UnusableQualTally::default();

        let result = parse_partial_from_record_buf(&buf, &flags, "good", &mut tally);

        let partial = result.expect("record with real quality must be kept");
        assert_eq!(partial.qual, vec![30u8 + 33; 10]);
        assert!(tally.missing.is_empty());
        assert!(tally.out_of_range.is_empty());
    }
}
