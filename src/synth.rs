//! Synthetic read generation from reference sequence + learned quality profile.
//!
//! Instead of cloning real reads (which produces exact duplicates flagged by dedup tools),
//! this module generates independent synthetic reads that have realistic quality scores
//! and correlated sequencing errors.
//!
//! The approach:
//! 1. Learn a `QualityProfile` from real reads: per-cycle quality distributions for R1 and R2
//! 2. Generate synthetic reads by sampling quality from the profile, reading reference bases,
//!    and introducing errors at the rate implied by the sampled quality score

use std::collections::HashMap;

use rand::rngs::StdRng;
use rand::Rng;

use crate::extract::reverse_complement;
use crate::haplotype::VariantHaplotype;
use crate::stats::FragmentDist;
use crate::types::{ReadPair, ReadPool};
use crate::reference::SharedReference;

/// Minimum observations in a (cycle, base) bin before we trust it.
/// Below this threshold, fall back to the cycle-only distribution.
const MIN_BASE_OBS: usize = 30;

/// Empirical per-cycle quality score distributions learned from real reads.
///
/// Two levels of conditioning:
/// 1. **Base-conditioned**: `(read_number, cycle, sequenced_base)` → quality distribution.
///    Captures base-specific effects like the Illumina GG quality dip.
/// 2. **Cycle-only fallback**: `(read_number, cycle)` → quality distribution.
///    Used when a base-conditioned bin has too few observations.
///
/// Both levels store sorted Vec<u8> for O(1) CDF sampling.
pub struct QualityProfile {
    /// Base-conditioned quality for read1.
    /// `r1_base_quals[cycle][base_idx]` = sorted Vec<u8> of Phred+33 values.
    /// base_idx: A=0, C=1, G=2, T=3.
    r1_base_quals: Vec<[Vec<u8>; 4]>,

    /// Base-conditioned quality for read2.
    r2_base_quals: Vec<[Vec<u8>; 4]>,

    /// Cycle-only fallback for read1.
    /// `r1_cycle_quals[cycle]` = sorted Vec<u8> of Phred+33 values (all bases pooled).
    r1_cycle_quals: Vec<Vec<u8>>,

    /// Cycle-only fallback for read2.
    r2_cycle_quals: Vec<Vec<u8>>,
}

impl QualityProfile {
    /// Learn quality profile from extracted read pairs.
    ///
    /// Uses the read's own sequence bases (seq1/seq2) as the conditioning context.
    /// These are in FASTQ orientation (matching the quality scores) and are ~99%
    /// correct, so they faithfully represent the base the sequencer was reading.
    pub fn from_read_pairs(pairs: &[ReadPair], read_length: usize) -> Self {
        let mut r1_base: Vec<[Vec<u8>; 4]> =
            (0..read_length).map(|_| Default::default()).collect();
        let mut r2_base: Vec<[Vec<u8>; 4]> =
            (0..read_length).map(|_| Default::default()).collect();
        let mut r1_cycle: Vec<Vec<u8>> = (0..read_length).map(|_| Vec::new()).collect();
        let mut r2_cycle: Vec<Vec<u8>> = (0..read_length).map(|_| Vec::new()).collect();

        for pair in pairs {
            let r1_len = read_length.min(pair.qual1.len()).min(pair.seq1.len());
            for c in 0..r1_len {
                let q = pair.qual1[c];
                r1_cycle[c].push(q);
                if let Some(bi) = base_index(pair.seq1[c]) {
                    r1_base[c][bi].push(q);
                }
            }

            let r2_len = read_length.min(pair.qual2.len()).min(pair.seq2.len());
            for c in 0..r2_len {
                let q = pair.qual2[c];
                r2_cycle[c].push(q);
                if let Some(bi) = base_index(pair.seq2[c]) {
                    r2_base[c][bi].push(q);
                }
            }
        }

        // Sort all distributions for CDF sampling.
        for cycle in r1_base.iter_mut() {
            for bin in cycle.iter_mut() {
                bin.sort_unstable();
            }
        }
        for cycle in r2_base.iter_mut() {
            for bin in cycle.iter_mut() {
                bin.sort_unstable();
            }
        }
        for v in r1_cycle.iter_mut() {
            v.sort_unstable();
        }
        for v in r2_cycle.iter_mut() {
            v.sort_unstable();
        }

        // Log summary.
        let r1_mean_start = mean_qual(&r1_cycle[0]);
        let r1_mean_mid = mean_qual(&r1_cycle[read_length / 2]);
        let r1_mean_end = mean_qual(&r1_cycle[read_length.saturating_sub(1)]);
        let r2_mean_start = mean_qual(&r2_cycle[0]);
        let r2_mean_end = mean_qual(&r2_cycle[read_length.saturating_sub(1)]);

        // Count how many (cycle, base) bins have enough data.
        let base_bins_total = read_length * 4 * 2; // R1 + R2
        let base_bins_ok = r1_base
            .iter()
            .chain(r2_base.iter())
            .flat_map(|cycle| cycle.iter())
            .filter(|bin| bin.len() >= MIN_BASE_OBS)
            .count();

        log::info!(
            "Quality profile: {} pairs, {} cycles. R1 mean Q: start={:.1} mid={:.1} end={:.1}, R2: start={:.1} end={:.1}. Base-conditioned bins: {}/{} usable",
            pairs.len(),
            read_length,
            r1_mean_start, r1_mean_mid, r1_mean_end,
            r2_mean_start, r2_mean_end,
            base_bins_ok, base_bins_total,
        );

        Self {
            r1_base_quals: r1_base,
            r2_base_quals: r2_base,
            r1_cycle_quals: r1_cycle,
            r2_cycle_quals: r2_cycle,
        }
    }

    /// Sample a quality score (Phred+33 ASCII) for a given read number, cycle,
    /// and sequenced base.
    ///
    /// Tries the base-conditioned distribution first. Falls back to cycle-only
    /// if the base bin has fewer than `MIN_BASE_OBS` observations.
    pub fn sample_quality(
        &self,
        read_num: u8,
        cycle: usize,
        base: u8,
        rng: &mut StdRng,
    ) -> u8 {
        let (base_quals, cycle_quals) = match read_num {
            1 => (&self.r1_base_quals, &self.r1_cycle_quals),
            _ => (&self.r2_base_quals, &self.r2_cycle_quals),
        };

        // Try base-conditioned first.
        if cycle < base_quals.len() {
            if let Some(bi) = base_index(base) {
                let bin = &base_quals[cycle][bi];
                if bin.len() >= MIN_BASE_OBS {
                    return bin[rng.gen_range(0..bin.len())];
                }
            }
        }

        // Fall back to cycle-only.
        if cycle < cycle_quals.len() && !cycle_quals[cycle].is_empty() {
            return cycle_quals[cycle][rng.gen_range(0..cycle_quals[cycle].len())];
        }

        b'!' + 20 // last resort: Q20
    }
}

/// Generates synthetic reads from reference sequence + learned quality profile.
///
/// When `haplotype_variants` is populated, synthetic reads carry the duplicated
/// haplotype's alleles at het SNP positions, producing correct BAF (2:1 ratio).
pub struct SynthReadGenerator<'a> {
    profile: QualityProfile,
    reference: &'a SharedReference,
    read_length: usize,
    /// Het SNP positions → allele on the duplicated haplotype.
    /// Synthetic reads substitute these bases instead of reference.
    haplotype_variants: HashMap<u64, u8>,
    /// Fraction of sequencing errors that are indels (vs substitutions).
    /// 0.0 = substitution-only (default), ~0.05 = typical Illumina.
    indel_error_rate: f64,
}

impl<'a> SynthReadGenerator<'a> {
    pub fn new(
        profile: QualityProfile,
        reference: &'a SharedReference,
        read_length: usize,
        indel_error_rate: f64,
    ) -> Self {
        Self {
            profile,
            reference,
            read_length,
            haplotype_variants: HashMap::new(),
            indel_error_rate,
        }
    }

    /// Read length this generator was configured for.
    pub fn read_length(&self) -> usize {
        self.read_length
    }

    /// Set the haplotype variant map for correct BAF in duplications.
    pub fn set_haplotype_variants(&mut self, variants: HashMap<u64, u8>) {
        if !variants.is_empty() {
            log::info!(
                "SynthReadGenerator: {} het SNP positions will carry haplotype alleles",
                variants.len(),
            );
        }
        self.haplotype_variants = variants;
    }

    /// Generate a single synthetic read at a reference position.
    ///
    /// Returns `(sequence, quality)` in forward-strand orientation relative to
    /// the reference. For read2 (reverse strand in BAM), the caller must
    /// reverse-complement the sequence and reverse the quality.
    ///
    /// `read_num`: 1 or 2 (controls which quality distribution is sampled).
    /// `reverse_cycles`: when true, sample quality from profile cycle `rl-1-c`
    ///   instead of `c`. Use this for R2 which will be reversed after generation,
    ///   so that the final FASTQ-oriented quality matches the learned profile.
    pub fn generate_read(
        &self,
        chrom: &str,
        ref_start: u64,
        read_num: u8,
        reverse_cycles: bool,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let rl = self.read_length;
        // Fetch extra ref bases in case indel errors shift our position.
        let fetch_extra = if self.indel_error_rate > 0.0 { 10 } else { 0 };
        let ref_end = ref_start + (rl + fetch_extra) as u64;

        let ref_seq = self
            .reference
            .fetch_sequence(chrom, ref_start, ref_end)
            .unwrap_or_else(|_| vec![b'N'; rl + fetch_extra]);

        let mut seq = Vec::with_capacity(rl);
        let mut qual = Vec::with_capacity(rl);
        let mut ref_idx = 0usize; // current position in ref_seq

        while seq.len() < rl && ref_idx < ref_seq.len() {
            let c = seq.len(); // cycle position in the read

            let ref_base = ref_seq[ref_idx].to_ascii_uppercase();

            let genomic_pos = ref_start + ref_idx as u64;
            let true_base = self
                .haplotype_variants
                .get(&genomic_pos)
                .copied()
                .unwrap_or(ref_base);

            let qual_cycle = if reverse_cycles { rl - 1 - c } else { c };
            let profile_base = if reverse_cycles {
                complement(true_base)
            } else {
                true_base
            };
            let q = self.profile.sample_quality(read_num, qual_cycle, profile_base, rng);

            if true_base == b'N' {
                seq.push(b'N');
                qual.push(q);
                ref_idx += 1;
                continue;
            }

            let phred = (q as f64 - 33.0).max(0.0);
            let p_err = 10.0_f64.powf(-phred / 10.0);

            if rng.gen::<f64>() < p_err {
                if self.indel_error_rate > 0.0 && rng.gen::<f64>() < self.indel_error_rate {
                    // Indel error: 50/50 insertion vs deletion.
                    if rng.gen::<bool>() {
                        // Insertion: add a random base without consuming ref.
                        seq.push(random_base(rng));
                        qual.push(q);
                        // Don't advance ref_idx — the ref base will be read next cycle.
                    } else {
                        // Deletion: skip this ref base entirely.
                        ref_idx += 1;
                        // Don't add to seq/qual — next iteration will read next ref base.
                    }
                } else {
                    // Substitution error.
                    seq.push(random_different_base(true_base, rng));
                    qual.push(q);
                    ref_idx += 1;
                }
            } else {
                seq.push(true_base);
                qual.push(q);
                ref_idx += 1;
            }
        }

        // Pad if we ran out of ref bases due to deletions.
        while seq.len() < rl {
            seq.push(b'N');
            qual.push(b'!' + 2); // Q2
        }

        // Truncate if insertions made it too long (shouldn't happen with while < rl, but safety).
        seq.truncate(rl);
        qual.truncate(rl);

        (seq, qual)
    }

    /// Generate a synthetic read pair for a fragment at a given position.
    ///
    /// R1 is forward-strand at `frag_start`. R2 is reverse-strand at
    /// `frag_start + frag_len - read_length`, stored in FASTQ orientation
    /// (reverse-complemented).
    pub fn generate_read_pair(
        &self,
        chrom: &str,
        frag_start: u64,
        frag_len: u64,
        name: &str,
        rng: &mut StdRng,
    ) -> Option<ReadPair> {
        let rl = self.read_length as u64;
        if frag_len < rl {
            return None;
        }

        let r2_start = frag_start + frag_len - rl;

        // Generate R1 (forward strand).
        let (seq1, qual1) = self.generate_read(chrom, frag_start, 1, false, rng);

        // Generate R2 (reverse strand) — generate forward then revcomp for FASTQ.
        // reverse_cycles=true so quality profile cycles align after qual2.reverse().
        let (mut seq2, mut qual2) = self.generate_read(chrom, r2_start, 2, true, rng);
        reverse_complement(&mut seq2);
        qual2.reverse();

        Some(ReadPair {
            name: name.to_string(),
            seq1,
            qual1,
            seq2,
            qual2,
            ref_start: frag_start,
            ref_end: frag_start + frag_len,
            insert_size: frag_len as i64,
            chrom: chrom.to_string(),
        })
    }

    /// Generate a synthetic depth-copy pair near an original read pair's position.
    ///
    /// Samples a new fragment length from `frag_dist` and adds small position
    /// jitter (±20bp) to avoid exact positional duplicates that dedup tools detect.
    pub fn generate_depth_pair(
        &self,
        chrom: &str,
        original: &ReadPair,
        name: &str,
        frag_dist: &FragmentDist,
        rng: &mut StdRng,
    ) -> Option<ReadPair> {
        let rl = self.read_length as i64;

        // Sample new fragment length.
        let new_frag = frag_dist.sample_in_range(rng, rl, 1500) as u64;

        // Position jitter: ±20bp to avoid exact positional duplicates.
        let jitter = rng.gen_range(-20i64..=20);
        let new_start = (original.ref_start as i64 + jitter).max(0) as u64;

        self.generate_read_pair(chrom, new_start, new_frag, name, rng)
    }

    /// Generate synthetic depth-copy pairs for a duplication region.
    ///
    /// For each original pair that overlaps [dup_start, dup_end) by at least 50%
    /// of its fragment length, either uses haplotype information (if available)
    /// or random selection at VAF rate to decide which pairs get a synthetic
    /// depth copy.
    ///
    /// When haplotype info is available:
    /// - Classified reads in `hap_set` → always copy (target haplotype).
    /// - Classified reads NOT in `hap_set` → never copy (other haplotype).
    /// - Unclassified reads (not in `classified_set`) → randomly copy at VAF rate.
    #[allow(clippy::too_many_arguments)]
    pub fn generate_dup_depth_copies(
        &self,
        pool: &ReadPool,
        dup_start: u64,
        dup_end: u64,
        vaf: f64,
        hap_set: &std::collections::HashSet<String>,
        classified_set: &std::collections::HashSet<String>,
        rng: &mut StdRng,
    ) -> Vec<ReadPair> {
        let use_hap = !hap_set.is_empty();
        let mut copies = Vec::new();

        for (i, pair) in pool.pairs.iter().enumerate() {
            // Include reads that overlap the DUP region by >= 50% of their
            // fragment length. The old filter (entirely inside) missed ~16%
            // of reads at the boundaries, leading to lower-than-expected depth.
            let overlap_start = pair.ref_start.max(dup_start);
            let overlap_end = pair.ref_end.min(dup_end);
            if overlap_end <= overlap_start {
                continue; // no overlap at all
            }
            let overlap = overlap_end - overlap_start;
            let frag_len = pair.ref_end.saturating_sub(pair.ref_start).max(1);
            if overlap * 2 < frag_len {
                continue; // less than 50% inside the DUP region
            }

            let should_copy = if use_hap {
                if classified_set.contains(&pair.name) {
                    if hap_set.contains(&pair.name) {
                        // Variant haplotype: copy at rate min(1, 2*vaf).
                        // At VAF=0.5 this copies all hap reads (50% of total).
                        // At VAF=0.3 this copies 60% of hap reads (30% of total).
                        rng.gen::<f64>() < (2.0 * vaf).min(1.0)
                    } else {
                        // Other haplotype: only copy when VAF > 0.5.
                        rng.gen::<f64>() < (2.0 * vaf - 1.0).max(0.0)
                    }
                } else {
                    rng.gen::<f64>() < vaf
                }
            } else {
                rng.gen::<f64>() < vaf
            };

            if should_copy {
                if let Some(synth) = self.generate_depth_pair(
                    &pair.chrom,
                    pair,
                    &format!("sim_dup_depth_{:06}", i),
                    &pool.frag_dist,
                    rng,
                ) {
                    copies.push(synth);
                }
            }
        }

        log::info!(
            "Generated {} synthetic depth copies (dup region {}-{})",
            copies.len(),
            dup_start,
            dup_end,
        );

        copies
    }

    /// Generate a synthetic read from an arbitrary sequence slice.
    ///
    /// Like `generate_read` but takes a pre-built sequence (e.g. from a
    /// `VariantHaplotype`) instead of fetching from the reference FASTA.
    /// Haplotype variants are already baked into the sequence.
    ///
    /// `read_num`: 1 or 2 (controls which quality distribution is sampled).
    /// `reverse_cycles`: when true, sample quality from profile cycle `rl-1-c`
    ///   instead of `c`. Use this for R2 (see `generate_read` docs).
    pub fn generate_read_from_seq(
        &self,
        seq: &[u8],
        read_num: u8,
        reverse_cycles: bool,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let rl = self.read_length.min(seq.len());
        let mut out_seq = Vec::with_capacity(rl);
        let mut out_qual = Vec::with_capacity(rl);
        let mut seq_idx = 0usize;

        while out_seq.len() < rl && seq_idx < seq.len() {
            let c = out_seq.len();
            let true_base = seq[seq_idx].to_ascii_uppercase();

            let qual_cycle = if reverse_cycles { rl - 1 - c } else { c };
            let profile_base = if reverse_cycles {
                complement(true_base)
            } else {
                true_base
            };
            let q = self.profile.sample_quality(read_num, qual_cycle, profile_base, rng);

            if true_base == b'N' {
                out_seq.push(b'N');
                out_qual.push(q);
                seq_idx += 1;
                continue;
            }

            let phred = (q as f64 - 33.0).max(0.0);
            let p_err = 10.0_f64.powf(-phred / 10.0);

            if rng.gen::<f64>() < p_err {
                if self.indel_error_rate > 0.0 && rng.gen::<f64>() < self.indel_error_rate {
                    if rng.gen::<bool>() {
                        // Insertion error.
                        out_seq.push(random_base(rng));
                        out_qual.push(q);
                    } else {
                        // Deletion error.
                        seq_idx += 1;
                    }
                } else {
                    out_seq.push(random_different_base(true_base, rng));
                    out_qual.push(q);
                    seq_idx += 1;
                }
            } else {
                out_seq.push(true_base);
                out_qual.push(q);
                seq_idx += 1;
            }
        }

        while out_seq.len() < rl {
            out_seq.push(b'N');
            out_qual.push(b'!' + 2);
        }

        out_seq.truncate(rl);
        out_qual.truncate(rl);

        (out_seq, out_qual)
    }

    /// Generate a synthetic read pair from a variant haplotype.
    ///
    /// R1 starts at `hap_frag_start` (forward), R2 at `hap_frag_start + frag_len - rl`
    /// (reverse-complement for FR orientation). Reference coordinates for the
    /// ReadPair are mapped back from the haplotype via `hap_to_ref`.
    pub fn generate_haplotype_read_pair(
        &self,
        haplotype: &VariantHaplotype,
        hap_frag_start: u64,
        frag_len: u64,
        name: &str,
        rng: &mut StdRng,
    ) -> Option<ReadPair> {
        let rl = self.read_length as u64;
        if frag_len < rl {
            return None;
        }

        let r2_hap_start = hap_frag_start + frag_len - rl;

        // Get sequences from the haplotype.
        let r1_seq = haplotype.get_sequence(hap_frag_start, self.read_length);
        let r2_seq = haplotype.get_sequence(r2_hap_start, self.read_length);

        if r1_seq.len() < self.read_length || r2_seq.len() < self.read_length {
            return None;
        }

        // Generate R1 (forward strand).
        let (seq1, qual1) = self.generate_read_from_seq(r1_seq, 1, false, rng);

        // Generate R2 (reverse strand): generate forward then revcomp.
        let (mut seq2, mut qual2) = self.generate_read_from_seq(r2_seq, 2, true, rng);
        reverse_complement(&mut seq2);
        qual2.reverse();

        // Map haplotype positions back to reference coordinates.
        let (chrom, r1_ref) = haplotype
            .hap_to_ref(hap_frag_start)
            .unwrap_or_else(|| (haplotype.primary_chrom().to_string(), 0));
        let (_, r2_ref) = haplotype
            .hap_to_ref(r2_hap_start + rl - 1)
            .unwrap_or_else(|| (haplotype.primary_chrom().to_string(), r1_ref + frag_len));

        // Normalize: ensure ref_start <= ref_end. For reads in inverted
        // segments, R1 maps to a higher ref position than R2 (reversed mapping).
        let ref_start = r1_ref.min(r2_ref);
        let ref_end = r1_ref.max(r2_ref) + 1; // exclusive

        Some(ReadPair {
            name: name.to_string(),
            seq1,
            qual1,
            seq2,
            qual2,
            ref_start,
            ref_end,
            insert_size: frag_len as i64,
            chrom,
        })
    }
}

/// Map a DNA base to an index: A=0, C=1, G=2, T=3. Returns None for N or other.
fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Complement of a DNA base: A↔T, C↔G.
fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => base, // N stays N
    }
}

/// Pick a random base (uniform among A, C, G, T).
fn random_base(rng: &mut StdRng) -> u8 {
    b"ACGT"[rng.gen_range(0..4)]
}

/// Pick a random base different from the given one (uniform among the other 3).
fn random_different_base(base: u8, rng: &mut StdRng) -> u8 {
    let others: &[u8] = match base {
        b'A' => b"CGT",
        b'C' => b"AGT",
        b'G' => b"ACT",
        b'T' => b"ACG",
        _ => b"ACGT",
    };
    others[rng.gen_range(0..others.len())]
}

/// Mean Phred quality (Q value, not ASCII) of a quality vector.
fn mean_qual(quals: &[u8]) -> f64 {
    if quals.is_empty() {
        return 0.0;
    }
    let sum: f64 = quals.iter().map(|&q| (q as f64 - 33.0).max(0.0)).sum();
    sum / quals.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn mock_read_pair(
        name: &str,
        qual1: Vec<u8>,
        qual2: Vec<u8>,
        ref_start: u64,
    ) -> ReadPair {
        let rl = qual1.len();
        ReadPair {
            name: name.to_string(),
            seq1: vec![b'A'; rl],
            qual1,
            seq2: vec![b'T'; rl],
            qual2,
            ref_start,
            ref_end: ref_start + rl as u64 * 2,
            insert_size: (rl * 2) as i64,
            chrom: "chr1".to_string(),
        }
    }

    #[test]
    fn test_quality_profile_from_pairs() {
        let rl = 10;
        // Create pairs with known quality patterns: Q30 at all cycles for R1,
        // Q20 at all cycles for R2.
        let q30 = vec![b'!' + 30; rl]; // Phred+33: Q30
        let q20 = vec![b'!' + 20; rl]; // Phred+33: Q20

        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| mock_read_pair(&format!("read_{}", i), q30.clone(), q20.clone(), i * 10))
            .collect();

        let profile = QualityProfile::from_read_pairs(&pairs, rl);

        // Cycle-only fallback should have all 100 observations per cycle.
        assert_eq!(profile.r1_cycle_quals.len(), rl);
        assert_eq!(profile.r2_cycle_quals.len(), rl);
        assert_eq!(profile.r1_cycle_quals[0].len(), 100);

        // All R1 cycle-only values should be Q30.
        assert!(profile.r1_cycle_quals[0].iter().all(|&q| q == b'!' + 30));
        // All R2 cycle-only values should be Q20.
        assert!(profile.r2_cycle_quals[0].iter().all(|&q| q == b'!' + 20));

        // Base-conditioned: mock_read_pair sets seq1=all-A, so A bin should have data,
        // other bins should be empty.
        assert_eq!(profile.r1_base_quals[0][0].len(), 100); // A bin
        assert_eq!(profile.r1_base_quals[0][1].len(), 0); // C bin
        assert_eq!(profile.r1_base_quals[0][2].len(), 0); // G bin
        assert_eq!(profile.r1_base_quals[0][3].len(), 0); // T bin

        // R2: seq2=all-T, so T bin should have data.
        assert_eq!(profile.r2_base_quals[0][3].len(), 100); // T bin
        assert_eq!(profile.r2_base_quals[0][0].len(), 0); // A bin
    }

    #[test]
    fn test_quality_sampling_distribution() {
        let rl = 10;
        // Mix of Q10 and Q30 at cycle 0 (50/50 split).
        let pairs: Vec<ReadPair> = (0..200)
            .map(|i| {
                let q = if i < 100 {
                    vec![b'!' + 10; rl]
                } else {
                    vec![b'!' + 30; rl]
                };
                mock_read_pair(&format!("read_{}", i), q.clone(), q, i * 10)
            })
            .collect();

        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let mut rng = StdRng::seed_from_u64(42);

        // Sample 1000 values at cycle 0 with base A (matching mock seq1).
        let samples: Vec<u8> = (0..1000)
            .map(|_| profile.sample_quality(1, 0, b'A', &mut rng))
            .collect();

        let q10_count = samples.iter().filter(|&&q| q == b'!' + 10).count();
        let q30_count = samples.iter().filter(|&&q| q == b'!' + 30).count();

        // Both should be roughly 50% (±10% tolerance).
        assert!(q10_count > 350, "too few Q10: {}", q10_count);
        assert!(q30_count > 350, "too few Q30: {}", q30_count);
    }

    #[test]
    fn test_base_conditioned_vs_fallback() {
        let rl = 10;
        // Create pairs where seq1 = all-A with Q30, so base-conditioned bin for A
        // is populated but bins for C, G, T are empty.
        let q30 = vec![b'!' + 30; rl];
        let q20 = vec![b'!' + 20; rl];

        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| mock_read_pair(&format!("read_{}", i), q30.clone(), q20.clone(), i * 10))
            .collect();

        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let mut rng = StdRng::seed_from_u64(42);

        // Querying with base=A should use the base-conditioned bin (all Q30).
        let q_a = profile.sample_quality(1, 0, b'A', &mut rng);
        assert_eq!(q_a, b'!' + 30);

        // Querying with base=C should fall back to cycle-only (also all Q30,
        // since all reads had Q30 regardless of base). The C bin has 0 obs.
        let q_c = profile.sample_quality(1, 0, b'C', &mut rng);
        assert_eq!(q_c, b'!' + 30); // cycle-only fallback, still Q30

        // Querying with base=N should fall back to cycle-only.
        let q_n = profile.sample_quality(1, 0, b'N', &mut rng);
        assert_eq!(q_n, b'!' + 30);
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'N'), b'N');
    }

    #[test]
    fn test_error_rate_matches_quality() {
        // At Q10, error rate should be ~10%.
        let phred = 10.0;
        let p_err = 10.0_f64.powf(-phred / 10.0);
        assert!((p_err - 0.1).abs() < 0.001);

        // At Q30, error rate should be ~0.1%.
        let phred = 30.0;
        let p_err = 10.0_f64.powf(-phred / 10.0);
        assert!((p_err - 0.001).abs() < 0.0001);
    }

    #[test]
    fn test_random_different_base() {
        let mut rng = StdRng::seed_from_u64(42);

        // Verify base is always different.
        for _ in 0..100 {
            let b = random_different_base(b'A', &mut rng);
            assert_ne!(b, b'A');
            assert!(b == b'C' || b == b'G' || b == b'T');
        }
    }

    #[test]
    fn test_mean_qual() {
        let quals = vec![b'!' + 30, b'!' + 30, b'!' + 30]; // all Q30
        assert!((mean_qual(&quals) - 30.0).abs() < 0.001);

        let quals = vec![b'!' + 10, b'!' + 30]; // mean = 20
        assert!((mean_qual(&quals) - 20.0).abs() < 0.001);
    }
}
