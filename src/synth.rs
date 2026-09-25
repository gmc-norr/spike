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
use crate::loh::SampleCopies;
use crate::reference::SharedReference;
use crate::stats::FragmentDist;
use crate::types::{ReadPair, ReadPool};

/// Minimum observations in a (cycle, base) bin before we trust it.
/// Below this threshold, fall back to the cycle-only distribution.
const MIN_BASE_OBS: usize = 30;

/// Number of bins for quantizing the previous quality score in the Markov model.
/// Bins: Q0-9 → 0, Q10-19 → 1, Q20-29 → 2, Q30+ → 3.
const PREV_Q_BINS: usize = 4;

/// Minimum observations in a Markov transition bin before we trust it.
const MIN_MARKOV_OBS: usize = 30;

/// Extra template bases fetched past a read's 3' end when indel errors are on,
/// so a deletion error is covered by real sequence instead of `N` padding (L1).
const INDEL_SLACK: usize = 10;

/// Quality reported for every `N` a synthetic read emits, whatever put it
/// there — a reference `N`, padding past a contig end, or padding after the
/// template ran out. An `N` is a no-call, and a real Illumina no-call is
/// always Q2; the learned profile knows nothing about `N` and would hand one
/// an ordinary score (often Q37) instead (L18).
const N_QUAL: u8 = b'!' + 2; // Q2

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

    /// Markov transition table for read1, conditioned on base.
    /// `r1_markov_base[cycle][base_idx][prev_q_bin]` = sorted Vec<u8> of Phred+33 values.
    r1_markov_base: Vec<[[Vec<u8>; PREV_Q_BINS]; 4]>,

    /// Markov transition table for read2, conditioned on base.
    r2_markov_base: Vec<[[Vec<u8>; PREV_Q_BINS]; 4]>,

    /// Markov transition table for read1, cycle-only (all bases pooled).
    /// `r1_markov_cycle[cycle][prev_q_bin]` = sorted Vec<u8>.
    r1_markov_cycle: Vec<[Vec<u8>; PREV_Q_BINS]>,

    /// Markov transition table for read2, cycle-only.
    r2_markov_cycle: Vec<[Vec<u8>; PREV_Q_BINS]>,
}

impl QualityProfile {
    /// Learn quality profile from extracted read pairs.
    ///
    /// Uses the read's own sequence bases (seq1/seq2) as the conditioning context.
    /// These are in FASTQ orientation (matching the quality scores) and are ~99%
    /// correct, so they faithfully represent the base the sequencer was reading.
    pub fn from_read_pairs(pairs: &[ReadPair], read_length: usize) -> Self {
        // Marginal accumulators (existing).
        let mut r1_base: Vec<[Vec<u8>; 4]> = (0..read_length).map(|_| Default::default()).collect();
        let mut r2_base: Vec<[Vec<u8>; 4]> = (0..read_length).map(|_| Default::default()).collect();
        let mut r1_cycle: Vec<Vec<u8>> = (0..read_length).map(|_| Vec::new()).collect();
        let mut r2_cycle: Vec<Vec<u8>> = (0..read_length).map(|_| Vec::new()).collect();

        // Markov transition accumulators.
        let mut r1_mkv_base: Vec<[[Vec<u8>; PREV_Q_BINS]; 4]> =
            (0..read_length).map(|_| Default::default()).collect();
        let mut r2_mkv_base: Vec<[[Vec<u8>; PREV_Q_BINS]; 4]> =
            (0..read_length).map(|_| Default::default()).collect();
        let mut r1_mkv_cycle: Vec<[Vec<u8>; PREV_Q_BINS]> =
            (0..read_length).map(|_| Default::default()).collect();
        let mut r2_mkv_cycle: Vec<[Vec<u8>; PREV_Q_BINS]> =
            (0..read_length).map(|_| Default::default()).collect();

        for pair in pairs {
            let r1_len = read_length.min(pair.qual1.len()).min(pair.seq1.len());
            for c in 0..r1_len {
                let q = pair.qual1[c];
                r1_cycle[c].push(q);
                if let Some(bi) = base_index(pair.seq1[c]) {
                    r1_base[c][bi].push(q);
                }
                // Markov: record transition from previous quality (cycle > 0).
                if c > 0 {
                    let pbin = prev_q_bin(pair.qual1[c - 1]);
                    r1_mkv_cycle[c][pbin].push(q);
                    if let Some(bi) = base_index(pair.seq1[c]) {
                        r1_mkv_base[c][bi][pbin].push(q);
                    }
                }
            }

            let r2_len = read_length.min(pair.qual2.len()).min(pair.seq2.len());
            for c in 0..r2_len {
                let q = pair.qual2[c];
                r2_cycle[c].push(q);
                if let Some(bi) = base_index(pair.seq2[c]) {
                    r2_base[c][bi].push(q);
                }
                if c > 0 {
                    let pbin = prev_q_bin(pair.qual2[c - 1]);
                    r2_mkv_cycle[c][pbin].push(q);
                    if let Some(bi) = base_index(pair.seq2[c]) {
                        r2_mkv_base[c][bi][pbin].push(q);
                    }
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
        for cycle in r1_mkv_base.iter_mut() {
            for base_bins in cycle.iter_mut() {
                for bin in base_bins.iter_mut() {
                    bin.sort_unstable();
                }
            }
        }
        for cycle in r2_mkv_base.iter_mut() {
            for base_bins in cycle.iter_mut() {
                for bin in base_bins.iter_mut() {
                    bin.sort_unstable();
                }
            }
        }
        for cycle in r1_mkv_cycle.iter_mut() {
            for bin in cycle.iter_mut() {
                bin.sort_unstable();
            }
        }
        for cycle in r2_mkv_cycle.iter_mut() {
            for bin in cycle.iter_mut() {
                bin.sort_unstable();
            }
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

        // Count Markov bin usability.
        let mkv_base_total = read_length * 4 * PREV_Q_BINS * 2;
        let mkv_base_ok = r1_mkv_base
            .iter()
            .chain(r2_mkv_base.iter())
            .flat_map(|cycle| cycle.iter().flat_map(|base| base.iter()))
            .filter(|bin| bin.len() >= MIN_MARKOV_OBS)
            .count();
        let mkv_cycle_total = read_length * PREV_Q_BINS * 2;
        let mkv_cycle_ok = r1_mkv_cycle
            .iter()
            .chain(r2_mkv_cycle.iter())
            .flat_map(|cycle| cycle.iter())
            .filter(|bin| bin.len() >= MIN_MARKOV_OBS)
            .count();

        log::info!(
            "Quality profile: {} pairs, {} cycles. R1 mean Q: start={:.1} mid={:.1} end={:.1}, R2: start={:.1} end={:.1}. \
             Base-conditioned bins: {}/{} usable. Markov bins: base {}/{}, cycle {}/{} usable",
            pairs.len(),
            read_length,
            r1_mean_start, r1_mean_mid, r1_mean_end,
            r2_mean_start, r2_mean_end,
            base_bins_ok, base_bins_total,
            mkv_base_ok, mkv_base_total,
            mkv_cycle_ok, mkv_cycle_total,
        );

        Self {
            r1_base_quals: r1_base,
            r2_base_quals: r2_base,
            r1_cycle_quals: r1_cycle,
            r2_cycle_quals: r2_cycle,
            r1_markov_base: r1_mkv_base,
            r2_markov_base: r2_mkv_base,
            r1_markov_cycle: r1_mkv_cycle,
            r2_markov_cycle: r2_mkv_cycle,
        }
    }

    /// Sample a quality score (Phred+33 ASCII) for a given read number, cycle,
    /// sequenced base, and optionally the previous quality score.
    ///
    /// Uses a 4-level fallback hierarchy:
    /// 1. Markov + base: `P(Q_i | cycle, base, prev_q_bin)` — best if enough data
    /// 2. Markov + cycle: `P(Q_i | cycle, prev_q_bin)` — drop base conditioning
    /// 3. Base-only: `P(Q_i | cycle, base)` — no Markov (cycle 0 or sparse bins)
    /// 4. Cycle-only: `P(Q_i | cycle)` — final fallback
    ///
    /// The returned byte is clamped to the valid Phred+33 range [b'!', b'~'] = [33, 126].
    ///
    /// This clamp is defense in depth, not the primary fix for M14 (missing
    /// donor qualities poisoning this model with byte 32): `extract.rs` now
    /// drops donor records with no stored quality, or with a raw score above
    /// Q93, before they ever reach this profile, and
    /// `fastq.rs::write_paired_fastq` refuses to write any out-of-range byte
    /// regardless of source. Kept here because nothing proves this sampler
    /// can never be fed bad data through some other path -- and because a
    /// poisoned donor is only ever *suppressed*, never written, so the writer
    /// would not catch it.
    pub fn sample_quality(
        &self,
        read_num: u8,
        cycle: usize,
        base: u8,
        prev_qual: Option<u8>,
        rng: &mut StdRng,
    ) -> u8 {
        // Phred+33: valid byte range is b'!' (Q0) to b'~' (Q93).
        const MIN_QUAL_BYTE: u8 = b'!';
        const MAX_QUAL_BYTE: u8 = b'!' + 93; // '~' = 126
        let q = self.sample_quality_inner(read_num, cycle, base, prev_qual, rng);
        q.clamp(MIN_QUAL_BYTE, MAX_QUAL_BYTE)
    }

    fn sample_quality_inner(
        &self,
        read_num: u8,
        cycle: usize,
        base: u8,
        prev_qual: Option<u8>,
        rng: &mut StdRng,
    ) -> u8 {
        let (base_quals, cycle_quals, mkv_base, mkv_cycle) = match read_num {
            1 => (
                &self.r1_base_quals,
                &self.r1_cycle_quals,
                &self.r1_markov_base,
                &self.r1_markov_cycle,
            ),
            _ => (
                &self.r2_base_quals,
                &self.r2_cycle_quals,
                &self.r2_markov_base,
                &self.r2_markov_cycle,
            ),
        };

        // If we have a previous quality, try Markov tables first.
        if let Some(pq) = prev_qual {
            let pbin = prev_q_bin(pq);

            // Level 1: Markov + base-conditioned.
            if cycle < mkv_base.len() {
                if let Some(bi) = base_index(base) {
                    let bin = &mkv_base[cycle][bi][pbin];
                    if bin.len() >= MIN_MARKOV_OBS {
                        return bin[rng.gen_range(0..bin.len())];
                    }
                }
            }

            // Level 2: Markov + cycle-only.
            if cycle < mkv_cycle.len() {
                let bin = &mkv_cycle[cycle][pbin];
                if bin.len() >= MIN_MARKOV_OBS {
                    return bin[rng.gen_range(0..bin.len())];
                }
            }
        }

        // Level 3: Base-conditioned marginal (no Markov).
        if cycle < base_quals.len() {
            if let Some(bi) = base_index(base) {
                let bin = &base_quals[cycle][bi];
                if bin.len() >= MIN_BASE_OBS {
                    return bin[rng.gen_range(0..bin.len())];
                }
            }
        }

        // Level 4: Cycle-only marginal.
        if cycle < cycle_quals.len() && !cycle_quals[cycle].is_empty() {
            return cycle_quals[cycle][rng.gen_range(0..cycle_quals[cycle].len())];
        }

        b'!' + 20 // last resort: Q20
    }
}

/// Generates synthetic reads from reference sequence + learned quality profile.
pub struct SynthReadGenerator<'a> {
    profile: QualityProfile,
    reference: &'a SharedReference,
    read_length: usize,
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
            indel_error_rate,
        }
    }

    /// Read length this generator was configured for.
    pub fn read_length(&self) -> usize {
        self.read_length
    }

    /// The reference this generator reads from.
    pub fn reference(&self) -> &'a SharedReference {
        self.reference
    }

    /// Generate a read from a template that is already in sequencing order.
    ///
    /// `template` is what the sequencer reads, 5'→3', with the sample's own
    /// alleles baked in; it should carry a few bases past `rl` so a deletion
    /// error is covered by real sequence instead of `N` padding (L1).
    /// Generating in sequencing order also keeps the quality Markov chain
    /// running with the cycle counter for both mates, not against it (L17).
    fn generate_from_template(
        &self,
        template: &[u8],
        rl: usize,
        read_num: u8,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut seq = Vec::with_capacity(rl);
        let mut qual = Vec::with_capacity(rl);
        let mut idx = 0usize; // current position in template
        let mut prev_qual: Option<u8> = None; // Markov chain state

        while seq.len() < rl && idx < template.len() {
            let c = seq.len(); // cycle position in the read
            let true_base = template[idx].to_ascii_uppercase();

            let q = self
                .profile
                .sample_quality(read_num, c, true_base, prev_qual, rng);

            if true_base == b'N' {
                // An `N` is a no-call and reports Q2, not the profile's `q`
                // (L18). `q` is still drawn, so an `N` costs the same single
                // random draw as any other base, and the Markov chain carries
                // the quality the read actually reports.
                seq.push(b'N');
                qual.push(N_QUAL);
                prev_qual = Some(N_QUAL);
                idx += 1;
                continue;
            }

            let phred = (q as f64 - 33.0).max(0.0);
            let p_err = 10.0_f64.powf(-phred / 10.0);

            if rng.gen::<f64>() < p_err {
                if self.indel_error_rate > 0.0 && rng.gen::<f64>() < self.indel_error_rate {
                    // Indel error: 50/50 insertion vs deletion.
                    if rng.gen::<bool>() {
                        // Insertion: add a random base without consuming template.
                        seq.push(random_base(rng));
                        qual.push(q);
                        prev_qual = Some(q);
                        // Don't advance idx — the template base is read next cycle.
                    } else {
                        // Deletion: skip this template base entirely.
                        idx += 1;
                        // Don't add to seq/qual — next iteration reads the next base.
                        // Don't update prev_qual — no quality was emitted.
                    }
                } else {
                    // Substitution error.
                    seq.push(random_different_base(true_base, rng));
                    qual.push(q);
                    prev_qual = Some(q);
                    idx += 1;
                }
            } else {
                seq.push(true_base);
                qual.push(q);
                prev_qual = Some(q);
                idx += 1;
            }
        }

        // Pad only if the template itself ran out (contig or haplotype end).
        while seq.len() < rl {
            seq.push(b'N');
            qual.push(N_QUAL);
        }

        // Truncate if insertions made it too long (shouldn't happen with while < rl, but safety).
        seq.truncate(rl);
        qual.truncate(rl);

        (seq, qual)
    }

    /// Generate a single synthetic read at a reference position.
    ///
    /// `alleles` maps reference position → the base of the sample copy the
    /// read comes from; other positions take the reference base.
    ///
    /// `ref_start` is the leftmost reference base the read covers either way.
    ///
    /// `read_num`: 1 or 2 (controls which quality distribution is sampled).
    /// `is_reverse`: when true the read comes off the reverse strand. Its
    ///   template is complemented and walked right to left, so generation runs
    ///   in sequencing order and the returned `(sequence, quality)` is already
    ///   in FASTQ orientation — the caller must not reverse it again. A
    ///   forward read is returned in reference orientation, which is the same
    ///   thing. Reads are `read_length` long either way.
    fn generate_read(
        &self,
        chrom: &str,
        ref_start: u64,
        alleles: &HashMap<u64, u8>,
        read_num: u8,
        is_reverse: bool,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let rl = self.read_length;
        // Fetch extra ref bases in case indel errors shift our position. They
        // go past the read's 3' end, which for a reverse read is to the left.
        let slack = if self.indel_error_rate > 0.0 { INDEL_SLACK as u64 } else { 0 };
        let (fetch_start, fetch_end) = if is_reverse {
            (ref_start.saturating_sub(slack), ref_start + rl as u64)
        } else {
            (ref_start, ref_start + rl as u64 + slack)
        };

        let ref_seq = self
            .reference
            .fetch_sequence(chrom, fetch_start, fetch_end)
            .unwrap_or_else(|_| vec![b'N'; (fetch_end - fetch_start) as usize]);

        // The sample's own allele where it has one, else the reference base.
        let base_at = |i: usize| {
            let ref_base = ref_seq[i].to_ascii_uppercase();
            alleles
                .get(&(fetch_start + i as u64))
                .copied()
                .unwrap_or(ref_base)
                .to_ascii_uppercase()
        };
        // Put the template in sequencing order: a reverse-strand read runs
        // right to left along the reference, complemented. `fetch_sequence`
        // clamps to the contig length, so a fetch that runs off the contig
        // end comes back short at its high end — which for a reverse read is
        // its 5' start. Represent that shortfall as `N` there, not as a
        // window shifted onto real bases from past the read's other end.
        let template: Vec<u8> = if is_reverse {
            let missing = (fetch_end - fetch_start).saturating_sub(ref_seq.len() as u64) as usize;
            std::iter::repeat_n(b'N', missing)
                .chain((0..ref_seq.len()).rev().map(|i| complement(base_at(i))))
                .collect()
        } else {
            (0..ref_seq.len()).map(base_at).collect()
        };

        self.generate_from_template(&template, rl, read_num, rng)
    }

    /// Generate a synthetic read pair for a fragment at a given position.
    ///
    /// The fragment is `[frag_start, frag_start + frag_len)` either way; a
    /// coin flip decides which mate comes off which end. F1R2: R1 forward at
    /// `frag_start`, R2 reverse at `frag_start + frag_len - read_length`.
    /// F2R1: the other way round. Real libraries are ~50/50 (M13), and a
    /// one-sided strand makes callers' read-orientation filters fire.
    /// The reverse mate is stored in FASTQ orientation (reverse-complemented).
    /// Both reads carry `alleles` (see `generate_read`).
    pub fn generate_read_pair(
        &self,
        chrom: &str,
        frag_start: u64,
        frag_len: u64,
        alleles: &HashMap<u64, u8>,
        name: &str,
        rng: &mut StdRng,
    ) -> Option<ReadPair> {
        let rl = self.read_length as u64;
        if frag_len < rl {
            return None;
        }

        let right_start = frag_start + frag_len - rl;

        // Draw the orientation before either read, so the stream stays in a
        // fixed order and the same --seed keeps giving the same output (M7).
        let r1_is_reverse = rng.gen::<bool>();
        // Either way the left end is read forward and the right end reverse —
        // the flip only decides which of the two is R1. read_num picks the
        // quality model, is_reverse handles the flip, so whichever mate
        // ends up as R1 keeps the R1 model.
        let (fwd_num, rev_num) = if r1_is_reverse { (2, 1) } else { (1, 2) };

        // Forward mate (left end of the fragment).
        let (fwd_seq, fwd_qual) =
            self.generate_read(chrom, frag_start, alleles, fwd_num, false, rng);

        // Reverse mate (right end) — generated in sequencing order, so it
        // already comes back in FASTQ orientation.
        let (rev_seq, rev_qual) =
            self.generate_read(chrom, right_start, alleles, rev_num, true, rng);

        let (seq1, qual1, seq2, qual2) = if r1_is_reverse {
            (rev_seq, rev_qual, fwd_seq, fwd_qual)
        } else {
            (fwd_seq, fwd_qual, rev_seq, rev_qual)
        };

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
        alleles: &HashMap<u64, u8>,
        name: &str,
        frag_dist: &FragmentDist,
        rng: &mut StdRng,
    ) -> Option<ReadPair> {
        let rl = self.read_length as i64;

        // Sample new fragment length.
        let new_frag = frag_dist.sample_in_range(rng, rl, crate::stats::MAX_FRAGMENT_LEN) as u64;

        // Position jitter: ±20bp to avoid exact positional duplicates.
        let jitter = rng.gen_range(-20i64..=20);
        let new_start = (original.ref_start as i64 + jitter).max(0) as u64;

        self.generate_read_pair(chrom, new_start, new_frag, alleles, name, rng)
    }

    /// Generate synthetic depth-copy pairs for a duplication region.
    ///
    /// For each original pair that overlaps [dup_start, dup_end) by at least 50%
    /// of its fragment length, decides by the pair's copy whether it gets a
    /// synthetic depth copy:
    /// - event copy → at rate min(1, 2·vaf),
    /// - other copy → at rate max(0, 2·vaf − 1),
    /// - unknown copy → at rate vaf.
    ///
    /// A depth copy carries the alleles of the copy it repeats. One of unknown
    /// copy comes from the event copy with the event copy's share of the
    /// added reads, min(1, 2·vaf) / (2·vaf).
    #[allow(clippy::too_many_arguments)]
    pub fn generate_dup_depth_copies(
        &self,
        pool: &ReadPool,
        dup_start: u64,
        dup_end: u64,
        vaf: f64,
        sample: &SampleCopies,
        name_prefix: &str,
        rng: &mut StdRng,
    ) -> Vec<ReadPair> {
        let p_event_copy = if vaf > 0.5 { 1.0 / (2.0 * vaf) } else { 1.0 };
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

            let copy = sample.read_copy.get(&pair.name).copied();
            if rng.gen::<f64>() < copy_rate(copy, vaf) {
                // At VAF <= 0.5 all added reads are the event copy's: skip
                // the draw so the random stream matches the simple case.
                let on_event_copy =
                    copy.unwrap_or_else(|| p_event_copy >= 1.0 || rng.gen::<f64>() < p_event_copy);
                let alleles = if on_event_copy { &sample.event_copy } else { &sample.other_copy };
                if let Some(synth) = self.generate_depth_pair(
                    &pair.chrom,
                    pair,
                    alleles,
                    &format!("{}_dup_depth_{:06}", name_prefix, i),
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
    /// `seq` is in haplotype orientation and should carry a few bases past the
    /// read's 3' end so deletion errors have real sequence to fall back on.
    /// The read is `min(read_length, seq.len())` bases long.
    ///
    /// `read_num`: 1 or 2 (controls which quality distribution is sampled).
    /// `is_reverse`: when true the read comes off the reverse strand — `seq`
    ///   is reverse-complemented first and the read is returned in FASTQ
    ///   orientation (see `generate_read` docs).
    fn generate_read_from_seq(
        &self,
        seq: &[u8],
        read_num: u8,
        is_reverse: bool,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let rl = self.read_length.min(seq.len());
        let mut template = seq.to_vec();
        if is_reverse {
            // Sequencing order: the reverse mate reads the other strand,
            // 3' end of `seq` first.
            reverse_complement(&mut template);
        }
        self.generate_from_template(&template, rl, read_num, rng)
    }

    /// Generate a synthetic read pair from a variant haplotype.
    ///
    /// A coin flip decides which end of the fragment each mate comes off, as
    /// in `generate_read_pair` (M13): the forward mate starts at
    /// `hap_frag_start`, the reverse mate at `hap_frag_start + frag_len - rl`
    /// and is reverse-complemented for FR orientation. Reference coordinates
    /// for the ReadPair are mapped back from the haplotype via `hap_to_ref`.
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

        let right_hap_start = hap_frag_start + frag_len - rl;

        // Get sequences from the haplotype, with extra bases past each read's
        // 3' end so deletion errors read real sequence, not `N` padding (L1).
        // For the reverse mate that end is to the left of `right_hap_start`.
        let slack = if self.indel_error_rate > 0.0 { INDEL_SLACK } else { 0 };
        let right_slack = slack.min(right_hap_start as usize);
        let left_seq = haplotype.get_sequence(hap_frag_start, self.read_length + slack);
        let right_seq = haplotype.get_sequence(
            right_hap_start - right_slack as u64,
            right_slack + self.read_length,
        );

        if left_seq.len() < self.read_length
            || right_seq.len() < right_slack + self.read_length
        {
            return None;
        }

        // Draw the orientation before either read, so the stream stays in a
        // fixed order and the same --seed keeps giving the same output (M7).
        let r1_is_reverse = rng.gen::<bool>();
        // read_num picks the quality model, is_reverse handles the flip,
        // so whichever mate ends up as R1 keeps the R1 model.
        let (fwd_num, rev_num) = if r1_is_reverse { (2, 1) } else { (1, 2) };

        // Forward mate (left end of the fragment).
        let (fwd_seq, fwd_qual) = self.generate_read_from_seq(left_seq, fwd_num, false, rng);

        // Reverse mate (right end): generated in sequencing order, so it
        // already comes back in FASTQ orientation.
        let (rev_seq, rev_qual) = self.generate_read_from_seq(right_seq, rev_num, true, rng);

        let (seq1, qual1, seq2, qual2) = if r1_is_reverse {
            (rev_seq, rev_qual, fwd_seq, fwd_qual)
        } else {
            (fwd_seq, fwd_qual, rev_seq, rev_qual)
        };

        // Map haplotype positions back to reference coordinates. An end inside
        // inserted sequence takes the nearest reference base: such pairs are
        // the one-end-anchored evidence for the insertion.
        // Both fragment ends are mapped, so the span does not depend on which
        // mate is R1.
        let (chrom, left_ref) = haplotype.hap_to_ref_nearest(hap_frag_start)?;
        let (_, right_ref) = haplotype.hap_to_ref_nearest(right_hap_start + rl - 1)?;

        // Normalize: ensure ref_start <= ref_end. For reads in inverted
        // segments the left end maps to a higher ref position (reversed mapping).
        let ref_start = left_ref.min(right_ref);
        let ref_end = left_ref.max(right_ref) + 1; // exclusive

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

/// Rate at which an original pair is replaced (or, for depth copies, repeated)
/// given its copy: event copy min(1, 2·vaf), other copy max(0, 2·vaf − 1),
/// unknown copy vaf. Averaged over both copies this is vaf.
pub fn copy_rate(copy: Option<bool>, vaf: f64) -> f64 {
    match copy {
        Some(true) => (2.0 * vaf).min(1.0),
        Some(false) => (2.0 * vaf - 1.0).max(0.0),
        None => vaf,
    }
}

/// Quantize a Phred+33 quality score into a previous-quality bin for the Markov model.
/// Bins: Q0-9 → 0, Q10-19 → 1, Q20-29 → 2, Q30+ → 3.
fn prev_q_bin(phred_plus_33: u8) -> usize {
    let phred = phred_plus_33.saturating_sub(b'!');
    match phred {
        0..=9 => 0,
        10..=19 => 1,
        20..=29 => 2,
        _ => 3,
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

    fn mock_read_pair(name: &str, qual1: Vec<u8>, qual2: Vec<u8>, ref_start: u64) -> ReadPair {
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
            .map(|_| profile.sample_quality(1, 0, b'A', None, &mut rng))
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
        let q_a = profile.sample_quality(1, 0, b'A', None, &mut rng);
        assert_eq!(q_a, b'!' + 30);

        // Querying with base=C should fall back to cycle-only (also all Q30,
        // since all reads had Q30 regardless of base). The C bin has 0 obs.
        let q_c = profile.sample_quality(1, 0, b'C', None, &mut rng);
        assert_eq!(q_c, b'!' + 30); // cycle-only fallback, still Q30

        // Querying with base=N should fall back to cycle-only.
        let q_n = profile.sample_quality(1, 0, b'N', None, &mut rng);
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

    // ── Markov quality model tests ──────────────────────────────────────

    #[test]
    fn test_prev_q_bin() {
        assert_eq!(prev_q_bin(b'!' + 0), 0); // Q0 → bin 0
        assert_eq!(prev_q_bin(b'!' + 9), 0); // Q9 → bin 0
        assert_eq!(prev_q_bin(b'!' + 10), 1); // Q10 → bin 1
        assert_eq!(prev_q_bin(b'!' + 19), 1); // Q19 → bin 1
        assert_eq!(prev_q_bin(b'!' + 20), 2); // Q20 → bin 2
        assert_eq!(prev_q_bin(b'!' + 29), 2); // Q29 → bin 2
        assert_eq!(prev_q_bin(b'!' + 30), 3); // Q30 → bin 3
        assert_eq!(prev_q_bin(b'!' + 40), 3); // Q40 → bin 3
    }

    #[test]
    fn test_markov_tables_populated() {
        let rl = 10;
        let q30 = vec![b'!' + 30; rl];
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| mock_read_pair(&format!("read_{}", i), q30.clone(), q30.clone(), i * 10))
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);

        // Cycle 0 has no previous quality → Markov bins should be empty.
        for pbin in 0..PREV_Q_BINS {
            assert_eq!(profile.r1_markov_base[0][0][pbin].len(), 0);
            assert_eq!(profile.r1_markov_cycle[0][pbin].len(), 0);
        }

        // Cycle 1+: prev_q_bin(Q30) = 3, base A (idx 0) should have 100 obs.
        assert_eq!(profile.r1_markov_base[1][0][3].len(), 100);
        assert_eq!(profile.r1_markov_cycle[1][3].len(), 100);

        // Other prev_q bins at cycle 1 should be empty (all reads had Q30).
        assert_eq!(profile.r1_markov_base[1][0][0].len(), 0);
        assert_eq!(profile.r1_markov_base[1][0][1].len(), 0);
        assert_eq!(profile.r1_markov_base[1][0][2].len(), 0);
    }

    #[test]
    fn test_markov_quality_correlation() {
        let rl = 10;
        // Half reads are all-Q10, half are all-Q30. The Markov model should
        // learn that Q30 follows Q30 and Q10 follows Q10.
        let pairs: Vec<ReadPair> = (0..200)
            .map(|i| {
                let q = if i < 100 {
                    vec![b'!' + 10; rl]
                } else {
                    vec![b'!' + 30; rl]
                };
                mock_read_pair(&format!("r_{}", i), q.clone(), q, i * 10)
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let mut rng = StdRng::seed_from_u64(42);

        // Sample Q at cycle 5, base A, with prev_qual=Q30.
        // Markov should strongly prefer Q30 (from the all-Q30 reads).
        let samples: Vec<u8> = (0..1000)
            .map(|_| profile.sample_quality(1, 5, b'A', Some(b'!' + 30), &mut rng))
            .collect();
        let q30_frac = samples.iter().filter(|&&q| q == b'!' + 30).count() as f64 / 1000.0;
        // Without Markov: ~50%. With Markov: ~100% (Q30→Q30 only from Q30 reads).
        assert!(
            q30_frac > 0.85,
            "Markov Q30|prev=Q30 should be >85%, got {:.1}%",
            q30_frac * 100.0
        );

        // Conversely, prev_qual=Q10 should strongly prefer Q10.
        let samples: Vec<u8> = (0..1000)
            .map(|_| profile.sample_quality(1, 5, b'A', Some(b'!' + 10), &mut rng))
            .collect();
        let q10_frac = samples.iter().filter(|&&q| q == b'!' + 10).count() as f64 / 1000.0;
        assert!(
            q10_frac > 0.85,
            "Markov Q10|prev=Q10 should be >85%, got {:.1}%",
            q10_frac * 100.0
        );
    }

    #[test]
    fn test_markov_cycle0_falls_back_to_marginal() {
        let rl = 10;
        // Mix of Q10 and Q30 reads.
        let pairs: Vec<ReadPair> = (0..200)
            .map(|i| {
                let q = if i < 100 {
                    vec![b'!' + 10; rl]
                } else {
                    vec![b'!' + 30; rl]
                };
                mock_read_pair(&format!("r_{}", i), q.clone(), q, i * 10)
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let mut rng = StdRng::seed_from_u64(42);

        // At cycle 0, prev_qual=None → should use marginal (50/50 Q10/Q30).
        let samples: Vec<u8> = (0..1000)
            .map(|_| profile.sample_quality(1, 0, b'A', None, &mut rng))
            .collect();
        let q10_count = samples.iter().filter(|&&q| q == b'!' + 10).count();
        let q30_count = samples.iter().filter(|&&q| q == b'!' + 30).count();
        assert!(q10_count > 350, "Cycle 0 should have ~50% Q10, got {}", q10_count);
        assert!(q30_count > 350, "Cycle 0 should have ~50% Q30, got {}", q30_count);
    }

    // ── Haplotype read pair generation tests ────────────────────────────

    use crate::haplotype::{HaplotypeSegment, SegmentOrigin, VariantHaplotype};
    use crate::reference::SharedReference;
    use std::collections::HashMap as StdHashMap;

    /// Build a mock deletion haplotype with known segments.
    fn mock_del_haplotype(flank: u64, del_size: u64) -> VariantHaplotype {
        let pattern = b"ACGT";
        let left_seq: Vec<u8> = (0..flank).map(|i| pattern[(i % 4) as usize]).collect();
        let right_seq: Vec<u8> = (0..flank)
            .map(|i| pattern[((del_size + flank + i) % 4) as usize])
            .collect();

        VariantHaplotype::from_segments(vec![
            HaplotypeSegment {
                sequence: left_seq,
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 0,
                    ref_end: flank,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: right_seq,
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: flank + del_size,
                    ref_end: flank + del_size + flank,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ])
    }

    fn mock_shared_ref() -> SharedReference {
        let pattern = b"ACGT";
        let seq: Vec<u8> = (0..100_000u64).map(|i| pattern[(i % 4) as usize]).collect();
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), seq);
        SharedReference::from_sequences(seqs)
    }

    fn mock_synth_gen(read_length: usize) -> (SharedReference, SynthReadGenerator<'static>) {
        let q30 = vec![b'!' + 30; read_length];
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| {
                mock_read_pair(
                    &format!("mock_{}", i),
                    q30.clone(),
                    q30.clone(),
                    i as u64 * 500,
                )
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, read_length);
        let reference = mock_shared_ref();
        // SAFETY: we leak the reference to get a 'static lifetime for tests
        let ref_static: &'static SharedReference = Box::leak(Box::new(reference));
        let gen = SynthReadGenerator::new(profile, ref_static, read_length, 0.0);
        // Return a dummy ref for keeping alive (already leaked)
        (mock_shared_ref(), gen)
    }

    #[test]
    fn test_hap_read_pair_allows_boundary_crossing_r1() {
        // DEL: flank=500, del=5000 → breakpoint at hap offset 500
        // Hap: [seg0: 0..500] [seg1: 500..1000]
        let hap = mock_del_haplotype(500, 5000);
        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // R1 spans [400..550] → crosses boundary at 500 → produces split-read evidence
        let result = gen.generate_haplotype_read_pair(&hap, 400, 450, "test", &mut rng);
        assert!(
            result.is_some(),
            "R1 crossing boundary should produce a split-read pair"
        );
    }

    #[test]
    fn test_hap_read_pair_discordant_and_split() {
        let hap = mock_del_haplotype(500, 5000);
        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // Fragment at hap_start=200, frag_len=450 → R1=[200..350], R2=[500..650]
        // R1 in seg0, R2 in seg1 → discordant pair
        let result = gen.generate_haplotype_read_pair(&hap, 200, 450, "test", &mut rng);
        assert!(
            result.is_some(),
            "R1 in seg0, R2 in seg1 should produce discordant pair"
        );

        // Fragment crossing boundary → split-read evidence
        let result = gen.generate_haplotype_read_pair(&hap, 250, 450, "test2", &mut rng);
        assert!(result.is_some(), "Boundary-crossing pair should succeed");

        let result = gen.generate_haplotype_read_pair(&hap, 210, 450, "test3", &mut rng);
        assert!(result.is_some(), "Discordant pair spanning deletion");
    }

    #[test]
    fn test_hap_read_pair_safe_zones_accepted() {
        let hap = mock_del_haplotype(500, 5000);
        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // Both arms entirely in left segment
        let result = gen.generate_haplotype_read_pair(&hap, 0, 400, "left", &mut rng);
        assert!(result.is_some(), "Both arms in left segment should work");

        // Both arms entirely in right segment
        let result = gen.generate_haplotype_read_pair(&hap, 500, 400, "right", &mut rng);
        assert!(result.is_some(), "Both arms in right segment should work");
    }

    #[test]
    fn test_discordant_pair_has_large_ref_span() {
        // DEL: flank=500, del=5000
        // Left segment: ref[0..500), Right segment: ref[5500..6000)
        // Breakpoint in hap at offset 500
        let hap = mock_del_haplotype(500, 5000);
        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // Discordant placement: R1 in seg0, R2 in seg1
        // hap_start=200, frag_len=450 → R1=[200..350], R2=[500..650]
        let pair = gen
            .generate_haplotype_read_pair(&hap, 200, 450, "disc", &mut rng)
            .expect("Discordant pair should be generated");

        // R1 maps to ref: 0 + 200 = 200
        // R2 maps to ref: 5500 + 150 - 1 = 5649 (r2_hap_start=500, in seg1 at offset 0 → ref 5500)
        // ref_span should be much larger than frag_len (450) due to deletion gap
        let ref_span = pair.ref_end - pair.ref_start;
        assert!(
            ref_span > 5000,
            "Discordant pair ref span ({}) should exceed deletion size (5000)",
            ref_span
        );

        // Verify coordinates: ref_start should be in left segment region
        assert!(
            pair.ref_start < 500,
            "R1 ref_start should be before deletion"
        );
        // ref_end should be in right segment region
        assert!(pair.ref_end > 5500, "R2 ref_end should be after deletion");
    }

    #[test]
    fn test_hap_read_pair_boundary_crossing_always_allowed() {
        // Small-variant-like haplotype: left flank + 1bp alt + right flank.
        // Reads crossing segment boundaries are always allowed (produces
        // split-read evidence when aligned by BWA-MEM).
        let hap = VariantHaplotype::from_segments(vec![
            HaplotypeSegment {
                sequence: vec![b'A'; 300],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1000,
                    ref_end: 1300,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: vec![b'T'; 1],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1300,
                    ref_end: 1301,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
            HaplotypeSegment {
                sequence: vec![b'C'; 300],
                origin: Some(SegmentOrigin {
                    chrom: "chr1".to_string(),
                    ref_start: 1301,
                    ref_end: 1601,
                    is_reverse: false,
                }),
                hap_offset: 0,
            },
        ]);

        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // R1 [250..400) crosses the left->alt/right boundary — should succeed
        let result = gen.generate_haplotype_read_pair(&hap, 250, 300, "sv", &mut rng);
        assert!(
            result.is_some(),
            "boundary-crossing reads should be allowed for split-read evidence"
        );
    }

    #[test]
    fn test_hap_read_pair_with_read1_inside_insertion_is_kept() {
        // ref [0,2000) | 1000 bp insertion | ref [2000,4000).
        let seg = |sequence: Vec<u8>, origin: Option<(u64, u64)>| HaplotypeSegment {
            sequence,
            origin: origin.map(|(ref_start, ref_end)| SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start,
                ref_end,
                is_reverse: false,
            }),
            hap_offset: 0,
        };
        let hap = VariantHaplotype::from_segments(vec![
            seg(vec![b'A'; 2000], Some((0, 2000))),
            seg(vec![b'T'; 1000], None),
            seg(vec![b'A'; 2000], Some((2000, 4000))),
        ]);
        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // R1 = hap [2800,2950) inside the insertion; R2 ends at hap 3199 = ref 2199.
        let pair = gen
            .generate_haplotype_read_pair(&hap, 2800, 400, "ins", &mut rng)
            .expect("a pair with one read in the insertion is real evidence");

        // R1's start maps to the nearest reference base: hap 3000 = ref 2000.
        assert_eq!((pair.chrom.as_str(), pair.ref_start, pair.ref_end), ("chr1", 2000, 2200));
    }

    // ── Fragment orientation (F1R2 vs F2R1) ─────────────────────────────

    /// A generator using `profile` over `seq` on chr1.
    fn mock_gen_with_profile(
        profile: QualityProfile,
        seq: Vec<u8>,
        read_length: usize,
        indel_error_rate: f64,
    ) -> SynthReadGenerator<'static> {
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), seq);
        let reference: &'static SharedReference =
            Box::leak(Box::new(SharedReference::from_sequences(seqs)));
        SynthReadGenerator::new(profile, reference, read_length, indel_error_rate)
    }

    /// A generator over `seq` on chr1, with every R1 base at Q`q1` and every
    /// R2 base at Q`q2`, and `indel_error_rate` of its errors indels.
    fn mock_gen_over_with_indels(
        seq: Vec<u8>,
        read_length: usize,
        q1: u8,
        q2: u8,
        indel_error_rate: f64,
    ) -> SynthReadGenerator<'static> {
        let qual1 = vec![b'!' + q1; read_length];
        let qual2 = vec![b'!' + q2; read_length];
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| {
                mock_read_pair(&format!("mock_{}", i), qual1.clone(), qual2.clone(), i as u64 * 500)
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, read_length);
        mock_gen_with_profile(profile, seq, read_length, indel_error_rate)
    }

    /// A generator over `seq` on chr1, with every R1 base at Q`q1` and every
    /// R2 base at Q`q2`, and no indel sequencing errors.
    fn mock_gen_over(seq: Vec<u8>, read_length: usize, q1: u8, q2: u8) -> SynthReadGenerator<'static> {
        mock_gen_over_with_indels(seq, read_length, q1, q2, 0.0)
    }

    /// A sequence with no palindromic structure, so a read off the left end of
    /// a fragment never equals the reverse complement of the right end.
    fn scrambled_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..len).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
    }

    /// True when at least 90% of `seq` is `base`. Over an all-A reference an
    /// R1 is all A forward and all T reverse, bar the odd sequencing error.
    fn mostly(seq: &[u8], base: u8) -> bool {
        seq.iter().filter(|&&b| b == base).count() * 10 >= seq.len() * 9
    }

    #[test]
    fn test_read_pair_orientation_splits_near_half() {
        // Real Illumina libraries are ~50/50 F1R2 / F2R1 (M13).
        let gen = mock_gen_over(vec![b'A'; 100_000], 150, 40, 40);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(1);
        let (mut fwd, mut rev) = (0usize, 0usize);
        for i in 0..400u64 {
            let pair = gen
                .generate_read_pair("chr1", 1000 + i * 100, 400, &no_alleles, "p", &mut rng)
                .unwrap();
            if mostly(&pair.seq1, b'A') {
                fwd += 1;
            } else if mostly(&pair.seq1, b'T') {
                rev += 1;
            }
        }
        assert_eq!(fwd + rev, 400, "every R1 should be one orientation or the other");
        assert!(fwd > 140 && rev > 140, "orientation split {} fwd / {} rev of 400", fwd, rev);
    }

    #[test]
    fn test_reversed_read_pair_covers_the_same_fragment() {
        // Flipping orientation must move which end each mate comes from, not
        // where the fragment sits: R1 reverse-complemented off the right end,
        // R2 forward off the left end.
        let rl = 150usize;
        let seq = scrambled_seq(50_000, 7);
        // Q93: no sequencing errors, but indel_error_rate > 0 turns on slack
        // (INDEL_SLACK past the 3' end) so the assertions below actually
        // depend on the slack landing on the correct side of the read.
        let gen = mock_gen_over_with_indels(seq.clone(), rl, 93, 93, 1.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(2);
        let (mut fwd, mut rev) = (0usize, 0usize);
        for i in 0..100u64 {
            let frag_start = 1000 + i * 37;
            let frag_len = 400u64;
            let pair = gen
                .generate_read_pair("chr1", frag_start, frag_len, &no_alleles, "p", &mut rng)
                .unwrap();
            let (s, e) = (frag_start as usize, (frag_start + frag_len) as usize);
            let left = seq[s..s + rl].to_vec();
            let mut right_rc = seq[e - rl..e].to_vec();
            reverse_complement(&mut right_rc);
            assert_ne!(left, right_rc, "the test sequence must not be palindromic");
            if pair.seq1 == left {
                assert_eq!(pair.seq2, right_rc, "F1R2 pair at {} has the wrong R2", frag_start);
                fwd += 1;
            } else if pair.seq1 == right_rc {
                assert_eq!(pair.seq2, left, "F2R1 pair at {} has the wrong R2", frag_start);
                rev += 1;
            } else {
                panic!("R1 at {} matches neither end of the fragment", frag_start);
            }
        }
        assert!(fwd > 0 && rev > 0, "orientation split {} fwd / {} rev of 100", fwd, rev);
    }

    #[test]
    fn test_read_pair_orientation_repeats_for_the_same_seed() {
        // The orientation draw comes out of the seeded stream in a fixed
        // order, so the same --seed still gives the same FASTQ (M7).
        let gen = mock_gen_over(vec![b'A'; 100_000], 150, 40, 40);
        let orientations = |seed: u64| -> Vec<bool> {
            let no_alleles = HashMap::new();
            let mut rng = StdRng::seed_from_u64(seed);
            (0..200u64)
                .map(|i| {
                    let pair = gen
                        .generate_read_pair("chr1", 1000 + i * 100, 400, &no_alleles, "p", &mut rng)
                        .unwrap();
                    mostly(&pair.seq1, b'T')
                })
                .collect()
        };

        let first = orientations(9);
        assert_eq!(first, orientations(9), "the same seed must give the same orientations");
        assert!(
            first.iter().any(|&r| r) && first.iter().any(|&r| !r),
            "both orientations should occur"
        );
    }

    #[test]
    fn test_reversed_r1_keeps_the_r1_quality_model() {
        // R1 at Q40, R2 at Q20: the mate that is R1 keeps the R1 model even
        // when it comes off the right end of the fragment.
        let gen = mock_gen_over(vec![b'A'; 100_000], 150, 40, 20);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(3);
        let mut reversed = 0usize;
        for i in 0..200u64 {
            let pair = gen
                .generate_read_pair("chr1", 1000 + i * 100, 400, &no_alleles, "p", &mut rng)
                .unwrap();
            assert!(pair.qual1.iter().all(|&q| q == b'!' + 40), "R1 lost the R1 quality model");
            assert!(pair.qual2.iter().all(|&q| q == b'!' + 20), "R2 lost the R2 quality model");
            if mostly(&pair.seq1, b'T') {
                reversed += 1;
            }
        }
        assert!(reversed > 0, "no pair came out in F2R1 orientation");
    }

    #[test]
    fn test_reversed_haplotype_pair_covers_the_same_fragment() {
        // Same property for the haplotype/tiling path.
        let rl = 150usize;
        let hap_seq = scrambled_seq(4000, 13);
        let hap = VariantHaplotype::from_segments(vec![HaplotypeSegment {
            sequence: hap_seq.clone(),
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: 0,
                ref_end: 4000,
                is_reverse: false,
            }),
            hap_offset: 0,
        }]);
        // This path reads the haplotype, not the reference; Q93: no
        // sequencing errors, but indel_error_rate > 0 turns on slack so the
        // assertions below depend on it landing on the correct side.
        let gen = mock_gen_over_with_indels(vec![b'A'; 1000], rl, 93, 93, 1.0);
        let mut rng = StdRng::seed_from_u64(4);
        let (mut fwd, mut rev) = (0usize, 0usize);
        for i in 0..100u64 {
            let hap_start = 100 + i * 17;
            let frag_len = 400u64;
            let pair = gen
                .generate_haplotype_read_pair(&hap, hap_start, frag_len, "h", &mut rng)
                .unwrap();
            let (s, e) = (hap_start as usize, (hap_start + frag_len) as usize);
            let left = hap_seq[s..s + rl].to_vec();
            let mut right_rc = hap_seq[e - rl..e].to_vec();
            reverse_complement(&mut right_rc);
            assert_ne!(left, right_rc, "the test haplotype must not be palindromic");
            if pair.seq1 == left {
                assert_eq!(pair.seq2, right_rc, "F1R2 pair at {} has the wrong R2", hap_start);
                fwd += 1;
            } else if pair.seq1 == right_rc {
                assert_eq!(pair.seq2, left, "F2R1 pair at {} has the wrong R2", hap_start);
                rev += 1;
            } else {
                panic!("R1 at hap {} matches neither end of the fragment", hap_start);
            }
            assert_eq!(
                (pair.ref_start, pair.ref_end),
                (hap_start, hap_start + frag_len),
                "the pair's reference span must not depend on orientation"
            );
        }
        assert!(fwd > 0 && rev > 0, "orientation split {} fwd / {} rev of 100", fwd, rev);
    }

    // ── Sequencing-order generation of the reverse mate (L1, L17) ───────

    /// Quality that starts at Q30 and, once it drops to Q10, never recovers.
    /// Only a Markov chain run in sequencing order reproduces that shape.
    fn decaying_qual(read_length: usize, rng: &mut StdRng) -> Vec<u8> {
        let mut qual = Vec::with_capacity(read_length);
        let mut q = b'!' + 30;
        for _ in 0..read_length {
            qual.push(q);
            if q > b'!' + 10 && rng.gen::<f64>() < 0.01 {
                q = b'!' + 10;
            }
        }
        qual
    }

    #[test]
    fn test_indel_deletion_errors_never_pad_a_read_with_n() {
        // A deletion sequencing error consumed a template base without
        // emitting one, so the read ran off the end of its template and the
        // shortfall was padded with `N` (L1) — at the 3' end of the forward
        // mate and, after the reverse-complement, at the 5' end of the
        // reverse one. Every base must come from real template sequence.
        let rl = 150usize;
        let hap_seq = scrambled_seq(4000, 21);
        let hap = VariantHaplotype::from_segments(vec![HaplotypeSegment {
            sequence: hap_seq,
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: 0,
                ref_end: 4000,
                is_reverse: false,
            }),
            hap_offset: 0,
        }]);
        // Q20 → 1% of bases mis-called, every error an indel: ~0.75 deletion
        // errors per read, so most reads used to need padding.
        let gen = mock_gen_over_with_indels(vec![b'A'; 1000], rl, 20, 20, 1.0);
        let mut rng = StdRng::seed_from_u64(5);
        let mut padded = 0usize;
        for i in 0..200u64 {
            let pair = gen
                .generate_haplotype_read_pair(&hap, 500 + i * 7, 400, "h", &mut rng)
                .unwrap();
            for seq in [&pair.seq1, &pair.seq2] {
                assert_eq!(seq.len(), rl, "a synthetic read came out short");
                if seq.contains(&b'N') {
                    padded += 1;
                }
            }
        }
        assert_eq!(padded, 0, "{} of 400 synthetic reads carry padded `N` bases", padded);
    }

    #[test]
    fn test_reverse_read_n_pads_past_the_contig_end_instead_of_shifting() {
        // SharedReference::fetch_sequence clamps `fetch_end` to the contig
        // length, so a reverse read whose claimed span overhangs the contig
        // end gets back a window shorter than requested. Reversing that
        // clamped-short window (instead of clamping the intended window
        // *before* reversing) put the contig's real last bases at the read's
        // 5' end — a full-length, N-free read silently shifted left of its
        // claimed coordinates, while ref_end still points past the contig.
        let rl = 150usize;
        let ref_start = 1000u64;
        let overhang = 5usize; // within 1..INDEL_SLACK
        let contig_len = ref_start as usize + rl - overhang;
        let seq = scrambled_seq(contig_len, 23);
        let gen = mock_gen_over_with_indels(seq.clone(), rl, 93, 93, 1.0); // Q93, slack on
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(17);

        let (read_seq, _qual) = gen.generate_read("chr1", ref_start, &no_alleles, 1, true, &mut rng);

        assert_eq!(read_seq.len(), rl, "a synthetic read came out short");
        let n_count = read_seq.iter().filter(|&&b| b == b'N').count();
        assert_eq!(
            n_count, overhang,
            "expected {} N bases for the part of the read past the contig end, got {}",
            overhang, n_count
        );
        assert!(
            read_seq[..overhang].iter().all(|&b| b == b'N'),
            "the missing bases belong at the read's 5' start (past the contig end), not scattered: {:?}",
            read_seq
        );

        // The rest of the read must be the real revcomp of the reference span
        // actually covered, not bases pulled from a window shifted to make up
        // the length.
        let mut expected_rc = seq[ref_start as usize..contig_len].to_vec();
        reverse_complement(&mut expected_rc);
        assert_eq!(
            &read_seq[overhang..],
            expected_rc.as_slice(),
            "real bases must come from the actually-covered span, not a shifted window"
        );
    }

    #[test]
    fn test_reference_n_bases_get_q2() {
        // A real sequencer reports a no-call as `N` at Q2. The learned
        // profile knows nothing about `N`, so it handed a reference `N` an
        // ordinary score (often Q37) and the read claimed a base it does not
        // have with high confidence (L18).
        let rl = 50usize;
        let mut seq = scrambled_seq(1000, 29);
        seq[520..540].fill(b'N');
        let gen = mock_gen_over(seq, rl, 37, 37);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(3);

        let (read_seq, read_qual) = gen.generate_read("chr1", 500, &no_alleles, 1, false, &mut rng);

        let n_cycles: Vec<usize> = (0..rl).filter(|&i| read_seq[i] == b'N').collect();
        assert_eq!(n_cycles.len(), 20, "expected the reference `N` run inside the read");
        for c in n_cycles {
            assert_eq!(
                read_qual[c],
                b'!' + 2,
                "reference `N` at cycle {} reported Q{}, not Q2",
                c,
                read_qual[c] - b'!'
            );
        }
    }

    #[test]
    fn test_contig_end_padding_n_gets_q2() {
        // Same rule for the `N` that stands in for bases past the contig end
        // on a reverse read: whatever puts an `N` in a read, it is a no-call
        // and reports Q2 (L18).
        let rl = 150usize;
        let ref_start = 1000u64;
        let overhang = 5usize;
        let contig_len = ref_start as usize + rl - overhang;
        let gen = mock_gen_over(scrambled_seq(contig_len, 23), rl, 37, 37);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(17);

        let (read_seq, read_qual) = gen.generate_read("chr1", ref_start, &no_alleles, 1, true, &mut rng);

        assert!(read_seq[..overhang].iter().all(|&b| b == b'N'), "expected padding at the 5' end");
        for (c, &q) in read_qual[..overhang].iter().enumerate() {
            assert_eq!(q, b'!' + 2, "contig-end padding `N` at cycle {} reported Q{}, not Q2", c, q - b'!');
        }
    }

    #[test]
    fn test_template_exhaustion_padding_n_gets_q2() {
        // The third way an `N` reaches a read: the template ran out, so the
        // tail is padded. Locks the third arm of the one rule so the three
        // cannot drift apart again (L18).
        let rl = 150usize;
        let ref_start = 1000u64;
        let short_by = 7usize;
        let contig_len = ref_start as usize + rl - short_by;
        let gen = mock_gen_over(scrambled_seq(contig_len, 31), rl, 37, 37);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(11);

        let (read_seq, read_qual) =
            gen.generate_read("chr1", ref_start, &no_alleles, 1, false, &mut rng);

        assert_eq!(read_seq.len(), rl, "a synthetic read came out short");
        assert!(read_seq[rl - short_by..].iter().all(|&b| b == b'N'), "expected a padded 3' tail");
        for (c, &q) in read_qual.iter().enumerate().skip(rl - short_by) {
            assert_eq!(q, b'!' + 2, "padding `N` at cycle {} reported Q{}, not Q2", c, q - b'!');
        }
    }

    #[test]
    fn test_reverse_read_from_seq_complements_lowercase_before_uppercasing() {
        // A reference FASTA's lowercase bases mark soft-masked repeats, a
        // legitimate input. `generate_read_from_seq` must complement first
        // and uppercase after (like `generate_read`), not the other way
        // round: complementing an un-uppercased 'a' the wrong way round
        // would silently turn it into 'A' instead of 'T'.
        let rl = 10usize;
        let gen = mock_gen_over(vec![b'A'; 1000], rl, 93, 93); // Q93: no errors
        let mut rng = StdRng::seed_from_u64(99);
        let seq = b"acgtacgtac".to_vec();

        let (read_seq, _qual) = gen.generate_read_from_seq(&seq, 1, true, &mut rng);

        // Hand-computed revcomp (not via `reverse_complement`, so the test
        // doesn't just check the function against itself): "acgtacgtac"
        // reversed is "catgcatgca", complemented is "GTACGTACGT".
        assert_eq!(read_seq, b"GTACGTACGT", "lowercase input must be complemented, not passed through unchanged");
    }

    #[test]
    fn test_reverse_mate_quality_decays_in_sequencing_order() {
        // The reverse mate was generated along the reference and reversed
        // afterwards, which ran its quality Markov chain backwards (L17):
        // each cycle was conditioned on the cycle after it, not before it.
        let rl = 150usize;
        let mut tr = StdRng::seed_from_u64(11);
        let pairs: Vec<ReadPair> = (0..2000)
            .map(|i| {
                let q1 = decaying_qual(rl, &mut tr);
                let q2 = decaying_qual(rl, &mut tr);
                mock_read_pair(&format!("t_{}", i), q1, q2, i as u64 * 500)
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let gen = mock_gen_with_profile(profile, vec![b'A'; 100_000], rl, 0.0);

        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(12);
        let half = rl / 2;
        let (mut fwd_head, mut fwd_tail) = (Vec::new(), Vec::new());
        let (mut rev_head, mut rev_tail) = (Vec::new(), Vec::new());
        for i in 0..200u64 {
            let pair = gen
                .generate_read_pair("chr1", 1000 + i * 100, 400, &no_alleles, "p", &mut rng)
                .unwrap();
            // Over an all-A reference the reverse mate reads as T.
            let ts = pair.seq1.iter().filter(|&&b| b == b'T').count();
            let as_ = pair.seq1.iter().filter(|&&b| b == b'A').count();
            let (rev_q, fwd_q) =
                if ts > as_ { (&pair.qual1, &pair.qual2) } else { (&pair.qual2, &pair.qual1) };
            // Cycle 0 has no Markov predecessor either way — start at cycle 1.
            fwd_head.extend_from_slice(&fwd_q[1..half]);
            fwd_tail.extend_from_slice(&fwd_q[half..]);
            rev_head.extend_from_slice(&rev_q[1..half]);
            rev_tail.extend_from_slice(&rev_q[half..]);
        }
        let (fh, ft) = (mean_qual(&fwd_head), mean_qual(&fwd_tail));
        let (rh, rt) = (mean_qual(&rev_head), mean_qual(&rev_tail));
        assert!(fh > ft + 3.0, "forward mate quality should decay: head {:.2} tail {:.2}", fh, ft);
        assert!(
            rh > rt + 3.0,
            "reverse mate quality should decay in sequencing order like the forward mate \
             (forward head {:.2} tail {:.2}): head {:.2} tail {:.2}",
            fh,
            ft,
            rh,
            rt
        );
    }
}
