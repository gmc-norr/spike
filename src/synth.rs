//! Synthetic read generation from reference sequence + learned quality profile.
//!
//! Instead of cloning real reads (which produces exact duplicates flagged by dedup tools),
//! this module generates independent synthetic reads that have realistic quality scores
//! and correlated sequencing errors.
//!
//! The approach:
//! 1. Learn a `QualityProfile` (`quality.rs`) from the donor reads: their
//!    qualities, and the errors they make, counted against the reference
//! 2. Generate synthetic reads by drawing each base's quality from the profile,
//!    reading the template base, and making an error at the rate the profile
//!    gives for that quality in that kind of read

use std::collections::HashMap;
use std::sync::Arc;

use rand::rngs::StdRng;
use rand::Rng;

use crate::extract::reverse_complement;
use crate::haplotype::VariantHaplotype;
use crate::loh::SampleCopies;
use crate::read_name::NameShape;
use crate::reference::SharedReference;
use crate::stats::FragmentDist;
use crate::types::{ReadPair, ReadPool};

/// Fewest sampled pairs the quality model is known to be good from: with
/// 50,000 pairs of HG002 35x it made 0.99% crashed reads, with 120,000 0.94%
/// (docs/superpowers/plans/2026-10-07-quality-model-v2.md).
const MIN_PROFILE_PAIRS: usize = 50_000;

/// The warning for a quality profile learned from `pairs` sampled pairs with
/// the context census `census`, or `None` when there are enough pairs.
pub(crate) fn thin_profile_warning(pairs: usize, census: &str) -> Option<String> {
    (pairs < MIN_PROFILE_PAIRS).then(|| {
        format!(
            "Quality profile learned from {} sampled pairs; it was measured good from {} \
             (HG002 35x), and with fewer its crashed reads and their errors are learned \
             from thinner counts. {}.",
            pairs, MIN_PROFILE_PAIRS, census
        )
    })
}

/// Extra template bases fetched past a read's 3' end when indel errors are on,
/// so a deletion error is covered by real sequence instead of `N` padding (L1).
const INDEL_SLACK: usize = 10;

/// Quality reported for every `N` a synthetic read emits, whatever put it
/// there — a reference `N`, padding past a contig end, or padding after the
/// template ran out. An `N` is a no-call, and a real Illumina no-call is
/// always Q2; the learned profile knows nothing about `N` and would hand one
/// an ordinary score (often Q37) instead (L18).
const N_QUAL: u8 = b'!' + 2; // Q2

pub use crate::quality::{ClassMix, QualityProfile};

/// Generates synthetic reads from reference sequence + learned quality profile.
pub struct SynthReadGenerator<'a> {
    /// The run's quality model, learned once from the startup sample.
    profile: Arc<QualityProfile>,
    /// The event's own read-class mix (`QualityProfile::class_mix`); `None`
    /// draws the sample's.
    class_mix: Option<ClassMix>,
    reference: &'a SharedReference,
    read_length: usize,
    /// Fraction of sequencing errors that are indels (vs substitutions).
    /// 0.0 = substitution-only (default), ~0.05 = typical Illumina.
    indel_error_rate: f64,
    /// Whether the input library's reads were adapter-trimmed before
    /// alignment. Then a read is `read_length` cycles, cut to its fragment
    /// when the fragment is shorter, and loses any 3' end that spells the
    /// adapter's start (see `adapter_suffix_len`). Otherwise every read is
    /// `read_length` long and no fragment is shorter than one read.
    adapter_trimmed: bool,
    /// The shape of the input's read names, which this generator's reads are
    /// named in (`read_name`).
    read_names: NameShape,
}

impl<'a> SynthReadGenerator<'a> {
    pub fn new(
        profile: impl Into<Arc<QualityProfile>>,
        reference: &'a SharedReference,
        read_length: usize,
        indel_error_rate: f64,
    ) -> Self {
        Self {
            profile: profile.into(),
            class_mix: None,
            reference,
            read_length,
            indel_error_rate,
            adapter_trimmed: false,
            read_names: NameShape::Other,
        }
    }

    /// This generator for a library whose reads were (`on`) or were not
    /// adapter-trimmed; `new` builds an untrimmed one.
    pub fn with_adapter_trim(mut self, on: bool) -> Self {
        self.adapter_trimmed = on;
        self
    }

    /// This generator drawing read classes with the event's own `mix`; `new`
    /// draws the sample's.
    pub fn with_class_mix(mut self, mix: ClassMix) -> Self {
        self.class_mix = Some(mix);
        self
    }

    /// This generator naming its reads in `shape`; `new` uses `Other`.
    pub fn with_read_names(mut self, shape: NameShape) -> Self {
        self.read_names = shape;
        self
    }

    /// The name of this generator's read `internal` (`NameShape::name`).
    pub fn read_name(&self, internal: &str) -> String {
        self.read_names.name(internal)
    }

    /// The shortest fragment this generator sequences: any, in a trimmed
    /// library, where a short fragment's reads are cut to it; one read's
    /// length otherwise.
    pub fn min_fragment_len(&self) -> i64 {
        if self.adapter_trimmed {
            1
        } else {
            self.read_length as i64
        }
    }

    /// How long each mate of a `frag_len` fragment is sequenced, before
    /// adapter trimming; `None` when no read can come off it.
    fn mate_length(&self, frag_len: u64) -> Option<u64> {
        let cycles = self.read_length as u64;
        let rl = if self.adapter_trimmed { cycles.min(frag_len) } else { cycles };
        (rl > 0 && frag_len >= rl).then_some(rl)
    }

    /// Cut a read's 3' end where it spells the adapter's start, as an
    /// adapter-trimmed library's were -- only when the fragment ran past both
    /// reads, so no adapter was read in full.
    fn trim_adapter_start(&self, frag_len: u64, seq: &mut Vec<u8>, qual: &mut Vec<u8>) {
        if self.adapter_trimmed && frag_len >= self.read_length as u64 {
            let keep = seq.len() - adapter_suffix_len(seq);
            seq.truncate(keep);
            qual.truncate(keep);
        }
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
    /// Generating in sequencing order also keeps the quality model's history
    /// running with the cycle counter for both mates, not against it (L17).
    ///
    /// `class` is the read's class from the pair's joint draw
    /// (`QualityProfile::draw_classes`); `None` draws one for this mate alone.
    fn generate_from_template(
        &self,
        template: &[u8],
        rl: usize,
        read_num: u8,
        class: Option<usize>,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut seq = Vec::with_capacity(rl);
        let mut qual = Vec::with_capacity(rl);
        let mut idx = 0usize; // current position in template
        let mut state = self.profile.start_read(read_num, class, rng);

        while seq.len() < rl && idx < template.len() {
            let c = seq.len(); // cycle position in the read
            let true_base = template[idx].to_ascii_uppercase();

            let q = self.profile.next_quality(read_num, &state, rl - c, rng);

            if true_base == b'N' {
                // An `N` is a no-call and reports Q2, not the profile's `q`
                // (L18). `q` is still drawn, so the quality draw stays one per
                // template base -- but the `p_err` draw below is skipped, so
                // an `N` consumes strictly fewer random numbers than a called
                // base and does not leave the stream unchanged.
                //
                // The history holds what the read *emitted* (the
                // indel-deletion branch below stays put for the same reason:
                // nothing was emitted), so Q2 enters it rather than the `q`
                // thrown away.
                seq.push(b'N');
                qual.push(N_QUAL);
                self.profile.emitted(&mut state, N_QUAL, b'N');
                self.profile.record_error(&mut state, false);
                idx += 1;
                continue;
            }

            let p_err = self.profile.error_rate(q, &state, rl - c);

            if rng.gen::<f64>() < p_err {
                if self.indel_error_rate > 0.0 && rng.gen::<f64>() < self.indel_error_rate {
                    // Indel error: 50/50 insertion vs deletion.
                    if rng.gen::<bool>() {
                        // Insertion: add a random base without consuming template.
                        seq.push(random_base(rng));
                        qual.push(q);
                        self.profile.emitted(&mut state, q, true_base);
                        self.profile.record_error(&mut state, true);
                        // Don't advance idx — the template base is read next cycle.
                    } else {
                        // Deletion: skip this template base entirely.
                        idx += 1;
                        // Don't add to seq/qual — next iteration reads the next base.
                        // Nothing was emitted, so the history -- qualities,
                        // run and errors -- stays as it is.
                    }
                } else {
                    // Substitution error.
                    seq.push(random_different_base(true_base, rng));
                    qual.push(q);
                    self.profile.emitted(&mut state, q, true_base);
                    self.profile.record_error(&mut state, true);
                    idx += 1;
                }
            } else {
                seq.push(true_base);
                qual.push(q);
                self.profile.emitted(&mut state, q, true_base);
                self.profile.record_error(&mut state, false);
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
    ///   thing. Reads are `rl` long either way.
    #[allow(clippy::too_many_arguments)]
    #[cfg(test)]
    fn generate_read(
        &self,
        chrom: &str,
        ref_start: u64,
        alleles: &HashMap<u64, u8>,
        read_num: u8,
        is_reverse: bool,
        rl: usize,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let template = self.template_at(chrom, ref_start, alleles, is_reverse, rl);
        self.generate_from_template(&template, rl, read_num, None, rng)
    }

    /// The template of a read at `ref_start` (see `generate_read`), in
    /// sequencing order, with any indel slack past its 3' end.
    fn template_at(
        &self,
        chrom: &str,
        ref_start: u64,
        alleles: &HashMap<u64, u8>,
        is_reverse: bool,
        rl: usize,
    ) -> Vec<u8> {
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
        if is_reverse {
            let missing = (fetch_end - fetch_start).saturating_sub(ref_seq.len() as u64) as usize;
            std::iter::repeat_n(b'N', missing)
                .chain((0..ref_seq.len()).rev().map(|i| complement(base_at(i))))
                .collect()
        } else {
            (0..ref_seq.len()).map(base_at).collect()
        }
    }

    /// Both mates' classes in one draw, so a poor pair is poor in both, given
    /// the runs of one base each mate's template holds (`fwd`, `rev`: the
    /// forward and reverse mate's, first `rl` bases) and the event's mix.
    fn draw_pair_classes(&self, fwd: &[u8], rev: &[u8], rl: usize, r1_is_reverse: bool, rng: &mut StdRng) -> (usize, usize) {
        let run = |t: &[u8]| crate::quality::template_run_bin(&t[..rl.min(t.len())]);
        let (r1, r2) = if r1_is_reverse { (run(rev), run(fwd)) } else { (run(fwd), run(rev)) };
        self.profile.draw_classes((r1, r2), self.class_mix.as_ref(), rng)
    }

    /// Generate a synthetic read pair for a fragment at a given position.
    ///
    /// The fragment is `[frag_start, frag_start + frag_len)` either way; a
    /// coin flip decides which mate comes off which end. F1R2: R1 forward at
    /// `frag_start`, R2 reverse at `frag_start + frag_len - rl`, where `rl`
    /// is `read_length` -- or, in an adapter-trimmed library, the fragment's
    /// length when that is shorter (see `mate_length`, `trim_adapter_start`).
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
        let rl = self.mate_length(frag_len)?;

        let right_start = frag_start + frag_len - rl;

        // Draw the orientation before either read, so the stream stays in a
        // fixed order and the same --seed keeps giving the same output (M7).
        let r1_is_reverse = rng.gen::<bool>();
        // Either way the left end is read forward and the right end reverse —
        // the flip only decides which of the two is R1. read_num picks the
        // quality model, is_reverse handles the flip, so whichever mate
        // ends up as R1 keeps the R1 model.
        let (fwd_num, rev_num) = if r1_is_reverse { (2, 1) } else { (1, 2) };
        let fwd_template = self.template_at(chrom, frag_start, alleles, false, rl as usize);
        let rev_template = self.template_at(chrom, right_start, alleles, true, rl as usize);
        let classes = self.draw_pair_classes(&fwd_template, &rev_template, rl as usize, r1_is_reverse, rng);
        let class_of = |num: u8| Some(if num == 1 { classes.0 } else { classes.1 });

        // Forward mate (left end of the fragment).
        let (mut fwd_seq, mut fwd_qual) =
            self.generate_from_template(&fwd_template, rl as usize, fwd_num, class_of(fwd_num), rng);

        // Reverse mate (right end) — generated in sequencing order, so it
        // already comes back in FASTQ orientation.
        let (mut rev_seq, mut rev_qual) =
            self.generate_from_template(&rev_template, rl as usize, rev_num, class_of(rev_num), rng);

        self.trim_adapter_start(frag_len, &mut fwd_seq, &mut fwd_qual);
        self.trim_adapter_start(frag_len, &mut rev_seq, &mut rev_qual);

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
            align: None,
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
        // Sample new fragment length.
        let new_frag = frag_dist.sample_in_range(
            rng,
            self.min_fragment_len(),
            crate::stats::MAX_FRAGMENT_LEN,
        ) as u64;

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
                    &self.read_name(&format!("{}_dup_depth_{:06}", name_prefix, i)),
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
    /// The read is `min(rl, seq.len())` bases long.
    ///
    /// `read_num`: 1 or 2 (controls which quality distribution is sampled).
    /// `is_reverse`: when true the read comes off the reverse strand — `seq`
    ///   is reverse-complemented first and the read is returned in FASTQ
    ///   orientation (see `generate_read` docs).
    #[cfg(test)]
    fn generate_read_from_seq(
        &self,
        seq: &[u8],
        rl: usize,
        read_num: u8,
        is_reverse: bool,
        rng: &mut StdRng,
    ) -> (Vec<u8>, Vec<u8>) {
        let rl = rl.min(seq.len());
        self.generate_from_template(&template_from_seq(seq, is_reverse), rl, read_num, None, rng)
    }

    /// Generate a synthetic read pair from a variant haplotype.
    ///
    /// A coin flip decides which end of the fragment each mate comes off, as
    /// in `generate_read_pair` (M13): the forward mate starts at
    /// `hap_frag_start`, the reverse mate at `hap_frag_start + frag_len - rl`
    /// and is reverse-complemented for FR orientation. Reference coordinates
    /// for the ReadPair are mapped back from the haplotype via `hap_to_ref`.
    /// Also returns where the pair lies on the haplotype ([`PairSpans`]).
    pub fn generate_haplotype_read_pair(
        &self,
        haplotype: &VariantHaplotype,
        hap_frag_start: u64,
        frag_len: u64,
        name: &str,
        rng: &mut StdRng,
    ) -> Option<(ReadPair, PairSpans)> {
        let rl = self.mate_length(frag_len)?;
        let rl_bases = rl as usize;

        let right_hap_start = hap_frag_start + frag_len - rl;

        // Get sequences from the haplotype, with extra bases past each read's
        // 3' end so deletion errors read real sequence, not `N` padding (L1).
        // For the reverse mate that end is to the left of `right_hap_start`.
        let slack = if self.indel_error_rate > 0.0 { INDEL_SLACK } else { 0 };
        let right_slack = slack.min(right_hap_start as usize);
        let left_seq = haplotype.get_sequence(hap_frag_start, rl_bases + slack);
        let right_seq = haplotype.get_sequence(
            right_hap_start - right_slack as u64,
            right_slack + rl_bases,
        );

        if left_seq.len() < rl_bases || right_seq.len() < right_slack + rl_bases {
            return None;
        }

        // Draw the orientation before either read, so the stream stays in a
        // fixed order and the same --seed keeps giving the same output (M7).
        let r1_is_reverse = rng.gen::<bool>();
        // read_num picks the quality model, is_reverse handles the flip,
        // so whichever mate ends up as R1 keeps the R1 model.
        let (fwd_num, rev_num) = if r1_is_reverse { (2, 1) } else { (1, 2) };
        let fwd_template = template_from_seq(left_seq, false);
        let rev_template = template_from_seq(right_seq, true);
        let classes = self.draw_pair_classes(&fwd_template, &rev_template, rl_bases, r1_is_reverse, rng);
        let class_of = |num: u8| Some(if num == 1 { classes.0 } else { classes.1 });

        // Forward mate (left end of the fragment).
        let (mut fwd_seq, mut fwd_qual) = self.generate_from_template(
            &fwd_template,
            rl_bases.min(fwd_template.len()),
            fwd_num,
            class_of(fwd_num),
            rng,
        );

        // Reverse mate (right end): generated in sequencing order, so it
        // already comes back in FASTQ orientation.
        let (mut rev_seq, mut rev_qual) = self.generate_from_template(
            &rev_template,
            rl_bases.min(rev_template.len()),
            rev_num,
            class_of(rev_num),
            rng,
        );

        self.trim_adapter_start(frag_len, &mut fwd_seq, &mut fwd_qual);
        self.trim_adapter_start(frag_len, &mut rev_seq, &mut rev_qual);

        // Each mate's 3' end was trimmed on its own, so each covers its own
        // length from its end of the fragment.
        let frag_end = hap_frag_start + frag_len;
        let spans = PairSpans {
            fragment: (hap_frag_start, frag_end),
            mates: [
                (hap_frag_start, hap_frag_start + fwd_seq.len() as u64),
                (frag_end - rev_seq.len() as u64, frag_end),
            ],
        };

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

        Some((
            ReadPair {
                name: name.to_string(),
                seq1,
                qual1,
                seq2,
                qual2,
                ref_start,
                ref_end,
                insert_size: frag_len as i64,
                chrom,
                align: None,
            },
            spans,
        ))
    }
}

/// A haplotype slice `seq` as its read meets it, in sequencing order: the
/// reverse mate reads the other strand, the 3' end of `seq` first.
fn template_from_seq(seq: &[u8], is_reverse: bool) -> Vec<u8> {
    let mut template = seq.to_vec();
    if is_reverse {
        reverse_complement(&mut template);
    }
    template
}

/// Where a planted pair lies on its haplotype, each `[start, end)`: the
/// fragment, and the bases each mate covers by its final length, the forward
/// mate first. A read with an indel error covers a base more or less.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PairSpans {
    pub fragment: (u64, u64),
    pub mates: [(u64, u64); 2],
}

/// The start of the TruSeq adapter, R1's and R2's alike. A library trimmed
/// before alignment loses any read end that spells a prefix of it: measured
/// on the HG001 and HG002 30x BAMs, 99.3-100% of 147-150 bp reads are
/// followed by exactly the bases they lost, and no 151 bp read ends in A
/// (docs/superpowers/plans/2026-09-28-read-length.md).
const ADAPTER_START: &[u8] = b"AGATCGGAAGAGC";

/// The longest end of `seq` that equals a prefix of `ADAPTER_START`: what
/// the trimmer cut from a read ending there. 0 when none does.
pub(crate) fn adapter_suffix_len(seq: &[u8]) -> usize {
    (1..=ADAPTER_START.len().min(seq.len()))
        .rev()
        .find(|&k| seq[seq.len() - k..] == ADAPTER_START[..k])
        .unwrap_or(0)
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
#[cfg(test)]
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
            align: None,
        }
    }


    fn on_threads<R: Send>(threads: usize, f: impl FnOnce() -> R + Send) -> R {
        rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(f)
    }



    #[test]
    fn test_quality_profile_from_pairs() {
        // R1 learns from the donors' R1 strings and R2 from their R2 strings.
        let rl = 10;
        let q30 = vec![b'!' + 30; rl];
        let q20 = vec![b'!' + 20; rl];
        let pairs: Vec<ReadPair> = (0..100)
            .map(|i| mock_read_pair(&format!("read_{}", i), q30.clone(), q20.clone(), i * 10))
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let gen = mock_gen_with_profile(profile, vec![b'A'; 10_000], rl, 0.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(5);
        for i in 0..50u64 {
            let p = gen.generate_read_pair("chr1", 100 + i * 20, 40, &no_alleles, "p", &mut rng).unwrap();
            assert!(p.qual1.iter().all(|&q| q == b'!' + 30), "R1 should be all Q30: {:?}", p.qual1);
            assert!(p.qual2.iter().all(|&q| q == b'!' + 20), "R2 should be all Q20: {:?}", p.qual2);
        }
    }

    #[test]
    fn test_quality_profile_is_the_same_learned_on_one_thread_or_many() {
        // The context counts are learned in chunks on the thread pool; they
        // are sums, so the thread count must not change the profile.
        let rl = 60usize;
        let mut tr = StdRng::seed_from_u64(77);
        let pairs = good_and_poor_pairs(rl, 5000, &mut tr);
        let one = on_threads(1, || QualityProfile::from_read_pairs(&pairs, rl));
        let many = on_threads(4, || QualityProfile::from_read_pairs(&pairs, rl));
        assert!(one == many, "the profile learned on 4 threads differs from the one on 1");
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

    // ── Quality model tests ──────────────────────────────────────

    #[test]
    fn test_prev_q_bin() {
        use crate::quality::prev_q_bin;
        assert_eq!(prev_q_bin(b'!' + 0), 0); // Q0 → bin 0
        assert_eq!(prev_q_bin(b'!' + 9), 0); // Q9 → bin 0
        assert_eq!(prev_q_bin(b'!' + 10), 1); // Q10 → bin 1
        assert_eq!(prev_q_bin(b'!' + 19), 1); // Q19 → bin 1
        assert_eq!(prev_q_bin(b'!' + 20), 2); // Q20 → bin 2
        assert_eq!(prev_q_bin(b'!' + 29), 2); // Q29 → bin 2
        assert_eq!(prev_q_bin(b'!' + 30), 3); // Q30 → bin 3
        assert_eq!(prev_q_bin(b'!' + 40), 3); // Q40 → bin 3
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
        let result = gen.generate_haplotype_read_pair(&hap, 400, 450, "test", &mut rng).map(|(pair, _)| pair);
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
        let result = gen.generate_haplotype_read_pair(&hap, 200, 450, "test", &mut rng).map(|(pair, _)| pair);
        assert!(
            result.is_some(),
            "R1 in seg0, R2 in seg1 should produce discordant pair"
        );

        // Fragment crossing boundary → split-read evidence
        let result = gen.generate_haplotype_read_pair(&hap, 250, 450, "test2", &mut rng).map(|(pair, _)| pair);
        assert!(result.is_some(), "Boundary-crossing pair should succeed");

        let result = gen.generate_haplotype_read_pair(&hap, 210, 450, "test3", &mut rng).map(|(pair, _)| pair);
        assert!(result.is_some(), "Discordant pair spanning deletion");
    }

    #[test]
    fn test_hap_read_pair_safe_zones_accepted() {
        let hap = mock_del_haplotype(500, 5000);
        let (_ref, gen) = mock_synth_gen(150);
        let mut rng = StdRng::seed_from_u64(42);

        // Both arms entirely in left segment
        let result = gen.generate_haplotype_read_pair(&hap, 0, 400, "left", &mut rng).map(|(pair, _)| pair);
        assert!(result.is_some(), "Both arms in left segment should work");

        // Both arms entirely in right segment
        let result = gen.generate_haplotype_read_pair(&hap, 500, 400, "right", &mut rng).map(|(pair, _)| pair);
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
            .generate_haplotype_read_pair(&hap, 200, 450, "disc", &mut rng).map(|(pair, _)| pair)
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
        let result = gen.generate_haplotype_read_pair(&hap, 250, 300, "sv", &mut rng).map(|(pair, _)| pair);
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
            .generate_haplotype_read_pair(&hap, 2800, 400, "ins", &mut rng).map(|(pair, _)| pair)
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
                .generate_haplotype_read_pair(&hap, hap_start, frag_len, "h", &mut rng).map(|(pair, _)| pair)
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

    // ── Read length: cycles, short fragments and adapter trimming ───────

    /// `mock_gen_over` at Q93 (no sequencing errors), for a library whose
    /// reads were adapter-trimmed.
    fn trimming_gen_over(seq: Vec<u8>, cycles: usize) -> SynthReadGenerator<'static> {
        mock_gen_over(seq, cycles, 93, 93).with_adapter_trim(true)
    }

    #[test]
    fn test_adapter_suffix_is_the_longest_read_end_that_spells_the_adapter_start() {
        assert_eq!(adapter_suffix_len(b"CCCCT"), 0);
        assert_eq!(adapter_suffix_len(b"CCCCA"), 1);
        assert_eq!(adapter_suffix_len(b"CCCAG"), 2);
        assert_eq!(adapter_suffix_len(b"CCAGA"), 3, "AGA, not just its last A");
        assert_eq!(adapter_suffix_len(b"CAGAT"), 4);
        assert_eq!(adapter_suffix_len(b"CCAGT"), 0);
        assert_eq!(adapter_suffix_len(b"GA"), 1);
        assert_eq!(adapter_suffix_len(b"TTAGATCGGAAGAGC"), 13);
        assert_eq!(adapter_suffix_len(b"AGATCGGAAGAGCA"), 1, "only a prefix of the adapter counts");
        assert_eq!(adapter_suffix_len(b""), 0);
    }

    #[test]
    fn test_min_fragment_is_one_base_in_a_trimmed_library_and_one_read_otherwise() {
        let gen = mock_gen_over(vec![b'A'; 1_000], 151, 30, 30);
        assert_eq!(gen.min_fragment_len(), 151);
        assert_eq!(gen.with_adapter_trim(true).min_fragment_len(), 1);
    }

    #[test]
    fn test_a_fragment_shorter_than_the_cycles_gives_reads_as_long_as_the_fragment() {
        let mut seq = scrambled_seq(20_000, 21);
        // The 100 bp fragment ends ...CA: its reads were cut at the adapter
        // itself, so that A stays.
        seq[5_098] = b'C';
        seq[5_099] = b'A';
        let gen = trimming_gen_over(seq.clone(), 151);
        let mut rng = StdRng::seed_from_u64(5);
        for frag_len in [60u64, 100, 150] {
            let pair = gen
                .generate_read_pair("chr1", 5_000, frag_len, &HashMap::new(), "p", &mut rng)
                .expect("a short fragment is sequenced");
            let fragment = seq[5_000..5_000 + frag_len as usize].to_vec();
            let mut fragment_rc = fragment.clone();
            reverse_complement(&mut fragment_rc);
            let mut got = [pair.seq1.clone(), pair.seq2.clone()];
            got.sort();
            let mut want = [fragment, fragment_rc];
            want.sort();
            assert_eq!(got, want, "a {} bp fragment is read whole from each end", frag_len);
            assert_eq!((pair.qual1.len(), pair.qual2.len()), (frag_len as usize, frag_len as usize));
            assert_eq!(
                (pair.ref_start, pair.ref_end, pair.insert_size),
                (5_000, 5_000 + frag_len, frag_len as i64)
            );
        }
        let untrimmed = mock_gen_over(seq, 151, 93, 93);
        assert!(
            untrimmed
                .generate_read_pair("chr1", 5_000, 100, &HashMap::new(), "p", &mut rng)
                .is_none(),
            "an untrimmed library has no fragment shorter than a read"
        );
    }

    #[test]
    fn test_a_read_loses_the_adapter_start_its_3_prime_end_spells() {
        let (s, frag_len, cycles) = (5_000usize, 400usize, 151usize);
        let e = s + frag_len;
        let mut seq = scrambled_seq(20_000, 22);
        // The forward read [s, s + 151) ends ...CA: it loses its A.
        seq[s + 149] = b'C';
        seq[s + 150] = b'A';
        // The reverse read is [e - 151, e) reverse-complemented, so its 3'
        // end is seq[e - 148..e - 152] backwards, complemented: ...CAGA, and
        // it loses AGA.
        seq[e - 151] = b'T';
        seq[e - 150] = b'C';
        seq[e - 149] = b'T';
        seq[e - 148] = b'G';
        let fwd = seq[s..s + 150].to_vec();
        let mut rev = seq[e - 151..e].to_vec();
        reverse_complement(&mut rev);
        rev.truncate(148);

        let gen = trimming_gen_over(seq.clone(), cycles);
        let mut rng = StdRng::seed_from_u64(6);
        let (mut f1r2, mut f2r1) = (0usize, 0usize);
        for _ in 0..20 {
            let pair = gen
                .generate_read_pair("chr1", s as u64, frag_len as u64, &HashMap::new(), "p", &mut rng)
                .unwrap();
            if pair.seq1 == fwd {
                assert_eq!(pair.seq2, rev);
                f1r2 += 1;
            } else {
                assert_eq!((&pair.seq1, &pair.seq2), (&rev, &fwd));
                f2r1 += 1;
            }
            assert_eq!((pair.qual1.len(), pair.qual2.len()), (pair.seq1.len(), pair.seq2.len()));
        }
        assert!(f1r2 > 0 && f2r1 > 0, "orientation split {} / {}", f1r2, f2r1);

        // A fragment exactly one read long: nothing past it was read, so its
        // reads are trimmed the same way. [8000, 8151) forward ends ...CA;
        // reverse-complemented it ends ...TG.
        seq[8_149] = b'C';
        seq[8_150] = b'A';
        seq[8_000] = b'C';
        seq[8_001] = b'A';
        let gen = trimming_gen_over(seq.clone(), cycles);
        let pair = gen
            .generate_read_pair("chr1", 8_000, 151, &HashMap::new(), "p", &mut rng)
            .unwrap();
        let mut whole_rc = seq[8_000..8_151].to_vec();
        reverse_complement(&mut whole_rc);
        let mut got = [pair.seq1.clone(), pair.seq2.clone()];
        got.sort();
        let mut want = [seq[8_000..8_150].to_vec(), whole_rc];
        want.sort();
        assert_eq!(got, want);

        let untrimmed = mock_gen_over(seq, cycles, 93, 93);
        let pair = untrimmed
            .generate_read_pair("chr1", s as u64, frag_len as u64, &HashMap::new(), "p", &mut rng)
            .unwrap();
        assert_eq!((pair.seq1.len(), pair.seq2.len()), (151, 151), "no trimming in an untrimmed library");
    }

    #[test]
    fn test_haplotype_pairs_are_sequenced_and_trimmed_the_same_way() {
        let (s, frag_len) = (1_000usize, 400usize);
        let e = s + frag_len;
        let mut hap_seq = scrambled_seq(4_000, 23);
        hap_seq[s + 149] = b'C'; // the forward read ends ...CA
        hap_seq[s + 150] = b'A';
        // The reverse read, [e - 151, e) reverse-complemented, ends ...TGC:
        // no adapter start ends in GC but the whole 13-mer, whose third-last
        // base is A.
        hap_seq[e - 151] = b'G';
        hap_seq[e - 150] = b'C';
        hap_seq[e - 149] = b'A';
        let hap = VariantHaplotype::from_segments(vec![HaplotypeSegment {
            sequence: hap_seq.clone(),
            origin: Some(SegmentOrigin {
                chrom: "chr1".to_string(),
                ref_start: 0,
                ref_end: 4_000,
                is_reverse: false,
            }),
            hap_offset: 0,
        }]);
        let fwd = hap_seq[s..s + 150].to_vec();
        let mut rev = hap_seq[e - 151..e].to_vec();
        reverse_complement(&mut rev);

        let gen = mock_gen_over(vec![b'A'; 1_000], 151, 93, 93).with_adapter_trim(true);
        let mut rng = StdRng::seed_from_u64(8);
        for _ in 0..10 {
            let pair = gen
                .generate_haplotype_read_pair(&hap, s as u64, frag_len as u64, "h", &mut rng).map(|(pair, _)| pair)
                .unwrap();
            let mut got = [pair.seq1.clone(), pair.seq2.clone()];
            got.sort();
            let mut want = [fwd.clone(), rev.clone()];
            want.sort();
            assert_eq!(got, want);
            assert_eq!((pair.ref_start, pair.ref_end), (s as u64, e as u64));
        }

        // Indel errors on (at Q93 none happen) fetch bases past each read's
        // 3' end; a short fragment's reads must still stop at the fragment.
        let gen = mock_gen_over_with_indels(vec![b'A'; 1_000], 151, 93, 93, 1.0).with_adapter_trim(true);
        let pair = gen.generate_haplotype_read_pair(&hap, 2_000, 90, "h", &mut rng).map(|(pair, _)| pair).unwrap();
        let fragment = hap_seq[2_000..2_090].to_vec();
        assert_eq!((pair.seq1.len(), pair.seq2.len()), (90, 90));
        assert!(pair.seq1 == fragment || pair.seq2 == fragment, "a short fragment is read whole");
        assert_eq!((pair.ref_start, pair.ref_end, pair.insert_size), (2_000, 2_090, 90));

        let untrimmed = mock_gen_over(vec![b'A'; 1_000], 151, 93, 93);
        assert!(untrimmed.generate_haplotype_read_pair(&hap, 2_000, 90, "h", &mut rng).map(|(pair, _)| pair).is_none());
    }

    #[test]
    fn test_depth_copies_draw_fragments_shorter_than_a_read_only_in_a_trimmed_library() {
        let seq = scrambled_seq(20_000, 24);
        let short = FragmentDist::from_stats(100.0, 0.0);
        let original = mock_read_pair("o", vec![b'!' + 93; 151], vec![b'!' + 93; 151], 5_000);
        let mut rng = StdRng::seed_from_u64(9);
        let pair = trimming_gen_over(seq.clone(), 151)
            .generate_depth_pair("chr1", &original, &HashMap::new(), "d", &short, &mut rng)
            .unwrap();
        assert_eq!((pair.insert_size, pair.seq1.len(), pair.seq2.len()), (100, 100, 100));
        let pair = mock_gen_over(seq, 151, 93, 93)
            .generate_depth_pair("chr1", &original, &HashMap::new(), "d", &short, &mut rng)
            .unwrap();
        assert_eq!(pair.insert_size, 151, "an untrimmed library clamps the fragment up to one read");
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
                .generate_haplotype_read_pair(&hap, 500 + i * 7, 400, "h", &mut rng).map(|(pair, _)| pair)
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

        let (read_seq, _qual) = gen.generate_read("chr1", ref_start, &no_alleles, 1, true, rl, &mut rng);

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

        let (read_seq, read_qual) = gen.generate_read("chr1", 500, &no_alleles, 1, false, rl, &mut rng);

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

        let (read_seq, read_qual) = gen.generate_read("chr1", ref_start, &no_alleles, 1, true, rl, &mut rng);

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
            gen.generate_read("chr1", ref_start, &no_alleles, 1, false, rl, &mut rng);

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

        let (read_seq, _qual) = gen.generate_read_from_seq(&seq, rl, 1, true, &mut rng);

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

    /// A read whose every base is Q11 with probability `p_low`, else Q37,
    /// independently: a good read at small `p_low`, a poor one at large.
    fn scattered_low_qual(read_length: usize, p_low: f64, rng: &mut StdRng) -> Vec<u8> {
        (0..read_length)
            .map(|_| if rng.gen::<f64>() < p_low { b'!' + 11 } else { b'!' + 37 })
            .collect()
    }

    /// The standard deviation of `reads`' mean qualities (Phred).
    fn read_mean_sd(reads: &[&Vec<u8>]) -> f64 {
        let means: Vec<f64> = reads.iter().map(|q| mean_phred(q)).collect();
        let m = means.iter().sum::<f64>() / means.len() as f64;
        (means.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (means.len() - 1) as f64).sqrt()
    }

    fn mean_phred(qual: &[u8]) -> f64 {
        qual.iter().map(|&q| (q - b'!') as f64).sum::<f64>() / qual.len() as f64
    }

    /// Half the donor pairs good (2% of bases low) and half poor (40% low),
    /// both mates alike: the read-to-read spread and the mate link of the
    /// sample, with nothing a one-base chain can see along a read.
    fn good_and_poor_pairs(rl: usize, n: usize, rng: &mut StdRng) -> Vec<ReadPair> {
        (0..n)
            .map(|i| {
                let p = if i % 2 == 0 { 0.02 } else { 0.4 };
                let (q1, q2) = (scattered_low_qual(rl, p, rng), scattered_low_qual(rl, p, rng));
                mock_read_pair(&format!("gp_{}", i), q1, q2, i as u64 * 500)
            })
            .collect()
    }

    #[test]
    fn test_generated_reads_keep_the_pools_read_to_read_spread() {
        // Real reads are good or poor as a whole: on HG002 35x the SD of each
        // read's mean quality is 2.08, and a first-order chain gives 0.52.
        let rl = 100usize;
        let mut tr = StdRng::seed_from_u64(31);
        let pairs = good_and_poor_pairs(rl, 4000, &mut tr);
        let pool_sd = read_mean_sd(&pairs.iter().flat_map(|p| [&p.qual1, &p.qual2]).collect::<Vec<_>>());
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let gen = mock_gen_with_profile(profile, vec![b'A'; 200_000], rl, 0.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(32);
        let made: Vec<ReadPair> = (0..2000u64)
            .map(|i| gen.generate_read_pair("chr1", 1000 + i * 50, 300, &no_alleles, "p", &mut rng).unwrap())
            .collect();
        let made_sd = read_mean_sd(&made.iter().flat_map(|p| [&p.qual1, &p.qual2]).collect::<Vec<_>>());
        assert!(
            (made_sd / pool_sd - 1.0).abs() <= 0.15,
            "generated read-mean SD {:.2} should be within 15% of the pool's {:.2}",
            made_sd,
            pool_sd
        );
    }

    #[test]
    fn test_the_two_mates_of_a_generated_pair_share_the_pools_mate_link() {
        // A pair's two reads come off one cluster pair: on HG002 35x the
        // correlation of R1's and R2's mean quality is 0.55.
        let rl = 100usize;
        let mut tr = StdRng::seed_from_u64(41);
        let pairs = good_and_poor_pairs(rl, 4000, &mut tr);
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let gen = mock_gen_with_profile(profile, vec![b'A'; 200_000], rl, 0.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(42);
        let (mut r1, mut r2) = (Vec::new(), Vec::new());
        for i in 0..2000u64 {
            let p = gen.generate_read_pair("chr1", 1000 + i * 50, 300, &no_alleles, "p", &mut rng).unwrap();
            r1.push(mean_phred(&p.qual1));
            r2.push(mean_phred(&p.qual2));
        }
        let n = r1.len() as f64;
        let (m1, m2) = (r1.iter().sum::<f64>() / n, r2.iter().sum::<f64>() / n);
        let cov: f64 = r1.iter().zip(&r2).map(|(a, b)| (a - m1) * (b - m2)).sum();
        let v1: f64 = r1.iter().map(|a| (a - m1).powi(2)).sum();
        let v2: f64 = r2.iter().map(|b| (b - m2).powi(2)).sum();
        let corr = cov / (v1 * v2).sqrt();
        assert!(corr > 0.5, "the mates' mean qualities should correlate above 0.5, got {:.3}", corr);
    }

    /// Donor pairs over `reference` (chr1): R1 forward at a random place and
    /// R2 reverse 200 bp on, both 100 bp. Good pairs have 10% of bases at Q11
    /// and err at 2% there; poor pairs 60% and 30%. Q37 bases never err.
    /// Two pairs in five are poor, so the 40% read-class cut falls between the
    /// two kinds rather than inside one: at half and half, one class holds both.
    fn donors_that_err_by_read(reference: &[u8], n: usize, rng: &mut StdRng) -> Vec<ReadPair> {
        use crate::types::MateAlignment;
        let rl = 100usize;
        let read = |start: usize, reverse: bool, p_low: f64, p_err: f64, rng: &mut StdRng| {
            let mut seq = reference[start..start + rl].to_vec();
            if reverse {
                reverse_complement(&mut seq);
            }
            let mut qual = Vec::with_capacity(rl);
            for b in seq.iter_mut() {
                if rng.gen::<f64>() < p_low {
                    qual.push(b'!' + 11);
                    if rng.gen::<f64>() < p_err {
                        *b = random_different_base(*b, rng);
                    }
                } else {
                    qual.push(b'!' + 37);
                }
            }
            (seq, qual)
        };
        (0..n)
            .map(|i| {
                let (p_low, p_err) = if i % 5 >= 2 { (0.1, 0.02) } else { (0.6, 0.30) };
                let start = rng.gen_range(100..reference.len() - 500);
                let (seq1, qual1) = read(start, false, p_low, p_err, rng);
                let (seq2, qual2) = read(start + 200, true, p_low, p_err, rng);
                ReadPair {
                    name: format!("d_{}", i),
                    seq1,
                    qual1,
                    seq2,
                    qual2,
                    ref_start: start as u64,
                    ref_end: (start + 300) as u64,
                    insert_size: 300,
                    chrom: "chr1".to_string(),
                    align: Some(Box::new([
                        MateAlignment { start: start as u64, reverse: false, cigar: vec![(b'M', rl as u32)] },
                        MateAlignment { start: (start + 200) as u64, reverse: true, cigar: vec![(b'M', rl as u32)] },
                    ])),
                }
            })
            .collect()
    }

    #[test]
    fn test_generated_errors_follow_the_donors_rate_by_kind_of_read() {
        // At the same quality, a poor real read errs far more than a good
        // one: aligned HG002 bases at Q25 err 0.07% in good reads and 0.95%
        // in poor ones. 10^(-Q/10) gives both the same rate.
        let rl = 100usize;
        // 0.5x, so almost no reference position has the 5 donor reads the
        // variant rule needs. At 3x, positions where one poor read erred among
        // 5-10 reads were taken for the sample's own variants (12% of donor
        // bases) and the poor reads' rate was learned as 0.24. At a real
        // sample's ~30x, one error is ~3% of a position's reads, under 10%.
        let reference = scrambled_seq(2_400_000, 61);
        let mut tr = StdRng::seed_from_u64(62);
        let donors = donors_that_err_by_read(&reference, 6000, &mut tr);
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), reference.clone());
        let shared = SharedReference::from_sequences(seqs);
        let profile = QualityProfile::from_donor_pairs(&donors, rl, &shared);
        let gen = mock_gen_with_profile(profile, reference.clone(), rl, 0.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(63);
        // [good, poor]: (wrong Q11 bases, Q11 bases)
        let mut tally = [(0u64, 0u64); 2];
        for i in 0..4000u64 {
            let start = 1000 + i * 90;
            let p = gen.generate_read_pair("chr1", start, 300, &no_alleles, "p", &mut rng).unwrap();
            let mut fwd = reference[start as usize..start as usize + rl].to_vec();
            fwd.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
            let mut rev = reference[start as usize + 200..start as usize + 300].to_vec();
            reverse_complement(&mut rev);
            for (seq, qual) in [(&p.seq1, &p.qual1), (&p.seq2, &p.qual2)] {
                let mism = |t: &[u8]| seq.iter().zip(t).filter(|(a, b)| a != b).count();
                let template = if mism(&fwd) < mism(&rev) { &fwd } else { &rev };
                let kind = (mean_phred(qual) < 30.0) as usize;
                for ((&b, &t), &q) in seq.iter().zip(template.iter()).zip(qual.iter()) {
                    if q == b'!' + 11 {
                        tally[kind].1 += 1;
                        tally[kind].0 += (b != t) as u64;
                    }
                }
            }
        }
        let rate = |k: usize| tally[k].0 as f64 / tally[k].1 as f64;
        assert!((rate(0) - 0.02).abs() < 0.05, "good reads' Q11 bases should err near 2%, got {:.3}", rate(0));
        assert!((rate(1) - 0.30).abs() < 0.05, "poor reads' Q11 bases should err near 30%, got {:.3}", rate(1));
    }

    #[test]
    fn test_generated_qualities_follow_the_pools_pattern_over_more_than_one_base() {
        // Donor qualities run Q37 Q37 Q11 Q11 Q37 Q37 ..., each read from a
        // random phase: the next quality is fixed by the last two, while the
        // last one alone leaves it a coin flip.
        let rl = 60usize;
        let mut tr = StdRng::seed_from_u64(91);
        let pattern = |phase: usize| -> Vec<u8> {
            (0..rl).map(|i| if (i + phase) % 4 < 2 { b'!' + 37 } else { b'!' + 11 }).collect()
        };
        let pairs: Vec<ReadPair> = (0..3000)
            .map(|i| mock_read_pair(&format!("pt_{}", i), pattern(tr.gen_range(0..4)), pattern(tr.gen_range(0..4)), i * 300))
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let gen = mock_gen_with_profile(profile, vec![b'A'; 400_000], rl, 0.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(92);
        let (mut kept, mut all) = (0usize, 0usize);
        for i in 0..1000u64 {
            let p = gen.generate_read_pair("chr1", 1000 + i * 200, 250, &no_alleles, "p", &mut rng).unwrap();
            for q in [&p.qual1, &p.qual2] {
                for i in 2..q.len() {
                    all += 1;
                    // Two alike are followed by a change, two different by a repeat.
                    kept += ((q[i] == q[i - 1]) == (q[i - 1] != q[i - 2])) as usize;
                }
            }
        }
        let share = kept as f64 / all as f64;
        assert!(share > 0.95, "only {:.3} of generated qualities follow the two-base pattern", share);
    }

    #[test]
    fn test_a_pool_of_many_quality_values_generates_only_those_values() {
        // Older Illumina BAMs keep ~40 quality values, not 4 bins (the
        // novoalign HG002 chr20 BAM has 31): the history bins them.
        let rl = 80usize;
        let mut tr = StdRng::seed_from_u64(71);
        let values: Vec<u8> = (2..42).map(|q| b'!' + q).collect();
        let pairs: Vec<ReadPair> = (0..1500)
            .map(|i| {
                let q1: Vec<u8> = (0..rl).map(|_| values[tr.gen_range(0..values.len())]).collect();
                let q2: Vec<u8> = (0..rl).map(|_| values[tr.gen_range(0..values.len())]).collect();
                mock_read_pair(&format!("v_{}", i), q1, q2, i as u64 * 300)
            })
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let gen = mock_gen_with_profile(profile, vec![b'A'; 300_000], rl, 0.0);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(72);
        let mut seen = std::collections::BTreeSet::new();
        for i in 0..500u64 {
            let p = gen.generate_read_pair("chr1", 1000 + i * 100, 250, &no_alleles, "p", &mut rng).unwrap();
            for &q in p.qual1.iter().chain(&p.qual2) {
                assert!(values.contains(&q), "Q{} is not one of the pool's values", q - b'!');
                seen.insert(q);
            }
        }
        assert!(seen.len() >= 35, "only {} of the pool's 40 values came out", seen.len());
    }

    #[test]
    fn test_an_n_reports_q2_and_the_read_goes_on() {
        // An `N` is a no-call at Q2 (L18); the read goes on to its full
        // length, drawing from the pool's own values after it.
        let rl = 60usize;
        let pairs: Vec<ReadPair> = (0..200)
            .map(|i| mock_read_pair(&format!("n_{}", i), vec![b'!' + 37; rl], vec![b'!' + 37; rl], i * 300))
            .collect();
        let profile = QualityProfile::from_read_pairs(&pairs, rl);
        let mut seq = scrambled_seq(2000, 81);
        seq[520] = b'N';
        let gen = mock_gen_with_profile(profile, seq, rl, 0.0);
        let mut rng = StdRng::seed_from_u64(82);
        let (read_seq, read_qual) = gen.generate_read("chr1", 500, &HashMap::new(), 1, false, rl, &mut rng);
        assert_eq!(read_seq.len(), rl);
        assert_eq!(read_seq[20], b'N');
        assert_eq!(read_qual[20], b'!' + 2, "the N reports Q2");
        assert!(read_qual[21..].iter().all(|&q| q == b'!' + 37), "after the N the read draws the pool's Q37 again");
    }

    // --- N7: how far a small pool's quality model drifts (a measurement) ---

    /// (R1, R2) quality strings, Phred+33.
    type QualPairs = Vec<(Vec<u8>, Vec<u8>)>;

    /// N7's three distances between two sets of quality strings: the mean
    /// over cycles and reads of |per-cycle mean Q difference|, |difference in
    /// the fraction under Q20|, and |difference in P(Q < 20 at c+1 given Q < 20
    /// at c)|.
    fn n7_distances(a: &QualPairs, b: &QualPairs, read_length: usize) -> [f64; 3] {
        struct Summary {
            cycle_sum: [Vec<f64>; 2],
            cycle_n: [Vec<f64>; 2],
            low: f64,
            all: f64,
            low_then_low: f64,
            low_then_any: f64,
        }
        let summarise = |set: &QualPairs| {
            let mut s = Summary {
                cycle_sum: [vec![0.0; read_length], vec![0.0; read_length]],
                cycle_n: [vec![0.0; read_length], vec![0.0; read_length]],
                low: 0.0,
                all: 0.0,
                low_then_low: 0.0,
                low_then_any: 0.0,
            };
            for (q1, q2) in set {
                for (r, q) in [q1, q2].into_iter().enumerate() {
                    let len = q.len().min(read_length);
                    for c in 0..len {
                        let phred = q[c].saturating_sub(33);
                        s.cycle_sum[r][c] += phred as f64;
                        s.cycle_n[r][c] += 1.0;
                        s.all += 1.0;
                        if phred < 20 {
                            s.low += 1.0;
                            if c + 1 < len {
                                s.low_then_any += 1.0;
                                if q[c + 1].saturating_sub(33) < 20 {
                                    s.low_then_low += 1.0;
                                }
                            }
                        }
                    }
                }
            }
            s
        };
        let (sa, sb) = (summarise(a), summarise(b));
        let (mut m1, mut cells) = (0.0, 0.0);
        for r in 0..2 {
            for c in 0..read_length {
                if sa.cycle_n[r][c] > 0.0 && sb.cycle_n[r][c] > 0.0 {
                    m1 += (sa.cycle_sum[r][c] / sa.cycle_n[r][c]
                        - sb.cycle_sum[r][c] / sb.cycle_n[r][c])
                        .abs();
                    cells += 1.0;
                }
            }
        }
        let frac = |s: &Summary| s.low / s.all;
        let persist = |s: &Summary| {
            if s.low_then_any > 0.0 {
                s.low_then_low / s.low_then_any
            } else {
                0.0
            }
        };
        [
            m1 / cells,
            (frac(&sa) - frac(&sb)).abs(),
            (persist(&sa) - persist(&sb)).abs(),
        ]
    }

    #[test]
    fn test_quality_profile_warns_below_the_measured_pool_size() {
        // Below 50,000 sampled pairs the model is learned from thinner counts
        // than it was measured good with: the run goes on, but it has to say
        // so, with the sample size, the size it needs, and the census.
        let census = "12 full contexts with at least 20 observations";
        let warning = thin_profile_warning(32, census).expect("32 pairs is under the measured size");
        for part in ["32", "50000", census] {
            assert!(warning.contains(part), "the warning must name {:?}: {}", part, warning);
        }
        assert!(thin_profile_warning(49_999, census).is_some(), "49,999 pairs is still under it");
        assert_eq!(thin_profile_warning(50_000, census), None, "50,000 pairs is enough");
    }

    #[test]
    fn test_n7_distances_measure_what_they_say() {
        // The ruler for N7's measurement, checked on inputs with known
        // answers before it is trusted on real ones.
        let q = |phreds: &[u8]| phreds.iter().map(|p| p + 33).collect::<Vec<u8>>();
        let a: QualPairs = vec![(q(&[30, 30, 10, 10]), q(&[30, 30, 30, 30])); 4];
        assert_eq!(n7_distances(&a, &a, 4), [0.0, 0.0, 0.0], "a set against itself");

        // Every quality 5 higher: M1 is 5. Q10 becomes Q15, still under Q20.
        let b: QualPairs = vec![(q(&[35, 35, 15, 15]), q(&[35, 35, 35, 35])); 4];
        assert_eq!(n7_distances(&a, &b, 4), [5.0, 0.0, 0.0]);

        // R1's last base lifted to Q30: under-Q20 falls from 2/8 to 1/8, and
        // a low base is followed by a low one 1/1 -> 0/1 of the time.
        let c: QualPairs = vec![(q(&[30, 30, 10, 30]), q(&[30, 30, 30, 30])); 4];
        let d = n7_distances(&a, &c, 4);
        assert!((d[0] - 20.0 / 8.0).abs() < 1e-12, "M1 {}", d[0]);
        assert!((d[1] - 1.0 / 8.0).abs() < 1e-12, "M2 {}", d[1]);
        assert!((d[2] - 1.0).abs() < 1e-12, "M3 {}", d[2]);
    }

    /// Draw R1 and R2 qualities from `profile` over each pair's own bases, as
    /// `generate_from_template` does: both classes in one draw, the history
    /// running along the read, `N` at Q2 included.
    fn n7_draw(profile: &QualityProfile, pairs: &[ReadPair], rng: &mut StdRng) -> QualPairs {
        fn draw(profile: &QualityProfile, read_num: u8, class: usize, bases: &[u8], rng: &mut StdRng) -> Vec<u8> {
            let mut state = profile.start_read(read_num, Some(class), rng);
            let len = bases.len();
            let mut out = Vec::with_capacity(len);
            for (c, &b) in bases.iter().enumerate() {
                let q = if b.eq_ignore_ascii_case(&b'N') {
                    N_QUAL
                } else {
                    profile.next_quality(read_num, &state, len - c, rng)
                };
                out.push(q);
                profile.emitted(&mut state, q, b);
            }
            out
        }
        pairs
            .iter()
            .map(|p| {
                let runs = (crate::quality::template_run_bin(&p.seq1), crate::quality::template_run_bin(&p.seq2));
                let (c1, c2) = profile.draw_classes(runs, None, rng);
                (draw(profile, 1, c1, &p.seq1, rng), draw(profile, 2, c2, &p.seq2, rng))
            })
            .collect()
    }

    /// N7's measurement, as docs/review/REVIEW.md's N7 plan locks it. Run by hand:
    /// `SPIKE_N7_BAM=<HG002 35x BAM> cargo test --release -- --ignored
    /// measure_n7_quality_drift --nocapture`. Prints one row per tolerance
    /// window and per repeat, the medians, and N*.
    #[test]
    #[ignore]
    fn measure_n7_quality_drift() {
        use crate::extract::extract_read_pairs;
        const READ_LENGTH: usize = 151;
        let bam = std::env::var("SPIKE_N7_BAM").expect("set SPIKE_N7_BAM to the HG002 35x BAM");

        let pool = extract_read_pairs(&bam, "chr20", 38_402_500, 38_432_500, 20, None)
            .unwrap()
            .pairs;
        let fnv = |name: &str| {
            name.bytes()
                .fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100_0000_01b3))
        };
        let (train, held): (Vec<ReadPair>, Vec<ReadPair>) =
            pool.into_iter().partition(|p| fnv(&p.name) % 2 == 0);
        let real: QualPairs = held.iter().map(|p| (p.qual1.clone(), p.qual2.clone())).collect();
        println!("pool\ttrain={}\theld_out={}", train.len(), held.len());

        let mut tol = [0.0f64; 3];
        for mb in [32u64, 33, 34, 35, 36, 37, 40, 41, 42, 43] {
            let start = mb * 1_000_000;
            let other: QualPairs = extract_read_pairs(&bam, "chr20", start, start + 30_000, 20, None)
                .unwrap()
                .pairs
                .iter()
                .map(|p| (p.qual1.clone(), p.qual2.clone()))
                .collect();
            let d = n7_distances(&other, &real, READ_LENGTH);
            println!("tolerance\t{}Mb\tpairs={}\t{:.4}\t{:.5}\t{:.5}", mb, other.len(), d[0], d[1], d[2]);
            for k in 0..3 {
                tol[k] = tol[k].max(d[k]);
            }
        }
        println!("T\t{:.4}\t{:.5}\t{:.5}", tol[0], tol[1], tol[2]);

        let sizes = [30, 60, 125, 250, 500, 1000, train.len()];
        let mut ok = Vec::new();
        for (si, &n) in sizes.iter().enumerate() {
            let mut per_metric: [Vec<f64>; 3] = Default::default();
            for rep in 0..20u64 {
                let mut rng = StdRng::seed_from_u64(1000 * si as u64 + rep);
                let subset: Vec<ReadPair> = rand::seq::index::sample(&mut rng, train.len(), n)
                    .iter()
                    .map(|i| train[i].clone())
                    .collect();
                let profile = QualityProfile::from_read_pairs(&subset, READ_LENGTH);
                let d = n7_distances(&n7_draw(&profile, &held, &mut rng), &real, READ_LENGTH);
                println!("size\t{}\t{}\t{:.4}\t{:.5}\t{:.5}", n, rep, d[0], d[1], d[2]);
                for k in 0..3 {
                    per_metric[k].push(d[k]);
                }
            }
            let median: Vec<f64> = per_metric
                .iter_mut()
                .map(|v| {
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    (v[9] + v[10]) / 2.0
                })
                .collect();
            println!("median\t{}\t{:.4}\t{:.5}\t{:.5}", n, median[0], median[1], median[2]);
            ok.push((0..3).all(|k| median[k] <= tol[k]));
        }
        // The smallest tested n inside the tolerance there and at every larger n.
        let n_star = (0..sizes.len())
            .find(|&i| ok[i..].iter().all(|&b| b))
            .map(|i| sizes[i]);
        println!("N*\t{:?}", n_star);
    }

    // --- K4: how many donor pairs the quality model needs (a measurement) ---

    /// K1's three numbers for `fake` strings against `real` ones: the ratio of
    /// their read-mean SDs, and two-proportion z for the share of perfect
    /// reads (every base at `real`'s top quality) and of crashed reads (at
    /// least 10 of the last 20 qualities below Q15).
    fn k1_metrics(fake: &QualPairs, real: &QualPairs) -> [f64; 3] {
        let reads = |set: &QualPairs| -> Vec<Vec<u8>> {
            set.iter().flat_map(|(a, b)| [a.clone(), b.clone()]).filter(|q| !q.is_empty()).collect()
        };
        let (f, r) = (reads(fake), reads(real));
        let top = *r.iter().flatten().max().unwrap();
        let sd = |v: &[Vec<u8>]| read_mean_sd(&v.iter().collect::<Vec<_>>());
        let share = |v: &[Vec<u8>], pred: &dyn Fn(&[u8]) -> bool| {
            (v.iter().filter(|q| pred(q)).count() as f64, v.len() as f64)
        };
        let z = |(k1, n1): (f64, f64), (k2, n2): (f64, f64)| {
            let p = (k1 + k2) / (n1 + n2);
            let se = (p * (1.0 - p) * (1.0 / n1 + 1.0 / n2)).sqrt();
            if se == 0.0 { 0.0 } else { (k1 / n1 - k2 / n2) / se }
        };
        let perfect = |q: &[u8]| q.iter().all(|&b| b == top);
        let crashed = |q: &[u8]| q.len() >= 20 && q[q.len() - 20..].iter().filter(|&&b| b < b'!' + 15).count() >= 10;
        [
            sd(&f) / sd(&r),
            z(share(&f, &perfect), share(&r, &perfect)),
            z(share(&f, &crashed), share(&r, &crashed)),
        ]
    }

    #[test]
    fn test_k1_metrics_measure_what_they_say() {
        let q = |v: &[u8]| v.iter().map(|p| p + 33).collect::<Vec<u8>>();
        let mut a: QualPairs = Vec::new();
        for i in 0..100 {
            let good = q(&[37; 30]);
            let poor = if i % 2 == 0 { q(&[37; 30]) } else { q(&[11; 30]) };
            a.push((good, poor));
        }
        let m = k1_metrics(&a, &a);
        assert!((m[0] - 1.0).abs() < 1e-12 && m[1] == 0.0 && m[2] == 0.0, "a set against itself: {:?}", m);
        // All average: the SD collapses, and there are no perfect or crashed reads.
        let b: QualPairs = (0..100).map(|_| (q(&[30; 30]), q(&[30; 30]))).collect();
        let m = k1_metrics(&b, &a);
        assert!(m[0] < 0.01 && m[1] < -3.0 && m[2] < -3.0, "an all-average set against a mixed one: {:?}", m);
    }

    /// K4 of docs/superpowers/plans/2026-10-07-quality-model.md. Run by hand:
    /// `SPIKE_N7_BAM=<HG002 35x BAM> cargo test --release -- --ignored
    /// measure_quality_pool_size --nocapture`. A 200 kb window (N7's 30 kb
    /// holds ~2,400 pairs, too few for a 5,000-pair pool), split by name into
    /// a train and a held-out half; for each pool size, 20 random pools from
    /// the train half generate the held-out pairs' strings, and K1's rule is
    /// applied to the medians: SD ratio 0.85-1.15, |z| < 3 for perfect and
    /// crashed reads.
    #[test]
    #[ignore]
    fn measure_quality_pool_size() {
        use crate::extract::extract_read_pairs;
        const READ_LENGTH: usize = 151;
        let bam = std::env::var("SPIKE_N7_BAM").expect("set SPIKE_N7_BAM to the HG002 35x BAM");
        let pool = extract_read_pairs(&bam, "chr20", 38_402_500, 38_602_500, 20, None).unwrap().pairs;
        let fnv = |name: &str| {
            name.bytes()
                .fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100_0000_01b3))
        };
        let (train, held): (Vec<ReadPair>, Vec<ReadPair>) = pool.into_iter().partition(|p| fnv(&p.name) % 2 == 0);
        let real: QualPairs = held.iter().map(|p| (p.qual1.clone(), p.qual2.clone())).collect();
        println!("pool\ttrain={}\theld_out={}", train.len(), held.len());
        let sizes = [500usize, 1000, 2000, 5000];
        let mut ok = Vec::new();
        for (si, &n) in sizes.iter().enumerate() {
            let n = n.min(train.len());
            let mut per_metric: [Vec<f64>; 3] = Default::default();
            for rep in 0..20u64 {
                let mut rng = StdRng::seed_from_u64(1000 * si as u64 + rep);
                let subset: Vec<ReadPair> = rand::seq::index::sample(&mut rng, train.len(), n)
                    .iter()
                    .map(|i| train[i].clone())
                    .collect();
                let profile = QualityProfile::from_read_pairs(&subset, READ_LENGTH);
                let m = k1_metrics(&n7_draw(&profile, &held, &mut rng), &real);
                println!("size\t{}\t{}\t{:.3}\t{:.2}\t{:.2}", n, rep, m[0], m[1], m[2]);
                for k in 0..3 {
                    per_metric[k].push(m[k]);
                }
            }
            let median: Vec<f64> = per_metric
                .iter_mut()
                .map(|v| {
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    (v[9] + v[10]) / 2.0
                })
                .collect();
            let pass = (0.85..=1.15).contains(&median[0]) && median[1].abs() < 3.0 && median[2].abs() < 3.0;
            println!("median\t{}\t{:.3}\t{:.2}\t{:.2}\t{}", n, median[0], median[1], median[2], if pass { "pass" } else { "fail" });
            ok.push(pass);
        }
        let n_star = (0..sizes.len()).find(|&i| ok[i..].iter().all(|&b| b)).map(|i| sizes[i]);
        println!("smallest passing size, and every larger one: {:?}", n_star);
    }

    // --- v2: runs, lows, error history, the event's class mix (plan 2026-10-07-quality-model-v2) ---

    /// A pair's two templates in sequencing order on `reference` (chr1): R1
    /// forward at `start`, R2 reverse at `start + 200`, both `rl` long.
    fn mate_templates(reference: &[u8], start: usize, rl: usize) -> (Vec<u8>, Vec<u8>) {
        let fwd = reference[start..start + rl].to_ascii_uppercase();
        let mut rev = reference[start + 200..start + 200 + rl].to_ascii_uppercase();
        reverse_complement(&mut rev);
        (fwd, rev)
    }

    /// The last cycle of the first run of 12+ of one base in `template`, if any.
    fn run12_end(template: &[u8]) -> Option<usize> {
        let mut len = 0;
        for i in 0..template.len() {
            len = if i > 0 && template[i] == template[i - 1] { len + 1 } else { 1 };
            if len >= 12 && (i + 1 == template.len() || template[i + 1] != template[i]) {
                return Some(i);
            }
        }
        None
    }

    /// Makes a donor mate's (called bases, qualities) from its template.
    type MakeMate<'a> = &'a dyn Fn(&[u8], &mut StdRng) -> (Vec<u8>, Vec<u8>);

    /// Donor pairs on `reference` (chr1) at random starts in `range`, R1
    /// forward and R2 reverse 200 bp on, `rl` long; `make(template, rng)`
    /// gives each mate's (called bases, qualities) in sequencing order.
    fn donors_on(
        reference: &[u8],
        range: std::ops::Range<usize>,
        n: usize,
        rl: usize,
        rng: &mut StdRng,
        make: MakeMate,
    ) -> Vec<ReadPair> {
        use crate::types::MateAlignment;
        (0..n)
            .map(|i| {
                let start = rng.gen_range(range.clone());
                let (t1, t2) = mate_templates(reference, start, rl);
                let (seq1, qual1) = make(&t1, rng);
                let (seq2, qual2) = make(&t2, rng);
                ReadPair {
                    name: format!("d_{}", i),
                    seq1,
                    qual1,
                    seq2,
                    qual2,
                    ref_start: start as u64,
                    ref_end: (start + 200 + rl) as u64,
                    insert_size: (200 + rl) as i64,
                    chrom: "chr1".to_string(),
                    align: Some(Box::new([
                        MateAlignment { start: start as u64, reverse: false, cigar: vec![(b'M', rl as u32)] },
                        MateAlignment { start: (start + 200) as u64, reverse: true, cigar: vec![(b'M', rl as u32)] },
                    ])),
                }
            })
            .collect()
    }

    /// Generated pairs at random starts in `range`, each mate with its template.
    fn made_on(
        gen: &SynthReadGenerator,
        reference: &[u8],
        range: std::ops::Range<usize>,
        n: usize,
        rl: usize,
        rng: &mut StdRng,
    ) -> Vec<(Vec<u8>, Vec<u8>, Vec<u8>)> {
        let no_alleles = HashMap::new();
        let mut out = Vec::new();
        for _ in 0..n {
            let start = rng.gen_range(range.clone());
            let p = gen
                .generate_read_pair("chr1", start as u64, (200 + rl) as u64, &no_alleles, "p", rng)
                .unwrap();
            let (fwd, rev) = mate_templates(reference, start, rl);
            for (seq, qual) in [(p.seq1, p.qual1), (p.seq2, p.qual2)] {
                let mism = |t: &[u8]| seq.iter().zip(t).filter(|(a, b)| a != b).count();
                let t = if mism(&fwd) <= mism(&rev) { fwd.clone() } else { rev.clone() };
                out.push((seq, qual, t));
            }
        }
        out
    }

    fn crashed(qual: &[u8]) -> bool {
        qual.len() >= 20 && qual[qual.len() - 20..].iter().filter(|&&q| q < b'!' + 15).count() >= 10
    }

    fn generator_on(reference: &[u8], donors: &[ReadPair], rl: usize) -> SynthReadGenerator<'static> {
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), reference.to_vec());
        let shared = SharedReference::from_sequences(seqs);
        let profile = QualityProfile::from_donor_pairs(donors, rl, &shared);
        mock_gen_with_profile(profile, reference.to_vec(), rl, 0.0)
    }

    /// `seq` with each base of `seq` wrong with chance `p_err`.
    fn with_errors(template: &[u8], qual: &[u8], p_err: &dyn Fn(usize, u8) -> f64, rng: &mut StdRng) -> Vec<u8> {
        template
            .iter()
            .zip(qual)
            .enumerate()
            .map(|(c, (&b, &q))| if rng.gen::<f64>() < p_err(c, q) { random_different_base(b, rng) } else { b })
            .collect()
    }

    #[test]
    fn test_reads_crash_after_a_run_of_one_base_as_the_donors_do() {
        // T1. A run of 12+ of one base sets off the rest of a read: on HG002
        // 35x, reads with one crash 15.0% of the time, against 0.46% without.
        let rl = 100usize;
        let mut reference = scrambled_seq(1_000_000, 81);
        for k in 0..450usize {
            let at = 1_000 + k * 2_000;
            reference[at - 1] = b'G';
            reference[at..at + 12].copy_from_slice(&[b'T'; 12]);
            reference[at + 12] = b'G';
        }
        let mut tr = StdRng::seed_from_u64(82);
        let donors = donors_on(&reference, 100..900_000, 8000, rl, &mut tr, &|t, rng| {
            let end = run12_end(t);
            let qual: Vec<u8> = (0..rl)
                .map(|c| {
                    let after = end.is_some_and(|e| c > e);
                    let p_low = if after { 0.7 } else { 0.02 };
                    if rng.gen::<f64>() < p_low { b'!' + 11 } else { b'!' + 37 }
                })
                .collect();
            (t.to_vec(), qual)
        });
        let early = |t: &[u8]| run12_end(t).is_some_and(|e| e <= 75);
        let donor_mates: Vec<(&Vec<u8>, &Vec<u8>, Vec<u8>)> = donors
            .iter()
            .flat_map(|p| {
                let (f, r) = mate_templates(&reference, p.ref_start as usize, rl);
                [(&p.seq1, &p.qual1, f), (&p.seq2, &p.qual2, r)]
            })
            .collect();
        let share = |v: &[&Vec<u8>]| v.iter().filter(|q| crashed(q)).count() as f64 / v.len().max(1) as f64;
        let d_run = share(&donor_mates.iter().filter(|m| early(&m.2)).map(|m| m.1).collect::<Vec<_>>());

        let gen = generator_on(&reference, &donors, rl);
        let mut rng = StdRng::seed_from_u64(83);
        let made = made_on(&gen, &reference, 100..900_000, 6000, rl, &mut rng);
        let g_run = share(&made.iter().filter(|m| early(&m.2)).map(|m| &m.1).collect::<Vec<_>>());
        let g_free = share(&made.iter().filter(|m| crate::quality::template_run_bin(&m.2) == 0).map(|m| &m.1).collect::<Vec<_>>());
        assert!(g_run >= 5.0 * g_free.max(0.001), "reads after a run crash {:.3}, run-free reads {:.3}: not 5x", g_run, g_free);
        assert!(g_run >= 0.7 * d_run, "reads after a run crash {:.3}, the donors' {:.3}", g_run, d_run);
        // And the crash starts after the run, not anywhere in the read.
        let (mut before, mut after) = ((0usize, 0usize), (0usize, 0usize));
        for (_, q, t) in &made {
            if let Some(e) = run12_end(t).filter(|e| (30..=60).contains(e)) {
                let low = |w: &[u8]| w.iter().filter(|&&x| x < b'!' + 15).count();
                before.0 += low(&q[e - 20..e]);
                before.1 += 20;
                after.0 += low(&q[e + 1..=e + 20]);
                after.1 += 20;
            }
        }
        let (b, a) = (before.0 as f64 / before.1 as f64, after.0 as f64 / after.1 as f64);
        assert!(a >= 3.0 * b, "low qualities after the run {:.3}, before it {:.3}", a, b);
    }

    #[test]
    fn test_runs_are_learned_from_the_template_not_from_a_dead_tails_calls() {
        // T2. A dead tail calls false runs of its own: on HG002 35x, 22% of
        // crashed reads on a run-free template show one. Learned from the
        // called bases, the model would tie poor reads to runs.
        let rl = 100usize;
        let mut reference = scrambled_seq(1_000_000, 84);
        for k in 0..40usize {
            let at = 900_000 + k * 2_000;
            reference[at - 1] = b'G';
            reference[at..at + 12].copy_from_slice(&[b'T'; 12]);
            reference[at + 12] = b'G';
        }
        let mut tr = StdRng::seed_from_u64(85);
        let donors = donors_on(&reference, 100..800_000, 8000, rl, &mut tr, &|t, rng| {
            if rng.gen::<f64>() < 0.2 {
                // A dead tail: Q11 at 80% over the last 30 cycles, all called A.
                let qual: Vec<u8> = (0..rl)
                    .map(|c| if c >= rl - 30 && rng.gen::<f64>() < 0.8 { b'!' + 11 } else { b'!' + 37 })
                    .collect();
                let mut seq = t.to_vec();
                seq[rl - 30..].fill(b'A');
                (seq, qual)
            } else {
                let qual: Vec<u8> = (0..rl).map(|_| if rng.gen::<f64>() < 0.02 { b'!' + 11 } else { b'!' + 37 }).collect();
                (t.to_vec(), qual)
            }
        });
        let pool_crash = donors.iter().flat_map(|p| [&p.qual1, &p.qual2]).filter(|q| crashed(q)).count() as f64
            / (2 * donors.len()) as f64;
        let gen = generator_on(&reference, &donors, rl);
        let mut rng = StdRng::seed_from_u64(86);
        let share = |v: &[(Vec<u8>, Vec<u8>, Vec<u8>)]| v.iter().filter(|m| crashed(&m.1)).count() as f64 / v.len() as f64;
        let free = made_on(&gen, &reference, 100..800_000, 4000, rl, &mut rng);
        let g_free = share(&free);
        assert!(
            (g_free / pool_crash - 1.0).abs() <= 0.3,
            "run-free templates crash {:.3}, the pool {:.3}",
            g_free,
            pool_crash
        );
        let with_run: Vec<_> = (0..40usize)
            .flat_map(|k| made_on(&gen, &reference, 900_000 + k * 2_000 - 40..900_000 + k * 2_000 - 39, 50, rl, &mut rng))
            .filter(|m| crate::quality::template_run_bin(&m.2) == 3)
            .collect();
        assert!(with_run.len() > 1000, "only {} reads over the runs", with_run.len());
        let g_run = share(&with_run);
        assert!(g_run <= 1.3 * pool_crash, "templates with a run crash {:.3}, the pool {:.3}", g_run, pool_crash);
    }

    #[test]
    fn test_errors_follow_a_dense_stretch_of_low_qualities() {
        // T3. Low qualities in a bad stretch are mixed with high ones (on
        // HG002 35x the run of lows ending at a read's last base has a median
        // of 1), and errors follow how dense they are, not how many in a row.
        // Every read has the same mean quality (one class) and one dense
        // stretch -- Q11 at every other cycle over 30 cycles, erring at 40% --
        // somewhere in cycles 20-70, with 10 more Q11s scattered elsewhere,
        // erring at 5%.
        let rl = 100usize;
        let reference = scrambled_seq(1_000_000, 87);
        let mut tr = StdRng::seed_from_u64(88);
        let donors = donors_on(&reference, 100..900_000, 8000, rl, &mut tr, &|t, rng| {
            let from = rng.gen_range(20..=40usize);
            let dense = |c: usize| c >= from && c < from + 30 && (c - from) & 1 == 0;
            let mut qual = vec![b'!' + 37; rl];
            (0..rl).filter(|&c| dense(c)).for_each(|c| qual[c] = b'!' + 11);
            let mut scattered = 0;
            while scattered < 10 {
                let c = rng.gen_range(0..rl);
                if !(from.saturating_sub(16)..from + 46).contains(&c) && qual[c] == b'!' + 37 {
                    qual[c] = b'!' + 11;
                    scattered += 1;
                }
            }
            let seq = with_errors(t, &qual, &|c, q| if q != b'!' + 11 { 0.0 } else if dense(c) { 0.40 } else { 0.05 }, rng);
            (seq, qual)
        });
        // Q11 bases with 6+ lows among the 16 cycles before them: inside a
        // dense stretch.
        let gen = generator_on(&reference, &donors, rl);
        let mut rng = StdRng::seed_from_u64(89);
        let made = made_on(&gen, &reference, 100..900_000, 4000, rl, &mut rng);
        let (mut wrong, mut all) = (0usize, 0usize);
        for (seq, q, t) in &made {
            for c in 16..rl {
                if q[c] == b'!' + 11 && q[c - 16..c].iter().filter(|&&x| x < b'!' + 15).count() >= 6 {
                    wrong += (seq[c] != t[c]) as usize;
                    all += 1;
                }
            }
        }
        assert!(all > 2000, "only {} Q11 bases in dense stretches", all);
        let rate = wrong as f64 / all as f64;
        assert!((rate - 0.40).abs() < 0.05, "Q11 bases in a dense stretch err {:.3}, not near 0.40", rate);
    }

    #[test]
    fn test_rare_errors_follow_the_lows_over_16_cycles_not_the_lows_in_a_row() {
        // T3b (added for mutant 4). Errors rare enough that the error history
        // stays near empty, but set by how many of the last 16 qualities were
        // low: Q11 bases with 6+ lows among the 16 cycles ending at them err
        // at 3%, others at 0.3%. Lows come at random inside a stretch, so the
        // run of lows in a row is short either way.
        let rl = 100usize;
        let reference = scrambled_seq(1_000_000, 95);
        let lows16 = |q: &[u8], c: usize| q[c.saturating_sub(15)..=c].iter().filter(|&&x| x < b'!' + 15).count();
        let mut tr = StdRng::seed_from_u64(96);
        let donors = donors_on(&reference, 100..900_000, 12_000, rl, &mut tr, &|t, rng| {
            let from = rng.gen_range(0..=60usize);
            let qual: Vec<u8> = (0..rl)
                .map(|c| if c >= from && c < from + 40 && rng.gen::<bool>() { b'!' + 11 } else { b'!' + 37 })
                .collect();
            let seq = with_errors(
                t,
                &qual,
                &|c, q| if q != b'!' + 11 { 0.0 } else if lows16(&qual, c) >= 6 { 0.03 } else { 0.003 },
                rng,
            );
            (seq, qual)
        });
        let gen = generator_on(&reference, &donors, rl);
        let mut rng = StdRng::seed_from_u64(97);
        let made = made_on(&gen, &reference, 100..900_000, 8000, rl, &mut rng);
        let (mut dense, mut sparse) = ((0usize, 0usize), (0usize, 0usize));
        for (seq, q, t) in &made {
            for c in 0..rl {
                if q[c] == b'!' + 11 {
                    let cell = if lows16(q, c) >= 6 { &mut dense } else { &mut sparse };
                    cell.0 += (seq[c] != t[c]) as usize;
                    cell.1 += 1;
                }
            }
        }
        let rate = |x: (usize, usize)| x.0 as f64 / x.1 as f64;
        assert!(dense.1 > 30_000 && sparse.1 > 10_000, "too few Q11 bases: {} dense, {} sparse", dense.1, sparse.1);
        assert!((rate(dense) / 0.03 - 1.0).abs() <= 0.3, "Q11 bases after 6+ lows err {:.4}, not near 0.03", rate(dense));
        assert!((rate(sparse) - 0.003).abs() <= 0.002, "Q11 bases after fewer lows err {:.4}, not near 0.003", rate(sparse));
    }

    #[test]
    fn test_reads_keep_the_donors_read_to_read_spread_of_errors() {
        // T4. Tails differ by read beyond their qualities: on HG002 35x the
        // per-read error rate of crashed tails has SD 0.160, against 0.108 if
        // every read erred alike. Here half the crashed reads err at 60% at
        // Q11 and half at 10%, with the same qualities.
        let rl = 100usize;
        let reference = scrambled_seq(1_000_000, 90);
        let mut tr = StdRng::seed_from_u64(91);
        let donors = donors_on(&reference, 100..900_000, 8000, rl, &mut tr, &|t, rng| {
            if rng.gen::<f64>() < 0.3 {
                let qual: Vec<u8> = (0..rl).map(|c| if c >= rl - 40 && c % 3 > 0 { b'!' + 11 } else { b'!' + 37 }).collect();
                let p = if rng.gen::<bool>() { 0.6 } else { 0.1 };
                let seq = with_errors(t, &qual, &|_, q| if q == b'!' + 11 { p } else { 0.0 }, rng);
                (seq, qual)
            } else {
                (t.to_vec(), vec![b'!' + 37; rl])
            }
        });
        let per_read = |seq: &[u8], q: &[u8], t: &[u8]| -> Option<f64> {
            let idx: Vec<usize> = (0..rl).filter(|&c| q[c] == b'!' + 11 && c >= rl - 40).collect();
            (idx.len() >= 10).then(|| idx.iter().filter(|&&c| seq[c] != t[c]).count() as f64 / idx.len() as f64)
        };
        let sd = |v: &[f64]| {
            let m = v.iter().sum::<f64>() / v.len() as f64;
            (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt()
        };
        let pool: Vec<f64> = donors
            .iter()
            .flat_map(|p| {
                let (f, r) = mate_templates(&reference, p.ref_start as usize, rl);
                [per_read(&p.seq1, &p.qual1, &f), per_read(&p.seq2, &p.qual2, &r)]
            })
            .flatten()
            .collect();
        let gen = generator_on(&reference, &donors, rl);
        let mut rng = StdRng::seed_from_u64(92);
        let made: Vec<f64> = made_on(&gen, &reference, 100..900_000, 4000, rl, &mut rng)
            .iter()
            .filter_map(|(s, q, t)| per_read(s, q, t))
            .collect();
        assert!(made.len() > 500, "only {} generated crashed reads", made.len());
        assert!(
            sd(&made) >= 0.8 * sd(&pool),
            "per-read error rate SD {:.3} against the pool's {:.3}",
            sd(&made),
            sd(&pool)
        );
    }

    #[test]
    fn test_an_events_reads_follow_its_own_class_mix() {
        // T5. Regions differ: on HG002 35x the crashed share runs 0.60-2.62%
        // over 20 blocks, and a block's own class mix halves the error of
        // predicting it. A sample of good and poor reads, and an event pool of
        // only poor ones.
        let rl = 100usize;
        let mut tr = StdRng::seed_from_u64(93);
        let sample = good_and_poor_pairs(rl, 4000, &mut tr);
        let pool: Vec<ReadPair> = good_and_poor_pairs(rl, 1000, &mut tr).into_iter().skip(1).step_by(2).collect();
        let profile = QualityProfile::from_read_pairs(&sample, rl);
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), vec![b'A'; 200_000]);
        let reference: &'static SharedReference = Box::leak(Box::new(SharedReference::from_sequences(seqs)));
        let mix = profile.class_mix(&pool, reference);
        let pool_class = pool.iter().flat_map(|p| [&p.qual1, &p.qual2]).map(|q| profile.class_of(q) as f64).sum::<f64>()
            / (2 * pool.len()) as f64;
        let gen = SynthReadGenerator::new(profile, reference, rl, 0.0).with_class_mix(mix);
        let no_alleles = HashMap::new();
        let mut rng = StdRng::seed_from_u64(94);
        let mut classes = Vec::new();
        for i in 0..2000u64 {
            let p = gen.generate_read_pair("chr1", 1000 + i * 50, 300, &no_alleles, "p", &mut rng).unwrap();
            classes.push(gen.profile.class_of(&p.qual1) as f64);
            classes.push(gen.profile.class_of(&p.qual2) as f64);
        }
        let made_class = classes.iter().sum::<f64>() / classes.len() as f64;
        assert!(
            (made_class - pool_class).abs() <= 1.0,
            "generated reads' mean class {:.2}, the event pool's {:.2}",
            made_class,
            pool_class
        );
    }

    /// K2b of docs/superpowers/plans/2026-10-07-quality-model-v2.md: the big
    /// held-out clip test, on this code. Run by hand (docs/analysis/quality-model-v2/k2b.sh):
    /// `SPIKE_K2B_BAM=<BAM> SPIKE_K2B_REF=<FASTA> SPIKE_K2B_OUT=<dir> SPIKE_K2B_SET=real|spike
    /// [SPIKE_K2B_SEED=21] cargo test --release -- --ignored measure_clip_share --nocapture`.
    /// Samples the BAM as a run does (`sample_input`), learns from the pairs
    /// whose name hashes even, and writes the odd ones with |TLEN| >= 151 and
    /// both mates full length to `<set>_R1.fq` / `_R2.fq`: as sequenced
    /// (`real`), or as spike makes them from the reference at each mate's own
    /// 5' position and strand, each block's pairs with that block's class mix
    /// (`spike`).
    #[test]
    #[ignore]
    fn measure_clip_share() {
        use std::io::Write;
        let env = |k: &str| std::env::var(k).unwrap_or_else(|_| panic!("set {}", k));
        let (bam, fasta, out, set) = (env("SPIKE_K2B_BAM"), env("SPIKE_K2B_REF"), env("SPIKE_K2B_OUT"), env("SPIKE_K2B_SET"));
        let seed: u64 = std::env::var("SPIKE_K2B_SEED").ok().and_then(|v| v.parse().ok()).unwrap_or(21);
        let rl = crate::bam_stats::compute_stats(&bam, 10_000, Some(&fasta)).unwrap().cycles;
        let fnv = |name: &str| name.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100_0000_01b3));
        let blocks = crate::quality::sample_input(&bam, &fasta, 20).unwrap();
        let (mut train, mut held) = (Vec::new(), Vec::new());
        for b in &blocks {
            let (tr, he): (Vec<ReadPair>, Vec<ReadPair>) = b.pairs.iter().cloned().partition(|p| fnv(&p.name) % 2 == 0);
            train.push(crate::quality::DonorBlock { pairs: tr, ..b.clone() });
            held.push(he);
        }
        let profile = QualityProfile::from_blocks(&train, rl);
        let mut f1 = std::io::BufWriter::new(std::fs::File::create(format!("{}/{}_R1.fq", out, set)).unwrap());
        let mut f2 = std::io::BufWriter::new(std::fs::File::create(format!("{}/{}_R2.fq", out, set)).unwrap());
        let put = |f: &mut std::io::BufWriter<std::fs::File>, name: &str, s: &[u8], q: &[u8]| {
            writeln!(f, "@{}\n{}\n+\n{}", name, String::from_utf8_lossy(s), String::from_utf8_lossy(q)).unwrap();
        };
        let none = SharedReference::from_sequences(StdHashMap::new());
        let gen = SynthReadGenerator::new(profile, &none, rl, 0.0);
        let mut rng = StdRng::seed_from_u64(seed);
        let mut written = 0usize;
        for (b, (tr, he)) in blocks.iter().zip(train.iter().zip(&held)) {
            let mix = gen.profile.class_mix(&tr.pairs, b);
            for p in he {
                let Some(al) = &p.align else { continue };
                if p.insert_size.abs() < rl as i64 || p.seq1.len() != rl || p.seq2.len() != rl {
                    continue;
                }
                written += 1;
                if set == "real" {
                    put(&mut f1, &p.name, &p.seq1, &p.qual1);
                    put(&mut f2, &p.name, &p.seq2, &p.qual2);
                    continue;
                }
                let t1 = crate::quality::template_of(b, &p.chrom, &al[0], rl);
                let t2 = crate::quality::template_of(b, &p.chrom, &al[1], rl);
                let runs = (crate::quality::template_run_bin(&t1), crate::quality::template_run_bin(&t2));
                let (c1, c2) = gen.profile.draw_classes(runs, Some(&mix), &mut rng);
                let (s1, q1) = gen.generate_from_template(&t1, rl, 1, Some(c1), &mut rng);
                let (s2, q2) = gen.generate_from_template(&t2, rl, 2, Some(c2), &mut rng);
                put(&mut f1, &p.name, &s1, &q1);
                put(&mut f2, &p.name, &s2, &q2);
            }
        }
        println!("set {} seed {}: {} held-out pairs written; learned from {} pairs in {} blocks", set, seed, written,
                 train.iter().map(|b| b.pairs.len()).sum::<usize>(), blocks.len());
    }
}
