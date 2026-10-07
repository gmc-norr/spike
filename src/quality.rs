//! Base qualities and sequencing errors, learned once from a sample of the
//! sample's own reads.
//!
//! **Where.** At startup spike reads 20 blocks of 50 kb spread over the parts
//! of the input its index says hold reads (`sample_input`), and learns
//! everything below from them once. One event's
//! donor pool (~2,800 pairs) was too thin for crashed reads and their errors
//! (docs/superpowers/plans/2026-10-07-quality-model.md, refuted).
//!
//! **Qualities.** fqzcomp's quality context (htscodecs `fqzcomp_qual`, as CRAM
//! uses it) turned into a generator. Each base's quality is drawn given:
//! - the last qualities of the read;
//! - how many cycles are left to its 3' end;
//! - whether its quality has already changed bin twice;
//! - the read's *class*: where its mean quality falls among the sample's;
//! - how many of its last 16 qualities were low;
//! - the longest run of one base the read has read so far, in its template.
//!
//! A run of one base sets off the rest of a read: after 12+ T's, the strand
//! reading into the run has 6.9-9.4x the low qualities and 14-22x the errors
//! for the rest of the read, and the other strand does not (HG002 35x,
//! docs/superpowers/plans/2026-10-07-quality-model-v2.md). The run is read
//! from the template -- the reference under a donor read, the haplotype under
//! spike's own -- because a dead tail calls false runs of its own.
//!
//! A pair's two classes are drawn together, from the sample's joint table,
//! weighted by each mate's run and by the event's own class mix
//! (`ClassMix`), so mates share state and a region keeps its share of poor
//! reads.
//!
//! **Errors.** Each base is wrong at the rate the sample's reads show for its
//! quality, the read's class, its recent low qualities, its distance from the
//! 3' end, its run, and the errors the read has already made, counted against
//! the reference (`block_errors`).

use std::collections::{HashMap, HashSet};
use std::hash::{BuildHasherDefault, Hasher};

use anyhow::{bail, Result};
use rand::rngs::StdRng;
use rand::Rng;
use rayon::prelude::*;

use crate::reference::SharedReference;
use crate::types::{MateAlignment, ReadPair};

/// Read classes (the "selector"), cut at these quantiles of the sample's read
/// mean qualities (both mates pooled). The low end is cut finely, because the
/// poor reads are few and differ most.
pub const READ_CLASSES: usize = 8;
const CLASS_QUANTILES: [f64; READ_CLASSES - 1] = [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8];

/// A context needs this many observations before its counts are drawn from;
/// below it, the draw backs off to a coarser context.
const MIN_CONTEXT_OBS: u32 = 20;

/// An error-table cell needs this many counted bases before its rate is used;
/// below it, the rate backs off to a coarser cell, then to 10^(-Q/10).
const MIN_ERROR_BASES: u64 = 200;

/// A quality below this (Phred) is "low".
const LOW_Q: u8 = 15;

/// Quality bins for the history of a large alphabet and the change flag:
/// Q0-9, Q10-19, Q20-29, Q30+.
pub const PREV_Q_BINS: usize = 4;

/// Bins for the low qualities among the last 16: 0, 1-2, 3-5, 6-9, 10+. Low
/// qualities in a bad tail are mixed with high ones to the read's end (the run
/// of lows ending at the last base has a median of 1), so they are counted
/// over a window, not in a row.
const LOWS_BINS: usize = 5;
/// Bins for the cycles left to a read's 3' end, the base itself counted:
/// 1-10, 11-30, 31-60, 61+.
const END_BINS: usize = 4;
/// Bins for the longest run of one base read so far (`run_bin`).
pub const RUN_BINS: usize = 4;
/// Bins for the errors among the read's last `ERR_WINDOW` bases: 0, 1, 2-3, 4+.
const ERR_BINS: usize = 4;
const ERR_WINDOW: u32 = 30;

/// Slots of quality history kept; an alphabet of more than 4 values keeps the
/// last quality exactly and the two before it in `PREV_Q_BINS` bins.
const HISTORY: usize = 5;
/// An empty history slot (the read's first cycles).
const EMPTY: u8 = u8::MAX;

/// Phred+33 byte drawn when the profile learned nothing at all.
const FALLBACK_QUAL: u8 = b'!' + 20;

/// The startup sample: this many blocks of this many bases.
pub const SAMPLE_BLOCKS: usize = 20;
pub const SAMPLE_BLOCK_LEN: u64 = 50_000;
/// Reference bases fetched past a block's reads on each side.
const BLOCK_REF_PAD: u64 = 1_000;
/// Reads of prior weight a run bin's class shares are shrunk by, toward the
/// sample's overall shares.
const RUN_PRIOR_READS: f64 = 8.0;
/// Fewest sampled pairs spike learns from, as for a donor pool.
const MIN_PAIRS_TO_LEARN: usize = 30;

/// A quality's bin in [`PREV_Q_BINS`]: Q0-9 → 0, Q10-19 → 1, Q20-29 → 2, Q30+ → 3.
pub fn prev_q_bin(phred_plus_33: u8) -> usize {
    match phred_plus_33.saturating_sub(33) {
        0..=9 => 0,
        10..=19 => 1,
        20..=29 => 2,
        _ => 3,
    }
}

fn lows_bin(recent_low: u16) -> usize {
    match recent_low.count_ones() {
        0 => 0,
        1..=2 => 1,
        3..=5 => 2,
        6..=9 => 3,
        _ => 4,
    }
}

/// A run of `len` bases `base`, binned by what it does to the cycles after it
/// (HG002 35x, low-quality share against the same cycles in general): A, G or
/// T 7-8 → 1 (about 1.2x), 9-11 → 2 (2.1-2.6x), 12+ → 3 (6.9-9.4x); C 5-6 → 2
/// (1.8x), C 7+ → 3 (5.9x and more); anything shorter, or `N`, → 0.
pub fn run_bin(base: u8, len: u32) -> u8 {
    match (base, len) {
        (b'C', 5..=6) => 2,
        (b'C', 7..) => 3,
        (b'A' | b'G' | b'T', 7..=8) => 1,
        (b'A' | b'G' | b'T', 9..=11) => 2,
        (b'A' | b'G' | b'T', 12..) => 3,
        _ => 0,
    }
}

/// The largest `run_bin` over `template` (sequencing order, either case).
pub fn template_run_bin(template: &[u8]) -> usize {
    let (mut best, mut base, mut len) = (0u8, b'N', 0u32);
    for &b in template {
        let b = b.to_ascii_uppercase();
        if b == base && b != b'N' {
            len += 1;
        } else {
            base = b;
            len = (b != b'N') as u32;
        }
        best = best.max(run_bin(base, len));
    }
    best as usize
}

fn err_bin(errors: u32) -> usize {
    match (errors & ((1 << ERR_WINDOW) - 1)).count_ones() {
        0 => 0,
        1 => 1,
        2..=3 => 2,
        _ => 3,
    }
}

fn end_bin(cycles_left: usize) -> usize {
    match cycles_left {
        0..=10 => 0,
        11..=30 => 1,
        31..=60 => 2,
        _ => 3,
    }
}

/// A hasher for the context keys: splitmix64's finaliser over the one `u64`
/// a key is. The keys are spike's own, so no protection against chosen keys
/// is needed, and the default hasher was most of the cost of learning.
#[derive(Default)]
struct KeyHasher(u64);

impl Hasher for KeyHasher {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.write_u64(self.0 ^ b as u64);
        }
    }

    fn write_u64(&mut self, x: u64) {
        let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        self.0 = z ^ (z >> 31);
    }
}

type ContextMap = HashMap<u64, Counts, BuildHasherDefault<KeyHasher>>;

/// One read being generated (or walked, when learning): its class and what
/// the context needs from the read so far.
#[derive(Debug, Clone)]
pub struct ReadState {
    class: usize,
    /// Symbols of the last qualities, most recent first; `EMPTY` before any.
    last: [u8; HISTORY],
    /// Times the quality changed bin (`prev_q_bin`) so far.
    changes: u32,
    /// Which of the last 16 qualities were below `LOW_Q` (bit 0 the last).
    recent_low: u16,
    /// The template's current run of one base, and the largest `run_bin` so far.
    run_base: u8,
    run_len: u32,
    run: u8,
    /// Which of the last bases were errors (bit 0 the last).
    errors: u32,
}

impl ReadState {
    fn new(class: usize) -> Self {
        Self {
            class,
            last: [EMPTY; HISTORY],
            changes: 0,
            recent_low: 0,
            run_base: b'N',
            run_len: 0,
            run: 0,
            errors: 0,
        }
    }
}

/// Counts per quality symbol, and their total.
#[derive(Debug, Clone, Default, PartialEq)]
struct Counts {
    by_symbol: Vec<u32>,
    total: u32,
}

/// Errors and counted bases in one error-table cell.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct ErrorCell {
    errors: u64,
    bases: u64,
}

impl ErrorCell {
    fn add(&mut self, other: ErrorCell) {
        self.errors += other.errors;
        self.bases += other.bases;
    }
}

/// The error table at its four levels, finest first.
#[derive(Debug, Clone, Default, PartialEq)]
struct ErrorTables {
    /// (symbol, class, lows, end, run, errors so far).
    full: Vec<ErrorCell>,
    /// (symbol, lows, run, errors so far).
    mid: Vec<ErrorCell>,
    /// (symbol, class).
    class: Vec<ErrorCell>,
    /// (symbol).
    symbol: Vec<ErrorCell>,
}

impl ErrorTables {
    fn new(alphabet: usize) -> Self {
        Self {
            full: vec![ErrorCell::default(); alphabet * READ_CLASSES * LOWS_BINS * END_BINS * RUN_BINS * ERR_BINS],
            mid: vec![ErrorCell::default(); alphabet * LOWS_BINS * RUN_BINS * ERR_BINS],
            class: vec![ErrorCell::default(); alphabet * READ_CLASSES],
            symbol: vec![ErrorCell::default(); alphabet],
        }
    }

    fn merge(mut self, other: ErrorTables) -> Self {
        for (a, b) in [
            (&mut self.full, &other.full),
            (&mut self.mid, &other.mid),
            (&mut self.class, &other.class),
            (&mut self.symbol, &other.symbol),
        ] {
            a.iter_mut().zip(b).for_each(|(x, y)| x.add(*y));
        }
        self
    }
}

/// How the sampled bases fared when counted for the error table.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ErrorCensus {
    /// Bases counted: aligned or in a clipped end judged a bad end.
    pub counted: u64,
    /// Bases at a reference position where the sample itself differs.
    pub variant_site: u64,
    /// Bases in a clipped end that is not the reference read badly.
    pub foreign_clip: u64,
    /// Inserted bases, `N`s, and bases off the fetched reference.
    pub other: u64,
    /// Pairs without an alignment (none in real runs).
    pub pairs_without_alignment: u64,
}

impl ErrorCensus {
    fn merge(mut self, o: ErrorCensus) -> Self {
        self.counted += o.counted;
        self.variant_site += o.variant_site;
        self.foreign_clip += o.foreign_clip;
        self.other += o.other;
        self.pairs_without_alignment += o.pairs_without_alignment;
        self
    }
}

/// Reference bases by window: the shared reference, or a sample block's own.
pub trait RefSource: Sync {
    /// Uppercase bases over 0-based `[start, end)` of `chrom`, clipped to what
    /// the source holds; `None` when it holds none of them.
    fn bases(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>>;
}

impl RefSource for SharedReference {
    fn bases(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        self.fetch_sequence(chrom, start, end)
            .ok()
            .filter(|s| !s.is_empty())
            .map(|s| s.to_ascii_uppercase())
    }
}

/// One block of the startup sample: its pairs, and the reference under them.
#[derive(Debug, Clone)]
pub struct DonorBlock {
    pub chrom: String,
    /// The reference bases `ref_seq` start here (0-based).
    pub ref_start: u64,
    pub ref_seq: Vec<u8>,
    pub pairs: Vec<ReadPair>,
}

impl DonorBlock {
    /// The uppercase reference base at 0-based `pos`, `None` off the block's
    /// reference or at an `N`.
    fn base(&self, pos: i64) -> Option<u8> {
        let i = pos - self.ref_start as i64;
        (i >= 0)
            .then(|| self.ref_seq.get(i as usize).map(u8::to_ascii_uppercase))
            .flatten()
            .filter(|&b| b != b'N')
    }
}

impl RefSource for DonorBlock {
    fn bases(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>> {
        if chrom != self.chrom || start < self.ref_start {
            return None;
        }
        let a = (start - self.ref_start) as usize;
        let b = ((end - self.ref_start) as usize).min(self.ref_seq.len());
        (a < b).then(|| self.ref_seq[a..b].to_ascii_uppercase())
    }
}

/// `n` reference bases from mate `m`'s 5' end in sequencing order: from its
/// alignment start less a leading soft clip, or -- reverse-complemented --
/// back from its alignment end plus a trailing one. Uppercase, padded with `N`
/// where the reference runs out. This is what the sequencer read, the
/// sample's own alleles aside, and what its run is read from.
pub fn template_of(reference: &dyn RefSource, chrom: &str, m: &MateAlignment, n: usize) -> Vec<u8> {
    let (start, end) = placed_span(m);
    let mut t = if !m.reverse {
        reference.bases(chrom, start, start + n as u64).unwrap_or_default()
    } else {
        let mut v = reference.bases(chrom, end.saturating_sub(n as u64), end).unwrap_or_default();
        // A reverse mate's 5' end is the window's right end: pad on the left
        // (its 3' end) when the window came back short.
        let short = n.saturating_sub(v.len());
        if short > 0 {
            let mut padded = vec![b'N'; short];
            padded.extend_from_slice(&v);
            v = padded;
        }
        crate::extract::reverse_complement(&mut v);
        v
    };
    t.resize(n, b'N');
    t
}

/// An event's own class mix: per mate, how much more (or less) often each
/// class occurs among the event's donor reads than the sample predicts from
/// their runs. 1 throughout means none.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassMix {
    r: [[f64; READ_CLASSES]; 2],
}

impl Default for ClassMix {
    fn default() -> Self {
        Self { r: [[1.0; READ_CLASSES]; 2] }
    }
}

/// The learned quality and error model.
#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub struct QualityProfile {
    /// The sample's quality bytes (Phred+33), sorted.
    alphabet: Vec<u8>,
    /// A byte's symbol: its index in `alphabet`, or the nearest one's.
    symbol_of: Vec<u8>,
    /// Whether the history keeps 5 exact qualities (an alphabet of ≤ 4).
    small: bool,
    /// Cycles left are binned as `min(7, cycles_left >> pshift)`.
    pshift: u32,
    /// Read-class cuts, on mean Phred quality.
    class_cuts: Vec<f64>,
    /// Joint (R1 class, R2 class) counts over sampled pairs, flattened.
    joint: Vec<u64>,
    /// Per mate, (run bin of the template, class) counts over sampled reads.
    class_run: Vec<u64>,
    /// Per mate, counts per context key (all backoff levels in one map).
    tables: [ContextMap; 2],
    /// The error table; empty when no errors were learned.
    errors: ErrorTables,
    census: Option<ErrorCensus>,
    n_pairs: usize,
    n_blocks: usize,
    /// `run_ratio(mate, run bin, class)` for all `2 * RUN_BINS *
    /// READ_CLASSES` cells, computed once in `learn` where `class_run`
    /// becomes final, and read by `draw_classes` instead of being recomputed
    /// for each of the 64 weights of each pair. `None` until it is filled;
    /// `run_ratio_at` panics on `None` rather than falling back to a neutral
    /// 1.0, which would silently change every draw.
    run_ratios: Option<Vec<f64>>,
}

impl QualityProfile {
    /// Learn the qualities from `pairs`, with no reference: runs are read from
    /// the called bases and every base errs at 10^(-Q/10). For tests.
    #[cfg(test)]
    pub fn from_read_pairs(pairs: &[ReadPair], read_length: usize) -> Self {
        let block = DonorBlock { chrom: String::new(), ref_start: 0, ref_seq: Vec::new(), pairs: pairs.to_vec() };
        Self::learn(std::slice::from_ref(&block), read_length, false)
    }

    /// Learn the qualities and the error table from donor `pairs`, against
    /// `reference`: one block per chromosome. For tests.
    #[cfg(test)]
    pub fn from_donor_pairs(pairs: &[ReadPair], read_length: usize, reference: &SharedReference) -> Self {
        let mut by_chrom: Vec<(String, Vec<ReadPair>)> = Vec::new();
        for p in pairs {
            match by_chrom.iter_mut().find(|(c, _)| *c == p.chrom) {
                Some((_, v)) => v.push(p.clone()),
                None => by_chrom.push((p.chrom.clone(), vec![p.clone()])),
            }
        }
        let blocks: Vec<DonorBlock> = by_chrom
            .into_iter()
            .map(|(chrom, pairs)| {
                let (lo, hi) = reads_span(&pairs);
                let start = lo.saturating_sub(BLOCK_REF_PAD);
                let ref_seq = reference.bases(&chrom, start, hi + BLOCK_REF_PAD).unwrap_or_default();
                DonorBlock { chrom, ref_start: start, ref_seq, pairs }
            })
            .collect();
        Self::learn(&blocks, read_length, true)
    }

    /// Learn the qualities and the error table from the startup sample.
    pub fn from_blocks(blocks: &[DonorBlock], read_length: usize) -> Self {
        Self::learn(blocks, read_length, true)
    }

    /// Sample the input (`sample_input`) and learn from it.
    pub fn from_input(alignment_path: &str, ref_path: &str, min_mapq: u8, read_length: usize) -> Result<Self> {
        let t0 = std::time::Instant::now();
        let blocks = sample_input(alignment_path, ref_path, min_mapq)?;
        let t_sample = t0.elapsed().as_secs_f64();
        let pairs: usize = blocks.iter().map(|b| b.pairs.len()).sum();
        if pairs < MIN_PAIRS_TO_LEARN {
            bail!(
                "the quality sample of {} holds {} read pair(s), fewer than the {} spike needs to \
                 learn the sample's base qualities and errors; check that the file is indexed and \
                 that --min-mapq is not filtering everything out",
                alignment_path,
                pairs,
                MIN_PAIRS_TO_LEARN
            );
        }
        let t1 = std::time::Instant::now();
        let profile = Self::from_blocks(&blocks, read_length);
        let first = blocks.first().map(|b| format!("{}:{}", b.chrom, b.ref_start)).unwrap_or_default();
        let last = blocks.last().map(|b| format!("{}:{}", b.chrom, b.ref_start)).unwrap_or_default();
        log::info!(
            "Quality sample: {} blocks of {} kb from {} to {}, {} pairs; read in {:.1} s, learned in {:.1} s",
            blocks.len(),
            SAMPLE_BLOCK_LEN / 1000,
            first,
            last,
            pairs,
            t_sample,
            t1.elapsed().as_secs_f64()
        );
        Ok(profile)
    }

    /// Learn everything, then fill the run-ratio table. The only way to a
    /// `QualityProfile`: `learn_counts` is private and called from here
    /// alone, so no profile can reach `draw_classes` with an unfilled table,
    /// the empty-alphabet profile included.
    fn learn(blocks: &[DonorBlock], read_length: usize, with_reference: bool) -> Self {
        let mut profile = Self::learn_counts(blocks, read_length, with_reference);
        profile.fill_run_ratios();
        profile
    }

    /// `run_ratio(m, h, c)` for every cell, from the final `class_run`
    /// counts. Called once, from `learn`.
    fn fill_run_ratios(&mut self) {
        let table: Vec<f64> = (0..2)
            .flat_map(|m| (0..RUN_BINS).flat_map(move |h| (0..READ_CLASSES).map(move |c| (m, h, c))))
            .map(|(m, h, c)| self.run_ratio(m, h, c))
            .collect();
        self.run_ratios = Some(table);
    }

    fn learn_counts(blocks: &[DonorBlock], read_length: usize, with_reference: bool) -> Self {
        let pairs = || blocks.iter().flat_map(|b| b.pairs.iter());
        let mut seen = [false; 256];
        for p in pairs() {
            for &q in p.qual1.iter().chain(&p.qual2) {
                seen[q as usize] = true;
            }
        }
        let alphabet: Vec<u8> = (0..=255u8).filter(|&b| seen[b as usize]).collect();
        let symbol_of: Vec<u8> = (0..=255u8)
            .map(|b| {
                alphabet
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, &a)| (a as i32 - b as i32).abs())
                    .map(|(i, _)| i as u8)
                    .unwrap_or(0)
            })
            .collect();
        let pshift = ((read_length.max(8) as f64 / 8.0).log2().round().max(0.0)) as u32;

        let mut profile = Self {
            small: alphabet.len() <= 4,
            alphabet,
            symbol_of,
            pshift,
            class_cuts: Vec::new(),
            joint: vec![0; READ_CLASSES * READ_CLASSES],
            class_run: vec![0; 2 * RUN_BINS * READ_CLASSES],
            tables: [ContextMap::default(), ContextMap::default()],
            errors: ErrorTables::default(),
            census: None,
            n_pairs: pairs().count(),
            n_blocks: blocks.len(),
            run_ratios: None,
        };
        if profile.alphabet.is_empty() {
            log::info!("Quality profile: no sampled qualities; every base is Q20");
            return profile;
        }

        // Read classes: cut the sample's read means at fixed quantiles.
        let mut means: Vec<f64> = pairs()
            .flat_map(|p| [&p.qual1, &p.qual2])
            .filter(|q| !q.is_empty())
            .map(|q| mean_phred(q))
            .collect();
        means.sort_by(|a, b| a.partial_cmp(b).unwrap());
        profile.class_cuts = CLASS_QUANTILES.iter().map(|&f| quantile(&means, f)).collect();

        // Each mate's template: the reference under it, or -- with no
        // reference -- its called bases.
        let templates: Vec<Vec<[Vec<u8>; 2]>> = blocks
            .par_iter()
            .map(|b| {
                b.pairs
                    .iter()
                    .map(|p| match (with_reference, &p.align) {
                        (true, Some(al)) => [
                            template_of(b, &p.chrom, &al[0], p.qual1.len()),
                            template_of(b, &p.chrom, &al[1], p.qual2.len()),
                        ],
                        _ => [p.seq1.to_ascii_uppercase(), p.seq2.to_ascii_uppercase()],
                    })
                    .collect()
            })
            .collect();

        for (b, ts) in blocks.iter().zip(&templates) {
            for (p, t) in b.pairs.iter().zip(ts) {
                let (c1, c2) = (profile.class_of(&p.qual1), profile.class_of(&p.qual2));
                profile.joint[c1 * READ_CLASSES + c2] += 1;
                for (m, c) in [c1, c2].into_iter().enumerate() {
                    profile.class_run[(m * RUN_BINS + template_run_bin(&t[m])) * READ_CLASSES + c] += 1;
                }
            }
        }

        // The context counts, per mate, in chunks on the thread pool; counts
        // are sums, so neither the chunking nor the join order changes them.
        let flat: Vec<(&ReadPair, &[Vec<u8>; 2])> = blocks
            .iter()
            .zip(&templates)
            .flat_map(|(b, ts)| b.pairs.iter().zip(ts.iter()))
            .collect();
        for mate in 0..2 {
            let table = flat
                .par_chunks(2048)
                .map(|chunk| {
                    let mut t = ContextMap::default();
                    for (p, tm) in chunk {
                        let qual = if mate == 0 { &p.qual1 } else { &p.qual2 };
                        profile.count_read(qual, &tm[mate], &mut t);
                    }
                    t
                })
                .reduce(ContextMap::default, |mut a, b| {
                    for (k, v) in b {
                        let e = a.entry(k).or_insert_with(|| Counts {
                            by_symbol: vec![0; v.by_symbol.len()],
                            total: 0,
                        });
                        for (x, y) in e.by_symbol.iter_mut().zip(&v.by_symbol) {
                            *x += y;
                        }
                        e.total += v.total;
                    }
                    a
                });
            profile.tables[mate] = table;
        }

        if with_reference {
            let a = profile.alphabet.len();
            let (errors, census) = blocks
                .par_iter()
                .zip(templates.par_iter())
                .map(|(b, ts)| profile.block_errors(b, ts))
                .reduce(
                    || (ErrorTables::new(a), ErrorCensus::default()),
                    |(e1, c1), (e2, c2)| (e1.merge(e2), c1.merge(c2)),
                );
            profile.errors = errors;
            profile.census = Some(census);
        }
        profile.log_summary();
        profile
    }

    /// The class of a read with quality string `qual`.
    pub(crate) fn class_of(&self, qual: &[u8]) -> usize {
        if qual.is_empty() {
            return READ_CLASSES / 2;
        }
        let m = mean_phred(qual);
        self.class_cuts.iter().filter(|&&c| c <= m).count()
    }

    fn phred(&self, symbol: u8) -> u8 {
        self.alphabet[symbol as usize].saturating_sub(33)
    }

    /// Add one quality to `state`, as the read emitted it, with the template
    /// base `base` it was read from.
    fn push(&self, state: &mut ReadState, qual: u8, base: u8) {
        let symbol = self.symbol_of[qual as usize];
        if state.last[0] != EMPTY
            && prev_q_bin(self.alphabet[state.last[0] as usize]) != prev_q_bin(self.alphabet[symbol as usize])
        {
            state.changes += 1;
        }
        state.last.rotate_right(1);
        state.last[0] = symbol;
        state.recent_low = (state.recent_low << 1) | (self.phred(symbol) < LOW_Q) as u16;
        let base = base.to_ascii_uppercase();
        if base == state.run_base && base != b'N' {
            state.run_len += 1;
        } else {
            state.run_base = base;
            state.run_len = (base != b'N') as u32;
        }
        state.run = state.run.max(run_bin(state.run_base, state.run_len));
    }

    /// The context keys for the next base, finest first: the full context
    /// (with the lows and the run), then (2-quality history, lows, run,
    /// position, class), (last quality, run, position, class), (last quality,
    /// run, position), (position), and the whole mate.
    fn keys(&self, state: &ReadState, cycles_left: usize) -> [u64; 6] {
        let a = self.alphabet.len() as u64;
        let slot = |k: usize| -> u64 {
            let s = state.last[k];
            if self.small || k == 0 {
                if s == EMPTY { if self.small { 4 } else { a } } else { s as u64 }
            } else if s == EMPTY {
                PREV_Q_BINS as u64
            } else {
                prev_q_bin(self.alphabet[s as usize]) as u64
            }
        };
        let base0 = if self.small { 5 } else { a + 1 };
        let deep = if self.small { HISTORY } else { 3 };
        let mut full = 0u64;
        let mut mult = 1u64;
        for k in 0..deep {
            full += slot(k) * mult;
            mult *= if k == 0 { base0 } else if self.small { 5 } else { PREV_Q_BINS as u64 + 1 };
        }
        let two = slot(0) + slot(1) * base0;
        let one = slot(0);
        let pos = (cycles_left >> self.pshift).min(7) as u64;
        let flag = (state.changes >= 2) as u64;
        let class = state.class as u64;
        let lows = lows_bin(state.recent_low) as u64;
        let run = state.run as u64;
        let tag = |level: u64| level << 58;
        [
            tag(0) | run << 36 | lows << 32 | full << 12 | pos << 8 | flag << 4 | class,
            tag(1) | run << 36 | lows << 32 | two << 12 | pos << 8 | class,
            tag(2) | run << 36 | one << 12 | pos << 8 | class,
            tag(3) | run << 36 | one << 12 | pos << 8,
            tag(4) | pos << 8,
            tag(5),
        ]
    }

    /// Count one sampled read's qualities into `table`, its template beside it.
    fn count_read(&self, qual: &[u8], template: &[u8], table: &mut ContextMap) {
        let mut state = ReadState::new(self.class_of(qual));
        let len = qual.len();
        for (c, &q) in qual.iter().enumerate() {
            let symbol = self.symbol_of[q as usize] as usize;
            for key in self.keys(&state, len - c) {
                let e = table
                    .entry(key)
                    .or_insert_with(|| Counts { by_symbol: vec![0; self.alphabet.len()], total: 0 });
                e.by_symbol[symbol] += 1;
                e.total += 1;
            }
            self.push(&mut state, q, template.get(c).copied().unwrap_or(b'N'));
        }
    }

    /// Draw a pair's two read classes: the sample's joint table, weighted per
    /// mate by how its template's run bin (`runs`) shifts the classes, and by
    /// the event's own class mix.
    pub fn draw_classes(&self, runs: (usize, usize), mix: Option<&ClassMix>, rng: &mut StdRng) -> (usize, usize) {
        let neutral = ClassMix::default();
        let mix = mix.unwrap_or(&neutral);
        let weights: Vec<f64> = (0..READ_CLASSES * READ_CLASSES)
            .map(|i| {
                let (c1, c2) = (i / READ_CLASSES, i % READ_CLASSES);
                self.joint[i] as f64
                    * self.run_ratio_at(0, runs.0, c1)
                    * self.run_ratio_at(1, runs.1, c2)
                    * mix.r[0][c1]
                    * mix.r[1][c2]
            })
            .collect();
        let total: f64 = weights.iter().sum();
        if total <= 0.0 {
            return (READ_CLASSES / 2, READ_CLASSES / 2);
        }
        let mut r = rng.gen::<f64>() * total;
        for (i, &w) in weights.iter().enumerate() {
            if r < w {
                return (i / READ_CLASSES, i % READ_CLASSES);
            }
            r -= w;
        }
        // Rounding left `r` at the very top: the last class with weight.
        let i = weights.iter().rposition(|&w| w > 0.0).unwrap_or(0);
        (i / READ_CLASSES, i % READ_CLASSES)
    }

    /// Mate `m`'s P(class `c`), over every run bin, smoothed.
    fn class_share(&self, m: usize, c: usize) -> f64 {
        let at = |hh: usize, cc: usize| self.class_run[(m * RUN_BINS + hh) * READ_CLASSES + cc] as f64;
        let n_c: f64 = (0..RUN_BINS).map(|hh| at(hh, c)).sum();
        let n: f64 = (0..RUN_BINS).flat_map(|hh| (0..READ_CLASSES).map(move |cc| (hh, cc))).map(|(hh, cc)| at(hh, cc)).sum();
        (n_c + 0.5) / (n + 0.5 * READ_CLASSES as f64)
    }

    /// Mate `m`'s P(class `c` | run bin `h`), from the sample, shrunk toward
    /// P(class) by `RUN_PRIOR_READS` reads, so a run bin the sample barely
    /// holds draws the sample's classes rather than every class alike.
    fn class_given_run(&self, m: usize, h: usize, c: usize) -> f64 {
        let row = &self.class_run[(m * RUN_BINS + h.min(RUN_BINS - 1)) * READ_CLASSES..][..READ_CLASSES];
        let n_h: u64 = row.iter().sum();
        (row[c] as f64 + RUN_PRIOR_READS * self.class_share(m, c)) / (n_h as f64 + RUN_PRIOR_READS)
    }

    /// The learned `run_ratio(m, h, c)`, read from the table, with the run
    /// bin clamped exactly as `class_given_run` clamps it.
    fn run_ratio_at(&self, m: usize, h: usize, c: usize) -> f64 {
        self.run_ratios.as_ref().expect("the run-ratio table is filled in QualityProfile::learn")
            [(m * RUN_BINS + h.min(RUN_BINS - 1)) * READ_CLASSES + c]
    }

    /// P(class | run bin) / P(class), for mate `m`.
    fn run_ratio(&self, m: usize, h: usize, c: usize) -> f64 {
        self.class_given_run(m, h, c) / self.class_share(m, c)
    }

    /// An event's class mix, from its donor `pairs` against `reference`: per
    /// mate and class, the pool's count over the count the sample predicts
    /// from the pool's runs, so a region's poor reads are kept and the part
    /// its runs already explain is not counted twice.
    pub fn class_mix(&self, pairs: &[ReadPair], reference: &dyn RefSource) -> ClassMix {
        let mut mix = ClassMix::default();
        if self.class_cuts.is_empty() || pairs.is_empty() {
            return mix;
        }
        let mut seen = [[0f64; READ_CLASSES]; 2];
        let mut runs = [[0f64; RUN_BINS]; 2];
        for p in pairs {
            for (m, (qual, seq)) in [(&p.qual1, &p.seq1), (&p.qual2, &p.seq2)].into_iter().enumerate() {
                if qual.is_empty() {
                    continue;
                }
                let template = match &p.align {
                    Some(al) => template_of(reference, &p.chrom, &al[m], qual.len()),
                    None => seq.to_ascii_uppercase(),
                };
                seen[m][self.class_of(qual)] += 1.0;
                runs[m][template_run_bin(&template)] += 1.0;
            }
        }
        for m in 0..2 {
            for (c, (r, n)) in mix.r[m].iter_mut().zip(seen[m]).enumerate() {
                let expected: f64 = (0..RUN_BINS).map(|h| runs[m][h] * self.class_given_run(m, h, c)).sum();
                *r = (n + 0.5) / (expected + 0.5);
            }
        }
        mix
    }

    /// A fresh read of class `class`, or of a class drawn from mate
    /// `read_num`'s own share when the caller has none (a read on its own).
    pub fn start_read(&self, read_num: u8, class: Option<usize>, rng: &mut StdRng) -> ReadState {
        let class = class.unwrap_or_else(|| {
            let (c1, c2) = self.draw_classes((0, 0), None, rng);
            if read_num == 1 { c1 } else { c2 }
        });
        ReadState::new(class)
    }

    /// Draw the next base's quality (Phred+33) for mate `read_num`, with
    /// `cycles_left` cycles left in the read, this base included.
    pub fn next_quality(&self, read_num: u8, state: &ReadState, cycles_left: usize, rng: &mut StdRng) -> u8 {
        if self.alphabet.is_empty() {
            return FALLBACK_QUAL;
        }
        let table = &self.tables[if read_num == 1 { 0 } else { 1 }];
        let keys = self.keys(state, cycles_left);
        let counts = keys
            .iter()
            .filter_map(|k| table.get(k))
            .find(|c| c.total >= MIN_CONTEXT_OBS)
            .or_else(|| table.get(&keys[5]).filter(|c| c.total > 0));
        let Some(counts) = counts else { return FALLBACK_QUAL };
        let mut r = rng.gen_range(0..counts.total);
        for (s, &n) in counts.by_symbol.iter().enumerate() {
            if r < n {
                return self.alphabet[s];
            }
            r -= n;
        }
        unreachable!("a draw below a context's total lands in it")
    }

    /// Record that the read emitted `qual` (an `N`'s Q2 included), reading
    /// template base `base`.
    pub fn emitted(&self, state: &mut ReadState, qual: u8, base: u8) {
        if !self.alphabet.is_empty() {
            self.push(state, qual, base);
        }
    }

    /// Record whether the base the read just emitted was an error.
    pub fn record_error(&self, state: &mut ReadState, error: bool) {
        state.errors = (state.errors << 1) | error as u32;
    }

    /// The chance that a base the read emits at quality `qual` (Phred+33),
    /// `cycles_left` from its 3' end, is wrong.
    pub fn error_rate(&self, qual: u8, state: &ReadState, cycles_left: usize) -> f64 {
        let nominal = 10f64.powf(-(qual.saturating_sub(33) as f64) / 10.0);
        if self.errors.symbol.is_empty() {
            return nominal;
        }
        let (full, mid, class, symbol) = self.error_cells(qual, state, cycles_left);
        [self.errors.full[full], self.errors.mid[mid], self.errors.class[class], self.errors.symbol[symbol]]
            .into_iter()
            .find(|c| c.bases >= MIN_ERROR_BASES)
            .map(|c| c.errors as f64 / c.bases as f64)
            .unwrap_or(nominal)
    }

    /// A base's cells at the error table's four levels, finest first.
    fn error_cells(&self, qual: u8, state: &ReadState, cycles_left: usize) -> (usize, usize, usize, usize) {
        let symbol = self.symbol_of[qual as usize] as usize;
        let low = qual.saturating_sub(33) < LOW_Q;
        let lows = lows_bin((state.recent_low << 1) | low as u16);
        let run = state.run as usize;
        let err = err_bin(state.errors);
        let full = ((((symbol * READ_CLASSES + state.class) * LOWS_BINS + lows) * END_BINS + end_bin(cycles_left))
            * RUN_BINS
            + run)
            * ERR_BINS
            + err;
        let mid = ((symbol * LOWS_BINS + lows) * RUN_BINS + run) * ERR_BINS + err;
        (full, mid, symbol * READ_CLASSES + state.class, symbol)
    }

    /// Count one sample block's errors against its reference.
    ///
    /// Every base of a mate is counted as an error or not, or left out:
    /// - an aligned base (M, =, X) is an error when it differs from the
    ///   reference; an `N` on either side is left out;
    /// - an inserted base is left out;
    /// - a soft-clipped end is placed where it would have aligned. An end of
    ///   1-4 bases, or one matching the reference at ≥ 50% of its placed
    ///   bases, is a bad end and counted base by base; any other (adapter,
    ///   chimeric, foreign sequence) is left out;
    /// - a reference position where ≥ 5 of the block's reads have an aligned
    ///   base and ≥ 10% of them differ is the sample's own variant, and left out.
    ///
    /// On the HG002 35x chr20 slice this rule gives the same bad-end soft-clip
    /// share as a full classification of every clip (1.19% against 1.18%).
    fn block_errors(&self, block: &DonorBlock, templates: &[[Vec<u8>; 2]]) -> (ErrorTables, ErrorCensus) {
        let mut tables = ErrorTables::new(self.alphabet.len());
        let mut census = ErrorCensus::default();
        let index = |pos: i64| -> Option<usize> {
            let i = pos - block.ref_start as i64;
            (i >= 0 && (i as usize) < block.ref_seq.len()).then_some(i as usize)
        };

        // The sample's own variant positions, from the aligned bases.
        let mut pileup = vec![(0u32, 0u32); block.ref_seq.len()];
        for p in &block.pairs {
            let Some(al) = &p.align else { continue };
            for (m, seq) in al.iter().zip([&p.seq1, &p.seq2]) {
                for b in walk(m, seq) {
                    if let Placed::Aligned(pos, base) = b {
                        if let (Some(i), Some(r)) = (index(pos), block.base(pos)) {
                            if base != b'N' {
                                pileup[i].0 += 1;
                                pileup[i].1 += (base != r) as u32;
                            }
                        }
                    }
                }
            }
        }
        let is_variant = |pos: i64| index(pos).is_some_and(|i| pileup[i].0 >= 5 && pileup[i].1 * 10 >= pileup[i].0);

        for (p, t) in block.pairs.iter().zip(templates) {
            let Some(al) = &p.align else {
                census.pairs_without_alignment += 1;
                continue;
            };
            for (mate, (m, (seq, qual))) in al.iter().zip([(&p.seq1, &p.qual1), (&p.seq2, &p.qual2)]).enumerate() {
                let placed = walk(m, seq);
                // Each placed base's verdict, in alignment order.
                let mut verdict: Vec<Option<bool>> = vec![None; placed.len()];
                let clip_ok = |range: std::ops::Range<usize>| -> bool {
                    if range.len() <= 4 {
                        return true;
                    }
                    let (mut same, mut compared) = (0u32, 0u32);
                    for b in &placed[range] {
                        if let Placed::Clipped(pos, base) = *b {
                            if let Some(r) = block.base(pos).filter(|_| base != b'N') {
                                compared += 1;
                                same += (base == r) as u32;
                            }
                        }
                    }
                    compared > 0 && same * 2 >= compared
                };
                let lead = placed.iter().take_while(|b| matches!(b, Placed::Clipped(..))).count();
                let trail = placed.iter().rev().take_while(|b| matches!(b, Placed::Clipped(..))).count();
                let lead_ok = clip_ok(0..lead);
                let trail_ok = clip_ok(placed.len() - trail..placed.len());
                for (i, b) in placed.iter().enumerate() {
                    let (pos, base, ok_clip) = match *b {
                        Placed::Aligned(pos, base) => (pos, base, true),
                        Placed::Clipped(pos, base) => (pos, base, if i < lead { lead_ok } else { trail_ok }),
                        Placed::Inserted => {
                            census.other += 1;
                            continue;
                        }
                    };
                    if !ok_clip {
                        census.foreign_clip += 1;
                        continue;
                    }
                    let Some(r) = block.base(pos).filter(|_| base != b'N') else {
                        census.other += 1;
                        continue;
                    };
                    if is_variant(pos) {
                        census.variant_site += 1;
                        continue;
                    }
                    verdict[i] = Some(base != r);
                    census.counted += 1;
                }
                // Into sequencing order, then counted with the read's state.
                if m.reverse {
                    verdict.reverse();
                }
                let template = &t[mate];
                let mut state = ReadState::new(self.class_of(qual));
                let len = qual.len();
                for (c, &q) in qual.iter().enumerate() {
                    let err = verdict.get(c).copied().flatten();
                    if let Some(err) = err {
                        let cell = ErrorCell { errors: err as u64, bases: 1 };
                        let (full, mid, class, symbol) = self.error_cells(q, &state, len - c);
                        tables.full[full].add(cell);
                        tables.mid[mid].add(cell);
                        tables.class[class].add(cell);
                        tables.symbol[symbol].add(cell);
                    }
                    self.push(&mut state, q, template.get(c).copied().unwrap_or(b'N'));
                    self.record_error(&mut state, err == Some(true));
                }
            }
        }
        (tables, census)
    }

    fn log_summary(&self) {
        let cuts: Vec<String> = self.class_cuts.iter().map(|c| format!("{:.1}", c)).collect();
        let full_used = self.tables[0]
            .iter()
            .chain(self.tables[1].iter())
            .filter(|(k, c)| *k >> 58 == 0 && c.total >= MIN_CONTEXT_OBS)
            .count();
        let errors = match self.census {
            Some(c) => format!(
                "Error table from {} sampled bases ({} at the sample's own variant sites, {} in clipped ends that are not the reference read badly, {} inserted, N or off the reference left out)",
                c.counted, c.variant_site, c.foreign_clip, c.other
            ),
            None => "No error table: errors at 10^(-Q/10)".to_string(),
        };
        let qualities: Vec<String> = self.alphabet.iter().map(|q| (q - 33).to_string()).collect();
        let census_line = format!("{} full contexts with at least {} observations", full_used, MIN_CONTEXT_OBS);
        log::info!(
            "Quality profile: {} pairs in {} block(s); qualities {}; read classes cut at mean Q {}; {}. {}",
            self.n_pairs,
            self.n_blocks,
            qualities.join(","),
            cuts.join(","),
            census_line,
            errors
        );
        if let Some(warning) = crate::synth::thin_profile_warning(self.n_pairs, &census_line) {
            log::warn!("{}", warning);
        }
    }
}

/// The reference span `pairs` are placed on, clipped ends included: (lowest
/// start, highest end); `(0, 0)` with no alignment.
fn reads_span(pairs: &[ReadPair]) -> (u64, u64) {
    let (lo, hi) = pairs
        .iter()
        .filter_map(|p| p.align.as_ref())
        .flat_map(|al| al.iter().map(placed_span))
        .fold((u64::MAX, 0), |(a, b), (s, e)| (a.min(s), b.max(e)));
    if lo == u64::MAX {
        (0, 0)
    } else {
        (lo, hi)
    }
}

/// The input's 16 kb windows that hold reads, from its index, in contig
/// order: a BAI's leaf bins that have chunks, or each CRAI slice's span.
/// Without either index, every contig's whole length, in 16 kb steps.
pub fn indexed_windows(alignment_path: &str, ref_path: &str) -> Result<Vec<(String, u64, u64)>> {
    let header = crate::bam_stats::read_header(alignment_path, Some(ref_path))?;
    let names: Vec<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
    if crate::extract::is_cram(alignment_path) {
        let crai = format!("{}.crai", alignment_path);
        if let Ok(index) = noodles::cram::crai::fs::read(&crai) {
            return Ok(crai_windows(&index, &names));
        }
    } else {
        let stem = alignment_path.strip_suffix(".bam").map(|s| format!("{}.bai", s));
        for path in std::iter::once(format!("{}.bai", alignment_path)).chain(stem) {
            if let Ok(index) = noodles::bam::bai::fs::read(&path) {
                return Ok(bai_windows(&index, &names));
            }
        }
    }
    log::warn!(
        "No .bai or .crai beside {}: the quality sample is spread over the contigs' lengths, \
         and blocks may fall where the file has no reads",
        alignment_path
    );
    let mut windows = Vec::new();
    for (name, rs) in header.reference_sequences() {
        let len = usize::from(rs.length()) as u64;
        windows.extend((0..len).step_by(1 << 14).map(|s| (name.to_string(), s, s + (1 << 14))));
    }
    Ok(windows)
}

/// `windows` on contigs `fasta` holds, in order, and the names of the other
/// contigs (sorted, once each).
fn windows_on_reference(windows: Vec<(String, u64, u64)>, fasta: &HashSet<String>) -> (Vec<(String, u64, u64)>, Vec<String>) {
    let (kept, other): (Vec<_>, Vec<_>) = windows.into_iter().partition(|w| fasta.contains(&w.0));
    let mut dropped: Vec<String> = other.into_iter().map(|w| w.0).collect();
    dropped.sort();
    dropped.dedup();
    (kept, dropped)
}

/// A BAI's leaf bins (16 kb, ids 4681-37448) that hold chunks, per contig.
fn bai_windows(index: &noodles::bam::bai::Index, names: &[String]) -> Vec<(String, u64, u64)> {
    const FIRST_LEAF: usize = 4681;
    const LEAVES: usize = 1 << 15;
    let mut out = Vec::new();
    for (rs, name) in index.reference_sequences().iter().zip(names) {
        let mut starts: Vec<u64> = rs
            .bins()
            .iter()
            .filter(|(id, bin)| (FIRST_LEAF..FIRST_LEAF + LEAVES).contains(*id) && !bin.chunks().is_empty())
            .map(|(id, _)| ((id - FIRST_LEAF) as u64) << 14)
            .collect();
        starts.sort_unstable();
        out.extend(starts.into_iter().map(|s| (name.clone(), s, s + (1 << 14))));
    }
    out
}

/// Each CRAI slice's reference span, in contig and position order.
fn crai_windows(index: &[noodles::cram::crai::Record], names: &[String]) -> Vec<(String, u64, u64)> {
    let mut spans: Vec<(usize, u64, u64)> = index
        .iter()
        .filter_map(|r| {
            let id = r.reference_sequence_id()?;
            let start = usize::from(r.alignment_start()?) as u64 - 1;
            Some((id, start, start + r.alignment_span() as u64))
        })
        .filter(|(id, _, _)| *id < names.len())
        .collect();
    spans.sort_unstable();
    spans.into_iter().map(|(id, s, e)| (names[id].clone(), s, e)).collect()
}

/// `n` blocks of `len` bases spread evenly over `windows` (in contig order):
/// block `k` starts at window `floor((k + 0.5) · W / n)`, moved to the end of
/// the block before it when they would overlap on the same contig.
pub fn choose_blocks(windows: &[(String, u64, u64)], n: usize, len: u64) -> Vec<(String, u64, u64)> {
    let w = windows.len();
    let mut out: Vec<(String, u64, u64)> = Vec::new();
    if w == 0 {
        return out;
    }
    for k in 0..n {
        let (chrom, start, _) = &windows[((2 * k + 1) * w) / (2 * n)];
        let mut start = *start;
        if let Some((prev, _, prev_end)) = out.last() {
            if prev == chrom && *prev_end > start {
                start = *prev_end;
            }
        }
        out.push((chrom.clone(), start, start + len));
    }
    out
}

/// The startup sample: `SAMPLE_BLOCKS` blocks of `SAMPLE_BLOCK_LEN` bases
/// placed by `choose_blocks` over the `indexed_windows` on contigs the FASTA
/// holds (an input with none is refused), their pairs extracted as
/// donor pools are (`--min-mapq`, the same filters), on the thread pool, and
/// deduplicated by name; each with the reference under its reads.
pub fn sample_input(alignment_path: &str, ref_path: &str, min_mapq: u8) -> Result<Vec<DonorBlock>> {
    // Only contigs the FASTA holds can be read against it: an input aligned
    // to a larger reference (decoys, HLA) lists others, and a block there
    // would abort the run.
    let fasta: HashSet<String> = crate::reference::fasta_contigs(ref_path)?.into_iter().map(|(name, _)| name).collect();
    let (windows, dropped) = windows_on_reference(indexed_windows(alignment_path, ref_path)?, &fasta);
    let examples = dropped.iter().take(3).cloned().collect::<Vec<_>>().join(", ");
    if windows.is_empty() && !dropped.is_empty() {
        bail!(
            "none of the contigs {} lists ({} of them, e.g. {}) are in the reference FASTA {}: \
             the input and the FASTA name their contigs differently (e.g. \"20\" against \"chr20\") \
             or come from different builds",
            alignment_path,
            dropped.len(),
            examples,
            ref_path
        );
    }
    if !dropped.is_empty() {
        log::info!(
            "Quality sample: {} contig(s) of {} are not in the reference FASTA and are left out (e.g. {})",
            dropped.len(),
            alignment_path,
            examples
        );
    }
    let placed = choose_blocks(&windows, SAMPLE_BLOCKS, SAMPLE_BLOCK_LEN);
    let mut blocks: Vec<DonorBlock> = placed
        .par_iter()
        .map(|(chrom, start, end)| -> Result<DonorBlock> {
            let pairs = crate::extract::extract_read_pairs(alignment_path, chrom, *start, *end, min_mapq, Some(ref_path))?.pairs;
            let (lo, hi) = reads_span(&pairs);
            let (ref_start, ref_seq) = if pairs.is_empty() {
                (*start, Vec::new())
            } else {
                crate::reference::fetch_window(ref_path, chrom, lo.saturating_sub(BLOCK_REF_PAD), hi + BLOCK_REF_PAD)?
            };
            Ok(DonorBlock { chrom: chrom.clone(), ref_start, ref_seq, pairs })
        })
        .collect::<Result<_>>()?;
    let mut seen: HashSet<String> = HashSet::new();
    for b in &mut blocks {
        b.pairs.retain(|p| seen.insert(p.name.clone()));
    }
    Ok(blocks)
}

/// One base of a mate, in alignment order, placed on the reference.
#[derive(Debug, Clone, Copy, PartialEq)]
enum Placed {
    /// An aligned base at reference position (0-based) and its letter.
    Aligned(i64, u8),
    /// A soft-clipped base, placed where it would have aligned.
    Clipped(i64, u8),
    /// An inserted base: no reference position.
    Inserted,
}

/// `seq` (sequencing order) of a mate aligned as `m`, base by base in
/// alignment order and placed on the reference.
fn walk(m: &MateAlignment, seq: &[u8]) -> Vec<Placed> {
    let mut s = seq.to_vec();
    if m.reverse {
        crate::extract::reverse_complement(&mut s);
    }
    let s: Vec<u8> = s.iter().map(u8::to_ascii_uppercase).collect();
    let mut out = Vec::with_capacity(s.len());
    let mut rp = m.start as i64;
    let mut qi = 0usize;
    let mut aligned_yet = false;
    for &(op, len) in &m.cigar {
        let len = len as usize;
        match op {
            b'S' => {
                let start = if aligned_yet { rp } else { rp - len as i64 };
                for k in 0..len {
                    if let Some(&b) = s.get(qi + k) {
                        out.push(Placed::Clipped(start + k as i64, b));
                    }
                }
                qi += len;
            }
            b'M' => {
                aligned_yet = true;
                for k in 0..len {
                    if let Some(&b) = s.get(qi + k) {
                        out.push(Placed::Aligned(rp + k as i64, b));
                    }
                }
                qi += len;
                rp += len as i64;
            }
            b'I' => {
                aligned_yet = true;
                out.extend(std::iter::repeat_n(Placed::Inserted, len.min(s.len().saturating_sub(qi))));
                qi += len;
            }
            b'D' => {
                aligned_yet = true;
                rp += len as i64;
            }
            _ => {}
        }
    }
    out
}

/// The reference span a mate's bases are placed on, clipped ends included.
fn placed_span(m: &MateAlignment) -> (u64, u64) {
    let lead: u64 = m.cigar.iter().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as u64).sum();
    let trail: u64 = m.cigar.iter().rev().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as u64).sum();
    let ref_len: u64 = m.cigar.iter().filter(|(op, _)| *op == b'M' || *op == b'D').map(|(_, l)| *l as u64).sum();
    (m.start.saturating_sub(lead), m.start + ref_len + trail)
}

fn mean_phred(qual: &[u8]) -> f64 {
    qual.iter().map(|&q| q.saturating_sub(33) as f64).sum::<f64>() / qual.len() as f64
}

/// The `f` quantile of sorted `v`, interpolated linearly (numpy's default).
fn quantile(v: &[f64], f: f64) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let x = f * (v.len() - 1) as f64;
    let (lo, hi) = (x.floor() as usize, x.ceil() as usize);
    v[lo] + (v[hi] - v[lo]) * (x - lo as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    fn mate(start: u64, reverse: bool, cigar: &[(u8, u32)]) -> MateAlignment {
        MateAlignment { start, reverse, cigar: cigar.to_vec() }
    }

    #[test]
    fn test_walk_places_clipped_ends_where_they_would_align() {
        let m = mate(100, false, &[(b'S', 2), (b'M', 3), (b'I', 1), (b'M', 2), (b'D', 1), (b'M', 1), (b'S', 2)]);
        let placed = walk(&m, b"ccAAAgTTCgg");
        assert_eq!(
            placed,
            vec![
                Placed::Clipped(98, b'C'),
                Placed::Clipped(99, b'C'),
                Placed::Aligned(100, b'A'),
                Placed::Aligned(101, b'A'),
                Placed::Aligned(102, b'A'),
                Placed::Inserted,
                Placed::Aligned(103, b'T'),
                Placed::Aligned(104, b'T'),
                Placed::Aligned(106, b'C'),
                Placed::Clipped(107, b'G'),
                Placed::Clipped(108, b'G'),
            ]
        );
        assert_eq!(placed_span(&m), (98, 109));
    }

    #[test]
    fn test_walk_flips_a_reverse_mate_back_to_alignment_order() {
        // Stored in sequencing order: the reverse complement of "AACGT".
        let m = mate(10, true, &[(b'M', 5)]);
        assert_eq!(
            walk(&m, b"ACGTT"),
            vec![
                Placed::Aligned(10, b'A'),
                Placed::Aligned(11, b'A'),
                Placed::Aligned(12, b'C'),
                Placed::Aligned(13, b'G'),
                Placed::Aligned(14, b'T'),
            ]
        );
    }

    #[test]
    fn test_quantile_interpolates_like_numpy() {
        let v = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(quantile(&v, 0.0), 1.0);
        assert_eq!(quantile(&v, 1.0), 4.0);
        assert!((quantile(&v, 0.5) - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_run_bins_follow_the_measured_effect() {
        assert_eq!(run_bin(b'T', 6), 0);
        assert_eq!(run_bin(b'T', 7), 1);
        assert_eq!(run_bin(b'A', 11), 2);
        assert_eq!(run_bin(b'G', 12), 3);
        assert_eq!(run_bin(b'C', 4), 0);
        assert_eq!(run_bin(b'C', 5), 2);
        assert_eq!(run_bin(b'C', 7), 3);
        assert_eq!(run_bin(b'N', 20), 0);
        assert_eq!(template_run_bin(b"ACGTttttttttttttACG"), 3, "12 T's in lowercase");
        assert_eq!(template_run_bin(b"ACGTTTTTTTTNTTTTACG"), 1, "an N breaks the run into 7 and 4");
    }

    // --- the run-ratio table (2026-10-07-quality-speed.md) ---

    /// A learned profile whose class/run counts are asymmetric: R1's class
    /// follows its template's run bin, R2's does not. Learned with no
    /// reference, so each mate's template is its own called bases and the run
    /// bin is set by the sequence here.
    fn asymmetric_profile() -> QualityProfile {
        let rl = 50usize;
        // ACGT repeating: the longest run is 1, so run bin 0.
        let plain = |n: usize| -> Vec<u8> { (0..rl).map(|c| b"ACGT"[(c + n) % 4]).collect() };
        // Exactly 12 T's, walled by A on both sides: run bin 3.
        let long_run = |n: usize| -> Vec<u8> {
            let mut s = plain(n);
            s[9] = b'A';
            s[10..22].fill(b'T');
            s[22] = b'A';
            s
        };
        // Exactly 7 T's, walled the same way: run bin 1.
        let short_run = |n: usize| -> Vec<u8> {
            let mut s = plain(n);
            s[9] = b'A';
            s[10..17].fill(b'T');
            s[17] = b'A';
            s
        };
        let pairs: Vec<ReadPair> = (0..400usize)
            .map(|i| {
                // R1: a long run means a poor read, no run means a good one.
                // The means are spread inside each half so every class is used.
                let (seq1, q1) = if i % 2 == 0 {
                    (long_run(i), 8 + (i / 2 % 10) as u8)
                } else {
                    (plain(i), 28 + (i / 2 % 10) as u8)
                };
                // R2: its run bin and its quality turn independently of each other.
                let seq2 = match i % 3 {
                    0 => plain(i),
                    1 => short_run(i),
                    _ => long_run(i),
                };
                let q2 = 5 + (i * 7 % 33) as u8;
                ReadPair {
                    name: format!("p{i}"),
                    seq1,
                    qual1: vec![b'!' + q1; rl],
                    seq2,
                    qual2: vec![b'!' + q2; rl],
                    ref_start: 0,
                    ref_end: rl as u64,
                    insert_size: rl as i64,
                    chrom: "chr1".to_string(),
                    align: None,
                }
            })
            .collect();
        QualityProfile::from_read_pairs(&pairs, rl)
    }

    /// The table must hold `run_ratio(m, h, c)` itself, bit for bit, in every
    /// one of its `2 * RUN_BINS * READ_CLASSES` cells.
    #[test]
    fn test_the_run_ratio_table_is_the_formula_cell_for_cell() {
        let p = asymmetric_profile();
        let cells = || {
            (0..2).flat_map(|m| (0..RUN_BINS).flat_map(move |h| (0..READ_CLASSES).map(move |c| (m, h, c))))
        };
        for (m, h, c) in cells() {
            assert_eq!(
                p.run_ratio_at(m, h, c).to_bits(),
                p.run_ratio(m, h, c).to_bits(),
                "cell (mate {m}, run bin {h}, class {c}): table {} against the formula {}",
                p.run_ratio_at(m, h, c),
                p.run_ratio(m, h, c)
            );
        }
        // And the table has teeth, so this test can go red.
        assert!(
            cells().any(|(m, h, c)| (p.run_ratio_at(m, h, c) - 1.0).abs() > 0.2),
            "every ratio is about 1: a neutral table would pass unseen"
        );
        assert!(
            (0..READ_CLASSES).any(|c| (p.run_ratio_at(0, 3, c) - p.run_ratio_at(1, 3, c)).abs() > 0.2),
            "the two mates' rows are alike: a swapped mate would pass unseen"
        );
        assert!(
            (0..READ_CLASSES).any(|c| (p.run_ratio_at(0, 3, c) - p.run_ratio_at(0, 0, c)).abs() > 0.2),
            "run bin 3 and run bin 0 are alike: a clamped run bin would pass unseen"
        );
    }

    /// `draw_classes` must read the table, not recompute the formula. Without
    /// this, delegating `run_ratio_at` straight back to `run_ratio` throws the
    /// whole speed fix away with every other test still green -- and only the
    /// one-off S1 timing would notice (whole-branch review, finding 2).
    #[test]
    fn test_draw_classes_reads_the_table_rather_than_recomputing_it() {
        let mut p = asymmetric_profile();
        let runs = (3usize, 0usize);
        let before: Vec<(usize, usize)> =
            (0..40u64).map(|s| p.draw_classes(runs, None, &mut StdRng::seed_from_u64(s))).collect();
        // One cell of mate 0's run bin 3 row made overwhelming. A `draw_classes`
        // that reads the table must now draw that class for mate 0 every time;
        // one that recomputes `run_ratio` cannot see the poke at all.
        let (poked, at) = (5usize, |m: usize, h: usize, c: usize| (m * RUN_BINS + h) * READ_CLASSES + c);
        p.run_ratios.as_mut().expect("the table is filled")[at(0, 3, poked)] = 1e12;
        let after: Vec<(usize, usize)> =
            (0..40u64).map(|s| p.draw_classes(runs, None, &mut StdRng::seed_from_u64(s))).collect();
        assert!(
            after.iter().all(|&(c1, _)| c1 == poked),
            "poking (mate 0, run bin 3, class {poked}) to 1e12 left mate 0 drawing {:?}: draw_classes is not reading the table",
            after.iter().map(|&(c1, _)| c1).collect::<HashSet<_>>()
        );
        assert!(
            before.iter().any(|&(c1, _)| c1 != poked),
            "mate 0 already drew class {poked} every time before the poke: the test proves nothing"
        );
    }

    /// `draw_classes` must pick what the weights built straight from
    /// `run_ratio` pick -- same product order, same one `f64` off the RNG.
    #[test]
    fn test_draw_classes_picks_what_the_run_ratio_formula_picks() {
        let p = asymmetric_profile();
        let neutral = ClassMix::default();
        let oracle = |runs: (usize, usize), rng: &mut StdRng| -> (usize, usize) {
            let weights: Vec<f64> = (0..READ_CLASSES * READ_CLASSES)
                .map(|i| {
                    let (c1, c2) = (i / READ_CLASSES, i % READ_CLASSES);
                    p.joint[i] as f64
                        * p.run_ratio(0, runs.0, c1)
                        * p.run_ratio(1, runs.1, c2)
                        * neutral.r[0][c1]
                        * neutral.r[1][c2]
                })
                .collect();
            let total: f64 = weights.iter().sum();
            assert!(total > 0.0, "the oracle's weights are all zero");
            let mut r = rng.gen::<f64>() * total;
            for (i, &w) in weights.iter().enumerate() {
                if r < w {
                    return (i / READ_CLASSES, i % READ_CLASSES);
                }
                r -= w;
            }
            let i = weights.iter().rposition(|&w| w > 0.0).unwrap_or(0);
            (i / READ_CLASSES, i % READ_CLASSES)
        };
        let mut checked = 0usize;
        let mut seen = HashSet::new();
        for runs in [(3usize, 0usize), (0, 3), (3, 1), (1, 3)] {
            for s in 0..200u64 {
                let got = p.draw_classes(runs, None, &mut StdRng::seed_from_u64(s));
                let want = oracle(runs, &mut StdRng::seed_from_u64(s));
                assert_eq!(got, want, "runs {runs:?}, seed {s}");
                seen.insert(got);
                checked += 1;
            }
        }
        assert_eq!(checked, 800);
        assert!(seen.len() > 4, "only {} distinct pairs drawn: too few to see a wrong weight", seen.len());
    }

    #[test]
    fn test_a_template_is_read_from_the_mates_5_prime_end_in_sequencing_order() {
        let mut seqs = std::collections::HashMap::new();
        seqs.insert("chr1".to_string(), b"AAAACCCCGGGGTTTT".to_vec());
        let r = SharedReference::from_sequences(seqs);
        // Forward, 2 bases soft-clipped at the start: the template starts 2 before the alignment.
        assert_eq!(template_of(&r, "chr1", &mate(4, false, &[(b'S', 2), (b'M', 4)]), 6), b"AACCCC".to_vec());
        // Reverse, aligned 8-12 with 2 clipped at the end: its 5' end is 14; read back, complemented.
        assert_eq!(template_of(&r, "chr1", &mate(8, true, &[(b'M', 4), (b'S', 2)]), 6), b"AACCCC".to_vec());
        // Off the contig: padded with N past the end.
        assert_eq!(template_of(&r, "chr1", &mate(12, false, &[(b'M', 4)]), 6), b"TTTTNN".to_vec());
    }

    use std::collections::HashMap as StdHashMap;

    fn reference_of(seq: &[u8]) -> SharedReference {
        let mut seqs = StdHashMap::new();
        seqs.insert("chr1".to_string(), seq.to_vec());
        SharedReference::from_sequences(seqs)
    }

    /// A donor pair on chr1 whose R1 is `r1` (already in sequencing order,
    /// forward) aligned as `cigar` at `start`, with an R2 that matches the
    /// reference exactly, 50 bp, forward, at `r2_start`.
    fn donor(name: &str, r1: &[u8], start: u64, cigar: &[(u8, u32)], reference: &[u8], r2_start: u64) -> ReadPair {
        let r2 = reference[r2_start as usize..r2_start as usize + 50].to_vec();
        ReadPair {
            name: name.to_string(),
            seq1: r1.to_vec(),
            qual1: vec![b'!' + 30; r1.len()],
            seq2: r2,
            qual2: vec![b'!' + 30; 50],
            ref_start: start,
            ref_end: r2_start + 50,
            insert_size: (r2_start + 50 - start) as i64,
            chrom: "chr1".to_string(),
            align: Some(Box::new([mate(start, false, cigar), mate(r2_start, false, &[(b'M', 50)])])),
        }
    }

    #[test]
    fn test_learn_errors_counts_bad_ends_and_leaves_out_adapters_variants_and_insertions() {
        let mut rng = StdRng::seed_from_u64(3);
        let reference: Vec<u8> = (0..4000).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect();
        let r = |a: usize, b: usize| reference[a..b].to_vec();
        let other = |b: u8| if b == b'A' { b'C' } else { b'A' };
        let mut pairs = Vec::new();

        // A bad end: 80 aligned bases, then 20 clipped bases that are the
        // reference read badly -- 6 of them wrong. 6 errors.
        let mut bad_end = r(100, 200);
        for k in [82usize, 85, 88, 91, 94, 97] {
            bad_end[k] = other(bad_end[k]);
        }
        pairs.push(donor("bad_end", &bad_end, 100, &[(b'M', 80), (b'S', 20)], &reference, 400));

        // An adapter: 80 aligned bases, then 20 clipped adapter bases.
        let mut adapter = r(600, 680);
        adapter.extend_from_slice(b"AGATCGGAAGAGCACACGTC");
        pairs.push(donor("adapter", &adapter, 600, &[(b'M', 80), (b'S', 20)], &reference, 900));

        // The sample's own variant: six reads all differ at position 1500.
        for i in 0..6u64 {
            let mut v = r(1450, 1550);
            v[50] = other(v[50]);
            pairs.push(donor(&format!("variant_{}", i), &v, 1450, &[(b'M', 100)], &reference, 1700 + i * 10));
        }

        // Two inserted bases.
        let mut ins = r(2000, 2050);
        ins.extend_from_slice(b"TT");
        ins.extend_from_slice(&reference[2050..2098]);
        pairs.push(donor("insertion", &ins, 2000, &[(b'M', 50), (b'I', 2), (b'M', 48)], &reference, 2300));

        let p = QualityProfile::from_donor_pairs(&pairs, 100, &reference_of(&reference));
        let census = p.census.expect("an error table was learned");
        let errors: u64 = p.errors.symbol.iter().map(|c| c.errors).sum();
        let bases: u64 = p.errors.symbol.iter().map(|c| c.bases).sum();
        assert_eq!(errors, 6, "only the bad end's 6 wrong bases are errors");
        assert_eq!(census.foreign_clip, 20, "the adapter's 20 clipped bases are left out");
        assert_eq!(census.variant_site, 6, "the six reads' bases at the variant are left out");
        assert_eq!(census.other, 2, "the two inserted bases are left out");
        assert_eq!(bases, census.counted, "every counted base is in the table");
        // 9 pairs: R1s 100 + 100 + 6 x 100 + 100 bases, R2s 9 x 50, less those left out.
        assert_eq!(census.counted, 100 + 80 + 600 - 6 + 98 + 9 * 50);
    }

    #[test]
    fn test_an_empty_profile_draws_q20_and_errs_at_the_nominal_rate() {
        let p = QualityProfile::from_read_pairs(&[], 100);
        let mut rng = StdRng::seed_from_u64(1);
        let s = p.start_read(1, None, &mut rng);
        assert_eq!(p.next_quality(1, &s, 50, &mut rng), b'!' + 20);
        assert!((p.error_rate(b'!' + 20, &s, 50) - 0.01).abs() < 1e-12);
    }

    // --- the startup sample's placement (T7, T8) ---

    /// T7: an input whose reads lie only in 30.0-31.2 Mb of a 60 Mb contig, as
    /// a slice does. Every block must land inside the reads, none overlapping.
    #[test]
    fn test_sample_blocks_land_inside_a_slice_and_do_not_overlap() {
        let windows: Vec<(String, u64, u64)> = (0..(1_200_000 / 16_384))
            .map(|k| ("chr20".to_string(), 30_000_000 + k * 16_384, 30_000_000 + (k + 1) * 16_384))
            .collect();
        let blocks = choose_blocks(&windows, SAMPLE_BLOCKS, SAMPLE_BLOCK_LEN);
        assert_eq!(blocks.len(), SAMPLE_BLOCKS);
        for (i, (c, s, e)) in blocks.iter().enumerate() {
            assert_eq!(c, "chr20");
            assert!(*s >= 30_000_000 && *e <= 31_250_000, "block {} at {}-{} is outside the reads", i, s, e);
            if i > 0 {
                assert!(blocks[i - 1].2 <= *s, "blocks {} and {} overlap", i - 1, i);
            }
        }
    }

    /// T7: blocks spread over every contig that has reads, in proportion to
    /// how much of it does.
    #[test]
    fn test_sample_blocks_spread_over_the_contigs_with_reads() {
        let mut windows = Vec::new();
        for (c, n) in [("chr1", 300u64), ("chr2", 100), ("chr3", 0)] {
            windows.extend((0..n).map(|k| (c.to_string(), k * 1_000_000, k * 1_000_000 + 16_384)));
        }
        let blocks = choose_blocks(&windows, 20, SAMPLE_BLOCK_LEN);
        let on = |c: &str| blocks.iter().filter(|b| b.0 == c).count();
        assert_eq!((on("chr1"), on("chr2"), on("chr3")), (15, 5, 0));
    }

    /// T7: the windows are read from a BAI's leaf bins that hold chunks, not
    /// from its other bins.
    #[test]
    fn test_bai_windows_are_the_leaf_bins_that_hold_reads() {
        let mut bai: Vec<u8> = Vec::new();
        bai.extend_from_slice(b"BAI\x01");
        bai.extend_from_slice(&2u32.to_le_bytes()); // two references
        let chunk = |b: &mut Vec<u8>| {
            b.extend_from_slice(&1u32.to_le_bytes());
            b.extend_from_slice(&(1u64 << 16).to_le_bytes());
            b.extend_from_slice(&(2u64 << 16).to_le_bytes());
        };
        // ref 0: leaf bins at 30 Mb and 30 Mb + 16 kb, plus a level-4 bin that is not a window
        let leaf = |pos: u64| 4681 + (pos >> 14) as u32;
        bai.extend_from_slice(&3u32.to_le_bytes());
        for id in [leaf(30_000_000), leaf(30_016_384), 585 + (30_000_000u64 >> 17) as u32] {
            bai.extend_from_slice(&id.to_le_bytes());
            chunk(&mut bai);
        }
        bai.extend_from_slice(&0u32.to_le_bytes()); // no linear index
        // ref 1: no bins
        bai.extend_from_slice(&0u32.to_le_bytes());
        bai.extend_from_slice(&0u32.to_le_bytes());
        let path = std::env::temp_dir().join(format!("spike_bai_windows_{}.bai", std::process::id()));
        std::fs::write(&path, bai).unwrap();
        let index = noodles::bam::bai::fs::read(&path).unwrap();
        std::fs::remove_file(&path).ok();
        let names = vec!["chrA".to_string(), "chrB".to_string()];
        let w = bai_windows(&index, &names);
        let leaf_start = |pos: u64| (pos >> 14) << 14;
        assert_eq!(
            w,
            vec![
                ("chrA".to_string(), leaf_start(30_000_000), leaf_start(30_000_000) + 16_384),
                ("chrA".to_string(), leaf_start(30_016_384), leaf_start(30_016_384) + 16_384),
            ]
        );
    }

    /// An input aligned to a larger reference than the FASTA (a decoy-aware
    /// CRAM with the no-alt FASTA) lists windows on contigs the FASTA lacks.
    /// A block there cannot be read, and the run aborted; such windows are
    /// left out, and every other window is kept in its order.
    #[test]
    fn test_sample_windows_on_contigs_the_reference_lacks_are_left_out() {
        let w = |c: &str, s: u64| (c.to_string(), s, s + 16_384);
        let windows = vec![w("chr1", 0), w("chrUn_JTFH01000277v1_decoy", 0), w("chr1", 16_384), w("chr2", 0), w("HLA-A*01:01:01:01", 0)];
        let fasta: HashSet<String> = ["chr1", "chr2", "chrM"].iter().map(|s| s.to_string()).collect();
        let (kept, dropped) = windows_on_reference(windows, &fasta);
        assert_eq!(kept, vec![w("chr1", 0), w("chr1", 16_384), w("chr2", 0)]);
        assert_eq!(dropped, vec!["HLA-A*01:01:01:01".to_string(), "chrUn_JTFH01000277v1_decoy".to_string()]);
    }

    /// A BAM in `dir` of 20 proper pairs on `contig` (200 kb), indexed by one
    /// 16 kb leaf bin so the sampler offers a window there, and a FASTA that
    /// holds only `chrA` (200 kb). Returns (BAM, FASTA).
    fn sample_fixture(dir: &std::path::Path, contig: &str) -> (String, String) {
        use noodles::sam::alignment::record::cigar::{op::Kind, Op};
        use noodles::sam::alignment::record::{Flags, MappingQuality};
        use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
        use noodles::sam::alignment::RecordBuf;
        std::fs::create_dir_all(dir).unwrap();
        let read: Vec<u8> = (0..100).map(|i| b"ACGT"[i % 4]).collect();
        let mut records = Vec::new();
        for i in 0..20usize {
            let start = 1_000 + i * 500;
            for first in [true, false] {
                let (pos, mate) = if first { (start, start + 200) } else { (start + 200, start) };
                records.push(
                    RecordBuf::builder()
                        .set_name(format!("p{i}"))
                        .set_flags(Flags::from(if first { 0x63u16 } else { 0x93u16 }))
                        .set_reference_sequence_id(0)
                        .set_alignment_start(noodles::core::Position::new(pos).unwrap())
                        .set_mapping_quality(MappingQuality::new(60).unwrap())
                        .set_cigar([Op::new(Kind::Match, 100)].into_iter().collect())
                        .set_mate_reference_sequence_id(0)
                        .set_mate_alignment_start(noodles::core::Position::new(mate).unwrap())
                        .set_template_length(if first { 300 } else { -300 })
                        .set_sequence(Sequence::from(read.clone()))
                        .set_quality_scores(QualityScores::from(vec![30u8; 100]))
                        .build(),
                );
            }
        }
        records.sort_by_key(|r| r.alignment_start());
        let bam = crate::extract::test_fixtures::write_one_contig_bam(&dir.join("sample.bam"), contig, 200_000, &records);
        // The fixture's .bai holds bin 0; make it the 16 kb leaf bin 4681 that
        // the sampler reads windows from.
        let bai_path = format!("{bam}.bai");
        let mut bai = std::fs::read(&bai_path).unwrap();
        bai[12..16].copy_from_slice(&4681u32.to_le_bytes());
        std::fs::write(&bai_path, bai).unwrap();
        let fasta = dir.join("ref.fa");
        std::fs::write(&fasta, format!(">chrA\n{}\n", "ACGT".repeat(50_000))).unwrap();
        std::fs::write(dir.join("ref.fa.fai"), "chrA\t200000\t6\t200000\t200001\n").unwrap();
        let windows = indexed_windows(&bam, fasta.to_str().unwrap()).unwrap();
        assert!(windows.iter().any(|w| w.0 == contig), "the fixture must offer a window on {contig}: {windows:?}");
        (bam, fasta.to_str().unwrap().to_string())
    }

    /// End to end: the sample learns from reads on a contig the FASTA holds.
    /// Guards the filter's wiring from the other side: one that drops every
    /// window would leave this sample empty.
    #[test]
    fn test_the_startup_sample_reads_the_contigs_the_fasta_holds() {
        let dir = std::env::temp_dir().join(format!("spike_sample_chra_{}", std::process::id()));
        let (bam, fasta) = sample_fixture(&dir, "chrA");
        let blocks = sample_input(&bam, &fasta, 20);
        std::fs::remove_dir_all(&dir).ok();
        let blocks = blocks.expect("the startup sample of a chrA BAM against a chrA FASTA");
        let pairs: usize = blocks.iter().map(|b| b.pairs.len()).sum();
        assert_eq!(pairs, 20, "every pair on chrA is sampled");
        assert!(blocks.iter().all(|b| b.chrom == "chrA"));
    }

    /// End to end: an input whose only reads sit on a contig the FASTA lacks
    /// (a decoy, or "20" against "chr20") is refused with the reason named.
    /// Before the fix a block landed there and the run aborted in
    /// `fetch_window`; with the windows dropped but nothing said, it aborted
    /// later blaming the index and --min-mapq.
    #[test]
    fn test_the_startup_sample_names_a_contig_the_fasta_lacks() {
        let dir = std::env::temp_dir().join(format!("spike_sample_decoy_{}", std::process::id()));
        let (bam, fasta) = sample_fixture(&dir, "chrUn_decoy");
        let got = sample_input(&bam, &fasta, 20);
        std::fs::remove_dir_all(&dir).ok();
        let err = format!("{:#}", got.expect_err("an input with no contig in the FASTA must be refused"));
        assert!(
            err.contains("none of") && err.contains("chrUn_decoy") && err.contains("name their contigs differently"),
            "the error must name the contig mismatch, got: {err}"
        );
    }

    /// T8: a CRAM's windows are its slices' spans, by contig and position.
    #[test]
    fn test_crai_windows_are_the_slices_spans() {
        let rec = |id: Option<usize>, start: usize, span: usize| {
            noodles::cram::crai::Record::new(id, noodles::core::Position::new(start), span, 0, 0, 0)
        };
        let index = vec![rec(Some(1), 501, 100), rec(Some(0), 30_000_001, 20_000), rec(None, 0, 0), rec(Some(0), 1, 10)];
        let names = vec!["chrA".to_string(), "chrB".to_string()];
        assert_eq!(
            crai_windows(&index, &names),
            vec![
                ("chrA".to_string(), 0, 10),
                ("chrA".to_string(), 30_000_000, 30_020_000),
                ("chrB".to_string(), 500, 600),
            ]
        );
    }
}
