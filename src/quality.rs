//! Base qualities and sequencing errors, learned from the sample's own reads.
//!
//! **Qualities.** fqzcomp's quality context (htscodecs `fqzcomp_qual`, as CRAM
//! uses it) turned into a generator. Each base's quality is drawn given:
//! - the last qualities of the read;
//! - how many cycles are left to its 3' end;
//! - whether its quality has already changed bin twice;
//! - and the read's *class*: where its mean quality falls among the pool's.
//!
//! A pair's two classes are drawn together, from the donor pairs' joint table,
//! so mates share state. A first-order chain (the model this replaces) forgets
//! a read's state within a few bases, so all its reads came out average: on
//! HG002 35x the SD of a read's mean quality was 0.52 against the sample's 2.08.
//!
//! **Errors.** Each base is wrong at the rate the donor reads show for its
//! quality, the read's class, the run of low qualities ending at it, and its
//! distance from the 3' end, counted against the reference (`learn_errors`).
//! A binned quality understates the errors of a read that is falling apart:
//! aligned bases at Q25 err 0.07% in good HG002 reads and 0.95% in poor ones.

use std::collections::HashMap;

use rand::rngs::StdRng;
use rand::Rng;
use rayon::prelude::*;

use crate::reference::SharedReference;
use crate::types::{MateAlignment, ReadPair};

/// Read classes (the "selector"), cut at these quantiles of the pool's read
/// mean qualities (both mates pooled). The low end is cut finely, because the
/// poor reads are few and differ most.
pub const READ_CLASSES: usize = 8;
const CLASS_QUANTILES: [f64; READ_CLASSES - 1] = [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8];

/// A context needs this many observations before its counts are drawn from;
/// below it, the draw backs off to a coarser context.
const MIN_CONTEXT_OBS: u32 = 20;

/// An error-table cell needs this many counted donor bases before its rate
/// is used; below it, the rate backs off to a coarser cell, then to 10^(-Q/10).
const MIN_ERROR_BASES: u64 = 200;

/// A quality below this (Phred) is "low" for the error table's low-quality run.
const LOW_Q: u8 = 15;

/// Quality bins for the history of a large alphabet and the change flag:
/// Q0-9, Q10-19, Q20-29, Q30+.
pub const PREV_Q_BINS: usize = 4;

/// Bins for the run of low qualities ending at a base: 0, 1, 2-3, 4-7, 8+.
const RUN_BINS: usize = 5;
/// Bins for the cycles left to a read's 3' end, the base itself counted:
/// 1-10, 11-30, 31-60, 61+.
const END_BINS: usize = 4;

/// Slots of quality history kept; an alphabet of more than 4 values keeps the
/// last quality exactly and the two before it in `PREV_Q_BINS` bins.
const HISTORY: usize = 5;
/// An empty history slot (the read's first cycles).
const EMPTY: u8 = u8::MAX;

/// Phred+33 byte drawn when the profile learned nothing at all.
const FALLBACK_QUAL: u8 = b'!' + 20;

/// A quality's bin in [`PREV_Q_BINS`]: Q0-9 → 0, Q10-19 → 1, Q20-29 → 2, Q30+ → 3.
pub fn prev_q_bin(phred_plus_33: u8) -> usize {
    match phred_plus_33.saturating_sub(33) {
        0..=9 => 0,
        10..=19 => 1,
        20..=29 => 2,
        _ => 3,
    }
}

fn run_bin(run: u32) -> usize {
    match run {
        0 => 0,
        1 => 1,
        2..=3 => 2,
        4..=7 => 3,
        _ => 4,
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

/// One read being generated (or walked, when learning): its class and what
/// the context needs from the qualities so far.
#[derive(Debug, Clone)]
pub struct ReadState {
    class: usize,
    /// Symbols of the last qualities, most recent first; `EMPTY` before any.
    last: [u8; HISTORY],
    /// Times the quality changed bin (`prev_q_bin`) so far.
    changes: u32,
    /// Qualities below `LOW_Q` in a row, ending at the last one.
    low_run: u32,
}

impl ReadState {
    fn new(class: usize) -> Self {
        Self { class, last: [EMPTY; HISTORY], changes: 0, low_run: 0 }
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

/// How the donor bases fared when counted for the error table.
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
    /// Donor pairs without an alignment (none in real runs).
    pub pairs_without_alignment: u64,
}

/// The learned quality and error model.
#[cfg_attr(test, derive(PartialEq))]
#[derive(Debug)]
pub struct QualityProfile {
    /// The pool's quality bytes (Phred+33), sorted.
    alphabet: Vec<u8>,
    /// A byte's symbol: its index in `alphabet`, or the nearest one's.
    symbol_of: Vec<u8>,
    /// Whether the history keeps 5 exact qualities (an alphabet of ≤ 4).
    small: bool,
    /// Cycles left are binned as `min(7, cycles_left >> pshift)`.
    pshift: u32,
    /// Read-class cuts, on mean Phred quality.
    class_cuts: Vec<f64>,
    /// Joint (R1 class, R2 class) counts over donor pairs, flattened.
    joint: Vec<u64>,
    /// Per mate, counts per context key (all backoff levels in one map).
    tables: [HashMap<u64, Counts>; 2],
    /// Error cells per (symbol, class, run bin, end bin), then per (symbol,
    /// class), then per symbol; empty when no errors were learned.
    errors_full: Vec<ErrorCell>,
    errors_class: Vec<ErrorCell>,
    errors_symbol: Vec<ErrorCell>,
    n_pairs: usize,
}

impl QualityProfile {
    /// Learn the qualities from `pairs`, with no error table: every base errs
    /// at 10^(-Q/10). For tests.
    #[cfg(test)]
    pub fn from_read_pairs(pairs: &[ReadPair], read_length: usize) -> Self {
        Self::learn(pairs, read_length, None)
    }

    /// Learn the qualities and the error table from donor `pairs`, counting
    /// their mismatches against `reference` (`learn_errors`).
    pub fn from_donor_pairs(pairs: &[ReadPair], read_length: usize, reference: &SharedReference) -> Self {
        Self::learn(pairs, read_length, Some(reference))
    }

    fn learn(pairs: &[ReadPair], read_length: usize, reference: Option<&SharedReference>) -> Self {
        let mut seen = [false; 256];
        for p in pairs {
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
            tables: [HashMap::new(), HashMap::new()],
            errors_full: Vec::new(),
            errors_class: Vec::new(),
            errors_symbol: Vec::new(),
            n_pairs: pairs.len(),
        };
        if profile.alphabet.is_empty() {
            log::info!("Quality profile: no donor qualities; every base is Q20");
            return profile;
        }

        // Read classes: cut the pool's read means at fixed quantiles.
        let mut means: Vec<f64> = pairs
            .iter()
            .flat_map(|p| [&p.qual1, &p.qual2])
            .filter(|q| !q.is_empty())
            .map(|q| mean_phred(q))
            .collect();
        means.sort_by(|a, b| a.partial_cmp(b).unwrap());
        profile.class_cuts = CLASS_QUANTILES.iter().map(|&f| quantile(&means, f)).collect();

        for p in pairs {
            let (c1, c2) = (profile.class_of(&p.qual1), profile.class_of(&p.qual2));
            profile.joint[c1 * READ_CLASSES + c2] += 1;
        }

        // The context counts, per mate, in chunks on the thread pool; counts
        // are sums, so neither the chunking nor the join order changes them.
        for mate in 0..2 {
            let table = pairs
                .par_chunks(2048)
                .map(|chunk| {
                    let mut t: HashMap<u64, Counts> = HashMap::new();
                    for p in chunk {
                        let qual = if mate == 0 { &p.qual1 } else { &p.qual2 };
                        profile.count_read(qual, &mut t);
                    }
                    t
                })
                .reduce(HashMap::new, |mut a, b| {
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

        let census = reference.map(|r| profile.learn_errors(pairs, r));
        profile.log_summary(census);
        profile
    }

    /// The class of a read with quality string `qual`.
    fn class_of(&self, qual: &[u8]) -> usize {
        if qual.is_empty() {
            return READ_CLASSES / 2;
        }
        let m = mean_phred(qual);
        self.class_cuts.iter().filter(|&&c| c <= m).count()
    }

    fn phred(&self, symbol: u8) -> u8 {
        self.alphabet[symbol as usize].saturating_sub(33)
    }

    /// Add one quality to `state`, as the read emitted it.
    fn push(&self, state: &mut ReadState, qual: u8) {
        let symbol = self.symbol_of[qual as usize];
        if state.last[0] != EMPTY
            && prev_q_bin(self.alphabet[state.last[0] as usize]) != prev_q_bin(self.alphabet[symbol as usize])
        {
            state.changes += 1;
        }
        state.last.rotate_right(1);
        state.last[0] = symbol;
        state.low_run = if self.phred(symbol) < LOW_Q { state.low_run + 1 } else { 0 };
    }

    /// The context keys for the next base, finest first: the full context,
    /// then (2-quality history, position, class), (last quality, position,
    /// class), (last quality, position), (position), and the whole mate.
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
        let tag = |level: u64| level << 58;
        [
            tag(0) | full << 8 | pos << 4 | flag << 3 | class,
            tag(1) | two << 8 | pos << 4 | class,
            tag(2) | one << 8 | pos << 4 | class,
            tag(3) | one << 8 | pos << 4,
            tag(4) | pos << 4,
            tag(5),
        ]
    }

    /// Count one donor read's qualities into `table`.
    fn count_read(&self, qual: &[u8], table: &mut HashMap<u64, Counts>) {
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
            self.push(&mut state, q);
        }
    }

    /// Draw a pair's two read classes, from the donor pairs' joint table.
    pub fn draw_classes(&self, rng: &mut StdRng) -> (usize, usize) {
        let total: u64 = self.joint.iter().sum();
        if total == 0 {
            return (READ_CLASSES / 2, READ_CLASSES / 2);
        }
        let mut r = rng.gen_range(0..total);
        for (i, &n) in self.joint.iter().enumerate() {
            if r < n {
                return (i / READ_CLASSES, i % READ_CLASSES);
            }
            r -= n;
        }
        unreachable!("a draw below the joint table's total lands in it")
    }

    /// A fresh read of class `class`, or of a class drawn from mate
    /// `read_num`'s own share when the caller has none (a read on its own).
    pub fn start_read(&self, read_num: u8, class: Option<usize>, rng: &mut StdRng) -> ReadState {
        let class = class.unwrap_or_else(|| {
            let (c1, c2) = self.draw_classes(rng);
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

    /// Record that the read emitted `qual` (an `N`'s Q2 included).
    pub fn emitted(&self, state: &mut ReadState, qual: u8) {
        if !self.alphabet.is_empty() {
            self.push(state, qual);
        }
    }

    /// The chance that a base the read emits at quality `qual` (Phred+33),
    /// `cycles_left` from its 3' end, is wrong.
    pub fn error_rate(&self, qual: u8, state: &ReadState, cycles_left: usize) -> f64 {
        let nominal = 10f64.powf(-(qual.saturating_sub(33) as f64) / 10.0);
        if self.errors_symbol.is_empty() {
            return nominal;
        }
        let symbol = self.symbol_of[qual as usize] as usize;
        let run = if qual.saturating_sub(33) < LOW_Q { state.low_run + 1 } else { 0 };
        let full = self.errors_full[self.error_cell(symbol, state.class, run_bin(run), end_bin(cycles_left))];
        let class = self.errors_class[symbol * READ_CLASSES + state.class];
        let alone = self.errors_symbol[symbol];
        [full, class, alone]
            .into_iter()
            .find(|c| c.bases >= MIN_ERROR_BASES)
            .map(|c| c.errors as f64 / c.bases as f64)
            .unwrap_or(nominal)
    }

    fn error_cell(&self, symbol: usize, class: usize, run: usize, end: usize) -> usize {
        ((symbol * READ_CLASSES + class) * RUN_BINS + run) * END_BINS + end
    }

    /// Count the donor reads' errors against `reference` into the error table.
    ///
    /// Every base of a donor mate is counted as an error or not, or left out:
    /// - an aligned base (M, =, X) is an error when it differs from the
    ///   reference; an `N` on either side is left out;
    /// - an inserted base is left out;
    /// - a soft-clipped end is placed where it would have aligned. An end of
    ///   1-4 bases, or one matching the reference at ≥ 50% of its placed
    ///   bases, is a bad end and counted base by base; any other (adapter,
    ///   chimeric, foreign sequence) is left out;
    /// - a reference position where ≥ 5 donor reads have an aligned base and
    ///   ≥ 10% of them differ is the sample's own variant, and left out.
    ///
    /// On the HG002 35x chr20 slice this rule gives the same bad-end soft-clip
    /// share as a full classification of every clip (1.19% against 1.18%).
    fn learn_errors(&mut self, pairs: &[ReadPair], reference: &SharedReference) -> ErrorCensus {
        let a = self.alphabet.len();
        self.errors_full = vec![ErrorCell::default(); a * READ_CLASSES * RUN_BINS * END_BINS];
        self.errors_class = vec![ErrorCell::default(); a * READ_CLASSES];
        self.errors_symbol = vec![ErrorCell::default(); a];
        let mut census = ErrorCensus::default();

        // The reference under each chromosome's donor reads, fetched once.
        let mut spans: HashMap<&str, (u64, u64)> = HashMap::new();
        for p in pairs {
            let Some(al) = &p.align else {
                census.pairs_without_alignment += 1;
                continue;
            };
            for m in al.iter() {
                let (s, e) = placed_span(m);
                let span = spans.entry(p.chrom.as_str()).or_insert((s, e));
                span.0 = span.0.min(s);
                span.1 = span.1.max(e);
            }
        }
        let refs: HashMap<&str, (u64, Vec<u8>)> = spans
            .into_iter()
            .filter_map(|(chrom, (s, e))| {
                let start = s.saturating_sub(1);
                reference
                    .fetch_sequence(chrom, start, e + 1)
                    .ok()
                    .map(|seq| (chrom, (start, seq.to_ascii_uppercase())))
            })
            .collect();
        let ref_base = |chrom: &str, pos: i64| -> Option<u8> {
            let (start, seq) = refs.get(chrom)?;
            let i = pos - *start as i64;
            (i >= 0).then(|| seq.get(i as usize).copied()).flatten().filter(|&b| b != b'N')
        };

        // The sample's own variant positions, from the aligned bases.
        let mut pileup: HashMap<(&str, i64), (u32, u32)> = HashMap::new();
        for p in pairs {
            let Some(al) = &p.align else { continue };
            for (m, seq) in al.iter().zip([&p.seq1, &p.seq2]) {
                for b in walk(m, seq) {
                    if let Placed::Aligned(pos, base) = b {
                        if base != b'N' {
                            if let Some(r) = ref_base(&p.chrom, pos) {
                                let e = pileup.entry((p.chrom.as_str(), pos)).or_insert((0, 0));
                                e.0 += 1;
                                e.1 += (base != r) as u32;
                            }
                        }
                    }
                }
            }
        }
        let is_variant = |chrom: &str, pos: i64| {
            pileup.get(&(chrom, pos)).is_some_and(|&(cov, mm)| cov >= 5 && mm * 10 >= cov)
        };

        for p in pairs {
            let Some(al) = &p.align else { continue };
            for (mate, (m, (seq, qual))) in
                al.iter().zip([(&p.seq1, &p.qual1), (&p.seq2, &p.qual2)]).enumerate()
            {
                let _ = mate;
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
                            if base != b'N' {
                                if let Some(r) = ref_base(&p.chrom, pos) {
                                    compared += 1;
                                    same += (base == r) as u32;
                                }
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
                    let Some(r) = ref_base(&p.chrom, pos).filter(|_| base != b'N') else {
                        census.other += 1;
                        continue;
                    };
                    if is_variant(&p.chrom, pos) {
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
                let mut state = ReadState::new(self.class_of(qual));
                let len = qual.len();
                for (c, &q) in qual.iter().enumerate() {
                    if let Some(Some(err)) = verdict.get(c) {
                        let symbol = self.symbol_of[q as usize] as usize;
                        let run = if q.saturating_sub(33) < LOW_Q { state.low_run + 1 } else { 0 };
                        let cell = ErrorCell { errors: *err as u64, bases: 1 };
                        let i = self.error_cell(symbol, state.class, run_bin(run), end_bin(len - c));
                        self.errors_full[i].add(cell);
                        self.errors_class[symbol * READ_CLASSES + state.class].add(cell);
                        self.errors_symbol[symbol].add(cell);
                    }
                    self.push(&mut state, q);
                }
            }
        }
        census
    }

    fn log_summary(&self, census: Option<ErrorCensus>) {
        let cuts: Vec<String> = self.class_cuts.iter().map(|c| format!("{:.1}", c)).collect();
        let full_used = self.tables[0]
            .iter()
            .chain(self.tables[1].iter())
            .filter(|(k, c)| *k >> 58 == 0 && c.total >= MIN_CONTEXT_OBS)
            .count();
        let errors = match census {
            Some(c) => format!(
                "Error table from {} donor bases ({} at the sample's own variant sites, {} in clipped ends that are not the reference read badly, {} inserted, N or off the reference left out)",
                c.counted, c.variant_site, c.foreign_clip, c.other
            ),
            None => "No error table: errors at 10^(-Q/10)".to_string(),
        };
        let qualities: Vec<String> = self.alphabet.iter().map(|q| (q - 33).to_string()).collect();
        let census_line = format!("{} full contexts with at least {} observations", full_used, MIN_CONTEXT_OBS);
        log::info!(
            "Quality profile: {} pairs; qualities {}; read classes cut at mean Q {}; {}. {}",
            self.n_pairs,
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

/// One base of a donor mate, in alignment order, placed on the reference.
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

        let mut p = QualityProfile::from_read_pairs(&pairs, 100);
        let census = p.learn_errors(&pairs, &reference_of(&reference));
        let errors: u64 = p.errors_symbol.iter().map(|c| c.errors).sum();
        let bases: u64 = p.errors_symbol.iter().map(|c| c.bases).sum();
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
}
