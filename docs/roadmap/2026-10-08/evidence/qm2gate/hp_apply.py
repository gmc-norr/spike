"""Throwaway Gate B change (2026-10-07): add to spike's quality model, switchable by SPIKE_QX bits,
  1  lows in the last 16 qualities (5 bins) in the quality context, levels 0-1
  2  the longest one-letter run read so far (sticky, 4 bins, C weighted) in the quality context, levels 0-3
  4  the error table's low-run bin replaced by the lows-in-last-16 bin
  8  the run bin added to the error table, with a (quality, lows/run bin, run) backoff level
and two ignored tests: zz_gate_large_sample (quality strings, held out) and zz_gate_errors
(predicted against observed errors on held-out pairs, given their real qualities).
Run from the qm worktree: python3 ../qm2gate/hp_apply.py; revert with git checkout -- src/."""
import re

def edit(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        assert s.count(old) == 1, (path, old[:80])
        s = s.replace(old, new)
    open(path, "w").write(s)

Q = "src/quality.rs"
edit(Q, [
("""fn end_bin(cycles_left: usize) -> usize {""",
"""/// Experiment switches (throwaway): SPIKE_QX bits, read once.
fn qx() -> u32 {
    static QX: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    *QX.get_or_init(|| std::env::var("SPIKE_QX").ok().and_then(|v| v.parse().ok()).unwrap_or(0))
}

/// Lows among the last 16 qualities: 0, 1-2, 3-5, 6-9, 10+.
fn lows_bin(recent_low: u16) -> usize {
    match recent_low.count_ones() {
        0 => 0,
        1..=2 => 1,
        3..=5 => 2,
        6..=9 => 3,
        _ => 4,
    }
}

/// A one-letter run's bin, by its measured effect on the cycles after it (HG002 35x):
/// A/G/T 7-8 -> 1 (~1.2x low quality), 9-11 -> 2 (~2.3x), 12+ -> 3 (~7-9x); C 5-6 -> 2, 7+ -> 3.
fn hp_bin(base: u8, len: u32) -> u8 {
    match (base, len) {
        (b'C', 5..=6) => 2,
        (b'C', 7..) => 3,
        (b'C', _) => 0,
        (b'A' | b'G' | b'T', 7..=8) => 1,
        (b'A' | b'G' | b'T', 9..=11) => 2,
        (b'A' | b'G' | b'T', 12..) => 3,
        _ => 0,
    }
}
const HP_BINS: usize = 4;

fn end_bin(cycles_left: usize) -> usize {"""),
("""    /// Qualities below `LOW_Q` in a row, ending at the last one.
    low_run: u32,
}

impl ReadState {
    fn new(class: usize) -> Self {
        Self { class, last: [EMPTY; HISTORY], changes: 0, low_run: 0 }
    }
}""",
"""    /// Qualities below `LOW_Q` in a row, ending at the last one.
    low_run: u32,
    /// Which of the last 16 qualities were below `LOW_Q` (bit 0 the last).
    recent_low: u16,
    /// The template's current one-letter run, and the largest run bin so far.
    run_base: u8,
    run_len: u32,
    hp: u8,
}

impl ReadState {
    fn new(class: usize) -> Self {
        Self { class, last: [EMPTY; HISTORY], changes: 0, low_run: 0, recent_low: 0, run_base: b'N', run_len: 0, hp: 0 }
    }
}"""),
("""    fn push(&self, state: &mut ReadState, qual: u8) {""", """    fn push(&self, state: &mut ReadState, qual: u8, base: u8) {
        let base = base.to_ascii_uppercase();
        if base == state.run_base && base != b'N' {
            state.run_len += 1;
        } else {
            state.run_base = base;
            state.run_len = (base != b'N') as u32;
        }
        state.hp = state.hp.max(hp_bin(state.run_base, state.run_len));"""),
("""        state.low_run = if self.phred(symbol) < LOW_Q { state.low_run + 1 } else { 0 };
    }""", """        let low = self.phred(symbol) < LOW_Q;
        state.low_run = if low { state.low_run + 1 } else { 0 };
        state.recent_low = (state.recent_low << 1) | low as u16;
    }"""),
("""        let class = state.class as u64;
        let tag = |level: u64| level << 58;
        [
            tag(0) | full << 8 | pos << 4 | flag << 3 | class,
            tag(1) | two << 8 | pos << 4 | class,
            tag(2) | one << 8 | pos << 4 | class,
            tag(3) | one << 8 | pos << 4,""", """        let class = state.class as u64;
        let lows = if qx() & 1 != 0 { lows_bin(state.recent_low) as u64 } else { 0 };
        let hp = if qx() & 2 != 0 { state.hp as u64 } else { 0 };
        let tag = |level: u64| level << 58;
        [
            tag(0) | hp << 36 | lows << 32 | full << 8 | pos << 4 | flag << 3 | class,
            tag(1) | hp << 36 | lows << 32 | two << 8 | pos << 4 | class,
            tag(2) | hp << 36 | one << 8 | pos << 4 | class,
            tag(3) | hp << 36 | one << 8 | pos << 4,"""),
("""    fn count_read(&self, qual: &[u8], table: &mut HashMap<u64, Counts>) {""",
 """    fn count_read(&self, qual: &[u8], seq: &[u8], table: &mut HashMap<u64, Counts>) {"""),
("""            self.push(&mut state, q);
        }
    }

    /// Draw a pair's""", """            self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
        }
    }

    /// Draw a pair's"""),
("""                        let qual = if mate == 0 { &p.qual1 } else { &p.qual2 };
                        profile.count_read(qual, &mut t);""", """                        let (qual, seq) = if mate == 0 { (&p.qual1, &p.seq1) } else { (&p.qual2, &p.seq2) };
                        profile.count_read(qual, seq, &mut t);"""),
("""    pub fn emitted(&self, state: &mut ReadState, qual: u8) {
        if !self.alphabet.is_empty() {
            self.push(state, qual);""", """    pub fn emitted(&self, state: &mut ReadState, qual: u8, base: u8) {
        if !self.alphabet.is_empty() {
            self.push(state, qual, base);"""),
("""        let symbol = self.symbol_of[qual as usize] as usize;
        let run = if qual.saturating_sub(33) < LOW_Q { state.low_run + 1 } else { 0 };
        let full = self.errors_full[self.error_cell(symbol, state.class, run_bin(run), end_bin(cycles_left))];
        let class = self.errors_class[symbol * READ_CLASSES + state.class];
        let alone = self.errors_symbol[symbol];
        [full, class, alone]""", """        let symbol = self.symbol_of[qual as usize] as usize;
        let (rb, hp) = self.error_bins(qual, state);
        let full = self.errors_full[self.error_cell(symbol, state.class, rb, end_bin(cycles_left), hp)];
        let mid = if qx() & 12 != 0 { self.errors_mid[(symbol * RUN_BINS + rb) * HP_BINS + hp] } else { ErrorCell::default() };
        let class = self.errors_class[symbol * READ_CLASSES + state.class];
        let alone = self.errors_symbol[symbol];
        [full, mid, class, alone]"""),
("""    fn error_cell(&self, symbol: usize, class: usize, run: usize, end: usize) -> usize {
        ((symbol * READ_CLASSES + class) * RUN_BINS + run) * END_BINS + end
    }""", """    fn error_cell(&self, symbol: usize, class: usize, run: usize, end: usize, hp: usize) -> usize {
        (((symbol * READ_CLASSES + class) * RUN_BINS + run) * END_BINS + end) * HP_BINS + hp
    }

    /// The error table's (low bin, run bin) for a base emitted at `qual` after `state`.
    fn error_bins(&self, qual: u8, state: &ReadState) -> (usize, usize) {
        let low = qual.saturating_sub(33) < LOW_Q;
        let rb = if qx() & 4 != 0 {
            lows_bin((state.recent_low << 1) | low as u16)
        } else {
            run_bin(if low { state.low_run + 1 } else { 0 })
        };
        let hp = if qx() & 8 != 0 { state.hp as usize } else { 0 };
        (rb, hp)
    }"""),
("""        self.errors_full = vec![ErrorCell::default(); a * READ_CLASSES * RUN_BINS * END_BINS];""",
 """        self.errors_full = vec![ErrorCell::default(); a * READ_CLASSES * RUN_BINS * END_BINS * HP_BINS];
        self.errors_mid = vec![ErrorCell::default(); a * RUN_BINS * HP_BINS];"""),
("""                        let symbol = self.symbol_of[q as usize] as usize;
                        let run = if q.saturating_sub(33) < LOW_Q { state.low_run + 1 } else { 0 };
                        let cell = ErrorCell { errors: *err as u64, bases: 1 };
                        let i = self.error_cell(symbol, state.class, run_bin(run), end_bin(len - c));
                        self.errors_full[i].add(cell);""", """                        let symbol = self.symbol_of[q as usize] as usize;
                        let (rb, hp) = self.error_bins(q, &state);
                        let cell = ErrorCell { errors: *err as u64, bases: 1 };
                        let i = self.error_cell(symbol, state.class, rb, end_bin(len - c), hp);
                        self.errors_full[i].add(cell);
                        self.errors_mid[(symbol * RUN_BINS + rb) * HP_BINS + hp].add(cell);"""),
("""                    self.push(&mut state, q);
                }
            }
        }
        census""", """                    self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
                }
                if let Some(v) = verdicts_out.as_mut() {
                    v.push(verdict.clone());
                }
            }
        }
        census"""),
("""    errors_full: Vec<ErrorCell>,
    errors_class""", """    errors_full: Vec<ErrorCell>,
    errors_mid: Vec<ErrorCell>,
    errors_class"""),
("""            errors_full: Vec::new(),
            errors_class""", """            errors_full: Vec::new(),
            errors_mid: Vec::new(),
            errors_class"""),
("""    fn learn_errors(&mut self, pairs: &[ReadPair], reference: &SharedReference) -> ErrorCensus {""",
 """    fn learn_errors(&mut self, pairs: &[ReadPair], reference: &SharedReference) -> ErrorCensus {
        self.learn_errors_v(pairs, reference, None)
    }

    /// Throwaway: the verdicts of `pairs` (R1, R2, ... in sequencing order) against `reference`,
    /// counted into a scratch table, for the held-out error check.
    pub fn verdicts_of(&self, pairs: &[ReadPair], reference: &SharedReference) -> Vec<Vec<Option<bool>>> {
        let mut scratch = Self { tables: [HashMap::new(), HashMap::new()], ..self.clone_shallow() };
        let mut out = Vec::new();
        scratch.learn_errors_v(pairs, reference, Some(&mut out));
        out
    }

    fn clone_shallow(&self) -> Self {
        Self {
            alphabet: self.alphabet.clone(),
            symbol_of: self.symbol_of.clone(),
            small: self.small,
            pshift: self.pshift,
            class_cuts: self.class_cuts.clone(),
            joint: self.joint.clone(),
            tables: [HashMap::new(), HashMap::new()],
            errors_full: Vec::new(),
            errors_mid: Vec::new(),
            errors_class: Vec::new(),
            errors_symbol: Vec::new(),
            n_pairs: self.n_pairs,
        }
    }

    /// Throwaway: walk a mate's real qualities and bases, giving each base's (state before it,
    /// cycles left) so a caller can ask `error_rate`.
    pub fn walk_states(&self, qual: &[u8], seq: &[u8]) -> Vec<ReadState> {
        let mut state = ReadState::new(self.class_of(qual));
        let mut out = Vec::with_capacity(qual.len());
        for (c, &q) in qual.iter().enumerate() {
            out.push(state.clone());
            self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
        }
        out
    }

    fn learn_errors_v(&mut self, pairs: &[ReadPair], reference: &SharedReference, mut verdicts_out: Option<&mut Vec<Vec<Option<bool>>>>) -> ErrorCensus {"""),
])
# the hp state is public to tests through a getter
edit(Q, [("""impl ReadState {
    fn new""", """impl ReadState {
    pub fn hp(&self) -> u8 {
        self.hp
    }

    fn new""")])

S = "src/synth.rs"
edit(S, [
("""                seq.push(b'N');
                qual.push(N_QUAL);
                self.profile.emitted(&mut state, N_QUAL);""", """                seq.push(b'N');
                qual.push(N_QUAL);
                self.profile.emitted(&mut state, N_QUAL, b'N');"""),
("""                        seq.push(random_base(rng));
                        qual.push(q);
                        self.profile.emitted(&mut state, q);""", """                        seq.push(random_base(rng));
                        qual.push(q);
                        self.profile.emitted(&mut state, q, true_base);"""),
("""                    seq.push(random_different_base(true_base, rng));
                    qual.push(q);
                    self.profile.emitted(&mut state, q);""", """                    seq.push(random_different_base(true_base, rng));
                    qual.push(q);
                    self.profile.emitted(&mut state, q, true_base);"""),
("""                seq.push(true_base);
                qual.push(q);
                self.profile.emitted(&mut state, q);""", """                seq.push(true_base);
                qual.push(q);
                self.profile.emitted(&mut state, q, true_base);"""),
("""                out.push(q);
                profile.emitted(&mut state, q);""", """                out.push(q);
                profile.emitted(&mut state, q, b);"""),
])
s = open(S).read()
end = s.rstrip().rfind("}")
tests = r'''
    /// Throwaway Gate B: 20 blocks of 100 kb (chr1-20 at 30% of their length), split by name hash.
    fn zz_gate_blocks(bam: &str) -> (Vec<ReadPair>, Vec<ReadPair>) {
        use crate::extract::extract_read_pairs;
        let lens: Vec<(&str, u64)> = vec![("chr1",248956422),("chr2",242193529),("chr3",198295559),("chr4",190214555),("chr5",181538259),("chr6",170805979),("chr7",159345973),("chr8",145138636),("chr9",138394717),("chr10",133797422),("chr11",135086622),("chr12",133275309),("chr13",114364328),("chr14",107043718),("chr15",101991189),("chr16",90338345),("chr17",83257441),("chr18",80373285),("chr19",58617616),("chr20",64444167)];
        let fnv = |name: &str| name.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ b as u64).wrapping_mul(0x100_0000_01b3));
        let (mut train, mut held) = (Vec::new(), Vec::new());
        for (c, l) in &lens {
            let s = (*l as f64 * 0.3) as u64;
            for p in extract_read_pairs(bam, c, s, s + 100_000, 20, None).unwrap().pairs {
                if fnv(&p.name) % 2 == 0 { train.push(p) } else { held.push(p) }
            }
        }
        (train, held)
    }

    #[test]
    #[ignore]
    fn zz_gate_large_sample() {
        let bam = std::env::var("SPIKE_N7_BAM").unwrap();
        let (train, held) = zz_gate_blocks(&bam);
        let real: QualPairs = held.iter().map(|p| (p.qual1.clone(), p.qual2.clone())).collect();
        println!("train {} held {}", train.len(), held.len());
        let crash = |set: &QualPairs| { let v: Vec<&Vec<u8>> = set.iter().flat_map(|(a,b)| [a,b]).collect(); 100.0 * v.iter().filter(|q| q.len() >= 20 && q[q.len()-20..].iter().filter(|&&x| x < b'!' + 15).count() >= 10).count() as f64 / v.len() as f64 };
        // crash share among held-out reads whose own bases hold a long run (bin 3), real and fake
        let profile = QualityProfile::from_read_pairs(&train, 151);
        let mut rng = StdRng::seed_from_u64(7);
        let fake = n7_draw(&profile, &held, &mut rng);
        let m = k1_metrics(&fake, &real);
        println!("QX {}: SD ratio {:.3} perfect z {:.2} crash z {:.2} (fake {:.2}% real {:.2}%)",
                 std::env::var("SPIKE_QX").unwrap_or_default(), m[0], m[1], m[2], crash(&fake), crash(&real));
        for want in 0..4u8 {
            let mut f = Vec::new();
            let mut r = Vec::new();
            for (i, p) in held.iter().enumerate() {
                for (k, (seq, q)) in [(&p.seq1, &p.qual1), (&p.seq2, &p.qual2)].into_iter().enumerate() {
                    let st = profile.walk_states(q, seq);
                    let hp = st.last().map(|s| s.hp()).unwrap_or(0);
                    if hp == want {
                        r.push((q.clone(), q.clone()));
                        let fq = if k == 0 { &fake[i].0 } else { &fake[i].1 };
                        f.push((fq.clone(), fq.clone()));
                    }
                }
            }
            println!("  reads whose own bases reach run bin {}: {} -- crashed fake {:.2}% real {:.2}%", want, r.len(), crash(&f), crash(&r));
        }
    }

    #[test]
    #[ignore]
    fn zz_gate_errors() {
        let bam = std::env::var("SPIKE_N7_BAM").unwrap();
        let fasta = std::env::var("SPIKE_REF").unwrap();
        let chroms: Vec<String> = (1..=20).map(|i| format!("chr{}", i)).collect();
        let chrom_refs: Vec<&str> = chroms.iter().map(|s| s.as_str()).collect();
        let reference = crate::reference::SharedReference::load(&fasta, &chrom_refs).unwrap();
        let (train, held) = zz_gate_blocks(&bam);
        let profile = QualityProfile::from_donor_pairs(&train, 151, &reference);
        let verdicts = profile.verdicts_of(&held, &reference);
        let mates: Vec<(&Vec<u8>, &Vec<u8>)> = held.iter().flat_map(|p| [(&p.qual1, &p.seq1), (&p.qual2, &p.seq2)]).collect();
        assert_eq!(verdicts.len(), mates.len());
        // groups: [all, crashed tails Q11, Q25, Q37, after run bin 3, good reads' last 40]
        let mut obs = [0f64; 6];
        let mut pred = [0f64; 6];
        let mut n = [0f64; 6];
        let mut per_read: Vec<(f64, f64, f64)> = Vec::new(); // crashed: observed, predicted, bases (Q11, last 40)
        let mut rng = StdRng::seed_from_u64(11);
        let mut per_read_sim: Vec<f64> = Vec::new();
        for ((q, seq), v) in mates.iter().zip(&verdicts) {
            let len = q.len();
            let crashed = len >= 20 && q[len - 20..].iter().filter(|&&x| x < b'!' + 15).count() >= 10;
            let states = profile.walk_states(q, seq);
            let (mut ro, mut rp, mut rn, mut rs) = (0f64, 0f64, 0f64, 0f64);
            for c in 0..len {
                let Some(Some(err)) = v.get(c) else { continue };
                let p = profile.error_rate(q[c], &states[c], len - c);
                let e = *err as u8 as f64;
                let phred = q[c] - 33;
                let tail = c + 40 >= len;
                let mut add = |g: usize| { obs[g] += e; pred[g] += p; n[g] += 1.0; };
                add(0);
                if crashed && tail {
                    match phred { 0..=14 => add(1), 15..=29 => add(2), _ => add(3) }
                }
                if states[c].hp() == 3 { add(4); }
                if !crashed && tail { add(5); }
                if crashed && tail && phred < 15 {
                    ro += e; rp += p; rn += 1.0; rs += (rng.gen::<f64>() < p) as u8 as f64;
                }
            }
            if rn >= 10.0 { per_read.push((ro / rn, rp / rn, rn)); per_read_sim.push(rs / rn); }
        }
        let names = ["all counted bases", "crashed tails, Q<15", "crashed tails, Q15-29", "crashed tails, Q30+", "after a run of bin 3", "not crashed, last 40"];
        println!("QX {}: held-out bases, observed error rate against the model's rate given the real qualities", std::env::var("SPIKE_QX").unwrap_or_default());
        for g in 0..6 {
            println!("  {:24} bases {:9}  observed {:.4}  predicted {:.4}  ratio {:.2}", names[g], n[g], obs[g] / n[g], pred[g] / n[g], pred[g] / obs[g]);
        }
        let sd = |v: &[f64]| { let m = v.iter().sum::<f64>() / v.len() as f64; (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt() };
        let o: Vec<f64> = per_read.iter().map(|x| x.0).collect();
        println!("  crashed reads (>= 10 Q<15 tail bases): {}; per-read observed rate SD {:.3}, simulated from the model {:.3}", o.len(), sd(&o), sd(&per_read_sim));
    }
'''
s = s[:end] + tests + s[end:]
open(S, "w").write(s)
print("applied")
