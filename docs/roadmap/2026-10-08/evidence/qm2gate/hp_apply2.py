"""Throwaway Gate B, step 2 (on top of hp_apply.py), more SPIKE_QX bits:
  16 the pair's two read classes drawn given each mate's longest one-letter run (the class is cut on
     the read's mean quality, which already holds the crash a run set off)
  32 the error table also keyed by the errors in the read's last 30 bases (0, 1, 2-3, 4+)
zz_gate_errors now also SIMULATES errors along the real qualities (errors feeding the history) and
compares the simulated tail rates and per-read spread with the observed ones.
Run from the qm worktree after hp_apply.py: python3 ../qm2gate/hp_apply2.py"""

def edit(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        assert s.count(old) == 1, (path, old[:80])
        s = s.replace(old, new)
    open(path, "w").write(s)

Q = "src/quality.rs"
edit(Q, [
("""const HP_BINS: usize = 4;""", """const HP_BINS: usize = 4;
/// Errors among the read's last 30 bases: 0, 1, 2-3, 4+.
const ERR_BINS: usize = 4;
fn err_bin(hist: u32) -> usize {
    match (hist & ((1 << 30) - 1)).count_ones() {
        0 => 0,
        1 => 1,
        2..=3 => 2,
        _ => 3,
    }
}"""),
("""    run_len: u32,
    hp: u8,
}""", """    run_len: u32,
    hp: u8,
    /// Which of the last 30 bases were errors (bit 0 the last).
    err_hist: u32,
}"""),
("""run_base: b'N', run_len: 0, hp: 0 }""", """run_base: b'N', run_len: 0, hp: 0, err_hist: 0 }"""),
("""    errors_mid: Vec<ErrorCell>,
    errors_class""", """    errors_mid: Vec<ErrorCell>,
    /// Per mate, (longest-run bin, class) counts over donor reads.
    class_hp: Vec<u64>,
    errors_class"""),
("""            errors_mid: Vec::new(),
            errors_class: Vec::new(),
            errors_symbol: Vec::new(),
            n_pairs: pairs.len(),""", """            errors_mid: Vec::new(),
            class_hp: vec![0; 2 * HP_BINS * READ_CLASSES],
            errors_class: Vec::new(),
            errors_symbol: Vec::new(),
            n_pairs: pairs.len(),"""),
("""            errors_mid: Vec::new(),
            errors_class: Vec::new(),
            errors_symbol: Vec::new(),
            n_pairs: self.n_pairs,""", """            errors_mid: Vec::new(),
            class_hp: self.class_hp.clone(),
            errors_class: Vec::new(),
            errors_symbol: Vec::new(),
            n_pairs: self.n_pairs,"""),
("""            profile.joint[c1 * READ_CLASSES + c2] += 1;
        }""", """            profile.joint[c1 * READ_CLASSES + c2] += 1;
            for (m, (seq, c)) in [(&p.seq1, c1), (&p.seq2, c2)].into_iter().enumerate() {
                profile.class_hp[(m * HP_BINS + read_hp(seq)) * READ_CLASSES + c] += 1;
            }
        }"""),
("""    /// A fresh read of class `class`,""", """    /// Throwaway: the pair's classes given each mate's longest-run bin (SPIKE_QX & 16),
    /// reweighting the joint by P(class | run bin) / P(class) per mate.
    pub fn draw_classes_given(&self, h1: usize, h2: usize, rng: &mut StdRng) -> (usize, usize) {
        if qx() & 16 == 0 {
            return self.draw_classes(rng);
        }
        let ratio = |m: usize, h: usize, c: usize| -> f64 {
            let at = |hh: usize, cc: usize| self.class_hp[(m * HP_BINS + hh) * READ_CLASSES + cc] as f64;
            let n_hc = at(h, c);
            let n_h: f64 = (0..READ_CLASSES).map(|cc| at(h, cc)).sum();
            let n_c: f64 = (0..HP_BINS).map(|hh| at(hh, c)).sum();
            let n: f64 = (0..HP_BINS).flat_map(|hh| (0..READ_CLASSES).map(move |cc| (hh, cc))).map(|(hh, cc)| at(hh, cc)).sum();
            ((n_hc + 0.5) / (n_h + 4.0)) / ((n_c + 0.5) / (n + 4.0))
        };
        let w: Vec<f64> = (0..READ_CLASSES * READ_CLASSES)
            .map(|i| self.joint[i] as f64 * ratio(0, h1, i / READ_CLASSES) * ratio(1, h2, i % READ_CLASSES))
            .collect();
        let total: f64 = w.iter().sum();
        if total <= 0.0 {
            return self.draw_classes(rng);
        }
        let mut r = rng.gen::<f64>() * total;
        for (i, &x) in w.iter().enumerate() {
            if r < x {
                return (i / READ_CLASSES, i % READ_CLASSES);
            }
            r -= x;
        }
        (READ_CLASSES / 2, READ_CLASSES / 2)
    }

    /// Throwaway: record whether the base just emitted was an error.
    pub fn record_error(&self, state: &mut ReadState, err: bool) {
        state.err_hist = (state.err_hist << 1) | err as u32;
    }

    /// Throwaway: walk a mate's real qualities and bases with errors feeding the history: the
    /// real ones (`rng` None; uncounted bases as no error) or simulated ones. Per base, the rate
    /// the model gave and the error (observed, or drawn).
    pub fn walk_errors(&self, qual: &[u8], seq: &[u8], verdict: &[Option<bool>], mut rng: Option<&mut StdRng>) -> Vec<(f64, bool)> {
        let mut state = ReadState::new(self.class_of(qual));
        let len = qual.len();
        let mut out = Vec::with_capacity(len);
        for (c, &q) in qual.iter().enumerate() {
            let p = self.error_rate(q, &state, len - c);
            let e = match rng.as_mut() {
                Some(r) => r.gen::<f64>() < p,
                None => verdict.get(c).copied().flatten().unwrap_or(false),
            };
            out.push((p, e));
            self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
            self.record_error(&mut state, e);
        }
        out
    }

    /// A fresh read of class `class`,"""),
("""        let full = self.errors_full[self.error_cell(symbol, state.class, rb, end_bin(cycles_left), hp)];
        let mid = if qx() & 12 != 0 { self.errors_mid[(symbol * RUN_BINS + rb) * HP_BINS + hp] } else { ErrorCell::default() };""",
"""        let eb = if qx() & 32 != 0 { err_bin(state.err_hist) } else { 0 };
        let full = self.errors_full[self.error_cell(symbol, state.class, rb, end_bin(cycles_left), hp) * ERR_BINS + eb];
        let mid = if qx() & 44 != 0 { self.errors_mid[((symbol * RUN_BINS + rb) * HP_BINS + hp) * ERR_BINS + eb] } else { ErrorCell::default() };"""),
("""        self.errors_full = vec![ErrorCell::default(); a * READ_CLASSES * RUN_BINS * END_BINS * HP_BINS];
        self.errors_mid = vec![ErrorCell::default(); a * RUN_BINS * HP_BINS];""",
"""        self.errors_full = vec![ErrorCell::default(); a * READ_CLASSES * RUN_BINS * END_BINS * HP_BINS * ERR_BINS];
        self.errors_mid = vec![ErrorCell::default(); a * RUN_BINS * HP_BINS * ERR_BINS];"""),
("""                        let i = self.error_cell(symbol, state.class, rb, end_bin(len - c), hp);
                        self.errors_full[i].add(cell);
                        self.errors_mid[(symbol * RUN_BINS + rb) * HP_BINS + hp].add(cell);""",
"""                        let eb = if qx() & 32 != 0 { err_bin(state.err_hist) } else { 0 };
                        let i = self.error_cell(symbol, state.class, rb, end_bin(len - c), hp) * ERR_BINS + eb;
                        self.errors_full[i].add(cell);
                        self.errors_mid[((symbol * RUN_BINS + rb) * HP_BINS + hp) * ERR_BINS + eb].add(cell);"""),
("""                    self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
                }
                if let Some(v) = verdicts_out.as_mut() {""", """                    self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
                    let e = matches!(verdict.get(c), Some(Some(true)));
                    self.record_error(&mut state, e);
                }
                if let Some(v) = verdicts_out.as_mut() {"""),
("""fn mean_phred(qual: &[u8]) -> f64 {""", """/// Throwaway: the largest one-letter-run bin over a read's bases.
pub fn read_hp(seq: &[u8]) -> usize {
    let (mut best, mut base, mut len) = (0u8, b'N', 0u32);
    for &b in seq {
        let b = b.to_ascii_uppercase();
        if b == base && b != b'N' {
            len += 1;
        } else {
            base = b;
            len = (b != b'N') as u32;
        }
        best = best.max(hp_bin(base, len));
    }
    best as usize
}

fn mean_phred(qual: &[u8]) -> f64 {"""),
])

S = "src/synth.rs"
s = open(S).read()
old = """                let (c1, c2) = profile.draw_classes(rng);
                (draw(profile, 1, c1, &p.seq1, rng), draw(profile, 2, c2, &p.seq2, rng))"""
assert s.count(old) == 1
s = s.replace(old, """                let (c1, c2) = profile.draw_classes_given(crate::quality::read_hp(&p.seq1), crate::quality::read_hp(&p.seq2), rng);
                (draw(profile, 1, c1, &p.seq1, rng), draw(profile, 2, c2, &p.seq2, rng))""")
start = s.index("    #[test]\n    #[ignore]\n    fn zz_gate_errors() {")
end = s.rstrip().rfind("}")
s = s[:start] + r'''    #[test]
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
        // groups: all, crashed tails Q<15, Q15-29, Q30+, not crashed last 40
        let names = ["all counted bases", "crashed tails, Q<15", "crashed tails, Q15-29", "crashed tails, Q30+", "not crashed, last 40"];
        let mut obs = [0f64; 5];
        let mut pred = [0f64; 5];
        let mut sim = [0f64; 5];
        let mut n = [0f64; 5];
        let (mut ro_all, mut rs_all) = (Vec::new(), Vec::new());
        let mut rng = StdRng::seed_from_u64(11);
        for ((q, seq), v) in mates.iter().zip(&verdicts) {
            let len = q.len();
            let crashed = len >= 20 && q[len - 20..].iter().filter(|&&x| x < b'!' + 15).count() >= 10;
            let real = profile.walk_errors(q, seq, v, None);
            let simulated = profile.walk_errors(q, seq, v, Some(&mut rng));
            let (mut ro, mut rs, mut rn) = (0f64, 0f64, 0f64);
            for c in 0..len {
                let Some(Some(err)) = v.get(c) else { continue };
                let e = *err as u8 as f64;
                let phred = q[c] - 33;
                let tail = c + 40 >= len;
                let groups: Vec<usize> = [Some(0), (crashed && tail).then(|| match phred { 0..=14 => 1, 15..=29 => 2, _ => 3 }), (!crashed && tail).then_some(4)]
                    .into_iter().flatten().collect();
                for g in groups {
                    obs[g] += e; pred[g] += real[c].0; sim[g] += simulated[c].1 as u8 as f64; n[g] += 1.0;
                }
                if crashed && tail && phred < 15 {
                    ro += e; rs += simulated[c].1 as u8 as f64; rn += 1.0;
                }
            }
            if rn >= 10.0 { ro_all.push(ro / rn); rs_all.push(rs / rn); }
        }
        println!("QX {}: held-out bases; observed errors, the model's rate given the real history, and errors simulated along the real qualities", std::env::var("SPIKE_QX").unwrap_or_default());
        for g in 0..5 {
            println!("  {:24} bases {:9}  observed {:.4}  predicted {:.4} ({:.2})  simulated {:.4} ({:.2})", names[g], n[g], obs[g] / n[g], pred[g] / n[g], pred[g] / obs[g], sim[g] / n[g], sim[g] / obs[g]);
        }
        let sd = |v: &[f64]| { let m = v.iter().sum::<f64>() / v.len() as f64; (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / v.len() as f64).sqrt() };
        println!("  crashed reads (>= 10 Q<15 tail bases): {}; per-read tail rate SD observed {:.3}, simulated {:.3}", ro_all.len(), sd(&ro_all), sd(&rs_all));
    }
''' + s[end:]
open(S, "w").write(s)
print("applied step 2")
