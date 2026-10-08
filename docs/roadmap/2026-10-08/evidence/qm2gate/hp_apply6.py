"""Throwaway Gate B, step 6 (on top of clip_all.patch): the two clip-test fixes.
  A (always on): with a reference, the run bin is learned from the reference under each donor mate
    (from its 5' end, sequencing order) -- quality context, class draw and error table -- not from
    its called bases.
  B (SPIKE_QX bit 128): slips. learn_slips counts, per base and run-length bin, how often a donor
    mate's alignment has an indel inside a template run (runs within 3 bp of a masked variant, or not
    fully aligned with a margin, left out); zz_generate reads a run one base short or long at that rate.
Run from the qm worktree: git apply ../qm2gate/clip_all.patch && python3 ../qm2gate/hp_apply6.py"""

def edit(path, pairs):
    s = open(path).read()
    for old, new in pairs:
        assert s.count(old) == 1, (path, old[:90])
        s = s.replace(old, new)
    open(path, "w").write(s)

Q = "src/quality.rs"
edit(Q, [
# --- fix A: templates in learn()
("""        for p in pairs {
            let (c1, c2) = (profile.class_of(&p.qual1), profile.class_of(&p.qual2));
            profile.joint[c1 * READ_CLASSES + c2] += 1;
            for (m, (seq, c)) in [(&p.seq1, c1), (&p.seq2, c2)].into_iter().enumerate() {
                profile.class_hp[(m * HP_BINS + read_hp(seq)) * READ_CLASSES + c] += 1;
            }
        }""", """        // Throwaway fix A: the bases the run bin is read from -- the reference under each mate when
        // there is one, else the called bases.
        let templ: Vec<[Vec<u8>; 2]> = pairs
            .iter()
            .map(|p| match (reference, &p.align) {
                (Some(r), Some(al)) => [
                    template_of(r, &p.chrom, &al[0], p.qual1.len()),
                    template_of(r, &p.chrom, &al[1], p.qual2.len()),
                ],
                _ => [p.seq1.clone(), p.seq2.clone()],
            })
            .collect();
        for (p, t) in pairs.iter().zip(&templ) {
            let (c1, c2) = (profile.class_of(&p.qual1), profile.class_of(&p.qual2));
            profile.joint[c1 * READ_CLASSES + c2] += 1;
            for (m, (seq, c)) in [(&t[0], c1), (&t[1], c2)].into_iter().enumerate() {
                profile.class_hp[(m * HP_BINS + read_hp(seq)) * READ_CLASSES + c] += 1;
            }
        }"""),
("""            let table = pairs
                .par_chunks(2048)
                .map(|chunk| {
                    let mut t: HashMap<u64, Counts> = HashMap::new();
                    for p in chunk {
                        let (qual, seq) = if mate == 0 { (&p.qual1, &p.seq1) } else { (&p.qual2, &p.seq2) };
                        profile.count_read(qual, seq, &mut t);
                    }
                    t
                })""", """            let table = pairs
                .par_chunks(2048)
                .zip(templ.par_chunks(2048))
                .map(|(chunk, tchunk)| {
                    let mut t: HashMap<u64, Counts> = HashMap::new();
                    for (p, tm) in chunk.iter().zip(tchunk) {
                        let qual = if mate == 0 { &p.qual1 } else { &p.qual2 };
                        profile.count_read(qual, &tm[mate], &mut t);
                    }
                    t
                })"""),
("""                let mut state = ReadState::new(self.class_of(qual));
                let len = qual.len();
                for (c, &q) in qual.iter().enumerate() {
                    if let Some(Some(err)) = verdict.get(c) {""", """                let mut state = ReadState::new(self.class_of(qual));
                let len = qual.len();
                let tmpl = template_of(reference, &p.chrom, m, len);
                for (c, &q) in qual.iter().enumerate() {
                    if let Some(Some(err)) = verdict.get(c) {"""),
("""                    self.push(&mut state, q, seq.get(c).copied().unwrap_or(b'N'));
                    let e = matches!(verdict.get(c), Some(Some(true)));""", """                    self.push(&mut state, q, tmpl[c]);
                    let e = matches!(verdict.get(c), Some(Some(true)));"""),
# --- fix B: slip table
("""            class_hp: vec![0; 2 * HP_BINS * READ_CLASSES],""", """            class_hp: vec![0; 2 * HP_BINS * READ_CLASSES],
            slips: Vec::new(),"""),
("""            class_hp: self.class_hp.clone(),""", """            class_hp: self.class_hp.clone(),
            slips: self.slips.clone(),"""),
("""    /// Per mate, (longest-run bin, class) counts over donor reads.
    class_hp: Vec<u64>,""", """    /// Per mate, (longest-run bin, class) counts over donor reads.
    class_hp: Vec<u64>,
    /// Throwaway fix B: per (base, run-length bin): slips, opportunities, deletions.
    slips: Vec<[u64; 3]>,"""),
("""    /// Throwaway: record whether the base just emitted was an error.""", """    /// Throwaway fix B: count slips in donor mates' alignments (see the module script).
    pub fn learn_slips(&mut self, pairs: &[ReadPair], reference: &SharedReference, mask: &HashMap<String, Vec<(i64, i64)>>) -> String {
        self.slips = vec![[0; 3]; 4 * SLIP_BINS];
        let mut masked = 0u64;
        let is_masked = |chrom: &str, lo: i64, hi: i64| -> bool {
            let Some(v) = mask.get(chrom) else { return false };
            let k = v.partition_point(|&(s, _)| s <= hi);
            v[..k].iter().rev().take_while(|&&(s, _)| s >= lo - 1000).any(|&(_, e)| e >= lo)
        };
        for p in pairs {
            let Some(al) = &p.align else { continue };
            for (m, qual) in al.iter().zip([&p.qual1, &p.qual2]) {
                let n = qual.len();
                let t = template_of(reference, &p.chrom, m, n);
                let lead: i64 = m.cigar.iter().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as i64).sum();
                let trail: i64 = m.cigar.iter().rev().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as i64).sum();
                let ref_len: i64 = m.cigar.iter().filter(|(op, _)| *op == b'M' || *op == b'D').map(|(_, l)| *l as i64).sum();
                let start = m.start as i64;
                let refpos = |i: usize| if !m.reverse { start - lead + i as i64 } else { start + ref_len + trail - 1 - i as i64 };
                let (mut dels, mut inss) = (Vec::new(), Vec::new());
                let mut rp = start;
                for &(op, l) in &m.cigar {
                    match op {
                        b'M' => rp += l as i64,
                        b'D' => {
                            dels.push((rp, rp + l as i64));
                            rp += l as i64;
                        }
                        b'I' => inss.push(rp),
                        _ => {}
                    }
                }
                let mut i = 0;
                while i < n {
                    let b = t[i];
                    let mut j = i;
                    while j + 1 < n && t[j + 1] == b {
                        j += 1;
                    }
                    let len = j - i + 1;
                    if let (Some(bi), Some(bin)) = (base_index(b), slip_bin(len)) {
                        let (x, y) = (refpos(i), refpos(j));
                        let (a, z) = (x.min(y), x.max(y));
                        if i >= 5 && j + 5 < n && a > start && z + 1 < start + ref_len {
                            if is_masked(&p.chrom, a - 3, z + 3) {
                                masked += 1;
                            } else {
                                let del = dels.iter().any(|&(s, e)| s <= z && e > a);
                                let ins = inss.iter().any(|&q| q >= a && q <= z + 1);
                                let c = &mut self.slips[bi * SLIP_BINS + bin];
                                c[1] += 1;
                                if del || ins {
                                    c[0] += 1;
                                    c[2] += del as u64;
                                }
                            }
                        }
                    }
                    i = j + 1;
                }
            }
        }
        let names = ["3-4", "5-6", "7-8", "9-11", "12-14", "15+"];
        let mut out = format!("slips (rate, runs, deletion share) per base; {} runs masked\\n", masked);
        for (bi, b) in "ACGT".chars().enumerate() {
            out += &format!("  {}:", b);
            for (k, nm) in names.iter().enumerate() {
                let c = self.slips[bi * SLIP_BINS + k];
                out += &format!("  {} {:.4} ({}, del {:.2})", nm, c[0] as f64 / c[1].max(1) as f64, c[1], c[2] as f64 / c[0].max(1) as f64);
            }
            out += "\\n";
        }
        out
    }

    /// Throwaway fix B: (slip rate, deletion share) for a template run, when SPIKE_QX & 128.
    pub fn slip_rate(&self, base: u8, len: usize) -> Option<(f64, f64)> {
        if qx() & 128 == 0 || self.slips.is_empty() {
            return None;
        }
        let (bi, bin) = (base_index(base)?, slip_bin(len)?);
        let mut c = self.slips[bi * SLIP_BINS + bin];
        if c[1] < 50 {
            c = (0..4).fold([0; 3], |a, k| {
                let x = self.slips[k * SLIP_BINS + bin];
                [a[0] + x[0], a[1] + x[1], a[2] + x[2]]
            });
        }
        if c[1] == 0 {
            return None;
        }
        Some((c[0] as f64 / c[1] as f64, if c[0] > 0 { c[2] as f64 / c[0] as f64 } else { 0.5 }))
    }

    /// Throwaway: record whether the base just emitted was an error."""),
("""/// Throwaway: the largest one-letter-run bin over a read's bases.""", """const SLIP_BINS: usize = 6;
fn slip_bin(len: usize) -> Option<usize> {
    match len {
        3..=4 => Some(0),
        5..=6 => Some(1),
        7..=8 => Some(2),
        9..=11 => Some(3),
        12..=14 => Some(4),
        15.. => Some(5),
        _ => None,
    }
}
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Throwaway: `n` reference bases from mate `m`'s 5' end, in sequencing order, uppercase, N-padded.
pub fn template_of(reference: &SharedReference, chrom: &str, m: &MateAlignment, n: usize) -> Vec<u8> {
    let lead: u64 = m.cigar.iter().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as u64).sum();
    let trail: u64 = m.cigar.iter().rev().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as u64).sum();
    let ref_len: u64 = m.cigar.iter().filter(|(op, _)| *op == b'M' || *op == b'D').map(|(_, l)| *l as u64).sum();
    let mut t = if !m.reverse {
        let s = m.start.saturating_sub(lead);
        reference.fetch_sequence(chrom, s, s + n as u64).unwrap_or_default()
    } else {
        let e = m.start + ref_len + trail;
        let mut v = reference.fetch_sequence(chrom, e.saturating_sub(n as u64), e).unwrap_or_default();
        crate::extract::reverse_complement(&mut v);
        v
    };
    t.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
    t.resize(n, b'N');
    t
}

/// Throwaway: the largest one-letter-run bin over a read's bases."""),
])

S = "src/synth.rs"
edit(S, [
("""    fn zz_generate(profile: &QualityProfile, template: &[u8], rl: usize, read_num: u8, class: usize, rng: &mut StdRng) -> (Vec<u8>, Vec<u8>) {
        let mut state = profile.start_read(read_num, Some(class), rng);
        let (mut seq, mut qual) = (Vec::with_capacity(rl), Vec::with_capacity(rl));
        for c in 0..rl {
            let tb = template[c];""", """    fn zz_generate(profile: &QualityProfile, template: &[u8], rl: usize, read_num: u8, class: usize, rng: &mut StdRng) -> (Vec<u8>, Vec<u8>) {
        // fix B: each run of one base is read one short (skip its last base) or one long (repeat it)
        let n = template.len();
        let (mut skip, mut extra) = (vec![false; n], vec![false; n]);
        let mut i = 0;
        while i < n {
            let mut j = i;
            while j + 1 < n && template[j + 1] == template[i] {
                j += 1;
            }
            if let Some((p, pdel)) = profile.slip_rate(template[i], j - i + 1) {
                if rng.gen::<f64>() < p {
                    if rng.gen::<f64>() < pdel { skip[j] = true } else { extra[j] = true }
                }
            }
            i = j + 1;
        }
        let mut state = profile.start_read(read_num, Some(class), rng);
        let (mut seq, mut qual) = (Vec::with_capacity(rl), Vec::with_capacity(rl));
        let mut idx = 0usize;
        while seq.len() < rl && idx < n {
            if skip[idx] {
                idx += 1;
                continue;
            }
            let c = seq.len();
            let tb = template[idx];
            if extra[idx] {
                extra[idx] = false;
            } else {
                idx += 1;
            }"""),
("""            profile.emitted(&mut state, q, tb);
            profile.record_error(&mut state, err);
        }
        (seq, qual)
    }""", """            profile.emitted(&mut state, q, tb);
            profile.record_error(&mut state, err);
        }
        seq.resize(rl, b'N');
        qual.resize(rl, N_QUAL);
        (seq, qual)
    }"""),
("""            let profile = QualityProfile::from_donor_pairs(&train, 151, &reference);
            let mut rng = StdRng::seed_from_u64(21);""", """            let mut profile = QualityProfile::from_donor_pairs(&train, 151, &reference);
            if let Ok(path) = std::env::var("SPIKE_GATE_MASK") {
                let mut mask: std::collections::HashMap<String, Vec<(i64, i64)>> = std::collections::HashMap::new();
                for line in std::fs::read_to_string(path).unwrap().lines() {
                    let f: Vec<&str> = line.split('\\t').collect();
                    mask.entry(f[0].to_string()).or_default().push((f[1].parse().unwrap(), f[2].parse().unwrap()));
                }
                mask.values_mut().for_each(|v| v.sort());
                print!("{}", profile.learn_slips(&train, &reference, &mask));
            }
            let mut rng = StdRng::seed_from_u64(21);"""),
])
print("applied step 6")
