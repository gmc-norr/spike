"""Throwaway Gate B, step 5 (on top of hp_all.patch = steps 1-4): zz_gate_clips writes FASTQ for the
clip test. SPIKE_GATE_SET=real writes the held-out pairs as sequenced; any other name makes fake reads
from the reference at each mate's own 5' position and strand with the model SPIKE_QX selects
(substitution errors only), trained on the train half. Output: $SPIKE_GATE_OUT/<set>_R1.fq, _R2.fq.
Run from the qm worktree: git apply ../qm2gate/hp_all.patch && python3 ../qm2gate/hp_apply5.py"""
S = "src/synth.rs"
s = open(S).read()
end = s.rstrip().rfind("}")
test = r'''
    /// Throwaway: a mate's template from the reference, 161 bases from its 5' end in sequencing order.
    fn zz_template(reference: &crate::reference::SharedReference, chrom: &str, m: &crate::types::MateAlignment) -> Vec<u8> {
        let lead: u64 = m.cigar.iter().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as u64).sum();
        let trail: u64 = m.cigar.iter().rev().take_while(|(op, _)| *op == b'S' || *op == b'H').filter(|(op, _)| *op == b'S').map(|(_, l)| *l as u64).sum();
        let ref_len: u64 = m.cigar.iter().filter(|(op, _)| *op == b'M' || *op == b'D').map(|(_, l)| *l as u64).sum();
        let mut t = if !m.reverse {
            let s = m.start.saturating_sub(lead);
            reference.fetch_sequence(chrom, s, s + 161).unwrap_or_default()
        } else {
            let e = m.start + ref_len + trail;
            let mut v = reference.fetch_sequence(chrom, e.saturating_sub(161), e).unwrap_or_default();
            crate::extract::reverse_complement(&mut v);
            v
        };
        t.iter_mut().for_each(|b| *b = b.to_ascii_uppercase());
        t.resize(161, b'N');
        t
    }

    /// Throwaway: one read from `template` as generate_from_template makes it, with the error history fed.
    fn zz_generate(profile: &QualityProfile, template: &[u8], rl: usize, read_num: u8, class: usize, rng: &mut StdRng) -> (Vec<u8>, Vec<u8>) {
        let mut state = profile.start_read(read_num, Some(class), rng);
        let (mut seq, mut qual) = (Vec::with_capacity(rl), Vec::with_capacity(rl));
        for c in 0..rl {
            let tb = template[c];
            let q = profile.next_quality(read_num, &state, rl - c, rng);
            if tb == b'N' {
                seq.push(b'N');
                qual.push(N_QUAL);
                profile.emitted(&mut state, N_QUAL, b'N');
                profile.record_error(&mut state, false);
                continue;
            }
            let err = rng.gen::<f64>() < profile.error_rate(q, &state, rl - c);
            seq.push(if err { random_different_base(tb, rng) } else { tb });
            qual.push(q);
            profile.emitted(&mut state, q, tb);
            profile.record_error(&mut state, err);
        }
        (seq, qual)
    }

    #[test]
    #[ignore]
    fn zz_gate_clips() {
        use std::io::Write;
        let bam = std::env::var("SPIKE_N7_BAM").unwrap();
        let fasta = std::env::var("SPIKE_REF").unwrap();
        let out = std::env::var("SPIKE_GATE_OUT").unwrap();
        let set = std::env::var("SPIKE_GATE_SET").unwrap();
        let chroms: Vec<String> = (1..=20).map(|i| format!("chr{}", i)).collect();
        let chrom_refs: Vec<&str> = chroms.iter().map(|s| s.as_str()).collect();
        let reference = crate::reference::SharedReference::load(&fasta, &chrom_refs).unwrap();
        let (train, held) = zz_gate_blocks(&bam);
        let held: Vec<ReadPair> = held
            .into_iter()
            .filter(|p| p.insert_size.abs() >= 151 && p.seq1.len() == 151 && p.seq2.len() == 151 && p.align.is_some())
            .collect();
        let mut f1 = std::io::BufWriter::new(std::fs::File::create(format!("{}/{}_R1.fq", out, set)).unwrap());
        let mut f2 = std::io::BufWriter::new(std::fs::File::create(format!("{}/{}_R2.fq", out, set)).unwrap());
        let mut put = |f: &mut std::io::BufWriter<std::fs::File>, name: &str, s: &[u8], q: &[u8]| {
            writeln!(f, "@{}\n{}\n+\n{}", name, String::from_utf8_lossy(s), String::from_utf8_lossy(q)).unwrap();
        };
        if set == "real" {
            for p in &held {
                put(&mut f1, &p.name, &p.seq1, &p.qual1);
                put(&mut f2, &p.name, &p.seq2, &p.qual2);
            }
        } else {
            let profile = QualityProfile::from_donor_pairs(&train, 151, &reference);
            let mut rng = StdRng::seed_from_u64(21);
            for p in &held {
                let al = p.align.as_ref().unwrap();
                let t1 = zz_template(&reference, &p.chrom, &al[0]);
                let t2 = zz_template(&reference, &p.chrom, &al[1]);
                let (c1, c2) = profile.draw_classes_given(crate::quality::read_hp(&t1[..151]), crate::quality::read_hp(&t2[..151]), &mut rng);
                let (s1, q1) = zz_generate(&profile, &t1, 151, 1, c1, &mut rng);
                let (s2, q2) = zz_generate(&profile, &t2, 151, 2, c2, &mut rng);
                put(&mut f1, &p.name, &s1, &q1);
                put(&mut f2, &p.name, &s2, &q2);
            }
        }
        println!("set {} pairs {} (train {})", set, held.len(), train.len());
    }
'''
s = s[:end] + test + s[end:]
open(S, "w").write(s)
print("applied step 5")
