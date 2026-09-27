#!/usr/bin/env python3
"""RF13 second attempt: `ins_planted` with the inserted bases required exactly.

The same row as rf13_planted.py -- the reads spike made for the event
(`evNNNN_hap_*` for truth ID `sim_ins_N`), over POS +/- 150, not secondary or
supplementary, any MAPQ or other flag; two 31-base junction probes from
`ref[..POS] + INS + ref[POS..]`; PASS at 1 or more carrying reads; not evaluable
(FAIL) when the ALT has no bases or the ID is not `sim_ins_N` -- with one change:

- **Carrying.** A 31-base window of a read carries a probe when every position of
  the probe that holds an **inserted** base matches exactly, and at most 2
  positions holding a **reference flank** base differ. The first attempt allowed
  2 substitutions anywhere, and a truth with the wrong 4 bases matched (RF13 in
  docs/review/REVIEW.md). The flank tolerance stays for the sample's own SNPs, which spike
  writes onto reference positions only (`apply_variants` skips a segment with no
  reference origin, `src/haplotype.rs:612`), and for sequencing errors.

And one change to the test: the wrong-letters control replaces **every** inserted
base with a different base, so a "wrong" truth can never equal the right one.

Usage:
  rf13b_planted.py one BAM REFERENCE TRUTH_VCF
  rf13b_planted.py k K_DIR REFERENCE
  rf13b_planted.py null BAM SITES REFERENCE
"""
import os
import random
import re
import subprocess
import sys

PAD = 150
FLANK = 15
K = 31
MAX_FLANK_MISMATCH = 2
ID_RE = re.compile(r"^sim_ins_(\d+)$")
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def fetch(ref, chrom, start0, end0):
    out = subprocess.run(
        ["samtools", "faidx", ref, f"{chrom}:{start0 + 1}-{end0}"],
        capture_output=True, text=True, check=True,
    ).stdout
    return "".join(out.split("\n")[1:]).upper()


def probes(ref, chrom, pos, ins):
    """(probe, is_inserted mask) pairs, both orientations."""
    L = len(ins)
    left = fetch(ref, chrom, pos - FLANK - 16, pos)
    right = fetch(ref, chrom, pos, pos + FLANK + 16)
    hap = left + ins + right
    mask = [False] * len(left) + [True] * L + [False] * len(right)
    o = len(left)
    out = []
    for a, b in ((o - FLANK, o + 16), (o + L - 16, o + L + FLANK)):
        p, m = hap[a:b], mask[a:b]
        out.append((p, m))
        out.append((revcomp(p), m[::-1]))
    uniq = []
    for pm in out:
        if pm not in uniq:
            uniq.append(pm)
    return uniq


def carries(seq, ps):
    seq = seq.upper()
    for p, m in ps:
        for i in range(len(seq) - K + 1):
            flank_mism = 0
            for j in range(K):
                if seq[i + j] != p[j]:
                    if m[j]:
                        break
                    flank_mism += 1
                    if flank_mism > MAX_FLANK_MISMATCH:
                        break
            else:
                return True
    return False


def carriers(bam, ref, chrom, pos, alt, rec_id, any_name=False):
    m = ID_RE.match(rec_id)
    if not m:
        return f"id {rec_id} is not sim_ins_N"
    if len(alt) <= 1 or alt.startswith("<"):
        return "alt has no bases"
    prefix = f"ev{int(m.group(1)):04}_hap_"
    ps = probes(ref, chrom, pos, alt[1:].upper())
    out = subprocess.run(
        ["samtools", "view", bam, f"{chrom}:{max(0, pos - PAD) + 1}-{pos + PAD}"],
        capture_output=True, text=True, check=True,
    ).stdout
    names = set()
    for line in out.splitlines():
        f = line.split("\t", 10)
        if (not any_name and not f[0].startswith(prefix)) or int(f[1]) & (0x100 | 0x800):
            continue
        if carries(f[9], ps):
            names.add(f[0])
    return names


def wrong_letters(alt, seed):
    """The ALT with every inserted base replaced by a different base."""
    rng = random.Random(seed)
    return alt[0] + "".join(rng.choice([b for b in "ACGT" if b != c]) for c in alt[1:].upper())


def ins_records(truth):
    for line in open(truth):
        if not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            if "SVTYPE=INS" in f[7]:
                yield f


def main(argv):
    mode = argv[0]
    if mode == "one":
        bam, ref, truth = argv[1:4]
        for f in ins_records(truth):
            c = carriers(bam, ref, f[0], int(f[1]), f[4], f[2])
            w = carriers(bam, ref, f[0], int(f[1]), wrong_letters(f[4], int(f[1])), f[2])
            print(f"{f[0]}:{f[1]} {f[2]} len {len(f[4]) - 1} -> "
                  f"{c if isinstance(c, str) else len(c)}; wrong letters -> "
                  f"{w if isinstance(w, str) else len(w)}")
        return
    if mode == "null":
        bam, sites, ref = argv[1:4]
        passed = n = 0
        for line in open(sites):
            name, chrom, pos, _ = line.rstrip("\n").split("\t")
            if name != "null":
                continue
            ins = "".join(random.Random(int(pos)).choice("ACGT") for _ in range(60))
            anchor = fetch(ref, chrom, int(pos) - 1, int(pos))
            c = carriers(bam, ref, chrom, int(pos), anchor + ins, "sim_ins_1")
            n += 1
            passed += not isinstance(c, str) and len(c) >= 1
        print(f"N5 null sites passing: {passed} of {n}")
        return
    kdir, ref = argv[1:3]
    pos_rows, neg_fail, unfiltered = [], [], []
    for run in sorted(os.listdir(kdir)):
        d = os.path.join(kdir, run)
        truth = os.path.join(d, "run", "truth.vcf")
        merged = os.path.join(d, "run", "merged.bam")
        if not (os.path.exists(truth) and os.path.exists(merged)):
            pos_rows.append((run, "did not reach validate"))
            continue
        f = next(ins_records(truth))
        chrom, pos, rec_id, alt = f[0], int(f[1]), f[2], f[4]
        c = carriers(merged, ref, chrom, pos, alt, rec_id)
        pos_rows.append((run, c if isinstance(c, str) else len(c)))
        n_id = int(ID_RE.match(rec_id).group(1)) if ID_RE.match(rec_id) else 0
        negs = {
            "N1 empty": (os.path.join(d, "slice.bam"), pos, alt, rec_id),
            "N2 wrong letters": (merged, pos, wrong_letters(alt, pos), rec_id),
            "N3 shifted +1000": (merged, pos + 1000, alt, rec_id),
            "N4 wrong id": (merged, pos, alt, f"sim_ins_{n_id + 1}"),
        }
        if run.startswith("h_"):
            u = carriers(os.path.join(d, "slice.bam"), ref, chrom, pos, alt, rec_id, any_name=True)
            unfiltered.append((run, len(u)))
        for label, (bam, p, a, i) in negs.items():
            c = carriers(bam, ref, chrom, p, a, i)
            n = 0 if isinstance(c, str) else len(c)
            if n >= 1:
                neg_fail.append((run, label, n))
    for run, c in pos_rows:
        print(f"  {run}: {c}")
    ran = [c for _, c in pos_rows if isinstance(c, int)]
    ok = sum(1 for c in ran if c >= 1)
    print(f"K+ runs reaching validate: {len(ran)} of {len(pos_rows)}")
    print(f"K+ carriers >= 1: {ok} of {len(ran)}  (min {min(ran) if ran else '-'}, "
          f"median {sorted(ran)[len(ran) // 2] if ran else '-'})")
    print(f"K- negatives with a carrier: {len(neg_fail)} of {4 * len(ran)} {neg_fail}")
    print(f"check of the check, hg sites' unspiked slice with the name filter off: {unfiltered}")


if __name__ == "__main__":
    main(sys.argv[1:])
