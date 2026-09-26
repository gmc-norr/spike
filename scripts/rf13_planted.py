#!/usr/bin/env python3
"""RF13 plan: the proposed `ins_planted` row, in Python, run before any Rust.

The row asks whether the reads spike made for an insertion are in the BAM at the
truth record's position, carrying its bases -- whatever the aligner made of them.

- **Whose reads.** The truth record's ID is `sim_ins_N`; spike names that event's
  reads `evNNNN_hap_*` (`format!("ev{:04}", N)` in `simulate.rs`, and truth.rs
  numbers IDs the same way). No other read counts, so no read of the sample's own
  can carry, whatever it spells.
- **Which records.** Those `samtools view` returns over POS +/- 150, that are not
  secondary or supplementary. Any MAPQ, and duplicate, QC-fail or unmapped flags
  do not exclude: how the aligner placed or scored the read is not the question.
- **Carrying.** Two junction probes of 31 bases from the event's own haplotype
  `ref[..POS] + INS + ref[POS..]` (0-based POS = the truth record's POS, where
  spike inserts): `hap[POS-15, POS+16)` and `hap[POS+L-16, POS+L+15)`. A read
  carries if any 31-base window of its sequence is within 2 substitutions of
  either probe, in either orientation. Two, because spike writes the sample's own
  SNPs onto the event copy and sequencing errors onto every read (case file, RF6).
- **Verdict.** PASS at 1 or more carrying reads: with the sample's reads excluded
  there is no background to rise above. Not evaluable (FAIL) when the ALT has no
  bases or the ID is not `sim_ins_N`.

Usage:
  rf13_planted.py one BAM REFERENCE TRUTH_VCF       carriers per INS record
  rf13_planted.py k K_DIR SITES REFERENCE           the kill test (K+ and N1-N4)
  rf13_planted.py null BAM SITES REFERENCE          N5 on the `null` sites
"""
import os
import random
import re
import subprocess
import sys

PAD = 150
FLANK = 15
K = 31
MAX_MISMATCH = 2
ID_RE = re.compile(r"^sim_ins_(\d+)$")
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def fetch(ref, chrom, start0, end0):
    """Reference bases [start0, end0), 0-based half-open, uppercase."""
    out = subprocess.run(
        ["samtools", "faidx", ref, f"{chrom}:{start0 + 1}-{end0}"],
        capture_output=True, text=True, check=True,
    ).stdout
    return "".join(out.split("\n")[1:]).upper()


def probes(ref, chrom, pos, ins):
    L = len(ins)
    left = fetch(ref, chrom, pos - FLANK - 16, pos)
    right = fetch(ref, chrom, pos, pos + FLANK + 16)
    hap = left + ins + right  # hap[i] is 0-based pos - len(left) + i
    o = len(left)
    ps = [hap[o - FLANK:o + 16], hap[o + L - 16:o + L + FLANK]]
    return sorted({p for p in ps} | {revcomp(p) for p in ps})


def carries(seq, ps):
    seq = seq.upper()
    for p in ps:
        for i in range(len(seq) - K + 1):
            w = seq[i:i + K]
            mism = 0
            for a, b in zip(w, p):
                if a != b:
                    mism += 1
                    if mism > MAX_MISMATCH:
                        break
            if mism <= MAX_MISMATCH:
                return True
    return False


def carriers(bam, ref, chrom, pos, alt, rec_id, any_name=False):
    """Carrying read names, or a string saying why the row is not evaluable.

    `any_name` drops the name filter: only the check of the check uses it.
    """
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


def ins_records(truth):
    for line in open(truth):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if "SVTYPE=INS" in f[7]:
            yield f


def main(argv):
    mode = argv[0]
    if mode == "one":
        bam, ref, truth = argv[1:4]
        for f in ins_records(truth):
            c = carriers(bam, ref, f[0], int(f[1]), f[4], f[2])
            print(f"{f[0]}:{f[1]} {f[2]} -> {c if isinstance(c, str) else len(c)}")
        return
    if mode == "null":
        bam, sites, ref = argv[1:4]
        passed = n = 0
        for line in open(sites):
            name, chrom, pos, _ = line.split("\t")
            if name != "null":
                continue
            ins = "".join(random.Random(int(pos)).choice("ACGT") for _ in range(60))
            anchor = fetch(ref, chrom, int(pos) - 1, int(pos))
            c = carriers(bam, ref, chrom, int(pos), anchor + ins, "sim_ins_1")
            n += 1
            passed += not isinstance(c, str) and len(c) >= 1
        print(f"N5 null sites passing: {passed} of {n}")
        return
    # mode == "k": every run directory under K_DIR.
    kdir, sites, ref = argv[1:4]
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
        rng = random.Random(pos)
        wrong = alt[0] + "".join(rng.choice("ACGT") for _ in range(len(alt) - 1))
        n_id = int(ID_RE.match(rec_id).group(1)) if ID_RE.match(rec_id) else 0
        negs = {
            "N1 empty": (os.path.join(d, "slice.bam"), pos, alt, rec_id),
            "N2 wrong letters": (merged, pos, wrong, rec_id),
            "N3 shifted +1000": (merged, pos + 1000, alt, rec_id),
            "N4 wrong id": (merged, pos, alt, f"sim_ins_{n_id + 1}"),
        }
        if run.startswith("h_"):
            c = carriers(os.path.join(d, "slice.bam"), ref, chrom, pos, alt, rec_id, any_name=True)
            unfiltered.append((run, len(c)))
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
