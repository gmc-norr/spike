#!/usr/bin/env python3
"""RF6 kill test K: does a junction k-mer find every correct deletion, and only those?

For each run under RUN_DIR (T2's C4 layout: RUN_DIR/events.txt and RUN_DIR/<n>/ with
slice.bam and run/{sim.bam,truth.vcf}), the junction is read off the run's own
truth.vcf the way validate reads a DEL (start = POS, end = END, both 0-based as
TruthEvent holds them):

    J = ref[start - 15, start) + ref[end, end + 16)      (31 bases)

A read carries the junction if its bases hold J or its reverse complement. Reads are
the distinct names of primary, mapped, non-duplicate, non-QC-fail, non-supplementary
records at MAPQ >= 20 overlapping 500 bp either side of either breakpoint --
validate's usable_alignment filter and split_reads' window.

    K1  run/sim.bam (the spiked reads):         >= 2 reads   on 40 of 40
    K2  slice.bam (the donor, never spiked):     0 reads     on 40 of 40
    K3  run/sim.bam, J built with END + 50:      0 reads     on 40 of 40
    K4  J or its reverse complement in the reference within 1000 bp of either
        breakpoint:                                          on 0 of 40

Every samtools call raises on a non-zero exit: a helper that returns stdout unchecked
will eventually return silence and be summed (.claude/judgment-gate-cases.md).

Usage: rf6_kill.py RUN_DIR REFERENCE
"""
import pathlib
import subprocess
import sys

LEFT, RIGHT, PAD, REF_PAD, MIN_MAPQ, FLOOR = 15, 16, 500, 1000, 20, 2
RUN, REF = pathlib.Path(sys.argv[1]), sys.argv[2]
COMP = str.maketrans("ACGTN", "TGCAN")


def sh(args):
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def fetch(chrom, start0, end0):
    """Reference bases at 0-based [start0, end0), upper case."""
    out = sh(["samtools", "faidx", REF, f"{chrom}:{start0 + 1}-{end0}"])
    return "".join(out.splitlines()[1:]).upper()


def revcomp(s):
    return s.translate(COMP)[::-1]


def carriers(bam, chrom, breakpoints, probes):
    names = set()
    for bp in breakpoints:
        region = f"{chrom}:{max(1, bp - PAD + 1)}-{bp + PAD}"
        out = sh(["samtools", "view", "-F", "0xF04", "-q", str(MIN_MAPQ), str(bam), region])
        for line in out.splitlines():
            fields = line.split("\t")
            if any(p in fields[9].upper() for p in probes):
                names.add(fields[0])
    return len(names)


def truth_del(path):
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
        assert info.get("SVTYPE") == "DEL", line
        return f[0], int(f[1]), int(info["END"])
    raise ValueError(f"no record in {path}")


events = (RUN / "events.txt").read_text().split()
rows = []
for n, spec in enumerate(events, start=1):
    d = RUN / str(n)
    chrom, start, end = truth_del(d / "run" / "truth.vcf")
    j = fetch(chrom, start - LEFT, start) + fetch(chrom, end, end + RIGHT)
    j_shift = fetch(chrom, start - LEFT, start) + fetch(chrom, end + 50, end + 50 + RIGHT)
    assert len(j) == LEFT + RIGHT and len(j_shift) == LEFT + RIGHT, (spec, j, j_shift)
    probes, probes_shift = (j, revcomp(j)), (j_shift, revcomp(j_shift))
    in_ref = any(
        p in fetch(chrom, max(0, bp - REF_PAD), bp + REF_PAD)
        for bp in (start, end) for p in probes
    )
    k1 = carriers(d / "run" / "sim.bam", chrom, (start, end), probes)
    k2 = carriers(d / "slice.bam", chrom, (start, end), probes)
    k3 = carriers(d / "run" / "sim.bam", chrom, (start, end), probes_shift)
    rows.append((n, spec, k1, k2, k3, in_ref))
    print(f"{n:3} {spec:28} J={j} sim={k1:3} donor={k2} shifted={k3} in_ref={in_ref}")

k1_ok = sum(1 for r in rows if r[2] >= FLOOR)
k2_ok = sum(1 for r in rows if r[3] == 0)
k3_ok = sum(1 for r in rows if r[4] == 0)
k4_hit = sum(1 for r in rows if r[5])
sims = sorted(r[2] for r in rows)
print(f"\nK1 sim >= {FLOOR}: {k1_ok} of {len(rows)}   (min {sims[0]}, median {sims[len(sims) // 2]})")
print(f"K2 donor == 0: {k2_ok} of {len(rows)}")
print(f"K3 shifted == 0: {k3_ok} of {len(rows)}")
print(f"K4 J in reference: {k4_hit} of {len(rows)}")
n = len(rows)
ok = k1_ok == n and k2_ok == n and k3_ok == n and k4_hit == 0
print(f"K: {'PASS' if ok else 'FAIL'}")
