#!/usr/bin/env python3
"""RF6 second plan, kill test K: a junction probe that tolerates the sample's own SNPs.

Reads a real_events.sh output directory run with MASTER's binary (so validate.json and
validate.json.exit are master's verdict) and KEEP_MERGED=1. For every event spike
accepted, the junction is read off the run's own truth.vcf the way validate reads a DEL
(start = POS, end = END, 0-based as TruthEvent holds them):

    J = ref[start - 15, start) + ref[end, end + 16)          (31 bases)

A read carries the junction if some 31-base window of its bases is within
MAX_MISMATCH (2) substitutions of J or of J's reverse complement. Reads are the
distinct names of primary, mapped, non-duplicate, non-QC-fail, non-supplementary
records at MAPQ >= 20 overlapping 500 bp either side of either breakpoint.

The guard: J (or its reverse complement) within GUARD_MISMATCH (4) substitutions of
any 31-base window of the reference within 1000 bp of either breakpoint makes the
row not evaluable.

    K1  every event master's validate passes (exit 0): guard silent AND >= 2 carriers
        in run/sim.bam
    K2  slice.bam (the donor, never spiked): 0 carriers -- except an event
        whose donor holds J *exactly* (0 mismatches) in >= 2 reads, which is the
        background carrying the same deletion; it is reported and left out
    K3  J built with END + 50: 0 carriers in run/sim.bam, every event
    K1b (reported, not judged) events whose master split_reads FAILs: how many
        have >= 2 carriers

Every samtools call raises on a non-zero exit (.claude/judgment-gate-cases.md).

Usage: rf6_kill2.py REFERENCE OUT_DIR [OUT_DIR ...]
"""
import json
import pathlib
import subprocess
import sys

LEFT, RIGHT, PAD, REF_PAD, MIN_MAPQ, FLOOR = 15, 16, 500, 1000, 20, 2
K = LEFT + RIGHT
MAX_MISMATCH, GUARD_MISMATCH = 2, 4
REF = sys.argv[1]
COMP = str.maketrans("ACGTN", "TGCAN")


def sh(args):
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def fetch(chrom, start0, end0):
    out = sh(["samtools", "faidx", REF, f"{chrom}:{start0 + 1}-{end0}"])
    return "".join(out.splitlines()[1:]).upper()


def revcomp(s):
    return s.translate(COMP)[::-1]


def best(seq, probes):
    """Fewest substitutions between any K-window of seq and any probe."""
    low = K + 1
    for i in range(len(seq) - K + 1):
        w = seq[i:i + K]
        for p in probes:
            d = sum(1 for a, b in zip(w, p) if a != b)
            if d < low:
                low = d
                if low == 0:
                    return 0
    return low


def carriers(bam, chrom, breakpoints, probes, max_mm):
    names = set()
    for bp in breakpoints:
        region = f"{chrom}:{max(1, bp - PAD + 1)}-{bp + PAD}"
        out = sh(["samtools", "view", "-F", "0xF04", "-q", str(MIN_MAPQ), str(bam), region])
        for line in out.splitlines():
            f = line.split("\t")
            if best(f[9].upper(), probes) <= max_mm:
                names.add(f[0])
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


def split_reads_row(report):
    for r in report["checks"]:
        if r["check"] == "split_reads":
            return r["pass"]
    raise ValueError("no split_reads row")


k1_n = k1_ok = k2_n = k2_ok = k3_n = k3_ok = k1b_n = k1b_ok = 0
k1_bad, k2_bad, k2_background, k3_bad = [], [], [], []
for out_dir in map(pathlib.Path, sys.argv[2:]):
    events = (out_dir / "events.txt").read_text().splitlines()
    for n, spec in enumerate(events, start=1):
        d = out_dir / str(n)
        if (d / "spike.exit").read_text().strip() != "0":
            print(f"{out_dir.name}/{n} {spec}: spike refused, not scored")
            continue
        master_exit = (d / "validate.json.exit").read_text().strip()
        split_pass = split_reads_row(json.loads((d / "validate.json").read_text()))
        chrom, start, end = truth_del(d / "run" / "truth.vcf")
        j = fetch(chrom, start - LEFT, start) + fetch(chrom, end, end + RIGHT)
        js = fetch(chrom, start - LEFT, start) + fetch(chrom, end + 50, end + 50 + RIGHT)
        assert len(j) == K and len(js) == K, (spec, j, js)
        probes, probes_s = (j, revcomp(j)), (js, revcomp(js))
        guard = min(best(fetch(chrom, max(0, bp - REF_PAD), bp + REF_PAD), probes)
                    for bp in (start, end))
        sim = carriers(d / "run" / "sim.bam", chrom, (start, end), probes, MAX_MISMATCH)
        donor = carriers(d / "slice.bam", chrom, (start, end), probes, MAX_MISMATCH)
        donor_exact = carriers(d / "slice.bam", chrom, (start, end), probes, 0)
        shifted = carriers(d / "run" / "sim.bam", chrom, (start, end), probes_s, MAX_MISMATCH)
        guard_fires = guard <= GUARD_MISMATCH
        tag = f"{out_dir.name}/{n}"
        print(f"{tag:14} {spec:30} master_exit={master_exit} split={'PASS' if split_pass else 'FAIL'} "
              f"guard_dist={guard} sim={sim} donor={donor} donor_exact={donor_exact} shifted={shifted}")
        if master_exit == "0":
            k1_n += 1
            if not guard_fires and sim >= FLOOR:
                k1_ok += 1
            else:
                k1_bad.append(tag)
        if donor_exact >= FLOOR:
            k2_background.append(tag)
        else:
            k2_n += 1
            if donor == 0:
                k2_ok += 1
            else:
                k2_bad.append(tag)
        k3_n += 1
        if shifted == 0:
            k3_ok += 1
        else:
            k3_bad.append(tag)
        if not split_pass:
            k1b_n += 1
            if not guard_fires and sim >= FLOOR:
                k1b_ok += 1

print(f"\nK1 master-passing events kept passing: {k1_ok} of {k1_n}  {k1_bad}")
print(f"K2 donor carriers == 0: {k2_ok} of {k2_n}  {k2_bad}; background holds J exactly: {k2_background}")
print(f"K3 shifted carriers == 0: {k3_ok} of {k3_n}  {k3_bad}")
print(f"K1b split_reads FAIL events the row would pass: {k1b_ok} of {k1b_n}")
ok = k1_ok == k1_n and k2_ok == k2_n and k3_ok == k3_n
print(f"K: {'PASS' if ok else 'FAIL'}")
