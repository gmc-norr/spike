#!/usr/bin/env python3
"""--edit-model origin, the physics test: a made-up genome with an exact twin.

chrT (80 kb): unique U1 [0,20000), S [20000,30000), unique U2 [30000,50000), an
exact copy of S [50000,60000), unique U3 [60000,80000). The event is a het 1 kb
deletion in the middle of the first S, so its 5 kb footprint lies inside S.

Truth: reads from the diploid genome carrying the deletion on one copy (two
seeds). Donor: reads from it without. spike runs on the donor three ways --
clean, --min-mapq 0, --edit-model origin -- each through align.sh and merge.sh.
Measured: depth over L and P (any MAPQ) over the unique windows' mean. The rule
is locked in docs/superpowers/plans/2026-09-27-edit-model-origin.md, Task 12.

Usage: origin_physics.py OUT_DIR SPIKE [THREADS]
"""
import os
import random
import subprocess
import sys

GENOME_SEED = 20260927
TRUTH_SEEDS = (1, 2)
DONOR_SEED = 3
SPIKE_SEED = 7
PER_COPY = 30
READ_LEN = 150
MARGIN = 0.05
DELETED = (24500, 25500)
L = (24600, 25400)
P = (54600, 55400)
BASE = [(5000, 15000), (35000, 45000), (65000, 75000)]
EVENT = "del:chrT:24500-25500"


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


def genome():
    rng = random.Random(GENOME_SEED)
    seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    u1, s, u2, u3 = seq(20000), seq(10000), seq(20000), seq(20000)
    return u1 + s + u2 + s + u3


def reads(out, name, copies, seed, ref, threads):
    """wgsim from `copies`, qualities to Q30, aligned; returns the BAM."""
    haps = f"{out}/{name}.haps.fa"
    write_fasta(haps, [(f"copy{i}", c) for i, c in enumerate(copies)])
    pairs = PER_COPY * sum(len(c) for c in copies) // (2 * READ_LEN)
    raw = [f"{out}/{name}.raw{i}.fq" for i in (1, 2)]
    sh(["wgsim", "-N", str(pairs), "-1", str(READ_LEN), "-2", str(READ_LEN), "-d", "400",
        "-s", "50", "-e", "0.001", "-r", "0", "-R", "0", "-X", "0", "-S", str(seed),
        haps, *raw], stdout=subprocess.DEVNULL)
    fq = [f"{out}/{name}.R{i}.fq" for i in (1, 2)]
    for src, dst in zip(raw, fq):
        with open(src) as fin, open(dst, "w") as fout:
            for i, line in enumerate(fin):
                line = line.rstrip("\n")
                if i % 4 == 0 and line[-2:] in ("/1", "/2"):
                    line = line[:-2]
                if i % 4 == 3:
                    line = "?" * len(line)
                fout.write(line + "\n")
    bam = f"{out}/{name}.bam"
    align = subprocess.Popen(["bwa-mem2", "mem", "-t", str(threads), "-R", r"@RG\tID:toy\tSM:TOY",
                              ref, *fq], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    sh(["samtools", "sort", "-@", str(threads), "-o", bam, "-"], stdin=align.stdout)
    if align.wait() != 0:
        raise SystemExit(f"bwa-mem2 failed for {name}")
    sh(["samtools", "index", bam])
    return bam


def spike(out, name, spike_bin, donor, ref, extra, threads):
    """spike -> align.sh -> merge.sh; returns (spike exit, merged BAM or None)."""
    d = f"{out}/{name}"
    with open(f"{d}.log", "w") as log:
        rc = subprocess.run([spike_bin, "--bam", donor, "--reference", ref, "--event", EVENT,
                             "--seed", str(SPIKE_SEED), "--allow-resistant", "-o", d, *extra],
                            stderr=log, stdout=log).returncode
        if rc != 0:
            return rc, None
        sh(["bash", f"{d}/align.sh", ref, str(threads)], stdout=log, stderr=log)
        sh(["bash", f"{d}/merge.sh", donor, ref, str(threads)], stdout=log, stderr=log)
    return 0, f"{d}/merged.bam"


def mean_depth(bam, start, end):
    out = subprocess.run(["samtools", "depth", "-a", "-r", f"chrT:{start + 1}-{end}", bam],
                         capture_output=True, text=True, check=True).stdout
    vals = [int(line.split("\t")[2]) for line in out.splitlines()]
    return sum(vals) / len(vals)


def ratios(bam):
    base = sum(mean_depth(bam, a, b) * (b - a) for a, b in BASE) / sum(b - a for a, b in BASE)
    return mean_depth(bam, *L) / base, mean_depth(bam, *P) / base


def check_counts(primary, mapq0, with_xa):
    """The footprint lies 2.5 kb inside an exact twin, so the aligner cannot
    place any read there: every primary read must be MAPQ 0 and list the twin
    in XA. Anything less and the donor lacks the physics this test is about.
    The first version raised only when MAPQ 0 reads lacked XA, so a donor in
    which the aligner resolved the twin (0 MAPQ 0) passed silently (PD-31)."""
    if primary == 0:
        raise SystemExit("harness broken: no primary read over the footprint")
    if mapq0 != primary:
        raise SystemExit(f"harness broken: {primary - mapq0} of {primary} primary reads over "
                         "the footprint have MAPQ above 0, so the aligner told the twins apart")
    if with_xa != mapq0:
        raise SystemExit(f"harness broken: {mapq0 - with_xa} of the {mapq0} MAPQ 0 reads carry no XA")


def harness_check(donor):
    out = subprocess.run(["samtools", "view", "-F", "0x904", donor, "chrT:22501-27500"],
                         capture_output=True, text=True, check=True).stdout.splitlines()
    mapq0 = [r for r in out if r.split("\t")[4] == "0"]
    with_xa = [r for r in mapq0 if "\tXA:Z:" in r]
    print(f"harness: {len(out)} primary records over the footprint, {len(mapq0)} MAPQ 0, "
          f"{len(with_xa)} of those with XA")
    check_counts(len(out), len(mapq0), len(with_xa))


def main(out, spike_bin, threads="8"):
    os.makedirs(out, exist_ok=True)
    g = genome()
    ref = f"{out}/ref.fa"
    write_fasta(ref, [("chrT", g)])
    sh(["samtools", "faidx", ref])
    sh(["bwa-mem2", "index", ref], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    deleted = g[:DELETED[0]] + g[DELETED[1]:]

    truths = [ratios(reads(out, f"truth{s}", [deleted, g], s, ref, threads)) for s in TRUTH_SEEDS]
    donor = reads(out, "donor", [g, g], DONOR_SEED, ref, threads)
    harness_check(donor)

    runs = {}
    for name, extra in [("clean", []), ("mapq0", ["--min-mapq", "0"]),
                        ("origin", ["--edit-model", "origin"])]:
        rc, merged = spike(out, name, spike_bin, donor, ref, extra, threads)
        runs[name] = (rc, ratios(merged) if merged else None)

    band = [(min(t[i] for t in truths) - MARGIN, max(t[i] for t in truths) + MARGIN) for i in (0, 1)]
    inside = lambda r: r is not None and all(band[i][0] <= r[i] <= band[i][1] for i in (0, 1))
    for s, t in zip(TRUTH_SEEDS, truths):
        print(f"truth{s}\texit -\tL {t[0]:.3f}\tP {t[1]:.3f}")
    for name, (rc, r) in runs.items():
        cells = f"L {r[0]:.3f}\tP {r[1]:.3f}" if r else "refused"
        print(f"{name}\texit {rc}\t{cells}\tinside {inside(r)}")
    print(f"band L [{band[0][0]:.3f}, {band[0][1]:.3f}]  P [{band[1][0]:.3f}, {band[1][1]:.3f}]")

    mapq0_rc, mapq0 = runs["mapq0"]
    clean_rc, clean = runs["clean"]
    origin_rc, origin = runs["origin"]
    with open(f"{out}/clean.log") as fh:
        clean_refused_as_expected = clean_rc != 0 and "has no donor coverage" in fh.read()
    # A control is evidence only when it ran as meant.
    if origin_rc != 0:
        print("NO VERDICT: origin failed to run -- fix it and run again")
    elif mapq0_rc != 0:
        print("NO VERDICT: the --min-mapq 0 control failed to run -- fix it and run again")
    elif clean_rc != 0 and not clean_refused_as_expected:
        print("NO VERDICT: the clean control failed, but not with its expected refusal")
    elif inside(mapq0) or (clean_rc == 0 and inside(clean)):
        print("VERDICT: INCONCLUSIVE -- the test cannot tell the models apart")
    elif inside(origin):
        print("VERDICT: SUPPORTED")
    else:
        print("VERDICT: REFUTED")


if __name__ == "__main__":
    main(*sys.argv[1:])
