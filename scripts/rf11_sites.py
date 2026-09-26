#!/usr/bin/env python3
"""RF11 plan: the seeded sites the kill test runs on, fixed before any count is seen.

Prints one line per site, `<set>\t<chrom>\t<pos>`, for three sets on one chromosome:

  null-rand  200 positions anywhere in the HG002 stvar benchmark
  null-rep   200 positions inside a RepeatMasker Simple_repeat or Low_complexity
             interval, also in the benchmark -- where aligners clip reads most
  k1         8 positions for the correct spike-ins, drawn like null-rand but with
             their own seed and at least 250 kb apart (slice_loop pads 100 kb)

Every site is at least 1,000 bp inside its benchmark interval and at least 1,000 bp
from every HG002 variant whose REF and ALT differ in length by 10 bp or more, so
nothing HG002 carries can give clipped or inserted reads within `ins_reads`'s
+/-100 bp. Null sites are at least 1,000 bp from each other.

`pos` is used as a VCF POS: it is what `check_ins_reads` takes as `event.start`.

Usage: rf11_sites.py BENCHMARK_BED HG002_VCF_GZ RMSK_BED [CHROM] > sites.tsv
"""
import bisect
import gzip
import random
import sys

CHROM = "chr20"
N_NULL = 200
N_K1 = 8
EDGE = 1_000
VARIANT_GAP = 1_000
NULL_GAP = 1_000
K1_GAP = 250_000
MIN_INDEL = 10
SEEDS = {"null-rand": 20260926, "null-rep": 20260927, "k1": 20260928}


def read_bed(path, chrom, keep=None):
    out = []
    with open(path) as bed:
        for line in bed:
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            if keep is not None and not keep(f):
                continue
            out.append((int(f[1]), int(f[2])))
    return sorted(out)


def indel_spans(vcf_path, chrom):
    """(start, end) of every variant whose REF and an ALT differ by >= MIN_INDEL bp."""
    spans = []
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[0] != chrom:
                continue
            pos, ref, alts = int(f[1]), f[3], f[4].split(",")
            if any(a.startswith("<") or abs(len(a) - len(ref)) >= MIN_INDEL for a in alts):
                spans.append((pos, pos + len(ref)))
    return sorted(spans)


def merged_exclusion(spans):
    """Each span grown by VARIANT_GAP on both sides, merged: sorted and disjoint."""
    out = []
    for s, e in sorted((s - VARIANT_GAP, e + VARIANT_GAP) for s, e in spans):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def near_indel(pos, excl, starts):
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and pos < excl[i][1]


def overlap(a, b):
    """Intersections of two sorted interval lists."""
    out, j = [], 0
    for s, e in a:
        while j < len(b) and b[j][1] <= s:
            j += 1
        k = j
        while k < len(b) and b[k][0] < e:
            lo, hi = max(s, b[k][0]), min(e, b[k][1])
            if hi > lo:
                out.append((lo, hi))
            k += 1
    return out


def draw(name, intervals, bench, n, gap, spans, starts):
    """n seeded positions inside `intervals` and >= EDGE inside a benchmark interval."""
    bench_starts = [s for s, _ in bench]
    weights = [e - s for s, e in intervals]
    rng = random.Random(SEEDS[name])
    chosen, tries = [], 0
    while len(chosen) < n:
        tries += 1
        if tries > 1_000_000:
            sys.exit(f"{name}: only {len(chosen)} of {n} sites after 1e6 tries")
        s, e = rng.choices(intervals, weights=weights)[0]
        pos = rng.randrange(s, e)
        b = bisect.bisect_right(bench_starts, pos) - 1
        if b < 0 or not (bench[b][0] + EDGE <= pos < bench[b][1] - EDGE):
            continue
        if near_indel(pos, spans, starts):
            continue
        if any(abs(pos - other) < gap for other in chosen):
            continue
        chosen.append(pos)
    return sorted(chosen)


def main(bench_path, vcf_path, rmsk_path, chrom=CHROM):
    bench = read_bed(bench_path, chrom)
    spans = merged_exclusion(indel_spans(vcf_path, chrom))
    starts = [s for s, _ in spans]
    rep = read_bed(
        rmsk_path,
        chrom,
        keep=lambda f: f[3].split("/")[1] in ("Simple_repeat", "Low_complexity"),
    )
    sets = {
        "null-rand": draw("null-rand", bench, bench, N_NULL, NULL_GAP, spans, starts),
        "null-rep": draw("null-rep", overlap(rep, bench), bench, N_NULL, NULL_GAP, spans, starts),
        "k1": draw("k1", bench, bench, N_K1, K1_GAP, spans, starts),
    }
    for name, positions in sets.items():
        for pos in positions:
            print(f"{name}\t{chrom}\t{pos}")


if __name__ == "__main__":
    main(*sys.argv[1:])
