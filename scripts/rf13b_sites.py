#!/usr/bin/env python3
"""RF13 second attempt: the kill test's sites, on fresh seeds.

Prints `<set>\t<chrom>\t<pos>\t<seq or ->` lines, drawn as rf13_sites.py draws
them but with new seeds, and with no site shared with rf11_sites.tsv or
rf13_sites.tsv:

  rand   6 random benchmark positions, isolated, >= 250 kb apart
  hg     12 real HG002 insertions with their own bases: 6 of 1-19 bp (the short
         end the first attempt never had real data for) and 6 of 50-300 bp
  null   200 random benchmark positions

Usage: rf13b_sites.py BENCHMARK_BED HG002_VCF_GZ RMSK_BED SEEN_SITES_TSV... > rf13b_sites.tsv
"""
import bisect
import gzip
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rf11_sites as base  # noqa: E402

CHROM = "chr20"
SEEDS = {"rand": 20261011, "hg": 20261012, "null": 20261013}


def real_insertions(vcf_path, bench, excl, lo, hi):
    bench_starts = [s for s, _ in bench]
    starts = [s for s, _ in excl]
    out = []
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != CHROM or "," in f[4]:
                continue
            pos, ref, alt = int(f[1]), f[3], f[4]
            if f[9].split(":")[0].replace("|", "/") not in ("0/1", "1/0", "1/1"):
                continue
            if not (len(ref) == 1 and alt[0] == ref and lo <= len(alt) - 1 <= hi):
                continue
            if set(alt[1:].upper()) - set("ACGT"):
                continue
            b = bisect.bisect_right(bench_starts, pos) - 1
            if b < 0 or not (bench[b][0] + base.EDGE <= pos < bench[b][1] - base.EDGE):
                continue
            # Isolated from every variant of >= 10 bp but itself. A 1-19 bp
            # insertion is below that size, so it has no interval of its own and
            # must fall in none.
            i = bisect.bisect_right(starts, pos) - 1
            inside = i >= 0 and pos < excl[i][1]
            own = (pos - base.VARIANT_GAP, pos + 1 + base.VARIANT_GAP)
            if inside and not (len(alt) - 1 >= base.MIN_INDEL and excl[i] == own):
                continue
            out.append((pos, alt[1:].upper()))
    return out


def main(bench_path, vcf_path, rmsk_path, *seen_paths):
    seen = {int(l.split("\t")[2]) for p in seen_paths for l in open(p)}
    bench = base.read_bed(bench_path, CHROM)
    excl = base.merged_exclusion(base.indel_spans(vcf_path, CHROM))
    starts = [s for s, _ in excl]
    base.SEEDS.update({"rand": SEEDS["rand"], "null": SEEDS["null"]})
    rand = base.draw("rand", bench, bench, 6, base.K1_GAP, excl, starts)
    null = base.draw("null", bench, bench, 200, base.NULL_GAP, excl, starts)
    rng = random.Random(SEEDS["hg"])
    short = [r for r in real_insertions(vcf_path, bench, excl, 1, 19) if r[0] not in seen]
    long_ = [r for r in real_insertions(vcf_path, bench, excl, 50, 300) if r[0] not in seen]
    hg = sorted(rng.sample(short, 6) + rng.sample(long_, 6))
    assert not (set(rand) | set(null) | {p for p, _ in hg}) & seen, "a site was seen before"
    for pos in rand:
        print(f"rand\t{CHROM}\t{pos}\t-")
    for pos, seq in hg:
        print(f"hg\t{CHROM}\t{pos}\t{seq}")
    for pos in null:
        print(f"null\t{CHROM}\t{pos}\t-")


if __name__ == "__main__":
    main(*sys.argv[1:])
