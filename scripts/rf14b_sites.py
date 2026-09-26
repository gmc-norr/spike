#!/usr/bin/env python3
"""RF14 second attempt: the kill test's deletion sites, on fresh seeds.

Prints `<set>\t<chrom>\t<pos>\t<end or ->` lines, drawn as rf14_sites.py draws
them but with new seeds, one more set, and no site within 2 kb of anything seen:

  rand   6 random benchmark positions for deletions of 50, 300, 1000 and 10000 bp
  hg     12 real HG002 deletions: 6 of 50-299 bp and 6 of 300 bp to 50 kb
  low    30 random benchmark positions for a 1000 bp deletion at VAF 0.1, where
         the first attempt saw 0-3 carriers and the chance-zero rate is unmeasured
  null   200 random benchmark positions

`rand` and `low` are one draw of 36 positions, >= 250 kb apart, split 6 and 30 in
draw order; a `low` site's 1000 bp end is as clean as a `rand` site's four ends
(the 10000 bp end is required of both, which only makes `low` stricter).

Usage: rf14b_sites.py BENCHMARK_BED HG002_VCF_GZ SEEN... > rf14b_sites.tsv
"""
import bisect
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rf11_sites as base  # noqa: E402
import rf14_sites as first  # noqa: E402

CHROM = first.CHROM
SEEDS = {"rand": 20261031, "hg": 20261032, "null": 20261033}
N_RAND, N_LOW = 6, 30


def draw_spans(bench, excl, starts, seen, n):
    """n positions whose start and every rf14 size's end are clean, in draw order."""
    bench_starts = [s for s, _ in bench]
    weights = [e - s for s, e in bench]
    rng = random.Random(SEEDS["rand"])
    chosen, tries = [], 0
    while len(chosen) < n:
        tries += 1
        if tries > 1_000_000:
            sys.exit(f"rand/low: only {len(chosen)} of {n} sites after 1e6 tries")
        s, e = rng.choices(bench, weights=weights)[0]
        pos = rng.randrange(s, e)
        points = [pos] + [pos + size for size in first.SIZES]
        blocks = {first.in_bench(p, bench, bench_starts) for p in points}
        if None in blocks or len(blocks) != 1:
            continue
        if any(base.near_indel(p, excl, starts) or first.near(p, seen, first.SEEN_GAP)
               for p in points):
            continue
        if any(abs(pos - other) < base.K1_GAP for other in chosen):
            continue
        chosen.append(pos)
    return chosen


def main(bench_path, vcf_path, *seen_paths):
    first.SEEDS.update(SEEDS)  # draw_null reads its seed from there
    seen = first.seen_positions(seen_paths)
    bench = base.read_bed(bench_path, CHROM)
    indels = base.indel_spans(vcf_path, CHROM)
    excl = base.merged_exclusion(indels)
    starts = [s for s, _ in excl]
    spans = draw_spans(bench, excl, starts, seen, N_RAND + N_LOW)
    rand, low = sorted(spans[:N_RAND]), sorted(spans[N_RAND:])
    rng = random.Random(SEEDS["hg"])
    fresh = lambda ds: [d for d in ds
                        if not (first.near(d[0], seen, first.SEEN_GAP)
                                or first.near(d[1], seen, first.SEEN_GAP))]
    short = fresh(first.real_deletions(vcf_path, bench, indels, 50, 299))
    long_ = fresh(first.real_deletions(vcf_path, bench, indels, 300, 50_000))
    hg = sorted(rng.sample(short, 6) + rng.sample(long_, 6))
    null = first.draw_null(bench, excl, starts, seen, spans)
    used = spans + [p for d in hg for p in d] + null
    assert not any(first.near(p, seen, first.SEEN_GAP) for p in used), "a site was seen before"
    print(f"# candidates: {len(short)} of 50-299 bp, {len(long_)} of 300 bp-50 kb; "
          f"{len(seen)} seen positions", file=sys.stderr)
    for pos in rand:
        print(f"rand\t{CHROM}\t{pos}\t-")
    for pos, end in hg:
        print(f"hg\t{CHROM}\t{pos}\t{end}")
    for pos in low:
        print(f"low\t{CHROM}\t{pos}\t-")
    for pos in null:
        print(f"null\t{CHROM}\t{pos}\t-")


if __name__ == "__main__":
    main(*sys.argv[1:])
