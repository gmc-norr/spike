#!/usr/bin/env python3
"""T7 control 2: 40 random 14 kb windows inside the HG002 SV benchmark on chr20.

The same shape `cr4_placements.py` draws, but 14 kb (a 10 kb span plus 2 kb of
flank on each side, the footprint T7 measures) and with the seed the T7 plan wrote
down before the draw, so the control cannot have been chosen to agree.

Prints one `--event`-shaped spec per line so `t7_footprint_variants.py` can read it;
the spans are 14 kb and the scanner is run with FLANK=0 on them.

Usage: t7_random_windows.py BENCHMARK_BED > windows.txt
"""
import random
import sys

CHROM = "chr20"
N = 40
LENGTH = 14_000
MIN_GAP = 20_000
SEED = 20260926        # locked in the T7 plan before the draw


def main(bed_path):
    intervals = []
    with open(bed_path) as bed:
        for line in bed:
            f = line.split("\t")
            if f[0] != CHROM:
                continue
            start, end = int(f[1]), int(f[2])
            if end - start >= LENGTH:
                intervals.append((start, end))
    weights = [end - start - LENGTH + 1 for start, end in intervals]
    rng = random.Random(SEED)
    chosen = []
    while len(chosen) < N:
        start, end = rng.choices(intervals, weights=weights)[0]
        pos = rng.randrange(start, end - LENGTH + 1)
        if all(abs(pos - other) >= LENGTH + MIN_GAP for other in chosen):
            chosen.append(pos)
    for pos in sorted(chosen):
        print(f"win:{CHROM}:{pos}-{pos + LENGTH}")


if __name__ == "__main__":
    main(sys.argv[1])
