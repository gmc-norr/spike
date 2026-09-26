#!/usr/bin/env python3
"""CR4 plan, criterion C4: seeded 10 kb deletions inside the HG002 SV benchmark.

Prints one `--event` spec per line: 40 deletions of 10,000 bp on one
chromosome (chr20 unless named), each wholly inside one interval of the
T2T-Q100 stvar benchmark BED, and at least 20,000 bp from every other one (the
footprint rule needs 7,000). Positions are drawn by interval length with a
fixed seed, so the list is fixed before any census is measured.

Usage: cr4_placements.py BENCHMARK_BED [CHROM] > events.txt
"""
import random
import sys

CHROM = "chr20"
N_EVENTS = 40
LENGTH = 10_000
MIN_GAP = 20_000
SEED = 20260925


def main(bed_path, chrom=CHROM):
    intervals = []
    with open(bed_path) as bed:
        for line in bed:
            fields = line.split("\t")
            if fields[0] != chrom:
                continue
            start, end = int(fields[1]), int(fields[2])
            if end - start >= LENGTH:
                intervals.append((start, end))
    weights = [end - start - LENGTH + 1 for start, end in intervals]

    rng = random.Random(SEED)
    chosen = []
    while len(chosen) < N_EVENTS:
        start, end = rng.choices(intervals, weights=weights)[0]
        pos = rng.randrange(start, end - LENGTH + 1)
        if all(abs(pos - other) >= LENGTH + MIN_GAP for other in chosen):
            chosen.append(pos)
    for pos in sorted(chosen):
        # spike keeps a coordinate DEL spec's numbers as given, as a 0-based
        # half-open span (test_parse_coordinate_deletion): END - START bases.
        print(f"del:{chrom}:{pos}-{pos + LENGTH}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
