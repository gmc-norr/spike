#!/usr/bin/env python3
"""RF14: the kill test's deletion sites, on fresh seeds, fixed before any count is seen.

Prints `<set>\t<chrom>\t<pos>\t<end or ->` lines:

  rand   6 random benchmark positions for spiked deletions of 50, 300, 1000 and
         10000 bp starting there. `pos` and each of the four ends are >= 1 kb
         inside one benchmark interval and >= 1 kb from every HG002 variant of
         >= 10 bp; the sites are >= 250 kb apart (slice_loop.sh pads 100 kb).
  hg     12 real HG002 deletions, spiked at their own POS and END: 6 of 50-299 bp
         and 6 of 300 bp to 50 kb. Pure (one-base ALT), single-ALT, het or hom-alt,
         both ends >= 1 kb inside one benchmark interval, and no other variant of
         >= 10 bp within 1 kb.
  null   200 random benchmark positions, >= 1 kb from every HG002 variant of >= 10 bp.

No site is within 2 kb of any position in SEEN files, which are read as:
`del:chrom:start-end` lines (both ends), VCF records (POS and INFO END), or
`<set>\t<chrom>\t<pos>...` TSV lines.

A spike deletion `del:chrom:S-E` is truth POS S and END E; a VCF deletion at POS p
with REF of n + 1 bases is `del:chrom:p-(p+n)`.

Usage: rf14_sites.py BENCHMARK_BED HG002_VCF_GZ SEEN... > rf14_sites.tsv
"""
import bisect
import gzip
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rf11_sites as base  # noqa: E402

CHROM = "chr20"
SEEDS = {"rand": 20261021, "hg": 20261022, "null": 20261023}
SIZES = (50, 300, 1000, 10000)
SEEN_GAP = 2000
DEL_RE = re.compile(r"del:(chr\w+):(\d+)-(\d+)")


def seen_positions(paths):
    out = set()
    for path in paths:
        opener = gzip.open if path.endswith(".gz") else open
        for line in opener(path, "rt"):
            m = DEL_RE.search(line)
            if m:
                out.update((int(m.group(2)), int(m.group(3))))
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 8 and f[1].isdigit():
                out.add(int(f[1]))
                e = re.search(r"(?:^|;)END=(\d+)", f[7])
                if e:
                    out.add(int(e.group(1)))
            elif len(f) >= 3 and f[2].isdigit():
                out.add(int(f[2]))
    return sorted(out)


def near(pos, sorted_list, gap):
    i = bisect.bisect_left(sorted_list, pos - gap)
    return i < len(sorted_list) and sorted_list[i] < pos + gap


def in_bench(pos, bench, bench_starts):
    b = bisect.bisect_right(bench_starts, pos) - 1
    if b < 0 or not (bench[b][0] + base.EDGE <= pos < bench[b][1] - base.EDGE):
        return None
    return b


def draw_rand(bench, excl, starts, seen):
    bench_starts = [s for s, _ in bench]
    weights = [e - s for s, e in bench]
    rng = random.Random(SEEDS["rand"])
    chosen, tries = [], 0
    while len(chosen) < 6:
        tries += 1
        if tries > 1_000_000:
            sys.exit(f"rand: only {len(chosen)} of 6 sites after 1e6 tries")
        s, e = rng.choices(bench, weights=weights)[0]
        pos = rng.randrange(s, e)
        points = [pos] + [pos + n for n in SIZES]
        blocks = {in_bench(p, bench, bench_starts) for p in points}
        if None in blocks or len(blocks) != 1:
            continue
        if any(base.near_indel(p, excl, starts) or near(p, seen, SEEN_GAP) for p in points):
            continue
        if any(abs(pos - other) < base.K1_GAP for other in chosen):
            continue
        chosen.append(pos)
    return sorted(chosen)


def draw_null(bench, excl, starts, seen, avoid):
    bench_starts = [s for s, _ in bench]
    weights = [e - s for s, e in bench]
    rng = random.Random(SEEDS["null"])
    chosen, tries = [], 0
    while len(chosen) < 200:
        tries += 1
        if tries > 1_000_000:
            sys.exit(f"null: only {len(chosen)} of 200 sites after 1e6 tries")
        s, e = rng.choices(bench, weights=weights)[0]
        pos = rng.randrange(s, e)
        # The null deletion is POS to POS + 1000, so both ends must be clean.
        points = (pos, pos + 1000)
        if any(in_bench(p, bench, bench_starts) is None for p in points):
            continue
        if any(base.near_indel(p, excl, starts) or near(p, seen, SEEN_GAP) for p in points):
            continue
        if any(abs(pos - other) < base.NULL_GAP for other in chosen + avoid):
            continue
        chosen.append(pos)
    return sorted(chosen)


def real_deletions(vcf_path, bench, indels, lo, hi):
    bench_starts = [s for s, _ in bench]
    istarts = [s for s, _ in indels]
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
            n = len(ref) - 1
            if not (len(alt) == 1 and ref[0] == alt and lo <= n <= hi):
                continue
            end = pos + n
            b = in_bench(pos, bench, bench_starts)
            if b is None or in_bench(end, bench, bench_starts) != b:
                continue
            # Isolated: the only variant of >= 10 bp within 1 kb is the event.
            i = bisect.bisect_right(istarts, end + base.VARIANT_GAP)
            others = [
                (s, e) for s, e in indels[max(0, i - 200):i]
                if e > pos - base.VARIANT_GAP and (s, e) != (pos, pos + len(ref))
            ]
            if others:
                continue
            out.append((pos, end))
    return out


def main(bench_path, vcf_path, *seen_paths):
    seen = seen_positions(seen_paths)
    bench = base.read_bed(bench_path, CHROM)
    indels = base.indel_spans(vcf_path, CHROM)
    excl = base.merged_exclusion(indels)
    starts = [s for s, _ in excl]
    rand = draw_rand(bench, excl, starts, seen)
    rng = random.Random(SEEDS["hg"])
    fresh = lambda ds: [d for d in ds if not (near(d[0], seen, SEEN_GAP) or near(d[1], seen, SEEN_GAP))]
    short = fresh(real_deletions(vcf_path, bench, indels, 50, 299))
    long_ = fresh(real_deletions(vcf_path, bench, indels, 300, 50_000))
    hg = sorted(rng.sample(short, 6) + rng.sample(long_, 6))
    null = draw_null(bench, excl, starts, seen, rand)
    used = rand + [p for d in hg for p in d] + null
    assert not any(near(p, seen, SEEN_GAP) for p in used), "a site was seen before"
    print(f"# candidates: {len(short)} of 50-299 bp, {len(long_)} of 300 bp-50 kb; "
          f"{len(seen)} seen positions", file=sys.stderr)
    for pos in rand:
        print(f"rand\t{CHROM}\t{pos}\t-")
    for pos, end in hg:
        print(f"hg\t{CHROM}\t{pos}\t{end}")
    for pos in null:
        print(f"null\t{CHROM}\t{pos}\t-")


if __name__ == "__main__":
    main(*sys.argv[1:])
