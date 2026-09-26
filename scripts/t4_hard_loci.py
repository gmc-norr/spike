#!/usr/bin/env python3
"""T4: find chr20's hardest-to-map 10 kb windows at ordinary depth.

Steps across chr20 in 100 kb strides and, at each stride, measures the 10 kb window
starting there over mapped, primary, non-duplicate, non-QC-fail records:

  low_share  the share with MAPQ < 20
  any_depth  mean read depth at any MAPQ

A window is eligible when its any_depth lies in [0.5, 2.0] times the median any_depth
over all scanned windows -- "ordinary total depth", so neither an empty window (chr20's
N runs) nor a pile-up can be picked. The bounds and the stride are the T4 plan's,
locked before any window was seen.

Prints every window as TSV, then the three eligible windows with the highest
low_share and the three with the lowest. Judges nothing else.

Usage: t4_hard_loci.py BAM CHROM_LEN [STRIDE] [WINDOW] [WORKERS]
"""
import concurrent.futures
import statistics
import subprocess
import sys

CHROM = "chr20"
EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400 | 0x800   # for `samtools view -F`
# `samtools depth`'s -G ADDS to a default filter-out list already holding UNMAP,
# SECONDARY, QCFAIL and DUP; its default min-MQ is 0. -J counts a CIGAR-deleted
# position as covered, matching validate's own count_depth_in_region.
DEPTH_FLAGS = "-a -Q 0 -q 0 -J -G SUPPLEMENTARY"
DEPTH_LO, DEPTH_HI = 0.5, 2.0


def sh(command):
    """Run a shell pipeline, refusing to return silence on a failure.

    An unrecognised samtools option prints usage to stderr, exits non-zero and
    writes nothing to stdout -- which a caller summing stdout reads as zero. That
    happened in T3's first census run. Fail loudly.
    """
    done = subprocess.run(command, shell=True, capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(f"command failed ({done.returncode}): {command}\n{done.stderr}")
    return done.stdout


def measure(args):
    bam, start, window = args
    region = f"{CHROM}:{start + 1}-{start + window}"
    allc = int(sh(f"samtools view -c -F {EXCLUDE} '{bam}' '{region}'").strip() or 0)
    hi = int(sh(f"samtools view -c -F {EXCLUDE} -q 20 '{bam}' '{region}'").strip() or 0)
    depth_out = sh(f"samtools depth {DEPTH_FLAGS} -r '{region}' '{bam}'")
    total = sum(int(l.split("\t")[2]) for l in depth_out.splitlines() if l)
    return {
        "start": start,
        "end": start + window,
        "n": allc,
        "low_share": (allc - hi) / allc if allc else float("nan"),
        "any_depth": total / window,
    }


def main(bam, chrom_len, stride="100000", window="10000", workers="8"):
    chrom_len, stride, window = int(chrom_len), int(stride), int(window)
    starts = list(range(0, chrom_len - window, stride))
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(workers)) as pool:
        rows = list(pool.map(measure, [(bam, s, window) for s in starts]))

    depths = [r["any_depth"] for r in rows]
    median = statistics.median(depths)
    lo, hi = DEPTH_LO * median, DEPTH_HI * median
    for r in rows:
        r["eligible"] = lo <= r["any_depth"] <= hi and r["n"] > 0

    print(f"# scanned {len(rows)} windows of {window} bp every {stride} bp on {CHROM}")
    print(f"# median any_depth {median:.2f}; eligibility [{lo:.2f}, {hi:.2f}]")
    eligible = [r for r in rows if r["eligible"]]
    print(f"# eligible {len(eligible)}; excluded {len(rows) - len(eligible)}")
    print("start\tend\tn_reads\tlow_share\tany_depth\teligible")
    for r in rows:
        print(f'{r["start"]}\t{r["end"]}\t{r["n"]}\t{r["low_share"]:.4f}\t'
              f'{r["any_depth"]:.2f}\t{int(r["eligible"])}')

    by_low = sorted(eligible, key=lambda r: r["low_share"])
    print()
    print("# THREE HARDEST (highest low_share among eligible)")
    for r in reversed(by_low[-3:]):
        print(f'HARD\t{CHROM}:{r["start"]}-{r["end"]}\tlow_share={r["low_share"]:.4f}\t'
              f'any_depth={r["any_depth"]:.2f}\tn={r["n"]}')
    print("# THREE EASIEST (lowest low_share among eligible) -- the control")
    for r in by_low[:3]:
        print(f'EASY\t{CHROM}:{r["start"]}-{r["end"]}\tlow_share={r["low_share"]:.4f}\t'
              f'any_depth={r["any_depth"]:.2f}\tn={r["n"]}')


if __name__ == "__main__":
    main(*sys.argv[1:6])
