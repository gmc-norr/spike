#!/usr/bin/env python3
"""Read evidence of deletions in a BAM, measured the same way for real and spiked ones.

Written apart from spike on purpose: `spike validate` would grade spike with
spike's own idea of a deletion. Locked in
docs/superpowers/plans/2026-09-28-transplant-round1.md.

Per deletion [start, end) (0-based, half-open), over primary, non-duplicate,
non-QC-fail, mapped reads:
  n_gap    reads with a CIGAR gap >= half the deletion, starting within 20 bp of start
  n_clip   reads with a soft clip >= 10 bp whose clip point is within 10 bp of start or end
  n_split  reads whose SA alignment joins the two ends, each within 20 bp
  n_disc   FR pairs, forward mate ending before start + 20 and starting within 1 kb
           left of it, reverse mate starting within 1 kb right of end, insert > bound
  n_any    read pairs (by name) showing any of the four
  flank_depth  mean depth over the 1 kb left of start and the 1 kb right of end
  J        n_any / flank_depth
  E1       mean depth over [start, end) / flank_depth, for deletions >= 300 bp

A spiked BAM is read without writing a merged file: the recipient BAM with the
names spike replaced skipped (--skip-names), plus spike's sim.bam (--extra-bam).

Usage: evidence.py --events EVENTS.tsv --bam BAM --bound N [--skip-names FILE]
                   [--extra-bam BAM ...] --out OUT.tsv
EVENTS.tsv: id, chrom, start, end (tab-separated, no header).
"""
import argparse
import sys

import pysam

FLANK = 1_000
GAP_SLACK = 20
CLIP_SLACK = 10
MIN_CLIP = 10
SPLIT_SLACK = 20
MIN_E1_LEN = 300
EXCLUDED = 0x4 | 0x100 | 0x200 | 0x400 | 0x800  # unmapped, secondary, QC fail, duplicate, supplementary
REF_OPS = (0, 2, 3, 7, 8)  # M, D, N, =, X
COLUMNS = ["n_gap", "n_clip", "n_split", "n_disc", "n_any", "flank_depth", "J", "E1"]


def ref_length(cigar):
    """Reference bases a CIGAR string spans."""
    n, total = "", 0
    for ch in cigar:
        if ch.isdigit():
            n += ch
        else:
            if ch in "MDN=X":
                total += int(n)
            n = ""
    return total


def reads(paths, chrom, lo, hi):
    """Counted reads overlapping [lo, hi), from every (BAM, names to skip) pair."""
    for path, skip in paths:
        with pysam.AlignmentFile(path) as bam:
            for r in bam.fetch(chrom, max(0, lo), hi):
                if r.flag & EXCLUDED or r.query_name in skip:
                    continue
                yield r


def depth(blocks_by_read, lo, hi):
    """Mean depth over [lo, hi) from each read's aligned blocks."""
    total = 0
    for blocks in blocks_by_read:
        for a, b in blocks:
            total += max(0, min(b, hi) - max(a, lo))
    return total / (hi - lo)


def is_gap(r, start, length):
    pos = r.reference_start
    for op, n in r.cigartuples:
        if op in (2, 3) and n >= 0.5 * length and abs(pos - start) <= GAP_SLACK:
            return True
        if op in REF_OPS:
            pos += n
    return False


def is_clip(r, start, end):
    ops = r.cigartuples
    points = []
    if ops[0][0] == 4 and ops[0][1] >= MIN_CLIP:
        points.append(r.reference_start)
    if ops[-1][0] == 4 and ops[-1][1] >= MIN_CLIP:
        points.append(r.reference_end)
    return any(abs(p - start) <= CLIP_SLACK or abs(p - end) <= CLIP_SLACK for p in points)


def is_split(r, chrom, start, end):
    if not r.has_tag("SA"):
        return False
    for entry in r.get_tag("SA").rstrip(";").split(";"):
        rname, pos, _strand, cigar = entry.split(",")[:4]
        if rname != chrom:
            continue
        sa_start = int(pos) - 1
        sa_end = sa_start + ref_length(cigar)
        if abs(r.reference_end - start) <= SPLIT_SLACK and abs(sa_start - end) <= SPLIT_SLACK:
            return True
        if abs(sa_end - start) <= SPLIT_SLACK and abs(r.reference_start - end) <= SPLIT_SLACK:
            return True
    return False


def is_discordant(r, start, end, bound):
    return (r.is_paired and not r.mate_is_unmapped and r.next_reference_id == r.reference_id
            and not r.is_reverse and r.mate_is_reverse
            and start - FLANK <= r.reference_start < start and r.reference_end <= start + GAP_SLACK
            and end - GAP_SLACK <= r.next_reference_start < end + FLANK
            and abs(r.template_length) > bound)


def measure(bam, chrom, start, end, bound, skip_names=frozenset(), extra_bams=()):
    """The evidence of the deletion [start, end) of chrom; see the module doc."""
    paths = [(bam, frozenset(skip_names))] + [(p, frozenset()) for p in extra_bams]
    length = end - start
    blocks, gap, clip, split, disc = [], set(), set(), set(), set()
    for r in reads(paths, chrom, start - FLANK - 200, end + FLANK + 200):
        blocks.append(r.get_blocks())
        name = r.query_name
        if is_gap(r, start, length):
            gap.add(name)
        if is_clip(r, start, end):
            clip.add(name)
        if is_split(r, chrom, start, end):
            split.add(name)
        if is_discordant(r, start, end, bound):
            disc.add(name)
    flank = (depth(blocks, start - FLANK, start) + depth(blocks, end, end + FLANK)) / 2
    any_ = gap | clip | split | disc
    row = {"n_gap": len(gap), "n_clip": len(clip), "n_split": len(split), "n_disc": len(disc),
           "n_any": len(any_), "flank_depth": flank}
    row["J"] = len(any_) / flank if flank > 0 else None
    row["E1"] = depth(blocks, start, end) / flank if length >= MIN_E1_LEN and flank > 0 else None
    return row


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--events", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--bound", type=int, required=True)
    ap.add_argument("--skip-names")
    ap.add_argument("--extra-bam", action="append", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    skip = set()
    if a.skip_names:
        with open(a.skip_names) as fh:
            skip = {line.strip() for line in fh if line.strip()}
    with open(a.events) as fh, open(a.out, "w") as out:
        out.write("\t".join(["id", "chrom", "start", "end"] + COLUMNS) + "\n")
        for line in fh:
            eid, chrom, start, end = line.rstrip("\n").split("\t")[:4]
            row = measure(a.bam, chrom, int(start), int(end), a.bound, skip, a.extra_bam)
            cells = [eid, chrom, start, end] + ["" if row[c] is None else
                                                (f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c]))
                                                for c in COLUMNS]
            out.write("\t".join(cells) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
