#!/usr/bin/env python3
"""N19: where do the indel failures N18 leaves sit, and what do their reads show?

Counts each site the way `spike validate` does after N18 (the pad rule, one
vote per fragment, reads reaching 10 bp past the repeat region, MAPQ >= 20),
re-using `n18_indel_flank.py`, then reports REVIEW.md's N19 readouts:
repeat context (R1), other gaps in voting reads (R2), what the MAPQ filter
removes (R3) and how far off the failing sites are (R4).

Usage:
    n19_indel_residual.py --bam B --ref R --indels I.vcf --check-log validate-debug.err
"""

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
import n18_indel_flank as n18  # noqa: E402

FLANK = 10


def usable_any_mapq(r):
    return not (
        r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail
    ) and r.query_name


def other_gap(r, lo, hi, kind, length, pos):
    """Whether the read has an I/D operation starting inside [lo, hi] other
    than one counted carrier operation (kind, length, within the pad)."""
    ref_pos = r.reference_start
    skipped_carrier = False
    for op, n in r.cigartuples:
        if op in (1, 2) and lo <= ref_pos <= hi:
            is_carrier_op = (
                ((op == 2 and kind == "D") or (op == 1 and kind == "I"))
                and n == length
                and abs(ref_pos - (pos + 1)) <= n18.INDEL_POS_PAD
            )
            if is_carrier_op and not skipped_carrier:
                skipped_carrier = True
            else:
                return True
        if op in (0, 2, 3, 7, 8):
            ref_pos += n
    return False


def measure(bam, fasta, sites):
    out = []
    for chrom, vcf_pos, ref, alt in sites:
        pos = vcf_pos - 1
        kind = "D" if len(ref) > len(alt) else "I"
        length = abs(len(ref) - len(alt))
        s, e = n18.repeat_region(fasta, chrom, pos, ref, alt)
        extension = (e - s) - (length if kind == "D" else 0)
        lo, hi = s - 1 - FLANK, e + FLANK
        reads = []
        for r in bam.fetch(chrom, pos, pos + len(ref) + 1):
            if not usable_any_mapq(r):
                continue
            v = n18.indel_vote(r, pos, len(ref), kind, length)
            if v is None or not n18.qualifies(r.reference_start, r.reference_end, s, e, FLANK, True):
                continue
            reads.append((r.query_name, v, r.mapping_quality >= n18.MIN_MAPQ, other_gap(r, lo, hi, kind, length, pos)))
        out.append((vcf_pos, kind, length, extension, reads))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--indels", required=True)
    ap.add_argument("--check-log", required=True)
    a = ap.parse_args()

    sites = measure(pysam.AlignmentFile(a.bam), pysam.FastaFile(a.ref), n18.load_sites(a.indels))

    # Ruler: the MAPQ >= 20 votes must be validate's own.
    spike = {}
    for line in open(a.check_log):
        m = re.search(r"allele_freq SNP \S+:(\d+)-\d+ \(unknown\): (\d+) carry, (\d+) span", line)
        if m:
            spike[int(m.group(1)) + 1] = (int(m.group(2)), int(m.group(3)))
    rows = []
    same = 0
    for vcf_pos, kind, length, extension, reads in sites:
        votes = n18.fragment_votes((n, v) for n, v, ok, _ in reads if ok)
        c, sp = votes.count("C"), votes.count("S")
        same += spike.get(vcf_pos) == (c, sp)
        rows.append((vcf_pos, kind, length, extension, reads, c, sp, n18.grade(c, c + sp)))
    print(f"ruler\tcounts match validate at {same}/{len(sites)} sites ({100 * same / len(sites):.2f}%)")
    if same / len(sites) < 0.99:
        print("RULER CHECK FAILED: nothing below is to be read")
        return 1

    def klass(ext):
        return "unique" if ext == 0 else ("short repeat" if ext < 10 else "long repeat")

    # R1
    per = defaultdict(Counter)
    for _, kind, length, ext, _, c, sp, g in rows:
        k = per[klass(ext)]
        k["sites"] += 1
        if g in ("pass", "fail"):
            k["evaluable"] += 1
            if g == "fail":
                k["fail"] += 1
                k["low" if c / (c + sp) < 0.5 else "high"] += 1
    for name in ("unique", "short repeat", "long repeat"):
        k = per[name]
        rate = 100 * k["fail"] / k["evaluable"] if k["evaluable"] else float("nan")
        print(f"R1\t{name}\tsites={k['sites']}\tevaluable={k['evaluable']}\tout_of_range={k['fail']} ({rate:.2f}%)\ttoo_low={k['low']}\ttoo_high={k['high']}")

    # R2
    tally = defaultdict(Counter)
    for *_, reads, c, sp, g in rows:
        if g not in ("pass", "fail"):
            continue
        for _, v, ok, other in reads:
            if ok:
                t = tally[g]
                t["reads"] += 1
                t["other"] += other
                t[f"{v}_reads"] += 1
                t[f"{v}_other"] += other
    for g in ("fail", "pass"):
        t = tally[g]
        print(
            f"R2\t{g}ing sites\tvoting reads={t['reads']}\twith another gap={t['other']} ({100 * t['other'] / t['reads']:.2f}%)"
            f"\tcarriers {100 * t['C_other'] / max(t['C_reads'], 1):.2f}%\tspanning {100 * t['S_other'] / max(t['S_reads'], 1):.2f}%"
        )

    # R3
    hi_c = hi_s = lo_c = lo_s = 0
    for *_, reads, c, sp, g in rows:
        hv = n18.fragment_votes((n, v) for n, v, ok, _ in reads if ok)
        lv = n18.fragment_votes((n, v) for n, v, ok, _ in reads if not ok)
        hi_c += hv.count("C"); hi_s += hv.count("S")
        lo_c += lv.count("C"); lo_s += lv.count("S")
    print(f"R3\tMAPQ>=20: carriers {hi_c} of {hi_c + hi_s} ({hi_c / (hi_c + hi_s):.3f})\tMAPQ<20: carriers {lo_c} of {lo_c + lo_s} ({lo_c / max(lo_c + lo_s, 1):.3f})")

    # R4
    bins = Counter()
    for *_, c, sp, g in rows:
        if g == "fail":
            f = c / (c + sp)
            bins["<0.2" if f < 0.2 else "0.2-0.35" if f < 0.35 else "0.35-0.65" if f <= 0.65 else "0.65-0.8" if f <= 0.8 else ">0.8"] += 1
    print("R4\t" + "\t".join(f"{b}={bins[b]}" for b in ("<0.2", "0.2-0.35", "0.35-0.65", "0.65-0.8", ">0.8")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
