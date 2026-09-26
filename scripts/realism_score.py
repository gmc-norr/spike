#!/usr/bin/env python3
"""Realism probe, scoring: today's checks on real HG002 variants, beside spike's.

Usage: realism_score.py PROBE_DIR BAM REFERENCE RF6_PIPELINE_RUNS RF11_K1_DIR RF11_SITES

PROBE_DIR holds realism_probe.py's outputs and `spike validate --json` on HG002's
own BAM for each truth VCF: real_ins.json, real_del.json and pipeline20.json (the
pipeline's 20 deletions, truth records taken from the RF6 runs). RF6_PIPELINE_RUNS
is the RF6 `real_events.sh` run over those 20 on NA18488 (one dir per event, each
with validate.json). RF11_K1_DIR and RF11_SITES are RF11's spike-in runs.
"""
import collections
import csv
import json
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf11_count import per_read  # noqa: E402


def fetch(ref_path, start, end):
    """chr20 bases [start, end], 1-based inclusive, uppercase."""
    out = subprocess.run(
        ["samtools", "faidx", ref_path, f"chr20:{start}-{end}"],
        capture_output=True, text=True, check=True,
    ).stdout
    return "".join(out.split("\n")[1:]).upper()


def rows(path):
    d = collections.defaultdict(dict)
    for c in json.load(open(path))["checks"]:
        m = re.search(r"chr20:(\d+)", c["event"])
        if m:
            d[m.group(1)][c["check"]] = c
    return d


def size_bin(n):
    n = int(n)
    return "20-29" if n < 30 else "30-39" if n < 40 else "40-49"


def main(probe, bam, ref_path, rf6_runs, k1_dir, k1_sites):
    events = {(r["class"], r["pos"]): r for r in csv.DictReader(open(f"{probe}/events.tsv"), delimiter="\t")}

    print("== Real HG002 insertions, 20-49 bp, in HG002's own BAM")
    ins = rows(f"{probe}/real_ins.json")
    tab = collections.defaultdict(lambda: collections.Counter())
    for pos, ch in ins.items():
        e = events[("ins", pos)]
        t = tab[(size_bin(e["len"]), e["zyg"])]
        t["n"] += 1
        t["ins_reads FAIL"] += not ch["ins_reads"]["pass"]
        s = ch.get("ins_sequence")
        if s is None:
            t["ins_sequence none"] += 1
        elif s["pass"]:
            t["ins_sequence PASS"] += 1
        elif "kmer in ref" in s["observed"] or "alt is" in s["observed"]:
            t["ins_sequence not evaluable"] += 1
        else:
            t["ins_sequence FAIL"] += 1
    for k in sorted(tab):
        print(f"  {k[0]} {k[1]}: {dict(tab[k])}")

    # How the aligner wrote the failing ones under 40 bp.
    alts = {}
    for line in open(f"{probe}/real_ins.vcf"):
        if not line.startswith("#"):
            f = line.split("\t")
            alts[f[1]] = len(f[4]) - 1
    fails = [p for p in ins if not ins[p]["ins_reads"]["pass"] and alts[p] < 40]
    survey = collections.Counter()
    for p in fails:
        reads = per_read(bam, "chr20", int(p))
        lengths = [i for i, _ in reads.values() if i > 0]
        survey["failing events under 40 bp"] += 1
        survey["with no I at all"] += not lengths
        survey["with an I in >= 2 reads, under 2 of them SVLEN long"] += len(lengths) >= 2
    print(f"  aligner's I operations at those that FAIL: {dict(survey)}")

    # Is the inserted sequence a copy of the reference right beside POS, or within
    # 2*SVLEN of it (a repeat expansion, left-normalised), or neither?
    seqs = {}
    for line in open(f"{probe}/real_ins.vcf"):
        if not line.startswith("#"):
            f = line.split("\t")
            seqs[f[1]] = f[4][1:].upper()
    copies = collections.Counter()
    for pos, s in seqs.items():
        p, n = int(pos), len(s)
        window = fetch(ref_path, p - 2 * n, p + 2 * n)
        beside = (window[2 * n + 1:3 * n + 1], window[n + 1:2 * n + 1])  # after, ending at POS
        kind = "exact copy beside POS" if s in beside else "copy within 2*SVLEN" if s in window else "not a nearby copy"
        copies[(kind, "n")] += 1
        copies[(kind, "ins_reads FAIL")] += not ins[pos]["ins_reads"]["pass"]
    for kind in ("exact copy beside POS", "copy within 2*SVLEN", "not a nearby copy"):
        print(f"  {kind}: ins_reads FAIL {copies[(kind, 'ins_reads FAIL')]} of {copies[(kind, 'n')]}")

    print("\n== Spike-ins at random unique sites (RF11 K1), same BAM, het")
    k1 = collections.Counter()
    for line in open(k1_sites):
        name, _, pos = line.split()
        if name != "k1":
            continue
        for length in (20, 30, 40, 45, 49):
            ch = rows(os.path.join(k1_dir, f"ins_{pos}_{length}", "validate.json"))[pos]
            b = size_bin(length)
            k1[(b, "n")] += 1
            k1[(b, "ins_reads FAIL")] += not ch["ins_reads"]["pass"]
    for b in ("20-29", "30-39", "40-49"):
        print(f"  {b}: ins_reads FAIL {k1[(b, 'ins_reads FAIL')]} of {k1[(b, 'n')]}")

    def fail_counts(d):
        c = collections.Counter()
        for ch in d.values():
            for name, v in ch.items():
                if not v["pass"]:
                    c[name] += 1
        return dict(sorted(c.items()))

    print("\n== Real HG002 deletions >= 500 bp, isolated, in HG002's own BAM")
    dels = rows(f"{probe}/real_del.json")
    print(f"  {len(dels)} events, FAIL per check: {fail_counts(dels)}")

    print("\n== The pipeline's 20 deletions: HG002's real reads vs spike's reads on NA18488 (RF6)")
    real = rows(f"{probe}/pipeline20.json")
    spike = collections.defaultdict(dict)
    for d in os.listdir(rf6_runs):
        path = os.path.join(rf6_runs, d, "validate.json")
        if os.path.exists(path):
            for pos, ch in rows(path).items():
                spike[pos].update(ch)
    print(f"  real  ({len(real)}): {fail_counts(real)}")
    print(f"  spike ({len(spike)}): {fail_counts(spike)}")
    for check in ("split_reads", "coverage_ratio"):
        same = sum(1 for p in real if p in spike and real[p][check]["pass"] == spike[p][check]["pass"])
        print(f"  {check}: same verdict on {same} of {len(real)}")


if __name__ == "__main__":
    main(*sys.argv[1:])
