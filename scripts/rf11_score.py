#!/usr/bin/env python3
"""RF11 plan: apply the locked decision rule to K1 and K2 (see "RF11" in docs/review/REVIEW.md).

Usage: rf11_score.py NULL_COUNTS_TSV K1_OUT_DIR SITES_TSV

NULL_COUNTS_TSV is `rf11_count.py sites` over the unedited BAM; K1_OUT_DIR is
`rf11_k1.sh`'s output. Prints every number the rule uses, then the verdict.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf11_count import count, per_read  # noqa: E402

LENGTHS = (20, 30, 40, 45, 49)
FLOORS = (20, 30, 40)
TODAY = 50
MIN_INS_READS = 2


def null_passes(rows, set_name, column):
    return sum(1 for r in rows if r["set"] == set_name and int(r[column]) >= MIN_INS_READS)


def main(null_path, k1_dir, sites_path):
    rows = list(csv.DictReader(open(null_path), delimiter="\t"))
    sets = ("null-rand", "null-rep")
    n_sites = {s: sum(1 for r in rows if r["set"] == s) for s in sets}

    # K2: null passes. Today's worst case is a 50 bp insertion (I or clip >= 50);
    # the fix at floor F is worst at an F bp insertion (I or clip >= F).
    print("K2, null sites passing (>= 2 reads):")
    baseline = {s: null_passes(rows, s, f"IS{TODAY}") for s in sets}
    for s in sets:
        today_short = {L: null_passes(rows, s, f"I{L}") for L in LENGTHS}
        print(f"  {s} ({n_sites[s]} sites): today at 50 bp {baseline[s]}; "
              f"today I-only at {today_short}")
    qualifying = []
    for F in FLOORS:
        fix = {s: null_passes(rows, s, f"IS{F}") for s in sets}
        ok = all(fix[s] <= max(2, baseline[s]) for s in sets)
        print(f"  floor {F}: fix passes {fix} vs allowed max(2, today at 50) -> "
              f"{'qualifies' if ok else 'does not qualify'}")
        if ok:
            qualifying.append(F)
    f_star = min(qualifying) if qualifying else None
    print(f"  F* = {f_star}")

    # K1: the correct spike-ins.
    k1 = []
    for line in open(sites_path):
        name, chrom, pos = line.split()
        if name != "k1":
            continue
        for L in LENGTHS:
            run = os.path.join(k1_dir, f"ins_{pos}_{L}")
            report = os.path.join(run, "validate.json")
            bam = os.path.join(run, "run", "merged.bam")
            if not (os.path.exists(report) and os.path.exists(bam)):
                k1.append({"pos": pos, "len": L, "ran": False})
                continue
            checks = json.load(open(report))["checks"]
            row = next(c for c in checks if c["check"] == "ins_reads")
            reads = per_read(bam, chrom, int(pos))
            k1.append({
                "pos": pos, "len": L, "ran": True,
                "observed": int(row["observed"]), "pass": row["pass"],
                "replica": count(reads, L, TODAY),
                "fix": {F: count(reads, L, F) for F in FLOORS},
            })

    ran = [r for r in k1 if r["ran"]]
    print(f"\nK1: {len(ran)} of {len(k1)} runs reached validate")
    for r in k1:
        if not r["ran"]:
            print(f"  ins:{r['pos']}:{r['len']}  did not reach validate")
    mismatches = [r for r in ran if r["replica"] != r["observed"]]
    print(f"  replica vs validate's observed: {len(ran) - len(mismatches)} of {len(ran)} equal")
    for r in mismatches:
        print(f"    MISMATCH ins:{r['pos']}:{r['len']} validate {r['observed']} replica {r['replica']}")
    print("  today's ins_reads, failing runs by length:")
    for L in LENGTHS:
        at = [r for r in ran if r["len"] == L]
        print(f"    {L} bp: {sum(1 for r in at if not r['pass'])} of {len(at)} fail   "
              f"observed {[r['observed'] for r in at]}")
    failing_today = sum(1 for r in ran if not r["pass"])
    if f_star is not None:
        eligible = [r for r in ran if r["len"] >= f_star]
        fixed = sum(1 for r in eligible if r["fix"][f_star] >= MIN_INS_READS)
        print(f"  fix at F*={f_star}: {fixed} of {len(eligible)} runs of >= {f_star} bp pass")
        for L in LENGTHS:
            at = [r for r in ran if r["len"] == L]
            print(f"    {L} bp: fix counts {[r['fix'][f_star] for r in at]}")

    # The locked rule.
    print("\nVerdict:")
    if len(ran) < 36:
        print(f"  INCONCLUSIVE: only {len(ran)} of 40 runs reached validate (36 needed)")
    elif mismatches:
        print("  INCONCLUSIVE: the replica does not reproduce validate's count")
    elif failing_today < 3:
        print(f"  INCONCLUSIVE: today's ins_reads fails {failing_today} of {len(ran)} (3 needed to act)")
    elif f_star is None:
        print("  REFUTED: no floor keeps the null sites within the allowance")
    elif fixed < 0.95 * len(eligible):
        print(f"  REFUTED: the fix passes {fixed} of {len(eligible)}, under 95%")
    else:
        print(f"  SUPPORTED: floor {f_star}; today fails {failing_today} of {len(ran)}, "
              f"the fix passes {fixed} of {len(eligible)}")


if __name__ == "__main__":
    main(*sys.argv[1:])
