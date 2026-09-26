#!/usr/bin/env python3
"""RF14 second attempt: score K as the plan locks it.

The row is rf14_planted.py's, unchanged. For every run under K_DIR (rf14b_k.sh's
layout) this reads the truth record, its SIM_RESIST, master's `split_reads` row
from the run's validate.json, and counts carriers for the own truth and for the
negatives N1-N4. Then it prints each locked criterion and its verdict:

  K+     VAF 0.5 runs (r_*, h_*): at least 32 of 36 reach validate, at least 24 of
         those are at SIM_RESIST <= 0.5, and all but at most 1 of those carry.
  LOW    VAF 0.1 runs (l_*) at SIM_RESIST <= 0.5 that reach validate, at least 20
         of them: the row fails fewer than master's split_reads does. The row's
         failure count and its 95% Wilson interval are reported.
  GUARD  no run anywhere where master's split_reads passes and the row has 0.
  K-     0 carriers on every judged negative of every run that reaches validate.

Runs above SIM_RESIST 0.5 -- events spike refuses by default (RF8) and
--allow-resistant forced through -- are judged only by GUARD and K-; their carriers
are reported.

Usage: rf14b_score.py K_DIR REFERENCE
"""
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rf14_planted as row  # noqa: E402

RESIST_LINE = 0.5


def wilson(x, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = x / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, centre - half), min(1.0, centre + half))


def master_split_pass(run_dir):
    rows = json.load(open(os.path.join(run_dir, "validate.json")))["checks"]
    r = [c for c in rows if c["check"] == "split_reads"]
    return bool(r) and r[0]["pass"], (r[0]["observed"] if r else "-")


def main(kdir, ref):
    runs = []
    judged_neg, neg_fail, blind, unfiltered = 0, [], [], []
    for name in sorted(os.listdir(kdir)):
        d = os.path.join(kdir, name)
        if not os.path.isdir(d):
            continue
        truth = os.path.join(d, "run", "truth.vcf")
        merged = os.path.join(d, "run", "merged.bam")
        if not (os.path.exists(truth) and os.path.exists(merged)):
            runs.append({"name": name, "reached": False})
            continue
        chrom, start, end, rec_id, f = next(row.del_records(truth))
        resist = float(re.search(r"SIM_RESIST=([0-9.]+)", f[7]).group(1))
        c = row.count(row.carriers(merged, ref, chrom, start, end, rec_id))
        sr_pass, sr_obs = master_split_pass(d)
        runs.append({"name": name, "reached": True, "resist": resist, "carriers": c,
                     "sr_pass": sr_pass, "sr_obs": sr_obs})
        n_id = int(row.ID_RE.match(rec_id).group(1))
        negs = {
            "N1 empty": (os.path.join(d, "slice.bam"), start, end, rec_id, False),
            "N2a END+50": (merged, start, end + 50, rec_id, True),
            "N2b START-50": (merged, start - 50, end, rec_id, True),
            "N3 moved +1000": (merged, start + 1000, end + 1000, rec_id, True),
            "N4 wrong id": (merged, start, end, f"sim_del_{n_id + 1}", False),
        }
        for label, (bam, s, e, i, can_blind) in negs.items():
            n = row.count(row.carriers(bam, ref, chrom, s, e, i))
            if can_blind and row.indistinguishable(ref, chrom, start, end, s, e):
                blind.append((name, label, n))
                continue
            judged_neg += 1
            if n >= 1:
                neg_fail.append((name, label, n))
        if name.startswith("h_"):
            u = row.carriers(os.path.join(d, "slice.bam"), ref, chrom, start, end, rec_id,
                             any_name=True)
            unfiltered.append((name, row.count(u)))

    for r in runs:
        if r["reached"]:
            print(f"  {r['name']}: carriers {r['carriers']}  resist {r['resist']:.3f}  "
                  f"master split_reads {'PASS' if r['sr_pass'] else 'FAIL'} {r['sr_obs']}")
        else:
            print(f"  {r['name']}: did not reach validate")

    half = [r for r in runs if r["name"][0] in "rh"]
    half_reached = [r for r in half if r["reached"]]
    half_judged = [r for r in half_reached if r["resist"] <= RESIST_LINE]
    half_miss = [r["name"] for r in half_judged if r["carriers"] < 1]
    k_plus = len(half_reached) >= 32 and len(half_judged) >= 24 and len(half_miss) <= 1
    print(f"K+ VAF 0.5: {len(half_reached)} of {len(half)} reach validate; "
          f"{len(half_judged)} at SIM_RESIST <= 0.5; carrying {len(half_judged) - len(half_miss)} "
          f"of {len(half_judged)} {half_miss} -> {'PASS' if k_plus else 'FAIL'}")

    low = [r for r in runs if r["name"].startswith("l_")]
    low_judged = [r for r in low if r["reached"] and r["resist"] <= RESIST_LINE]
    row_fail = sum(1 for r in low_judged if r["carriers"] < 1)
    sr_fail = sum(1 for r in low_judged if not r["sr_pass"])
    lo, hi = wilson(row_fail, len(low_judged))
    if len(low_judged) < 20:
        low_verdict = "INCONCLUSIVE"
    else:
        low_verdict = "PASS" if row_fail < sr_fail else "FAIL"
    print(f"LOW VAF 0.1: {len(low_judged)} of {len(low)} judged; row fails {row_fail} "
          f"(95% {lo:.3f}-{hi:.3f}), master split_reads fails {sr_fail} -> {low_verdict}")
    print(f"    carriers: {sorted(r['carriers'] for r in low_judged)}")

    guard = [r["name"] for r in runs if r["reached"] and r["sr_pass"] and r["carriers"] < 1]
    print(f"GUARD master split_reads PASS but row 0: {len(guard)} {guard} -> "
          f"{'PASS' if not guard else 'FAIL'}")

    high = [(r["name"], r["resist"], r["carriers"]) for r in runs
            if r["reached"] and r["resist"] > RESIST_LINE]
    print(f"above SIM_RESIST 0.5 (reported): {len(high)} {high}")
    print(f"K- judged negatives with a carrier: {len(neg_fail)} of {judged_neg} {neg_fail} -> "
          f"{'PASS' if not neg_fail else 'FAIL'}")
    print(f"indistinguishable negatives (not judged): {len(blind)} {blind}")
    print(f"check of the check, h_ sites' unspiked slice with the name filter off: {unfiltered}")


if __name__ == "__main__":
    main(*sys.argv[1:])
