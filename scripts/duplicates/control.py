"""C1 of the duplicates check (docs/superpowers/plans/2026-10-04-duplicates.md): a duplicate set whose
kept pair is removed has a member that is no longer flagged duplicate after marking again.

Usage: control.py pick  --bam TAGGED_BAM --region chr:start-end --seed 7 --sets OUT.tsv --drop OUT.txt
       control.py count --bam TAGGED_BAM --region chr:start-end --sets SETS.tsv
"""
import argparse
import random
import sys

import pysam

import measure

N_SETS, PASS_REMOVED, PASS_KEPT = 200, 95, 5


def pick_sets(bam, chrom, start, end, seed, n):
    """`n` duplicate sets (Picard's DI tag) with exactly one non-duplicate pair, at least one
    duplicate pair, and only proper primary pairs: [{"rep": name, "dups": [names]}]."""
    sets = {}
    for r in bam.fetch(chrom, start, end):
        if r.flag & (0x100 | 0x800) or not r.has_tag("DI"):
            continue
        s = sets.setdefault(r.get_tag("DI"), {"rep": set(), "dups": set(), "bad": False})
        if not r.is_proper_pair or r.is_unmapped:
            s["bad"] = True
        (s["dups"] if r.is_duplicate else s["rep"]).add(r.query_name)
    usable = sorted(di for di, s in sets.items() if not s["bad"] and len(s["rep"]) == 1 and s["dups"])
    if len(usable) < n:
        sys.exit(f"control: {len(usable)} usable duplicate sets, {n} wanted")
    return [{"rep": next(iter(sets[di]["rep"])), "dups": sorted(sets[di]["dups"])}
            for di in random.Random(seed).sample(usable, n)]


def freed(bam, chrom, start, end, sets):
    """How many sets have a duplicate member with a primary record that is not flagged duplicate."""
    unflagged = {r.query_name for r in bam.fetch(chrom, start, end)
                 if not (r.flag & (0x4 | 0x100 | 0x800 | 0x400))}
    return sum(any(d in unflagged for d in s["dups"]) for s in sets)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("pick", "count"))
    for a in ("bam", "region", "sets"):
        ap.add_argument(f"--{a}", required=True)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--drop")
    a = ap.parse_args()
    chrom, start, end = measure.region(a.region)
    bam = pysam.AlignmentFile(a.bam)
    if a.mode == "pick":
        sets = pick_sets(bam, chrom, start - 1, end, a.seed, N_SETS)
        with open(a.sets, "w") as f:
            for i, s in enumerate(sets):
                f.write(f"{'removed' if i < N_SETS // 2 else 'kept'}\t{s['rep']}\t{','.join(s['dups'])}\n")
        with open(a.drop, "w") as f:
            f.write("".join(f"{s['rep']}\n" for s in sets[:N_SETS // 2]))
        print(f"control: picked {len(sets)} sets; dropping the kept pair of the first {N_SETS // 2}", flush=True)
        return
    groups = {"removed": [], "kept": []}
    for line in open(a.sets):
        tag, rep, dups = line.rstrip("\n").split("\t")
        groups[tag].append({"rep": rep, "dups": dups.split(",")})
    removed = freed(bam, chrom, start - 1, end, groups["removed"])
    kept = freed(bam, chrom, start - 1, end, groups["kept"])
    ok = removed >= PASS_REMOVED and kept <= PASS_KEPT
    print(f"freed sets: {removed} of {len(groups['removed'])} with the kept pair removed (>= {PASS_REMOVED}), "
          f"{kept} of {len(groups['kept'])} untouched (<= {PASS_KEPT}); C1: {'PASS' if ok else 'FAIL'}", flush=True)


if __name__ == "__main__":
    main()
