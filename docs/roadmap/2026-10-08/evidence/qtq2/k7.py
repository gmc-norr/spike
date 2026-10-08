"""K7 of docs/superpowers/plans/2026-10-07-quality-model-v2.md. For each round trip (het, m4, hom, snv):
allowance = the most allele_freq rows whose pass/fail differs between any two of master's seeds
42, 2, 3, 4, 5; pass = every row passing on master seed 42 passes on v2 seed 42, except allele_freq rows
that fail on v2, up to the allowance. Rows: `spike validate` text, fixed columns.
Usage: k7.py (from qtq2)"""
import itertools, os
def rows(path):
    out = {}
    for l in open(path):
        if len(l) < 98 or l.startswith(("Event", "---")) or not l[97:].strip():
            continue
        key = (l[:36].strip(), l[36:55].strip())
        out[key] = l[97:].strip().startswith("PASS")
    return out
seeds = [42, 2, 3, 4, 5]
verdict = True
for run in ["het", "m4", "hom", "snv"]:
    m = {s: rows(f"master_s{s}/{run}.validate.txt") for s in seeds}
    v = rows(f"v2_s42/{run}.validate.txt")
    af = lambda k: k[1] == "allele_freq"
    allow = max(sum(1 for k in m[a] if af(k) and k in m[b] and m[a][k] != m[b][k]) for a, b in itertools.combinations(seeds, 2))
    lost = [k for k in m[42] if m[42][k] and not v.get(k, False)]
    lost_af = [k for k in lost if af(k)]
    lost_other = [k for k in lost if not af(k)]
    ok = not lost_other and len(lost_af) <= allow
    verdict &= ok
    print(f"{run}: master s42 {sum(m[42].values())}/{len(m[42])}, v2 {sum(v.values())}/{len(v)}; "
          f"master seeds pass {[sum(m[s].values()) for s in seeds]}; allele_freq allowance {allow}; "
          f"lost allele_freq {len(lost_af)}, lost other {len(lost_other)} -> {'PASS' if ok else 'FAIL'}")
    for k in lost:
        print("   lost:", k)
print("K7:", "PASS" if verdict else "FAIL")
