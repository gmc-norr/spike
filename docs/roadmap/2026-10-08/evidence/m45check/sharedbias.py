import csv, sys, random
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
d = defaultdict(dict)
for r in rows:
    if r["value"] and r["set"] == "shared":
        d[(r["group"], r["which"], r["metric"])][(r["chrom"], r["pos"], r["ref"], r["alt"])] = float(r["value"])
rng = random.Random(7)
def pct(xs, p):
    xs = sorted(xs); k = (len(xs)-1)*p/100; lo = int(k); hi = min(lo+1, len(xs)-1)
    return xs[lo] + (xs[hi]-xs[lo])*(k-lo)
for (g, w, m), a in sorted(d.items()):
    if w != "HG001": continue
    b = d[(g, "HG002", m)]
    diffs = [a[e] - b[e] for e in a if e in b]
    bs = []
    for _ in range(2000):
        s = [rng.choice(diffs) for _ in diffs]; bs.append(sum(s)/len(s))
    print(g, m, len(diffs), "mean(HG001-HG002)", round(sum(diffs)/len(diffs), 4), "90%CI", round(pct(bs,5),4), round(pct(bs,95),4))
