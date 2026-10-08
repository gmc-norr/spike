import csv, sys
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
g = defaultdict(list)
for r in rows:
    if r["value"] and r["v"]:
        g[(r["set"], r["group"], r["which"], r["metric"])].append((float(r["value"]), float(r["v"])))
for k in sorted(g):
    xs = g[k]; n = len(xs)
    if n < 3: continue
    m = sum(x for x, _ in xs)/n
    var = sum((x-m)**2 for x, _ in xs)/(n-1)
    mv = sum(v for _, v in xs)/n
    if k[1] in ("SNV","DEL5-19","INS20-49","DUP50-299") and k[3] in ("A","J"):
        print(*k, n, round(var/mv, 2))
