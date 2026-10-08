import csv, sys
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
d = defaultdict(dict)
for r in rows:
    if r["value"]:
        d[(r["set"], r["group"], r["which"], r["metric"])][(r["chrom"], r["pos"], r["ref"], r["alt"])] = float(r["value"])
for (s, g, w, m), fk in sorted(d.items()):
    if w != "fake_normal": continue
    real = d[(s, g, "real", m)]
    diffs = [fk[e] - real[e] for e in fk if e in real]
    mr = sum(real[e] for e in fk if e in real)/len(diffs)
    print(s, g, m, len(diffs), "mean_d", round(sum(diffs)/len(diffs), 4), "mean_real", round(mr, 4))
