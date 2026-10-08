import csv, sys, random
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
d = defaultdict(dict)
for r in rows:
    if r["value"]:
        d[(r["set"], r["group"], r["which"], r["metric"])][(r["chrom"], r["pos"], r["ref"], r["alt"])] = float(r["value"])
def pct(xs, p):
    xs = sorted(xs); k = (len(xs)-1)*p/100; lo = int(k); hi = min(lo+1, len(xs)-1)
    return xs[lo] + (xs[hi]-xs[lo])*(k-lo)
rng = random.Random(3)
for s in ("forward", "reverse"):
    for (s2, g, w, m), fk in sorted(d.items()):
        if s2 != s or w != "fake_normal": continue
        real = d[(s, g, "real", m)]
        fr = [(fk[e], real[e]) for e in fk if e in real]
        # forward B1 shift per unit half-dose (forward controls, as round 3 does for reverse)
        fb = d[("forward", g, "fake_B1", m)]; fn = d[("forward", g, "fake_normal", m)]
        sh = [fb[e] - fn[e] for e in fb if e in fn]
        ests = []
        for _ in range(2000):
            a = [rng.choice(fr) for _ in fr]; b = [rng.choice(sh) for _ in sh]
            md = sum(x - y for x, y in a)/len(a); ms = sum(b)/len(b)
            ests.append(1 + 0.5*md/ms * -1 if False else 1 - 0.5*md/ms)
        md = sum(x - y for x, y in fr)/len(fr); ms = sum(sh)/len(sh)
        # effective dose: 1 + 0.5 * md / (-ms)  (ms negative = loss per half dose)
        eff = 1 + 0.5*md/(-ms)
        e2 = [1 + (x - 1)*-1 for x in ests]  # fix sign below
        boots = []
        rng2 = random.Random(5)
        for _ in range(2000):
            a = [rng2.choice(fr) for _ in fr]; b = [rng2.choice(sh) for _ in sh]
            mdb = sum(x - y for x, y in a)/len(a); msb = sum(b)/len(b)
            boots.append(1 + 0.5*mdb/(-msb))
        print(s, g, m, "mean_d", round(md, 4), "B1 shift", round(ms, 4), "eff_dose", round(eff, 3), "90%CI", round(pct(boots, 5), 3), round(pct(boots, 95), 3))
