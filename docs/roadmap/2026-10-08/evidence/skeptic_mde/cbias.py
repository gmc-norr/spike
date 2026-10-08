import csv, sys, random, statistics as st
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
D = defaultdict(dict)
for r in rows:
    if r["value"] == "" : continue
    k = (r["chrom"], r["pos"], r["ref"], r["alt"])
    D[(r["set"], r["group"], r["metric"])].setdefault(k, {})[r["which"]] = float(r["value"])
rng = random.Random(11)
def ci(d):
    bs = sorted(st.mean(rng.choice(d) for _ in d) for _ in range(2000)); return bs[50], bs[1949]
def pct(xs,p):
    xs=sorted(xs); return xs[int(p*(len(xs)-1))]
print("group metric | shared mean(HG002-HG001) [CI] | fwd d (HG002bg fake - HG001 real) | rev d (HG001bg fake - HG002 real) | fwd d - c | rev d + c")
for (s,g,m), ev in sorted(D.items()):
    if s != "shared": continue
    c = [x["HG002"]-x["HG001"] for x in ev.values() if "HG001" in x and "HG002" in x]
    f = [x["fake_normal"]-x["real"] for x in D[("forward",g,m)].values() if "fake_normal" in x and "real" in x]
    r = [x["fake_normal"]-x["real"] for x in D[("reverse",g,m)].values() if "fake_normal" in x and "real" in x]
    lo,hi = ci(c)
    print(f"{g} {m} | {st.mean(c):+.4f} [{lo:+.4f},{hi:+.4f}] n={len(c)} | {st.mean(f):+.4f} | {st.mean(r):+.4f} | {st.mean(f)-st.mean(c):+.4f} | {st.mean(r)+st.mean(c):+.4f}")
# across-site variance / mean v (R31 check)
print()
D2 = defaultdict(lambda: ([],[]))
for r in rows:
    if r["value"]=="" or r["v"]=="": continue
    if r["which"] in ("real","HG001","HG002"):
        D2[(r["set"],r["group"],r["metric"],r["which"])][0].append(float(r["value"])); D2[(r["set"],r["group"],r["metric"],r["which"])][1].append(float(r["v"]))
for k,(x,v) in sorted(D2.items()):
    if k[1] in ("SNV","DEL5-19","INS20-49","DUP50-299"):
        print(k, f"var_across/mean_v4 = {st.pvariance(x)/(4*st.mean(v)) if st.mean(v)>0 else float('nan'):.2f}", f"var_across/mean_v = {st.pvariance(x)/st.mean(v):.2f}")
