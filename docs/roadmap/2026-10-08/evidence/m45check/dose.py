import csv, random, sys
from collections import defaultdict
path = sys.argv[1]
rows = list(csv.DictReader(open(path), delimiter="\t"))
def f(x): return None if x in ("", None) else float(x)
data = defaultdict(dict)  # (set, group, which, metric) -> {event: (value, v)}
for r in rows:
    ev = (r["kind"], r["chrom"], r["pos"], r["ref"], r["alt"])
    data[(r["set"], r["group"], r["which"], r["metric"])][ev] = (f(r["value"]), f(r["v"]))
def items(fake, real):
    out = []
    for e, (x, vx) in fake.items():
        if e in real:
            y, vy = real[e]
            if None not in (x, vx, y, vy):
                out.append((e, x - y, vx + vy))
    return out
def R(it):
    n = sum(v for _, v in it); return sum(d*d for d, _ in it)/n if n > 0 else None
def pct(xs, p):
    xs = sorted(xs); k = (len(xs)-1)*p/100; lo = int(k); hi = min(lo+1, len(xs)-1)
    return xs[lo] + (xs[hi]-xs[lo])*(k-lo)
def boot(fi, ci, draws=2000, seed=3):
    rng = random.Random(seed); out = []
    for _ in range(draws):
        rf = R([rng.choice(fi) for _ in fi]); rc = R([rng.choice(ci) for _ in ci])
        if rf is not None and rc: out.append(rf/rc)
    return pct(out, 5), pct(out, 95)
cells = [("SNV","A"),("DEL5-19","A"),("DEL5-19","E"),("INS5-19","A"),("INS5-19","E"),("DEL1-4","A"),("DEL1-4","E"),("DEL20-49","A"),("DEL20-49","E"),("INS1-4","A"),("INS1-4","E"),("INS20-49","A"),("INS20-49","E"),("DUP50-299","J")]
print("group metric vaf mode rho rho5 rho95 status")
for g, m in cells:
    real = data[("forward", g, "real", m)]
    fn = data[("forward", g, "fake_normal", m)]
    fb = data[("forward", g, "fake_B1", m)]
    sh1 = data[("shared", g, "HG001", m)]; sh2 = data[("shared", g, "HG002", m)]
    ci = [(d, v) for _, d, v in items(sh1, sh2)]
    n_it = {e: (d, v) for e, d, v in items(fn, real)}
    b_it = {e: (d, v) for e, d, v in items(fb, real)}
    common = [e for e in n_it if e in b_it]
    meanshift = sum(b_it[e][0]-n_it[e][0] for e in common)/len(common)
    for vaf in (0.5, 0.45, 0.40, 0.35, 0.30, 0.25):
        t = (0.5 - vaf)/0.25
        for mode in ("mean", "event"):
            if mode == "mean":
                fi = [(n_it[e][0] + t*meanshift, n_it[e][1]) for e in common]
            else:
                fi = [(n_it[e][0] + t*(b_it[e][0]-n_it[e][0]), n_it[e][1]) for e in common]
            rho = R(fi)/R(ci); lo, hi = boot(fi, ci, 1000)
            st = "pass" if hi <= 2.25 else ("fail" if lo > 2.25 else "unsure")
            print(g, m, vaf, mode, round(rho,2), round(lo,2), round(hi,2), st)
    print(g, m, "meanshift", round(meanshift,4))
