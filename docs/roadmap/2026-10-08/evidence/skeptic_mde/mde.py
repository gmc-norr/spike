import csv, sys, random, statistics as st
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
D = defaultdict(dict)  # (set,group,metric) -> {key: {which:(val,v)}}
for r in rows:
    if r["value"] == "" or r["v"] == "": continue
    k = (r["chrom"], r["pos"], r["ref"], r["alt"])
    D[(r["set"], r["group"], r["metric"])].setdefault(k, {})[r["which"]] = (float(r["value"]), float(r["v"]))
def R(pairs):
    num = sum((a-b)**2 for a,b,va,vb in pairs); den = sum(va+vb for a,b,va,vb in pairs)
    return num/den if den>0 else float('nan')
def pct(xs, p):
    xs = sorted(xs); i = p*(len(xs)-1); lo=int(i); hi=min(lo+1,len(xs)-1); return xs[lo]+(xs[hi]-xs[lo])*(i-lo)
rng = random.Random(7)
print("set group metric | n_f n_c | mean_d(normal) [95%CI] mean_real | rel bias % [CI] | B1 mean_d | dose-equiv of normal bias (VAF) | interp: VAF where rho=2.25 ; VAF where P(rho5>2.25)>=0.8 (sim)")
for (s,g,m), ev in sorted(D.items()):
    if s=="shared": continue
    sh = D[("shared", g, m)]
    c = [(x["HG001"][0], x["HG002"][0], x["HG001"][1], x["HG002"][1]) for x in sh.values() if "HG001" in x and "HG002" in x]
    fn = [(x["fake_normal"][0], x["real"][0], x["fake_normal"][1], x["real"][1]) for x in ev.values() if "fake_normal" in x and "real" in x]
    fb = [(x["fake_B1"][0], x["real"][0], x["fake_B1"][1], x["real"][1]) for x in ev.values() if "fake_B1" in x and "real" in x]
    d = [a-b for a,b,_,_ in fn]; mr = st.mean(b for a,b,_,_ in fn)
    bs = []
    for _ in range(2000):
        smp = [rng.choice(d) for _ in d]; bs.append(st.mean(smp))
    md = st.mean(d); lo, hi = pct(bs,.025), pct(bs,.975)
    line = f"{s} {g} {m} | {len(fn)} {len(c)} | {md:+.4f} [{lo:+.4f},{hi:+.4f}] {mr:.3f} | {100*md/mr:+.1f}% [{100*lo/mr:+.1f},{100*hi/mr:+.1f}]"
    if fb:
        db = st.mean(a-b for a,b,_,_ in fb)
        # linear dose model: fake value at dose t (t=0 normal, t=1 B1 half dose): x_t = x_n + t*(x_b - x_n)
        # simulate rho at t via pairing the same events, adding noise: use per-event mixture
        Rc = R(c)
        def sim_rho(t, draws=400):
            out=[]
            for _ in range(draws):
                fi = [rng.randrange(len(fn)) for _ in fn]; ci = [rng.randrange(len(c)) for _ in c]
                num=den=0.0
                key = list(ev.values())
                pairs=[]
                for i in fi:
                    pass
                out.append(None)
            return out
        # analytic point: R_f(t) = R_n + t^2*(R_b - R_n)  (bias^2 scaling, noise unchanged)
        Rn, Rb = R(fn), R(fb)
        tstar = ((2.25*Rc - Rn)/(Rb-Rn))**0.5 if Rb>Rn and 2.25*Rc>Rn else float('nan')
        line += f" | B1 d {db:+.4f} | normal bias = {md/db:+.2f} of half-dose -> VAF {0.5-0.25*md/db:.3f} | rho=2.25 at t={tstar:.2f} (VAF {0.5-0.25*tstar:.3f}); Rn {Rn:.3f} Rb {Rb:.3f} Rc {Rc:.3f}"
    print(line)
