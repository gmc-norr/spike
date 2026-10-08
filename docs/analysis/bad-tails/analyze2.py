"""What happens in the bad tails? Step 3: is how wrong a tail is a property of the read that its
qualities do not show? Uses ALL crashed reads that have no clip of another cause (so clip or not;
grouping by clip would select on the errors themselves). Tail = last 40 cycles.
Usage: analyze2.py IN.pkl"""
import sys, pickle, collections, numpy as np

R = pickle.load(open(sys.argv[1], "rb"))
other = [r for r in R if r["crashed"] and r["anyclip"] and not r["badclip"]]
print("crashed reads with a clip of another cause (left out):", len(other), collections.Counter(r["cause"] for r in other).most_common())
cr = [r for r in R if r["crashed"] and (r["badclip"] or not r["anyclip"])]
print("crashed reads used:", len(cr), " of which bad-end clipped:", sum(r["badclip"] for r in cr))
T = 40

def tail(r, sel=None):
    m = r["counted"].copy(); m[:-T] = False
    if sel is not None:
        m &= sel
    e = (r["called"] != r["ref"]) & m
    return e.sum(), m.sum()

print("\n8a. HOW WRONG IS EACH CRASHED READ'S TAIL? (Q11 bases in the last 40 cycles)")
k = np.array([tail(r, r["q"] == 11) for r in cr], float)
keep = k[:, 1] >= 10
k = k[keep]; rs = [r for r, ok in zip(cr, keep) if ok]
rate = k[:, 0] / k[:, 1]; p = k[:, 0].sum() / k[:, 1].sum()
print(f"  reads with >= 10 such bases: {len(k)}; pooled error rate {p:.3f}")
hist = np.histogram(rate, bins=[0, .1, .2, .3, .4, .5, .6, .7, 1.01])[0]
print("  reads by tail error rate: " + "  ".join(f"{a:.1f}-{b:.1f}:{100*h/len(k):4.1f}%" for a, b, h in zip([0, .1, .2, .3, .4, .5, .6, .7], [.1, .2, .3, .4, .5, .6, .7, 1], hist)))
# the same if every read had the pooled rate (binomial, same base counts)
rng = np.random.default_rng(1)
sim = rng.binomial(k[:, 1].astype(int), p) / k[:, 1]
hs = np.histogram(sim, bins=[0, .1, .2, .3, .4, .5, .6, .7, 1.01])[0]
print("  if every read erred alike:  " + "  ".join(f"{a:.1f}-{b:.1f}:{100*h/len(k):4.1f}%" for a, b, h in zip([0, .1, .2, .3, .4, .5, .6, .7], [.1, .2, .3, .4, .5, .6, .7, 1], hs)))
print(f"  spread of the per-read rate: SD {rate.std():.3f} against {sim.std():.3f} if every read erred alike")

print("\n8b. DOES ANYTHING THE QUALITIES SHOW PREDICT IT? (mean tail Q11 error rate per bin; SD within the bin; binomial SD)")
def by(label, vals, edges):
    vals = np.asarray(vals)
    print(f"  by {label}:")
    for a, b in zip(edges[:-1], edges[1:]):
        m = (vals >= a) & (vals < b)
        if m.sum() < 30:
            continue
        pk = k[m, 0].sum() / k[m, 1].sum()
        sim = rng.binomial(k[m, 1].astype(int), pk) / k[m, 1]
        print(f"    [{a:g}, {b:g}): reads {m.sum():5d}  rate {pk:.3f}  SD {rate[m].std():.3f}  (binomial {sim.std():.3f})")
by("Q11 count in the last 40", [(r["q"][-T:] == 11).sum() for r in rs], [0, 15, 20, 25, 30, 41])
by("read mean quality", [r["q"].mean() for r in rs], [0, 26, 28, 30, 32, 40])
by("Q37 count in the last 40", [(r["q"][-T:] == 37).sum() for r in rs], [0, 5, 10, 15, 41])
by("first cycle of the last 40 with Q11 (from start)", [151 - T + int(np.argmax(r["q"][-T:] == 11)) for r in rs], [0, 115, 120, 125, 152])
by("Q11 count in the first 100 cycles", [(r["q"][:100] == 11).sum() for r in rs], [0, 5, 10, 20, 40, 101])
mates = collections.defaultdict(dict)
for r in R:
    mates[r["name"]][r["r1"]] = r
mc = [int(mates[r["name"]].get(not r["r1"], {"crashed": None})["crashed"] or 0) if (not r["r1"]) in mates[r["name"]] else -1 for r in rs]
by("mate crashed (1) or not (0); -1 mate not kept", mc, [-1, 0, 1, 2])

print("\n8c. IS IT THE SAME READ ALL ALONG? (crashed reads: Q11 error rate in cycles 111-130 against 131-150, per read)")
a = np.array([tail(r, (r["q"] == 11) & (np.arange(151) < 131)) for r in cr], float)
b = np.array([tail(r, (r["q"] == 11) & (np.arange(151) >= 131)) for r in cr], float)
ok = (a[:, 1] >= 5) & (b[:, 1] >= 5)
ra, rb = a[ok, 0] / a[ok, 1], b[ok, 0] / b[ok, 1]
print(f"  reads {ok.sum()}: correlation of the two halves' rates {np.corrcoef(ra, rb)[0,1]:.2f}")
mt = []
for r in cr:
    o = mates[r["name"]].get(not r["r1"])
    if o is None or not o["crashed"] or not r["r1"]:
        continue
    x, y = tail(r, r["q"] == 11), tail(o, o["q"] == 11)
    if x[1] >= 10 and y[1] >= 10:
        mt.append((x[0] / x[1], y[0] / y[1]))
mt = np.array(mt)
if len(mt) > 10:
    print(f"  pairs where both mates crashed and were kept: {len(mt)}; correlation of R1's and R2's tail rates {np.corrcoef(mt[:,0], mt[:,1])[0,1]:.2f}")

print("\n8d. IN A WRONG TAIL, ARE EVEN THE HIGH QUALITIES WRONG? (last 40 cycles, by the read's Q11 tail rate)")
for lo_, hi_ in ((0, .2), (.2, .4), (.4, 1.01)):
    sel = [r for r, x in zip(rs, rate) if lo_ <= x < hi_]
    out = []
    for q in (11, 25, 37):
        e = n = 0
        for r in sel:
            x = tail(r, r["q"] == q); e += x[0]; n += x[1]
        out.append(f"Q{q} {e/max(1,n):.3f} ({n})")
    print(f"  reads whose Q11 tail rate is {lo_:.1f}-{hi_:.1f} ({len(sel)}): " + "  ".join(out))
