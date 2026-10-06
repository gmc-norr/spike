"""E3: quality-string models fitted on half the real read pairs, scored on the other half.
M1 first-order chain (cycle, mate, previous-Q bin) -- spike today, without the base term
M2 copy a whole real pair's two quality strings -- the parked quality-tails branch
M3 hidden read class: 4 classes by pair mean quality; one chain per class; both mates share the class
"""
import pickle, zlib, numpy as np
R = pickle.load(open("reads.pkl", "rb"))
pairs = {}
for d in R: pairs.setdefault(d["name"], {})[d["mate"]] = d["q"]
names = sorted(n for n, p in pairs.items() if len(p) == 2)
test = [n for n in names if zlib.crc32(n.encode()) % 2]; train = [n for n in names if not zlib.crc32(n.encode()) % 2]
def arr(ns, m): return np.array([pairs[n][m] for n in ns], dtype=np.int16)
T1, T2 = arr(train, 1), arr(train, 2); H1, H2 = arr(test, 1), arr(test, 2)
vals = np.array([2, 11, 25, 37]); rng = np.random.default_rng(3)
pb = lambda q: np.minimum(q // 10, 3)
def fit_chain(X):
    n = len(X); p0 = np.array([(X[:, 0] == v).mean() for v in vals]); T = np.zeros((151, 4, 4))
    for c in range(1, 151):
        marg = np.array([(X[:, c] == v).mean() for v in vals])
        for b in range(4):
            s = pb(X[:, c - 1]) == b; cnt = np.array([(X[s, c] == v).sum() for v in vals])
            T[c, b] = cnt / cnt.sum() if cnt.sum() >= 30 else marg
    return p0, T
def run_chain(model, n):
    p0, T = model; out = np.zeros((n, 151), np.int16); out[:, 0] = vals[rng.choice(4, n, p=p0)]
    for c in range(1, 151):
        cdf = np.cumsum(T[c, pb(out[:, c - 1])], axis=1); u = rng.random(n)
        out[:, c] = vals[(u[:, None] > cdf).sum(axis=1).clip(0, 3)]
    return out
n = len(test)
# M1
m1 = (run_chain(fit_chain(T1), n), run_chain(fit_chain(T2), n))
# M2
idx = rng.integers(0, len(train), n); m2 = (T1[idx], T2[idx])
# M3
pm = (T1.mean(axis=1) + T2.mean(axis=1)) / 2; cuts = np.quantile(pm, [0.25, 0.5, 0.75]); cls = np.searchsorted(cuts, pm)
draw = rng.integers(0, len(train), n); dcls = cls[draw]
o1 = np.zeros((n, 151), np.int16); o2 = np.zeros((n, 151), np.int16)
for k in range(4):
    s = dcls == k
    o1[s] = run_chain(fit_chain(T1[cls == k]), s.sum()); o2[s] = run_chain(fit_chain(T2[cls == k]), s.sum())
m3 = (o1, o2)
def score(X1, X2):
    X = np.vstack([X1, X2]); rm = X.mean(axis=1)
    crash = ((X[:, -20:] < 15).sum(axis=1) >= 10).mean() * 100
    both = (((X1[:, -20:] < 15).sum(1) >= 10) & ((X2[:, -20:] < 15).sum(1) >= 10)).mean() * 100
    prev = pb(X[:, :-1]) == 1; p11 = (X[:, 1:][prev] == 11).mean() * 100
    return [rm.std(ddof=1), 100 * (X.min(1) == 37).mean(), 100 * (rm < 33).mean(), crash, both,
            np.corrcoef(X1.mean(1), X2.mean(1))[0, 1], p11, X.std(axis=0, ddof=1).mean()]
hdr = ["SD read mean", "all Q37 %", "mean<Q33 %", "crashed %", "both mates crashed %", "R1-R2 corr", "P(Q11|prev low) %", "per-cycle SD"]
rows = [("held-out real", score(H1, H2)), ("M1 chain (spike now)", score(*m1)), ("M2 copy real strings", score(*m2)), ("M3 4-class chain", score(*m3))]
print(f"train pairs {len(train)}, test pairs {n}")
print(f"{'':22s}" + "".join(f"{h:>21s}" for h in hdr))
for lab, s in rows: print(f"{lab:22s}" + "".join(f"{v:21.3f}" for v in s))
