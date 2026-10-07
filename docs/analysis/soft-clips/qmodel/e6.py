"""E6: an fqzcomp-style context model used as a generator.
Context: cycle bin (8), previous two qualities, delta (number of quality changes so far,
binned like a running read state), mate; optional base. Falls back to a smaller context
when the full one has < 20 training observations. Variant 'link': read 2 also sees read 1's
mean-quality class (4 bins) so mates can share a state.
Errors: per (Q, delta bin), from the training reads' own mismatches, clipped bases included."""
import pickle, zlib, gzip, collections, numpy as np, pysam
R = pickle.load(open("reads.pkl", "rb"))
rng = np.random.default_rng(11)
vals = np.array([2, 11, 25, 37]); qi = {2: 0, 11: 1, 25: 2, 37: 3}
code = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3, ord("N"): 4}
def dbin(d): return d if d < 4 else (4 if d < 9 else (5 if d < 17 else 6))
pairs = {}
for d in R: pairs.setdefault(d["name"], {})[d["mate"]] = d
names = sorted(n for n, p in pairs.items() if len(p) == 2)
test = [n for n in names if zlib.crc32(n.encode()) % 2]; train = [n for n in names if not zlib.crc32(n.encode()) % 2]
def r1class(q): return int(np.searchsorted([33, 35, 36.5], q.mean()))
def walk(q, base=None):
    """yield (i, ctx-without-base, q index) along one read"""
    delta = 0; out = []
    for i in range(len(q)):
        q1 = qi[int(q[i-1])] if i > 0 else 4; q2 = qi[int(q[i-2])] if i > 1 else 4
        if i > 1 and q[i-1] != q[i-2]: delta += 1
        out.append((i // 8, q1, q2, dbin(delta)))
    return out
C = collections.defaultdict(lambda: np.zeros(4)); E = collections.defaultdict(lambda: np.zeros(2))
for n in train:
    p = pairs[n]; link = r1class(p[1]["q"])
    for m in (1, 2):
        d = p[m]; ctxs = walk(d["q"])
        for i, c in enumerate(ctxs):
            s = qi[int(d["q"][i])]; b = code[int(d["base"][i])]; L = link if m == 2 else 4
            for key in (("full",) + c + (m, L, b), ("nob",) + c + (m, L), ("nolink",) + c + (m,), ("nolinkb",) + c + (m, b),
                        ("pq",) + c[:3] + (m,), ("p",) + (c[0], m)):
                C[key][s] += 1
            if d["ok"][i]: E[(int(d["q"][i]), c[3])][int(d["mm"][i])] += 1
def draw(keys):
    for k in keys:
        c = C.get(k)
        if c is not None and c.sum() >= 20: return rng.choice(4, p=c / c.sum())
    raise RuntimeError
def gen(n, link, use_base, templates=None):
    out1 = np.zeros((n, 151), np.int16); out2 = np.zeros((n, 151), np.int16)
    for j in range(n):
        L = 4
        for m, out in ((1, out1), (2, out2)):
            q = out[j]; delta = 0
            for i in range(151):
                q1 = qi[int(q[i-1])] if i > 0 else 4; q2 = qi[int(q[i-2])] if i > 1 else 4
                if i > 1 and q[i-1] != q[i-2]: delta += 1
                c = (i // 8, q1, q2, dbin(delta)); LL = L if (m == 2 and link) else 4
                b = code[ord(templates[j][m-1][i])] if use_base else None
                keys = ([("full",) + c + (m, LL, b)] if use_base and link else []) + \
                       ([("nolinkb",) + c + (m, b)] if use_base and not link else []) + \
                       ([("nob",) + c + (m, LL)] if link else []) + [("nolink",) + c + (m,), ("pq",) + c[:3] + (m,), ("p",) + (c[0], m)]
                q[i] = vals[draw(keys)]
            if m == 1: L = r1class(q)
    return out1, out2
def score(X1, X2):
    X = np.vstack([X1, X2]); rm = X.mean(axis=1)
    return [rm.std(ddof=1), 100 * (X.min(1) == 37).mean(), 100 * (rm < 33).mean(),
            100 * ((X[:, -20:] < 15).sum(1) >= 10).mean(), np.corrcoef(X1.mean(1), X2.mean(1))[0, 1], X.std(axis=0, ddof=1).mean()]
H1 = np.array([pairs[n][1]["q"] for n in test]); H2 = np.array([pairs[n][2]["q"] for n in test])
n = 3000
hdr = ["SD read mean", "all Q37 %", "mean<Q33 %", "crashed %", "R1-R2 corr", "per-cycle SD"]
print(f"{'':28s}" + "".join(f"{h:>14s}" for h in hdr))
print(f"{'held-out real':28s}" + "".join(f"{v:14.3f}" for v in score(H1, H2)))
G = {}
for lab, link in (("fqz-style", False), ("fqz-style + mate link", True)):
    G[lab] = gen(n, link, False)
    print(f"{lab:28s}" + "".join(f"{v:14.3f}" for v in score(*G[lab])))
pickle.dump(dict(E={k: v for k, v in E.items()}, G=G), open("e6.pkl", "wb"))
print("error rate % by (Q, delta bin), train reads incl. clipped bases:")
for qv in (11, 25, 37):
    print(f"  Q{qv}: " + "  ".join(f"d{b}:{100*E[(qv,b)][1]/max(1,E[(qv,b)].sum()):.3f}" for b in range(7)))
