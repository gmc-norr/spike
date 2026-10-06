"""E5: how much do the bases tell about quality, next to previous qualities and read state?
Held-out cross-entropy (bits per quality) of count models with +0.5 smoothing, trained on the
train half of the read pairs, scored on the test half (same split as e3.py)."""
import pickle, zlib, numpy as np, collections
R = pickle.load(open("reads.pkl", "rb"))
code = {ord("A"): 0, ord("C"): 1, ord("G"): 2, ord("T"): 3, ord("N"): 4}
qi = {2: 0, 11: 1, 25: 2, 37: 3}
def feats(d):
    q = [qi[int(x)] for x in d["q"]]; b = [code[int(x)] for x in d["base"]]; n = len(q)
    rows = []; delta = 0
    for i in range(n):
        q1 = q[i - 1] if i > 0 else 4; q2 = q[i - 2] if i > 1 else 4
        if i > 1 and q[i - 1] != q[i - 2]: delta += 1
        dbin = min(delta, 3) if delta < 4 else (4 if delta < 9 else (5 if delta < 17 else 6))
        bb = lambda j: b[j] if 0 <= j < n else 5
        rows.append((q[i], dict(pos=i // 8, q1=q1, q2=q2, delta=dbin, mate=d["mate"],
                                b0=bb(i), bm1=bb(i - 1), bp1=bb(i + 1), bm2=bb(i - 2), bp2=bb(i + 2))))
    return rows
tr, te = [], []
for d in R:
    (te if zlib.crc32(d["name"].encode()) % 2 else tr).extend(feats(d))
print("train symbols", len(tr), "test symbols", len(te))
models = [
    ("position", ["pos"]),
    ("+ previous Q", ["pos", "q1"]),
    ("+ previous Q + base (spike's level 1)", ["pos", "q1", "b0"]),
    ("+ 2 previous Q", ["pos", "q1", "q2"]),
    ("+ 2 previous Q + read state (delta)", ["pos", "q1", "q2", "delta"]),
    ("  ... + mate", ["pos", "q1", "q2", "delta", "mate"]),
    ("  ... + base", ["pos", "q1", "q2", "delta", "b0"]),
    ("  ... + 3 bases (b-1, b, b+1)", ["pos", "q1", "q2", "delta", "bm1", "b0", "bp1"]),
    ("  ... + 5 bases", ["pos", "q1", "q2", "delta", "bm2", "bm1", "b0", "bp1", "bp2"]),
    ("  ... + 2 bases before (b-2, b-1)", ["pos", "q1", "q2", "delta", "bm2", "bm1"]),
]
for name, keys in models:
    C = collections.defaultdict(lambda: np.zeros(4))
    for s, f in tr: C[tuple(f[k] for k in keys)][s] += 1
    bits = 0.0
    for s, f in te:
        c = C.get(tuple(f[k] for k in keys)); c = c if c is not None else np.zeros(4)
        bits -= np.log2((c[s] + 0.5) / (c.sum() + 2.0))
    print(f"{name:42s} {bits / len(te):.4f} bits/quality   contexts {len(C)}")
