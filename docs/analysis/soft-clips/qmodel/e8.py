"""E8: can a model free of a compressor's limits beat fqzcomp's context?
Held-out bits per quality on the window test reads; trained on the slice outside the windows.
Every model pays for its per-read selector (R1: -log2 p(sel), R2: -log2 p(sel2 | sel1)).
Big contexts back off to smaller ones (hierarchical Dirichlet: p = (c + a*p_parent)/(n + a))."""
import os
import pysam, pickle, zlib, numpy as np
wins = [(p - 5500, p + 5500) for p in (38550000 + i * 60000 for i in range(25))]
qlut = np.zeros(64, np.int64); qlut[[2, 11, 25, 37]] = [0, 1, 2, 3]; vals = np.array([2, 11, 25, 37])
lut = np.full(256, 4, np.int64); lut[[65, 67, 71, 84]] = [0, 1, 2, 3]
comp = bytes.maketrans(b"ACGTN", b"TGCAN")
reads = {}
for r in pysam.AlignmentFile(os.environ["SLICE"]):
    if r.flag & 0xF0C or r.query_length != 151 or any(a <= r.reference_start <= b for a, b in wins): continue
    q = np.array(r.query_qualities); s = r.query_sequence.encode()
    if r.is_reverse: q = q[::-1]; s = s.translate(comp)[::-1]
    reads.setdefault(r.query_name, {})[1 if r.is_read1 else 2] = (qlut[q], lut[np.frombuffer(s, np.uint8)])
pr = [p for p in reads.values() if len(p) == 2]
TR = [(np.array([p[m][0] for p in pr]), np.array([p[m][1] for p in pr])) for m in (1, 2)]
R = pickle.load(open("reads.pkl", "rb")); P = {}
for d in R: P.setdefault(d["name"], {})[d["mate"]] = d
names = sorted(n for n, p in P.items() if len(p) == 2)
test = [n for n in names if zlib.crc32(n.encode()) % 2]
TE = [(np.array([qlut[P[n][m]["q"]] for n in test]), np.array([lut[P[n][m]["base"]] for n in test])) for m in (1, 2)]
NS = 8
cuts = np.quantile(np.concatenate([vals[TR[0][0]].mean(1), vals[TR[1][0]].mean(1)]), [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8])
sel = lambda Q: np.searchsorted(cuts, vals[Q].mean(1), side="right")
def feats(Q, B, S, mate):
    n = len(Q); F = {}
    sh = lambda X, k, f: np.concatenate([np.full((n, k), f), X[:, :-k]], 1)
    shl = lambda X, k, f: np.concatenate([X[:, k:], np.full((n, k), f)], 1)
    h = np.zeros(Q.shape, np.int64)
    for k in range(1, 9): h = h * 5 + sh(Q, k, 4)
    F["h8"] = h; F["h5"] = h // 5 ** 3; F["h2"] = h // 5 ** 6; F["h1"] = h // 5 ** 7
    i = np.arange(151); F["pos3"] = np.broadcast_to(np.minimum(7, (151 - i) >> 4), Q.shape); F["pos19"] = np.broadcast_to(i // 8, Q.shape)
    ch = np.zeros(Q.shape, np.int64); ch[:, 2:] = Q[:, 1:-1] != Q[:, :-2]; d = np.cumsum(ch, 1)
    F["ch2"] = (d >= 2).astype(np.int64); F["d7"] = np.where(d < 4, d, np.where(d < 9, 4, np.where(d < 17, 5, 6)))
    low = (sh(Q, 1, 3) <= 1).astype(np.int64)
    F["lowrun"] = np.minimum(8, np.cumsum(low, 1) - np.maximum.accumulate(np.where(low == 0, np.cumsum(low, 1), 0), 1))
    F["sel"] = np.broadcast_to(S[:, None], Q.shape); F["mate"] = np.full(Q.shape, mate)
    F["b0"] = B; F["bp1"] = shl(B, 1, 5); F["bm1"] = sh(B, 1, 5)
    return F
SIZES = dict(h8=5**8, h5=5**5, h2=25, h1=5, pos3=8, pos19=19, ch2=2, d7=7, lowrun=9, sel=NS, mate=3, b0=6, bp1=6, bm1=6)
def cid(F, keys):
    c = np.zeros(F["h1"].shape, np.int64)
    for k in keys: c = c * SIZES[k] + F[k]
    return c
FTR = [feats(TR[m][0], TR[m][1], sel(TR[m][0]), m + 1) for m in (0, 1)]
FTE = [feats(TE[m][0], TE[m][1], sel(TE[m][0]), m + 1) for m in (0, 1)]
YTR = np.concatenate([TR[0][0].ravel(), TR[1][0].ravel()]); YTE = np.concatenate([TE[0][0].ravel(), TE[1][0].ravel()])
def level(keys):
    ctr = np.concatenate([cid(FTR[0], keys).ravel(), cid(FTR[1], keys).ravel()])
    cte = np.concatenate([cid(FTE[0], keys).ravel(), cid(FTE[1], keys).ravel()])
    u, inv = np.unique(ctr, return_inverse=True)
    cnt = np.bincount(inv * 4 + YTR, minlength=len(u) * 4).reshape(len(u), 4).astype(float)
    j = np.searchsorted(u, cte); hit = (j < len(u)) & (u[np.minimum(j, len(u) - 1)] == cte)
    out = np.zeros((len(cte), 4)); out[hit] = cnt[j[hit]]
    return out
def chain(levels, a=4.0):
    p = None
    for keys in levels:
        c = level(keys)
        p = (c + 0.5) / (c.sum(1, keepdims=True) + 2.0) if p is None else (c + a * p) / (c.sum(1, keepdims=True) + a)
    return -np.log2(p[np.arange(len(YTE)), YTE]).mean()
# selector cost, bits per quality
S1, S2 = sel(TR[0][0]), sel(TR[1][0]); s1, s2 = sel(TE[0][0]), sel(TE[1][0])
p1 = np.bincount(S1, minlength=NS) + 0.5; p1 /= p1.sum()
p21 = np.bincount(S1 * NS + S2, minlength=NS * NS).reshape(NS, NS) + 0.5; p21 /= p21.sum(1, keepdims=True)
selcost = (-np.log2(p1[s1]).sum() - np.log2(p21[s1, s2]).sum()) / len(YTE)
print(f"selector cost {selcost:.4f} bits/quality; train pairs {len(pr)}, test pairs {len(test)}")
base0 = ["mate", "pos19", "h1"]
fqz = ["mate", "h5", "pos3", "ch2", "sel"]
M = [
    ("fqzcomp context (no backoff)", [fqz], True),
    ("fqzcomp context, backoff", [base0, fqz], True),
    ("fqzcomp + base (fqzcomp5-like), backoff", [base0, fqz, fqz + ["b0"]], True),
    ("bigger: 8 prev Q, 19 pos, delta, low-run, sel", [base0, fqz, ["mate", "h8", "pos19", "d7", "lowrun", "sel"]], True),
    ("bigger + bases (b-1, b, b+1)", [base0, fqz, ["mate", "h8", "pos19", "d7", "lowrun", "sel"], ["mate", "h8", "pos19", "d7", "lowrun", "sel", "bm1", "b0", "bp1"]], True),
]
for name, levels, uses_sel in M:
    b = chain(levels) + (selcost if uses_sel else 0)
    print(f"{name:48s} {b:.4f} bits/quality (incl. selector)")
