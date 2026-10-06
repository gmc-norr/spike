"""E7: fqzcomp's quality context (htscodecs fqzcomp_qual, strat 0, as derived for NovaSeq 4-bin:
last 5 qualities, position from the 3' end min(7, remaining>>4), 'had >= 2 changes' bit, per-read
selector = read-mean class) used as a GENERATOR. Trained on the slice outside the windows; mates'
selectors drawn jointly. Errors: Bernoulli per (Q, selector) from the training reads' own
mismatches, clipped bases included. Scored on the held-out window reads; then aligned (e4 templates)."""
import os
import pysam, pickle, zlib, gzip, numpy as np
rng = np.random.default_rng(5)
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
Q1 = np.array([p[1][0] for p in pr]); Q2 = np.array([p[2][0] for p in pr])
B1 = np.array([p[1][1] for p in pr]); B2 = np.array([p[2][1] for p in pr])
print("training pairs", len(pr))
NS = 8
allmean = np.concatenate([vals[Q1].mean(1), vals[Q2].mean(1)])
cuts = np.quantile(allmean, [0.01, 0.03, 0.08, 0.2, 0.4, 0.6, 0.8])
sel = lambda Q: np.searchsorted(cuts, vals[Q].mean(1), side="right")
def ctxs(Q, S, B=None):
    n = len(Q); hist = np.zeros(Q.shape, np.int64); h = np.zeros(n, np.int64); ch = np.zeros(n, np.int64)
    C = np.zeros(Q.shape, np.int64)
    for i in range(151):
        pos = min(7, (151 - i) >> 4)
        c = ((h * 8 + pos) * 2 + (ch >= 2)) * NS + S
        if B is not None: c = c * 5 + B[:, i]
        C[:, i] = c
        if i > 0: ch += Q[:, i] != Q[:, i - 1]
        h = ((h << 2) + Q[:, i]) & 1023
    return C
def table(Qs, Ss, Bs, base):
    size = 1024 * 8 * 2 * NS * (5 if base else 1)
    T = np.zeros((2, size, 4))
    for m in (0, 1):
        C = ctxs(Qs[m], Ss[m], Bs[m] if base else None)
        T[m] = np.bincount((C * 4 + Qs[m]).ravel(), minlength=size * 4).reshape(size, 4)
    return T
S1, S2 = sel(Q1), sel(Q2)
joint = np.bincount(S1 * NS + S2, minlength=NS * NS).astype(float); joint /= joint.sum()
T = table((Q1, Q2), (S1, S2), (B1, B2), False)
Tb = table((Q1, Q2), (S1, S2), (B1, B2), True)
# held-out bits/quality on the window test reads
R = pickle.load(open("reads.pkl", "rb")); P = {}
for d in R: P.setdefault(d["name"], {})[d["mate"]] = d
names = sorted(n for n, p in P.items() if len(p) == 2)
test = [n for n in names if zlib.crc32(n.encode()) % 2]; train = [n for n in names if not zlib.crc32(n.encode()) % 2]
H = [np.array([qlut[P[n][m]["q"]] for n in test]) for m in (1, 2)]
HB = [np.array([lut[P[n][m]["base"]] for n in test]) for m in (1, 2)]
for lab, TT, base in (("fqzcomp context", T, False), ("fqzcomp context + base (fqzcomp5-like)", Tb, True)):
    bits = []
    for m in (0, 1):
        C = ctxs(H[m], sel(H[m]), HB[m] if base else None); cnt = TT[m][C]
        bits.append(-np.log2((np.take_along_axis(cnt, H[m][..., None], -1)[..., 0] + 0.5) / (cnt.sum(-1) + 2)))
    print(f"held-out bits/quality, {lab}: {np.concatenate(bits).mean():.4f}")
def generate(n):
    js = rng.choice(NS * NS, n, p=joint); S = (js // NS, js % NS); out = []
    for m in (0, 1):
        Q = np.zeros((n, 151), np.int64); h = np.zeros(n, np.int64); ch = np.zeros(n, np.int64)
        for i in range(151):
            pos = min(7, (151 - i) >> 4); c = ((h * 8 + pos) * 2 + (ch >= 2)) * NS + S[m]
            cnt = T[m][c] + 0.01; cdf = np.cumsum(cnt / cnt.sum(1, keepdims=True), 1)
            Q[:, i] = (rng.random(n)[:, None] > cdf).sum(1).clip(0, 3)
            if i > 0: ch += Q[:, i] != Q[:, i - 1]
            h = ((h << 2) + Q[:, i]) & 1023
        out.append(vals[Q])
    return out, S
def score(X1, X2):
    X = np.vstack([X1, X2]); rm = X.mean(axis=1)
    return [rm.std(ddof=1), 100 * (X.min(1) == 37).mean(), 100 * (rm < 33).mean(),
            100 * ((X[:, -20:] < 15).sum(1) >= 10).mean(), np.corrcoef(X1.mean(1), X2.mean(1))[0, 1], X.std(axis=0, ddof=1).mean()]
hdr = ["SD read mean", "all Q37 %", "mean<Q33 %", "crashed %", "R1-R2 corr", "per-cycle SD"]
print(f"{'':22s}" + "".join(f"{h:>14s}" for h in hdr))
print(f"{'held-out real':22s}" + "".join(f"{v:14.3f}" for v in score(vals[H[0]], vals[H[1]])))
(G1, G2), GS = generate(len(test))
print(f"{'fqzcomp generator':22s}" + "".join(f"{v:14.3f}" for v in score(G1, G2)))
# errors per (Q, selector) from the window train reads, clipped bases included
E = np.zeros((4, NS, 2))
for n in train:
    for m in (1, 2):
        d = P[n][m]; s = sel(qlut[d["q"]][None])[0]; k = d["ok"]
        np.add.at(E, (qlut[d["q"][k]], s, d["mm"][k].astype(int)), 1)
rate = E[..., 1] / np.maximum(1, E.sum(-1))
print("error % by Q (rows Q2,Q11,Q25,Q37) x selector class (worst..best):")
print(np.round(100 * rate[1:], 3))
# build reads on e4's templates (reference at the held-out pairs' own positions), align, count clips
fa = pysam.FastaFile(os.environ["REF"])
def template(d):
    p = d["pos"]; dif = np.diff(p)
    if (p < 0).any() or not (np.all(dif == 1) or np.all(dif == -1)): return None
    t = fa.fetch("chr20", int(p.min()), int(p.max()) + 1).upper()
    return t if dif[0] == 1 else t.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
tests = [(n, template(P[n][1]), template(P[n][2])) for n in test]
tests = [t for t in tests if t[1] and t[2] and "N" not in t[1] + t[2]]
(G1, G2), (S1g, S2g) = generate(len(tests))
ACGT = np.frombuffer(b"ACGT", np.uint8)
def apply(t, q, s):
    b = np.frombuffer(t.encode(), np.uint8).copy(); e = rng.random(151) < rate[qlut[q], s]
    for i in np.flatnonzero(e): b[i] = rng.choice([x for x in ACGT if x != b[i]])
    b[q == 2] = ord("N"); return b.tobytes().decode()
with gzip.open("e7_R1.fq.gz", "wt") as f1, gzip.open("e7_R2.fq.gz", "wt") as f2:
    for i, (n, t1, t2) in enumerate(tests):
        f1.write(f"@fqz_{i}\n{apply(t1, G1[i], S1g[i])}\n+\n{''.join(chr(int(x)+33) for x in G1[i])}\n")
        f2.write(f"@fqz_{i}\n{apply(t2, G2[i], S2g[i])}\n+\n{''.join(chr(int(x)+33) for x in G2[i])}\n")
print("pairs written", len(tests))
