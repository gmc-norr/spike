"""E5b: E5 with far more training data. Train: every read in the chr20 slice outside the 25
windows (+-3 kb); test: the real reads in the windows (reads.pkl). Bits per quality, held out."""
import os
import pysam, pickle, numpy as np
wins = [(p - 5500, p + 5500) for p in (38550000 + i * 60000 for i in range(25))]
lut = np.full(256, 4, np.int64); lut[[65, 67, 71, 84]] = [0, 1, 2, 3]
qlut = np.full(64, 0, np.int64); qlut[[2, 11, 25, 37]] = [0, 1, 2, 3]
comp = bytes.maketrans(b"ACGTN", b"TGCAN")
Q, B, M = [], [], []; seen = set()
for r in pysam.AlignmentFile(os.environ["SLICE"]):
    if r.flag & 0xF0C or r.query_length != 151: continue
    if any(a <= r.reference_start <= b for a, b in wins): continue
    key = (r.query_name, r.is_read1)
    if key in seen: continue
    seen.add(key)
    q = np.array(r.query_qualities); s = r.query_sequence.encode()
    if r.is_reverse: q = q[::-1]; s = s.translate(comp)[::-1]
    Q.append(qlut[q]); B.append(lut[np.frombuffer(s, np.uint8)]); M.append(1 if r.is_read1 else 2)
Rt = pickle.load(open("reads.pkl", "rb"))
Qt = np.array([qlut[d["q"]] for d in Rt]); Bt = np.array([lut[d["base"]] for d in Rt]); Mt = np.array([d["mate"] for d in Rt])
Q = np.array(Q); B = np.array(B); M = np.array(M)
print("train reads", len(Q), "test reads", len(Qt))
def ctx(Q, B, M, keys):
    n = len(Q); P = np.broadcast_to(np.arange(151) // 8, Q.shape)
    sh = lambda X, k, fill: np.concatenate([np.full((n, k), fill), X[:, :-k]], axis=1)
    shl = lambda X, k, fill: np.concatenate([X[:, k:], np.full((n, k), fill)], axis=1)
    q1 = sh(Q, 1, 4); q2 = sh(Q, 2, 4)
    ch = np.zeros_like(Q); ch[:, 2:] = (Q[:, 1:-1] != Q[:, :-2]); d = np.cumsum(ch, axis=1)
    db = np.where(d < 4, d, np.where(d < 9, 4, np.where(d < 17, 5, 6)))
    F = dict(pos=(P, 19), q1=(q1, 5), q2=(q2, 5), delta=(db, 7), mate=(np.broadcast_to(M[:, None], Q.shape), 3),
             b0=(B, 6), bm1=(sh(B, 1, 5), 6), bp1=(shl(B, 1, 5), 6), bm2=(sh(B, 2, 5), 6), bp2=(shl(B, 2, 5), 6),
             bm3=(sh(B, 3, 5), 6))
    c = np.zeros(Q.shape, np.int64); size = 1
    for k in keys:
        v, s = F[k]; c = c * s + v; size *= s
    return c, size
models = [
    ("position + 2 previous Q + delta (fqzcomp-like)", ["pos", "q1", "q2", "delta"]),
    ("  ... + base", ["pos", "q1", "q2", "delta", "b0"]),
    ("  ... + 3 bases (b-1, b, b+1)", ["pos", "q1", "q2", "delta", "bm1", "b0", "bp1"]),
    ("  ... + 3 bases before (b-3..b-1)", ["pos", "q1", "q2", "delta", "bm3", "bm2", "bm1"]),
    ("  ... + 5 bases", ["pos", "q1", "q2", "delta", "bm2", "bm1", "b0", "bp1", "bp2"]),
]
for name, keys in models:
    c, size = ctx(Q, B, M, keys)
    cnt = np.bincount((c * 4 + Q).ravel(), minlength=size * 4).reshape(size, 4).astype(float)
    ct, _ = ctx(Qt, Bt, Mt, keys)
    p = (cnt[ct, Qt] + 0.5) / (cnt[ct].sum(axis=-1) + 2.0)
    print(f"{name:48s} {-np.log2(p).mean():.4f} bits/quality   contexts used {(cnt.sum(1) > 0).sum()}")
