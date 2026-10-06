"""E4: do soft clips come back? Reads are rebuilt from the reference at held-out real pairs'
own positions (so only the quality/error model differs), then aligned with the run's bwa-mem2.
 real   the held-out real reads themselves, realigned the same way
 chain  spike-like first-order chain qualities + Phred errors (10^(-Q/10))
 copyQ  a training pair's two quality strings + Phred errors        (parked quality-tails)
 copyQE a training pair's two quality strings + that pair's own error positions
"""
import os
import pickle, zlib, gzip, numpy as np, pysam
exec(open("e3.py").read().split("n = len(test)")[0])          # pairs, train/test, fit_chain, run_chain
R = pickle.load(open("reads.pkl", "rb"))
NE = pickle.load(open("../clips/nonerr_mask.pkl", "rb"))
for d in R:
    m = NE.get((d["name"], d["mate"] == 1))
    if m is not None: d["ok"] = d["ok"] & ~m
D = {(d["name"], d["mate"]): d for d in R}
fa = pysam.FastaFile(os.environ["REF"])
comp = str.maketrans("ACGTN", "TGCAN")
def template(d):
    p = d["pos"]; dif = np.diff(p)
    if (p < 0).any() or not (np.all(dif == 1) or np.all(dif == -1)): return None
    lo, hi = p.min(), p.max(); s = fa.fetch("chr20", int(lo), int(hi) + 1).upper()
    return s if dif[0] == 1 else s.translate(comp)[::-1]
tests = []
for n in test:
    t1, t2 = template(D[(n, 1)]), template(D[(n, 2)])
    if t1 and t2 and "N" not in t1 + t2: tests.append((n, t1, t2))
print("test pairs with indel-free templates:", len(tests), "of", len(test))
ACGT = np.frombuffer(b"ACGT", np.uint8)
def errs(seq, q, mask):
    s = np.frombuffer(seq.encode(), np.uint8).copy()
    for i in np.flatnonzero(mask):
        s[i] = rng.choice([b for b in ACGT if b != s[i]])
    s[q == 2] = ord("N")
    return s.tobytes().decode()
phred = lambda q: rng.random(len(q)) < 10 ** (-q / 10)
nt = len(tests); ch1 = run_chain(fit_chain(T1), nt); ch2 = run_chain(fit_chain(T2), nt)
idx = rng.integers(0, len(train), nt)
out = {m: [] for m in ("copyQEm",)}
for i, (n, t1, t2) in enumerate(tests):
    d1, d2 = D[(n, 1)], D[(n, 2)]
    s1, s2 = D[(train[idx[i]], 1)], D[(train[idx[i]], 2)]
    out["copyQEm"].append((errs(t1, s1["q"], s1["mm"] & s1["ok"]), s1["q"], errs(t2, s2["q"], s2["mm"] & s2["ok"]), s2["q"]))
with gzip.open("e4m_R1.fq.gz", "wt") as f1, gzip.open("e4m_R2.fq.gz", "wt") as f2:
    for m, rows in out.items():
        for i, (a, qa, b, qb) in enumerate(rows):
            f1.write(f"@{m}_{i}\n{a}\n+\n{''.join(chr(int(x)+33) for x in qa)}\n")
            f2.write(f"@{m}_{i}\n{b}\n+\n{''.join(chr(int(x)+33) for x in qb)}\n")
print("written", {m: len(v) for m, v in out.items()})
