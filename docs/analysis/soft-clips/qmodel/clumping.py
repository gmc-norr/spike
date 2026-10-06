"""Clumping check: share of reads with >= 3 errors in their last 20 cycles.
real: held-out test reads' own mismatches (clipped bases included, variant sites masked);
generated: mismatches against the template each was built from."""
import os
import pickle, zlib, gzip, numpy as np, pysam
R = pickle.load(open("reads.pkl", "rb")); P = {}
for d in R: P.setdefault(d["name"], {})[d["mate"]] = d
names = sorted(n for n, p in P.items() if len(p) == 2)
test = [n for n in names if zlib.crc32(n.encode()) % 2]
fa = pysam.FastaFile(os.environ["REF"])
def template(d):
    p = d["pos"]; dif = np.diff(p)
    if (p < 0).any() or not (np.all(dif == 1) or np.all(dif == -1)): return None
    t = fa.fetch("chr20", int(p.min()), int(p.max()) + 1).upper()
    return t if dif[0] == 1 else t.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
tests = [(n, template(P[n][1]), template(P[n][2])) for n in test]
tests = [t for t in tests if t[1] and t[2] and "N" not in t[1] + t[2]]
def tail_counts(mm): return mm[:, -20:].sum(1)
real = np.array([P[n][m]["mm"] & P[n][m]["ok"] for n, _, _ in tests for m in (1, 2)])
def gen(path_prefix, tag):
    out = []
    for m, f in ((1, f"{path_prefix}_R1.fq.gz" if path_prefix != "e4" else "R1.fq.gz"), (2, f"{path_prefix}_R2.fq.gz" if path_prefix != "e4" else "R2.fq.gz")):
        rows = {}
        with gzip.open(f, "rt") as fh:
            while True:
                h = fh.readline()
                if not h: break
                s = fh.readline().strip(); fh.readline(); fh.readline()
                name = h[1:].strip()
                if name.startswith(tag + "_"): rows[int(name.split("_")[1])] = s
        for i, (n, t1, t2) in enumerate(tests):
            t = t1 if m == 1 else t2; s = rows[i]
            out.append((i, m, np.array([a != b and a != "N" for a, b in zip(s, t)])))
    out.sort(key=lambda x: (x[0], x[1])); return np.array([x[2] for x in out])
sets = {"real": real, "copy quality + error spots": gen("e4", "copyQE"),
        "fqzcomp + per-base errors (e7)": gen("e7", "fqz"), "fqzcomp + richer per-base errors (e9)": gen("e9", "fqzE")}
for k, X in sets.items():
    t = tail_counts(X)
    print(f"{k:40s} reads {len(X)}  >=3 errors in last 20: {100*(t>=3).mean():.2f}%   >=5: {100*(t>=5).mean():.2f}%   mean errors/read {X.sum(1).mean():.3f}")
