import os
import pysam, math, gzip, pickle, zlib, numpy as np
def wilson(k,n,z=1.96):
    c=(k+z*z/2)/(n+z*z); h=z*math.sqrt(k*(n-k)/n+z*z/4)/(n+z*z); return 100*(c-h),100*(c+h)
for s in ("e4m", "e9m"):
    k=n=0
    for r in pysam.AlignmentFile(f"{s}.bam"):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        n+=1; k+=any(op==4 for op,_ in r.cigartuples)
    lo,hi=wilson(k,n); print(f"{s}: soft-clipped {100*k/n:.2f}% ({lo:.2f}-{hi:.2f}) n={k}/{n}")
# clumping vs templates
R = pickle.load(open("reads.pkl", "rb")); P = {}
for d in R: P.setdefault(d["name"], {})[d["mate"]] = d
names = sorted(n for n, p in P.items() if len(p) == 2); test = [n for n in names if zlib.crc32(n.encode()) % 2]
fa = pysam.FastaFile(os.environ["REF"])
def template(d):
    p = d["pos"]; dif = np.diff(p)
    if (p < 0).any() or not (np.all(dif == 1) or np.all(dif == -1)): return None
    t = fa.fetch("chr20", int(p.min()), int(p.max()) + 1).upper()
    return t if dif[0] == 1 else t.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
tests = [(n, template(P[n][1]), template(P[n][2])) for n in test]; tests = [t for t in tests if t[1] and t[2] and "N" not in t[1] + t[2]]
for pre, tag in (("e4m", "copyQEm"), ("e9m", "fqzM")):
    X = []
    for m in (1, 2):
        rows = {}
        with gzip.open(f"{pre}_R{m}.fq.gz", "rt") as fh:
            while True:
                h = fh.readline()
                if not h: break
                s = fh.readline().strip(); fh.readline(); fh.readline(); rows[int(h[1:].strip().split("_")[1])] = s
        for i, (n, t1, t2) in enumerate(tests):
            t = t1 if m == 1 else t2; X.append(np.array([a != b and a != "N" for a, b in zip(rows[i], t)]))
    X = np.array(X); t = X[:, -20:].sum(1)
    print(f"{pre}: >=5 errors in last 20 {100*(t>=5).mean():.2f}%  >=3 {100*(t>=3).mean():.2f}%  errors/read {X.sum(1).mean():.3f}")
