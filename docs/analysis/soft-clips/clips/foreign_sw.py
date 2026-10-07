"""The 'rest' foreign clips: do they align near the clip site once insertions and deletions are
allowed? Smith-Waterman (match +1, mismatch -1, gap -2) of each clip against the reference
from 300 bp before to 300 bp after its boundary. 'Local' = aligns over >= 80% of the clip at
>= 80% identity. Also tried: the reverse complement (an inverted piece)."""
import os, pickle, collections, numpy as np, pysam
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; F = pickle.load(open("foreign_cls.pkl", "rb"))
fa = pysam.FastaFile(os.environ["REF"]); comp = str.maketrans("ACGTN", "TGCAN")
def sw(q, t):
    best = (0, 0, 0)  # score, matches, aligned length
    H = [[0] * (len(t) + 1) for _ in range(len(q) + 1)]; M = [[0] * (len(t) + 1) for _ in range(len(q) + 1)]; Lg = [[0] * (len(t) + 1) for _ in range(len(q) + 1)]
    for i in range(1, len(q) + 1):
        for j in range(1, len(t) + 1):
            d = H[i-1][j-1] + (1 if q[i-1] == t[j-1] else -1); u = H[i-1][j] - 2; l = H[i][j-1] - 2
            h = max(0, d, u, l); H[i][j] = h
            if h == 0: M[i][j] = Lg[i][j] = 0
            elif h == d: M[i][j] = M[i-1][j-1] + (q[i-1] == t[j-1]); Lg[i][j] = Lg[i-1][j-1] + 1
            elif h == u: M[i][j] = M[i-1][j]; Lg[i][j] = Lg[i-1][j] + 1
            else: M[i][j] = M[i][j-1]; Lg[i][j] = Lg[i][j-1] + 1
            if h > best[0]: best = (h, M[i][j], Lg[i][j])
    return best
res = collections.Counter(); ex = []; per = {}
for i, k in F:
    if k != "rest": continue
    c = C[i]; q = c["clip_ref_orient"][:120]
    t = fa.fetch("chr20", max(0, c["boundary"] - 300), c["boundary"] + 300).upper()
    for lab, qq in (("local", q), ("local-inverted", q.translate(comp)[::-1])):
        s, m, L = sw(qq, t)
        if L >= 0.8 * len(qq) and m >= 0.8 * L: res[lab] += 1; per[i] = lab; break
    else:
        res["not local"] += 1; ex.append(c); per[i] = "not local"
n = sum(res.values())
print(f"'rest' foreign clips {n}: " + ", ".join(f"{k} {v} ({100*v/n:.0f}%)" for k, v in res.most_common()))
L = [c["L"] for c in ex]; print(f"  not local: median length {np.median(L):.0f}, mean Q {np.mean([c['meanq'] for c in ex]):.1f}, 3' {100*np.mean([c['three'] for c in ex]):.0f}%")
pickle.dump(per, open("foreign_sw.pkl", "wb"))
