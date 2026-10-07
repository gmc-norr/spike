"""Re-score the earlier test with causes separated. Real test reads (the 5,190 held-out pairs of e4)
are labelled by their own clips' causes (original alignment). 'Error-type' clips = badend + hiQerr + polyG."""
import os
import pickle, zlib, collections, numpy as np, gzip, pysam
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; cls = pickle.load(open("cls2.pkl", "rb"))
cause = collections.defaultdict(set)
for c, k in zip(C, cls): cause[(c["name"], c["r1"])].add(k)
R = pickle.load(open("../qmodel/reads.pkl", "rb")); P = {}
for d in R: P.setdefault(d["name"], {})[d["mate"]] = d
names = sorted(n for n, p in P.items() if len(p) == 2)
test = [n for n in names if zlib.crc32(n.encode()) % 2]
fa = pysam.FastaFile(os.environ["REF"])
def ok_template(d):
    p = d["pos"]; dif = np.diff(p)
    return not ((p < 0).any() or not (np.all(dif == 1) or np.all(dif == -1)))
tests = [n for n in test if ok_template(P[n][1]) and ok_template(P[n][2])]
keys = [(n, m == 1) for n in tests for m in (1, 2)]
N = len(keys)
cc = collections.Counter()
for k in keys:
    for x in cause.get(k, ()): cc[x] += 1
anyclip = sum(1 for k in keys if k in cause)
errtype = {"badend", "hiQerr", "polyG"}
nonerr = {"adapter", "site", "foreign", "chimera"}
only_err = sum(1 for k in keys if cause.get(k) and cause[k] <= errtype)
print(f"held-out real reads {N} (original alignment): soft-clipped {100*anyclip/N:.2f}%")
print("  by cause: " + ", ".join(f"{k} {100*v/N:.2f}%" for k, v in cc.most_common()))
print(f"  clipped ONLY for error-type causes (bad end, high-Q mismatches, polyG): {100*only_err/N:.2f}%")
# clumping in real reads, excluding reads with any non-error clip
mm = {(d["name"], d["mate"] == 1): d["mm"] & d["ok"] for d in R}
keep = [k for k in keys if not (cause.get(k, set()) & nonerr)]
X = np.array([mm[k] for k in keep]); t = X[:, -20:].sum(1)
print(f"real reads without adapter/site/foreign/chimera clips: {len(keep)}; >=5 errors in last 20: {100*(t>=5).mean():.2f}%  >=3: {100*(t>=3).mean():.2f}%  errors/read {X.sum(1).mean():.3f}")
