"""E1: how much does quality depend on the base? E2: error rate by Q, by read class, clipped bases included."""
import pickle, numpy as np
R = pickle.load(open("reads.pkl", "rb"))
Q = np.vstack([d["q"] for d in R]); B = np.vstack([d["base"] for d in R])
MM = np.vstack([d["mm"] for d in R]); OK = np.vstack([d["ok"] for d in R]); CL = np.vstack([d["clip"] for d in R])
print("E1  quality by called base (sequencing order), all real reads")
for b in b"ACGT":
    s = B == b
    print(f"  {chr(b)}: n={s.sum():8d}  mean Q {Q[s].mean():.2f}  share below Q37 {100*(Q[s]<37).mean():.2f}%")
# previous base effect: share <37 by (prev base, base)
low = (Q < 37)
print("  share below Q37 by previous base -> this base (top 4 and bottom 2 of 16)")
rows = []
for pb in b"ACGT":
    for b in b"ACGT":
        s = np.zeros_like(low); s[:, 1:] = (B[:, :-1] == pb) & (B[:, 1:] == b)
        rows.append((100 * low[s].mean(), chr(pb) + chr(b), s.sum()))
rows.sort(reverse=True)
for r in rows[:4] + rows[-2:]: print(f"    {r[1]}: {r[0]:.2f}%  (n={r[2]})")
print("  overall share below Q37: %.2f%%" % (100 * low.mean()))

rm = Q.mean(axis=1)
poor = rm < 33
crash = (Q[:, -20:] < 15).sum(axis=1) >= 10
print(f"\nE2  error rate by Q (mismatch vs reference, variant sites masked); reads: good {(~poor).sum()}, poor (mean<Q33) {poor.sum()}, crashed {crash.sum()}")
print("  Q   nominal  | all reads: aligned only, +clipped | good reads +clipped | poor reads +clipped | crashed reads +clipped")
for qv in (2, 11, 25, 37):
    nom = 10 ** (-qv / 10) * 100
    def rate(sel, incl_clip=True):
        s = (Q == qv) & OK & sel[:, None]
        if not incl_clip: s &= ~CL
        return 100 * MM[s].mean() if s.sum() else float("nan"), int(s.sum())
    a0 = rate(np.ones(len(R), bool), False); a1 = rate(np.ones(len(R), bool))
    g = rate(~poor); p = rate(poor); c = rate(crash)
    print(f"  Q{qv:<2} {nom:7.3f}% | {a0[0]:7.3f}% {a1[0]:7.3f}% (n={a1[1]}) | {g[0]:7.3f}% (n={g[1]}) | {p[0]:7.3f}% (n={p[1]}) | {c[0]:7.3f}% (n={c[1]})")
clipped_reads = CL.any(axis=1)
print(f"\n  soft-clipped reads: all {100*clipped_reads.mean():.2f}%  good {100*clipped_reads[~poor].mean():.2f}%  poor {100*clipped_reads[poor].mean():.2f}%  crashed {100*clipped_reads[crash].mean():.2f}%")
