"""What happens in the bad tails? Step 5: how much do runs of one letter (homopolymers) in the
template set off quality loss, on an unbiased sample (collect.py ... all)?
For each run of one base in the reference under the read (sequencing direction) that ends at cycle e:
the low-quality (< Q15) share and the error rate in the 20 cycles after it, against what reads have
at those same cycles in general (observed / expected). And per read: P(crash) by the longest run.
Usage: analyze4.py ALL.pkl"""
import sys, pickle, collections, numpy as np

R = pickle.load(open(sys.argv[1], "rb"))
B = "ACGT"
n = len(R)
Q = np.array([r["q"] for r in R]); LOW = Q < 15
E = np.array([(r["called"] != r["ref"]) & r["counted"] for r in R]); C = np.array([r["counted"] for r in R])
base_low = LOW.mean(0)                         # per cycle, all reads
base_err = E.sum(0) / np.maximum(1, C.sum(0))
crashed = np.array([r["crashed"] for r in R])
print(f"reads {n}; crashed {100*crashed.mean():.2f}%")

def runs(t):
    """(base, length, last cycle) of each run of one base in the template t (codes 0-3)."""
    out = []; i = 0
    while i < len(t):
        j = i
        while j + 1 < len(t) and t[j + 1] == t[i]:
            j += 1
        if t[i] < 4:
            out.append((int(t[i]), j - i + 1, j))
        i = j + 1
    return out

W = 20
LB = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 11), (12, 151)]
acc = {(b, lb): np.zeros(6) for b in range(4) for lb in LB}  # obs low, exp low, obs err, exp err, n err bases, windows
by_cycle = {(lb, cb): np.zeros(3) for lb in LB for cb in ((10, 60), (61, 100), (101, 125))}  # obs low, exp low, windows
longest = []
for k, r in enumerate(R):
    best = 0
    for b, L, e in runs(r["ref"]):
        if e > 125 or e < 10:
            continue
        best = max(best, L)
        w = slice(e + 1, e + 1 + W)
        lb = next(x for x in LB if x[0] <= L <= x[1])
        a = acc[(b, lb)]
        a[0] += LOW[k, w].sum(); a[1] += base_low[w].sum()
        cm = C[k, w]
        a[2] += E[k, w].sum(); a[3] += (base_err[w] * cm).sum(); a[4] += cm.sum(); a[5] += 1
        cb = next(x for x in ((10, 60), (61, 100), (101, 125)) if x[0] <= e <= x[1])
        bc = by_cycle[(lb, cb)]; bc[0] += LOW[k, w].sum(); bc[1] += base_low[w].sum(); bc[2] += 1
    longest.append(best)
longest = np.array(longest)

print(f"\n10a. LOW-QUALITY SHARE IN THE {W} CYCLES AFTER A RUN, observed / expected at those cycles (runs ending at cycles 10-125)")
print("      run length:      " + "".join(f"{f'{a}-{b}' if b < 151 else f'{a}+':>16s}" for a, b in LB))
for b in range(4):
    print(f"      {B[b]} (seq. dir.)     " + "".join(f"{acc[(b, lb)][0]/max(1,acc[(b, lb)][1]):>9.2f} ({int(acc[(b, lb)][5]):>5d})" for lb in LB))
print(f"\n10b. ERROR RATE IN THE {W} CYCLES AFTER A RUN, observed / expected")
for b in range(4):
    print(f"      {B[b]} (seq. dir.)     " + "".join(f"{acc[(b, lb)][2]/max(1,acc[(b, lb)][3]):>9.2f} ({int(acc[(b, lb)][4]):>6d})"[:16].rjust(16) for lb in LB))
print("\n10c. LOW SHARE AFTER A RUN (any base), observed / expected, by the cycle where the run ends")
for cb in ((10, 60), (61, 100), (101, 125)):
    print(f"      ends at cycles {cb[0]:3d}-{cb[1]:3d}  " + "".join(f"{by_cycle[(lb, cb)][0]/max(1,by_cycle[(lb, cb)][1]):>9.2f} ({int(by_cycle[(lb, cb)][2]):>5d})" for lb in LB))

print("\n10d. P(READ CRASHES) BY THE LONGEST RUN IN ITS TEMPLATE ENDING AT CYCLES 10-125")
for a, b in ((0, 4), (5, 6), (7, 8), (9, 11), (12, 151)):
    m = (longest >= a) & (longest <= b)
    print(f"      longest run {a:2d}-{b if b < 151 else '':<3}: reads {m.sum():6d} ({100*m.mean():5.1f}%)  crashed {100*crashed[m].mean():5.2f}%")
m7 = longest >= 7
pc, pc_no = crashed.mean(), crashed[~m7].mean()
print(f"      crashed reads that hold a run >= 7: {100*crashed[m7].sum()/crashed.sum():.1f}% of all crashed reads "
      f"(share of reads with one: {100*m7.mean():.1f}%); crash rate without such a run {100*pc_no:.2f}% vs overall {100*pc:.2f}% "
      f"-> {100*(1-pc_no/pc):.0f}% of crashes are set off by runs >= 7")
