"""What happens in the bad tails? Step 6: once a run of one base sets a read off, does it recover?
Low-quality share observed / expected, in windows after the run's end, for runs ending at cycles
10-90 (so 60 cycles follow). Usage: analyze5.py ALL.pkl"""
import sys, pickle, numpy as np
exec(open("tails/analyze4.py").read().split("W = 20")[0].split('"""', 2)[2])   # loads R, LOW, base_low, runs()
WINS = ((1, 5), (6, 10), (11, 20), (21, 40), (41, 60))
for lab, lo_, hi_ in (("runs 7-8", 7, 8), ("runs 9-11", 9, 11), ("runs 12+", 12, 151)):
    acc = np.zeros((len(WINS), 2)); nw = 0
    for k, r in enumerate(R):
        for b, L, e in runs(r["ref"]):
            if not (10 <= e <= 90 and lo_ <= L <= hi_):
                continue
            nw += 1
            for i, (a, z) in enumerate(WINS):
                w = slice(e + a, e + z + 1)
                acc[i, 0] += LOW[k, w].sum(); acc[i, 1] += base_low[w].sum()
    print(f"  {lab:9s} ({nw:4d} runs): " + "  ".join(f"cycles +{a}-{z}: {acc[i,0]/acc[i,1]:5.2f}" for i, (a, z) in enumerate(WINS)))
