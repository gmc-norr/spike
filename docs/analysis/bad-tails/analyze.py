"""What happens in the bad tails? Step 2: the measurements, on collect.py's output.
Groups: good (no clip, not crashed; 1 in 10), crash-noclip (crashed, no clip of any cause),
crash-clip (crashed, bad-end 3' clip), clip-nocrash (bad-end 3' clip, not crashed).
Usage: analyze.py IN.pkl"""
import sys, pickle, collections, numpy as np

R = pickle.load(open(sys.argv[1], "rb"))
B = "ACGTN"
groups = {
    "good": [r for r in R if r["good"]],
    "crash-noclip": [r for r in R if r["crashed"] and not r["anyclip"]],
    "crash-clip": [r for r in R if r["crashed"] and r["badclip"]],
    "clip-nocrash": [r for r in R if r["badclip"] and not r["crashed"]],
}
alpha = sorted({int(x) for r in R for x in np.unique(r["q"])})
print("quality values:", alpha)
for g, rs in groups.items():
    print(f"  {g:13s} {len(rs):6d} reads")

print("\n1. SHAPE OF THE TAIL")
for g in ("crash-noclip", "crash-clip", "clip-nocrash"):
    rs = groups[g]
    low_end = np.mean([r["q"][-1] < 15 for r in rs])
    # a read's low suffix: the run of < Q15 ending at the last base
    suf = []
    step = []
    for r in rs:
        lowm = r["q"] < 15
        k = 0
        while k < len(lowm) and lowm[-1 - k]:
            k += 1
        suf.append(k)
        first = np.argmax(lowm) if lowm.any() else len(lowm)
        step.append(lowm[first:].all() if lowm.any() else False)  # every base after the first low one is low
    suf = np.array(suf)
    print(f"  {g:13s} last base < Q15: {100*low_end:5.1f}%   low run at the end, bases: median {np.median(suf):.0f} "
          f"(quartiles {np.percentile(suf,25):.0f}-{np.percentile(suf,75):.0f})   once low, low to the end: {100*np.mean(step):5.1f}%")
for g in ("crash-clip", "clip-nocrash"):
    L = np.array([r["clip3"] for r in groups[g]])
    print(f"  {g:13s} 3' clip length: median {np.median(L):.0f} (quartiles {np.percentile(L,25):.0f}-{np.percentile(L,75):.0f}); "
          f"1-4 bp {100*np.mean(L<=4):.0f}%")
# within the last 40 cycles of crashed reads: share of each quality value, by cycles from the end
print("  share at each quality, crashed reads (clip or not), by cycles from the 3' end:")
cr = groups["crash-noclip"] + groups["crash-clip"]
for lo_, hi_ in ((1, 1), (2, 5), (6, 10), (11, 20), (21, 40), (41, 80)):
    vals = np.concatenate([r["q"][-hi_:len(r["q"]) - lo_ + 1] for r in cr])
    print(f"    {lo_:2d}-{hi_:2d} from end: " + "  ".join(f"Q{a}:{100*np.mean(vals==a):5.1f}%" for a in alpha))

print("\n2. ERROR RATE BY QUALITY (counted bases; mismatches / bases)")
def rate(rs, sel):
    mm = n = 0
    for r in rs:
        m = r["counted"] & sel(r)
        n += m.sum()
        mm += (m & (r["called"] != r["ref"])).sum()
    return mm, n
print(f"  {'group':13s} {'part':8s} " + " ".join(f"{'Q'+str(a):>16s}" for a in alpha))
for g, rs in groups.items():
    for part, f in (("aligned", lambda r: ~r["clipped"]), ("clipped", lambda r: r["clipped"])):
        if g in ("good", "crash-noclip") and part == "clipped":
            continue
        cells = []
        for a in alpha:
            mm, n = rate(rs, lambda r, a=a, f=f: f(r) & (r["q"] == a))
            cells.append(f"{(mm/n if n else float('nan')):.3f} ({n:>7d})")
        print(f"  {g:13s} {part:8s} " + " ".join(f"{c:>16s}" for c in cells))

print("\n3. WHICH BASE IS CALLED (sequencing orientation)")
def spectrum(rs, sel):
    called = collections.Counter(); refb = collections.Counter(); mmcalled = collections.Counter(); pair = collections.Counter()
    for r in rs:
        m = r["counted"] & sel(r)
        for c, rb in zip(r["called"][m], r["ref"][m]):
            called[c] += 1; refb[rb] += 1
            if c != rb:
                mmcalled[c] += 1; pair[(rb, c)] += 1
    return called, refb, mmcalled, pair
for label, rs, sel in (
    ("good reads, Q37 bases", groups["good"], lambda r: r["q"] == 37),
    ("good reads, Q<15 bases", groups["good"], lambda r: r["q"] < 15),
    ("crashed tails, Q11", cr, lambda r: (r["q"] == 11) & (np.arange(151) >= 111)),
    ("crashed tails, Q25", cr, lambda r: (r["q"] == 25) & (np.arange(151) >= 111)),
    ("crashed tails, Q37", cr, lambda r: (r["q"] == 37) & (np.arange(151) >= 111)),
    ("bad-end clipped bases", groups["crash-clip"] + groups["clip-nocrash"], lambda r: r["clipped"] & (np.arange(151) >= 75)),
):
    called, refb, mmc, pair = spectrum(rs, sel)
    nc, nr, nm = sum(called.values()), sum(refb.values()), sum(mmc.values())
    if not nc: continue
    print(f"  {label:24s} bases {nc:7d}  called " + " ".join(f"{B[i]}:{100*called[i]/nc:4.1f}" for i in range(4)) +
          "   reference " + " ".join(f"{B[i]}:{100*refb[i]/nr:4.1f}" for i in range(4)) +
          f"   wrong calls ({nm}) " + " ".join(f"{B[i]}:{100*mmc[i]/max(1,nm):4.1f}" for i in range(4)))
    top = pair.most_common(4)
    print(f"  {'':24s} top ref>called: " + ", ".join(f"{B[a]}>{B[b]} {100*v/max(1,nm):.0f}%" for (a, b), v in top))

print("\n4. OUT OF STEP? (wrong calls that equal the next or previous reference base in sequencing order)")
for label, rs, sel in (
    ("good reads, Q37", groups["good"], lambda r: r["q"] == 37),
    ("good reads, Q11", groups["good"], lambda r: r["q"] == 11),
    ("crashed tails, Q11", cr, lambda r: (r["q"] == 11) & (np.arange(151) >= 111)),
    ("crashed tails, Q25", cr, lambda r: (r["q"] == 25) & (np.arange(151) >= 111)),
    ("crashed tails, Q37", cr, lambda r: (r["q"] == 37) & (np.arange(151) >= 111)),
    ("bad-end clipped bases", groups["crash-clip"] + groups["clip-nocrash"], lambda r: r["clipped"] & (np.arange(151) >= 75)),
):
    nxt = prv = n = 0
    for r in rs:
        m = r["counted"] & sel(r) & (r["called"] != r["ref"]) & (r["rnext"] != r["ref"]) & (r["rprev"] != r["ref"]) & (r["rnext"] < 4) & (r["rprev"] < 4)
        n += m.sum(); nxt += (m & (r["called"] == r["rnext"])).sum(); prv += (m & (r["called"] == r["rprev"])).sum()
    print(f"  {label:24s} wrong calls {n:6d}: = next base {100*nxt/max(1,n):4.1f}%   = previous base {100*prv/max(1,n):4.1f}%   (chance about 33% each)")
# per clip of >= 8 bases: identity at offset 0, +1, -1
idn = []
for r in groups["crash-clip"] + groups["clip-nocrash"]:
    if r["clip3"] < 8:
        continue
    m = r["clipped"] & r["counted"] & (np.arange(151) >= 75) & (r["rnext"] < 4) & (r["rprev"] < 4)
    if m.sum() < 8:
        continue
    c = r["called"][m]
    idn.append((np.mean(c == r["ref"][m]), np.mean(c == r["rnext"][m]), np.mean(c == r["rprev"][m])))
idn = np.array(idn)
if len(idn):
    print(f"  clips >= 8 bp ({len(idn)}): identity in step {np.mean(idn[:,0]):.2f}, one ahead {np.mean(idn[:,1]):.2f}, one behind {np.mean(idn[:,2]):.2f}; "
          f"clips better one ahead or behind than in step by 0.2+: {100*np.mean(np.maximum(idn[:,1], idn[:,2]) > idn[:,0] + 0.2):.1f}%")

print("\n5. DO ERRORS COME IN RUNS? (crashed tails, last 40 cycles, counted bases with Q < 15)")
for g in ("crash-noclip", "crash-clip"):
    after_mm = [0, 0]; after_ok = [0, 0]
    for r in groups[g]:
        e = (r["called"] != r["ref"]); c = r["counted"]; low = r["q"] < 15
        for i in range(111, 150):
            if c[i] and c[i + 1] and low[i + 1]:
                tgt = after_mm if e[i] else after_ok
                tgt[0] += e[i + 1]; tgt[1] += 1
    print(f"  {g:13s} P(wrong | base before wrong) {after_mm[0]/max(1,after_mm[1]):.3f} ({after_mm[1]})   "
          f"P(wrong | base before right) {after_ok[0]/max(1,after_ok[1]):.3f} ({after_ok[1]})")

print("\n6. ERROR RATE AT Q < 15 BY HOW FAR INTO THE LOW RUN (crashed reads; the low run ending at the base)")
bins = ((1, 1), (2, 3), (4, 7), (8, 15), (16, 31), (32, 151))
for g in ("crash-noclip", "crash-clip"):
    acc = {b: [0, 0] for b in bins}
    for r in groups[g]:
        run = 0
        for i in range(151):
            run = run + 1 if r["q"][i] < 15 else 0
            if run and r["counted"][i]:
                for b in bins:
                    if b[0] <= run <= b[1]:
                        acc[b][0] += r["called"][i] != r["ref"][i]; acc[b][1] += 1
    print(f"  {g:13s} " + "  ".join(f"{a}-{b}: {acc[(a,b)][0]/max(1,acc[(a,b)][1]):.3f}" for a, b in bins))

print("\n7. CRASHED READS: WHAT SEPARATES CLIP FROM NO CLIP")
for g in ("crash-noclip", "crash-clip"):
    rs = groups[g]
    mm30 = np.array([((r["called"] != r["ref"]) & r["counted"])[-30:].sum() for r in rs])
    low30 = np.array([(r["q"][-30:] < 15).sum() for r in rs])
    lastmm = []
    for r in rs:
        e = np.where(((r["called"] != r["ref"]) & r["counted"]))[0]
        lastmm.append(150 - e[-1] if len(e) else 151)
    lastmm = np.array(lastmm)
    print(f"  {g:13s} wrong in last 30: median {np.median(mm30):.0f} (quartiles {np.percentile(mm30,25):.0f}-{np.percentile(mm30,75):.0f})   "
          f"low in last 30: median {np.median(low30):.0f}   last wrong base, cycles from end: median {np.median(lastmm):.0f}")
    print(f"  {'':13s} wrong in last 30, share of reads: " + " ".join(f"{k}:{100*np.mean(mm30==k):4.1f}" for k in range(0, 9)) + f" 9+:{100*np.mean(mm30>=9):4.1f}")
