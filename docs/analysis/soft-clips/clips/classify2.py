"""Why does each real soft clip exist? (v2: both TruSeq adapters with a 0-3 bp offset; degraded ends by identity.)
 adapter   3' clip holding the adapter start AGATCGGAAGAGC (shared by the R1 and R2 adapters) at offset 0-3
           (<= 1 mismatch, >= 5 bases compared); or a 1-5 bp 3' clip on a fragment shorter than the read
 site      >= 3 other reads clip at the same boundary (+-2 bp): the sample or the site, not the read
 chimera   clipped bases differ from the reference here (identity < 0.5) and the read has SA, or the clip
           (>= 20 bp) realigns elsewhere with MAPQ >= 20
 foreign   identity < 0.5, maps nowhere (repeats, unplaced sequence)
 badend    identity >= 0.5 and mean clipped quality < 30: reference sequence read badly
 hiQerr    identity >= 0.5 and mean quality >= 30: reference sequence with a few high-quality mismatches
 polyG     >= 80% G over >= 5 bases (checked before badend)"""
import pickle, collections, pysam, numpy as np
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; total = D["total"]
AD = "AGATCGGAAGAGC"
def adapter(c):
    if not c["three"]: return False
    s = c["seq_seqorient"]
    for k in range(0, 4):
        t = s[k:]; n = min(len(t), len(AD))
        if n >= 5 and sum(a != b for a, b in zip(t[:n], AD[:n])) <= (1 if n >= 10 else 0): return True
    return c["L"] <= 5 and 0 < c["tlen"] < 151
bcount = collections.Counter((c["side"], c["boundary"]) for c in C)
site = lambda c: sum(bcount[(c["side"], c["boundary"] + d)] for d in range(-2, 3)) - 1 >= 3
remap = {}
for r in pysam.AlignmentFile("clips_ge20.bam"):
    if not (r.is_secondary or r.is_supplementary or r.is_unmapped): remap[int(r.query_name)] = r.mapping_quality
order = ["adapter", "site", "chimera", "foreign", "polyG", "badend", "hiQerr"]
cls = []
for i, c in enumerate(C):
    if adapter(c): k = "adapter"
    elif site(c): k = "site"
    elif c["L"] >= 5 and c["seq_seqorient"].count("G") / c["L"] >= 0.8: k = "polyG"
    elif c["ident"] < 0.5: k = "chimera" if (c["sa"] or remap.get(i, -1) >= 20) else "foreign"
    elif c["meanq"] < 30: k = "badend"
    else: k = "hiQerr"
    cls.append(k)
pickle.dump(cls, open("cls2.pkl", "wb"))
# per read: a read's cause = its first clip's cause in this priority
prio = {k: i for i, k in enumerate(order)}; rc = {}
for c, k in zip(C, cls):
    key = (c["name"], c["r1"]); rc[key] = k if key not in rc or prio[k] < prio[rc[key]] else rc[key]
rcnt = collections.Counter(rc.values())
print(f"primary reads {total}; soft-clipped reads {len(rc)} ({100*len(rc)/total:.2f}%); clips {len(C)}")
print(f"{'cause':8s} {'reads':>6s} {'% of clipped reads':>19s} {'% of all reads':>15s} {'median L':>9s} {'mean Q':>7s} {'ident':>6s} {'3prime':>7s}")
for k in order:
    idx = [i for i, x in enumerate(cls) if x == k]
    L = np.median([C[i]["L"] for i in idx]); mq = np.mean([C[i]["meanq"] for i in idx]); idn = np.mean([C[i]["ident"] for i in idx]); th = np.mean([C[i]["three"] for i in idx])
    print(f"{k:8s} {rcnt[k]:6d} {100*rcnt[k]/len(rc):18.1f}% {100*rcnt[k]/total:14.3f}% {L:9.0f} {mq:7.1f} {idn:6.2f} {100*th:6.0f}%")
fr = collections.Counter()
for c in C:
    pass
