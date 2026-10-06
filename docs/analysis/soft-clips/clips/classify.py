"""Why does each real soft clip exist? Checks in priority order:
 adapter  3' clip that starts with the TruSeq adapter (AGATCGGAAGAGC; >= 6 bases, <= 1 mismatch per 10)
 site     >= 3 other reads clip at the same boundary (+-2 bp): a property of the sample/site, not the read
 chimera  SA tag, or the clipped part (>= 20 bp) realigns elsewhere with MAPQ >= 20
 polyG    >= 80% G (NovaSeq no-signal), >= 5 bases
 lowQ     mean quality of the clipped bases < 20
 nearref  the clipped bases match the reference where they would align (>= 75%): a few errors/variants at the end
 other    high-quality non-reference sequence, explained by none of the above"""
import pickle, collections, pysam, numpy as np
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; total = D["total"]
AD = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
def adapter(c):
    if not c["three"] or c["L"] < 6: return False
    s = c["seq_seqorient"]; n = min(len(s), len(AD)); mm = sum(a != b for a, b in zip(s[:n], AD[:n]))
    return mm <= n // 10
bcount = collections.Counter((c["side"], c["boundary"]) for c in C)
def site(c):
    return sum(bcount[(c["side"], c["boundary"] + d)] for d in range(-2, 3)) - 1 >= 3
remap = {}
for r in pysam.AlignmentFile("clips_ge20.bam"):
    if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
    remap[int(r.query_name)] = r.mapping_quality
order = ["adapter", "site", "chimera", "polyG", "lowQ", "nearref", "other"]
cls = []
for i, c in enumerate(C):
    if adapter(c): k = "adapter"
    elif site(c): k = "site"
    elif c["sa"] or remap.get(i, -1) >= 20: k = "chimera"
    elif c["L"] >= 5 and c["seq_seqorient"].count("G") / c["L"] >= 0.8: k = "polyG"
    elif c["meanq"] < 20: k = "lowQ"
    elif c["ident"] >= 0.75: k = "nearref"
    else: k = "other"
    cls.append(k)
cnt = collections.Counter(cls)
print(f"reads {total}, clips {len(C)}")
print(f"{'cause':9s} {'clips':>6s} {'share':>7s} {'3prime':>7s} {'median L':>9s} {'mean Q':>7s} {'ident':>6s} {'tlen<151':>9s}")
for k in order:
    idx = [i for i, x in enumerate(cls) if x == k]
    if not idx: print(f"{k:9s} {0:6d}"); continue
    L = np.median([C[i]["L"] for i in idx]); mq = np.mean([C[i]["meanq"] for i in idx]); idn = np.mean([C[i]["ident"] for i in idx])
    th = np.mean([C[i]["three"] for i in idx]); short = np.mean([0 < C[i]["tlen"] < 151 for i in idx])
    print(f"{k:9s} {len(idx):6d} {100*len(idx)/len(C):6.1f}% {100*th:6.0f}% {L:9.0f} {mq:7.1f} {idn:6.2f} {100*short:8.0f}%")
# by length bin
print("\nshare of clips by cause within length bins:")
bins = [(1, 4), (5, 19), (20, 151)]
for a, b in bins:
    idx = [i for i, c in enumerate(C) if a <= c["L"] <= b]; cc = collections.Counter(cls[i] for i in idx)
    print(f"  {a}-{b} bp (n={len(idx)}): " + ", ".join(f"{k} {100*cc[k]/len(idx):.0f}%" for k in order if cc[k]))
pickle.dump(cls, open("cls.pkl", "wb"))
