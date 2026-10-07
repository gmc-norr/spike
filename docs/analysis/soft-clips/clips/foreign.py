"""What are the 'foreign' soft clips (clipped bases < 50% identical to the reference where they
would align, not at a shared spot, not mapping elsewhere with MAPQ >= 20)?
Per clip, in priority order:
 short      1-4 bases: too few to tell foreign from a few errors
 adapter2   an Illumina adapter/primer 12-mer anywhere in the clip, or the read's own sequence
            spells the adapter start where its fragment ends (fragment shorter than the read)
 chimericpair  the mate maps to another chromosome or > 1,500 bp away: a chimeric fragment
 samplevar  within 150 bp of a >= 10 bp HG002 variant in the T2T-Q100 benchmark: the sample's own sequence
 repeatmap  the clip (>= 20 bp) aligns somewhere with MAPQ < 20: repeat sequence
 lowcomp    one base >= 60% of the clip, or a 1-3 bp unit repeated over >= 80% of it
 rest       none of the above"""
import os, pickle, collections, pysam, numpy as np
SLICE = os.environ["SLICE"]; Q100 = os.environ["Q100_VCF"]
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; cls = pickle.load(open("cls2.pkl", "rb"))
idx = [i for i, k in enumerate(cls) if k == "foreign"]
need = {(C[i]["name"], C[i]["r1"]) for i in idx}
info = {}
comp = str.maketrans("ACGTN", "TGCAN")
for r in pysam.AlignmentFile(SLICE):
    if r.flag & 0xF0C: continue
    k = (r.query_name, r.is_read1)
    if k in need:
        s = r.query_sequence
        info[k] = dict(mchrom=r.next_reference_name, mpos=r.next_reference_start, pos=r.reference_start,
                       seq=s.translate(comp)[::-1] if r.is_reverse else s, tlen=abs(r.template_length))
ADS = {"TruSeq R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", "TruSeq R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
       "P5": "AATGATACGGCGACCACCGAGATCTACAC", "P7": "CAAGCAGAAGACGGCATACGAGAT", "Nextera": "CTGTCTCTTATACACATCT"}
kmers = set()
for a in ADS.values():
    for s in (a, a.translate(comp)[::-1]):
        kmers |= {s[j:j + 12] for j in range(len(s) - 11)}
def adapter_anywhere(s): return any(s[j:j + 12] in kmers for j in range(len(s) - 11))
def adapter_at_fragment_end(c, inf):
    t = inf["tlen"]
    if not (0 < t < 151): return False
    tail = inf["seq"][t:t + 13]
    n = len(tail)
    return n >= 5 and sum(a != b for a, b in zip(tail, "AGATCGGAAGAGC"[:n])) <= (1 if n >= 10 else 0)
# HG002 Q100 variants >= 10 bp on chr20 in the slice
sv = []
for rec in pysam.VariantFile(Q100).fetch("chr20", 38400000, 40300000):
    L = max(len(rec.ref), max(len(a) for a in rec.alts if not a.startswith("<")) if rec.alts else 0)
    if abs(len(rec.ref) - min(len(a) for a in rec.alts)) >= 10 or L >= 10: sv.append(rec.pos)
sv = np.array(sorted(sv))
def near_sv(b): j = np.searchsorted(sv, b); return any(abs(sv[k] - b) <= 150 for k in (j - 1, j) if 0 <= k < len(sv))
remap_any = {}
for r in pysam.AlignmentFile("clips_ge20.bam"):
    if not (r.is_secondary or r.is_supplementary or r.is_unmapped): remap_any[int(r.query_name)] = r.mapping_quality
def lowcomp(s):
    if max(s.count(b) for b in "ACGT") / len(s) >= 0.6: return True
    for u in (2, 3):
        unit = s[:u]
        if sum(s[j:j + u] == unit for j in range(0, len(s) - u + 1, u)) * u >= 0.8 * len(s): return True
    return False
order = ["short", "adapter2", "chimericpair", "samplevar", "repeatmap", "lowcomp", "rest"]
out = []
for i in idx:
    c = C[i]; inf = info[(c["name"], c["r1"])]; s = c["seq_seqorient"]
    if c["L"] <= 4: k = "short"
    elif adapter_anywhere(s) or (c["three"] and adapter_at_fragment_end(c, inf)): k = "adapter2"
    elif inf["mchrom"] != "chr20" or abs(inf["mpos"] - inf["pos"]) > 1500: k = "chimericpair"
    elif near_sv(c["boundary"]): k = "samplevar"
    elif i in remap_any: k = "repeatmap"
    elif lowcomp(s): k = "lowcomp"
    else: k = "rest"
    out.append((i, k))
pickle.dump(out, open("foreign_cls.pkl", "wb"))
cnt = collections.Counter(k for _, k in out); n = len(out)
print(f"foreign clips {n}; HG002 Q100 variants >= 10 bp in the slice: {len(sv)}")
print(f"{'kind':13s} {'clips':>6s} {'share':>6s} {'median L':>9s} {'mean Q':>7s} {'3prime':>7s}")
for k in order:
    ii = [i for i, x in out if x == k]
    if not ii: continue
    print(f"{k:13s} {len(ii):6d} {100*len(ii)/n:5.1f}% {np.median([C[i]['L'] for i in ii]):9.0f} {np.mean([C[i]['meanq'] for i in ii]):7.1f} {100*np.mean([C[i]['three'] for i in ii]):6.0f}%")
rest = [C[i] for i, x in out if x == "rest"]
print("\nrest, examples (L, end, mean Q, clip):")
for c in rest[:12]: print(f"  {c['L']:3d} {'3p' if c['three'] else '5p'} {c['meanq']:5.1f} {c['seq_seqorient'][:60]}")
