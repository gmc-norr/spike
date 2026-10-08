"""Reported part of K1 (v2 plan): the crash share by the template's longest-run bin. From k1.py: spike's reads against the real reads
in the 25 windows of the SNV run (each 5 kb, centred on an SNV), primary, non-duplicate, both
mates present or not. Pass: read-mean SD ratio 0.85-1.15; |z| < 3 for the share of perfect reads
(every base at the real reads' top quality) and of crashed reads (>= 10 of the last 20 < Q15).
Usage: k1.py MERGED_BAM"""
import sys, math, numpy as np, pysam
def z(k1, n1, k2, n2):
    p = (k1 + k2) / (n1 + n2); se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2)); return (k1 / n1 - k2 / n2) / se if se else 0.0
regs = [("chr20", p - 2500, p + 2500) for p in (38550000 + i * 60000 for i in range(25))]
reads = {"spike": {}, "real": {}}
REF = pysam.FastaFile("/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
comp_t = str.maketrans('ACGTN', 'TGCAN')
def hp_bin(b, n):
    if b == 'C': return 3 if n >= 7 else 2 if n >= 5 else 0
    return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
def read_hp(t):
    best, prev, n = 0, '', 0
    for b in t:
        n = n + 1 if b == prev and b != 'N' else (b != 'N'); prev = b; best = max(best, hp_bin(b, n))
    return best
bins = {"spike": {}, "real": {}}
bam = pysam.AlignmentFile(sys.argv[1])
for c, a, b in regs:
    for r in bam.fetch(c, a, b):
        if r.flag & 0xF0C: continue
        kind = "spike" if r.query_name.startswith("SPIKE_") else "real"
        q = np.array(r.query_qualities)
        if r.is_reverse: q = q[::-1]
        reads[kind].setdefault(r.query_name, {})[1 if r.is_read1 else 2] = q
        cig = r.cigartuples; lead = cig[0][1] if cig[0][0] == 4 else 0; trail = cig[-1][1] if cig[-1][0] == 4 else 0
        t = REF.fetch(c, max(0, r.reference_start - lead), r.reference_end + trail).upper()
        h = read_hp(t.translate(comp_t)[::-1] if r.is_reverse else t)
        acc = bins[kind].setdefault(h, [0, 0]); acc[0] += 1; acc[1] += int(len(q) >= 20 and (q[-20:] < 15).sum() >= 10)
out = {}
top = max(int(q.max()) for d in reads["real"].values() for q in d.values())
for kind, d in reads.items():
    allq = [q for m in d.values() for q in m.values()]
    means = np.array([q.mean() for q in allq])
    perfect = sum(int((q == top).all()) for q in allq)
    crash = sum(int((q[-20:] < 15).sum() >= 10) for q in allq if len(q) >= 20)
    pairs = [(m[1].mean(), m[2].mean()) for m in d.values() if 1 in m and 2 in m]
    corr = np.corrcoef(np.array(pairs).T)[0, 1]
    out[kind] = dict(n=len(allq), sd=means.std(ddof=1), perfect=perfect, crash=crash, corr=corr, lt33=(means < 33).mean())
s, r = out["spike"], out["real"]
ratio = s["sd"] / r["sd"]; zp = z(s["perfect"], s["n"], r["perfect"], r["n"]); zc = z(s["crash"], s["n"], r["crash"], r["n"])
print(f"top quality Q{top}; reads: spike {s['n']}, real {r['n']}")
print(f"read-mean SD: spike {s['sd']:.3f}, real {r['sd']:.3f}, ratio {ratio:.3f} -> {'PASS' if 0.85 <= ratio <= 1.15 else 'FAIL'}")
print(f"perfect reads: spike {s['perfect']}/{s['n']} = {100*s['perfect']/s['n']:.2f}%, real {r['perfect']}/{r['n']} = {100*r['perfect']/r['n']:.2f}%, z {zp:+.2f} -> {'PASS' if abs(zp) < 3 else 'FAIL'}")
print(f"crashed reads: spike {s['crash']}/{s['n']} = {100*s['crash']/s['n']:.2f}%, real {r['crash']}/{r['n']} = {100*r['crash']/r['n']:.2f}%, z {zc:+.2f} -> {'PASS' if abs(zc) < 3 else 'FAIL'}")
print(f"reported: R1-R2 correlation spike {s['corr']:.3f}, real {r['corr']:.3f}; mean < Q33 spike {100*s['lt33']:.2f}%, real {100*r['lt33']:.2f}%")
print("K1:", "PASS" if 0.85 <= ratio <= 1.15 and abs(zp) < 3 and abs(zc) < 3 else "FAIL")
for kind in ("spike", "real"):
    print(f"crash by run bin, {kind}: " + "  ".join(f"{h}: {100*bins[kind].get(h,[0,0])[1]/max(1,bins[kind].get(h,[0,0])[0]):.2f}% ({bins[kind].get(h,[0,0])[0]})" for h in range(4)))
