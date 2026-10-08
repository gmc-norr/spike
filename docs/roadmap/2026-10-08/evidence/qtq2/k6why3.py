"""Crash share by template run bin: spike's remade reads vs the donors in K6's windows."""
import pysam, numpy as np
REF = "/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
BAM = "/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.GRCh38.chr20.bam"
WIN = [("chr20", 38547500, 38552500), ("chr20", 38895000, 38915000)]
fa = pysam.FastaFile(REF)
comp = str.maketrans("ACGTN", "TGCAN")
def run_bin(b, n):
    if b == "C": return 3 if n >= 7 else 2 if n >= 5 else 0
    if b in "AGT": return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
    return 0
def tbin(r):
    L = r.query_length
    if r.is_reverse:
        e = r.reference_end; t = fa.fetch(r.reference_name, max(0, e - L), e).upper().translate(comp)[::-1]
    else:
        s = r.reference_start; t = fa.fetch(r.reference_name, s, s + L).upper()
    best, base, n = 0, "N", 0
    for b in t:
        if b == base and b != "N": n += 1
        else: base, n = b, int(b != "N")
        best = max(best, run_bin(base, n))
    return best
def crash(r):
    q = np.array(r.query_qualities); q = q[::-1] if r.is_reverse else q
    return int((q[-20:] < 15).sum() >= 10)
def table(name, reads):
    t = np.zeros((2, 4, 2), int)
    for r in reads:
        t[0 if r.is_read1 else 1, tbin(r), 0] += 1; t[0 if r.is_read1 else 1, tbin(r), 1] += crash(r)
    for m in range(2):
        cells = "  ".join(f"{100*t[m,h,1]/t[m,h,0]:5.1f}% ({t[m,h,0]:4d})" if t[m,h,0] else "    -        " for h in range(4))
        print(f"{name:22s} R{m+1}: {cells}")
names = set(l.strip() for l in open("k6/replaced_reads.txt"))
src = pysam.AlignmentFile(BAM)
def win(f, keep):
    seen = set()
    for c, a, b in WIN:
        for r in f.fetch(c, a, b):
            if r.flag & 0xF0C or r.is_unmapped: continue
            k = (r.query_name, r.is_read1)
            if k in seen or not keep(r): continue
            seen.add(k); yield r
print("crash share by run bin 0 / 1 / 2 / 3 (reads)")
table("donors in windows", win(src, lambda r: True))
table("spike remade", win(pysam.AlignmentFile("k6/merged.bam"), lambda r: r.query_name.startswith("SPIKE_")))
