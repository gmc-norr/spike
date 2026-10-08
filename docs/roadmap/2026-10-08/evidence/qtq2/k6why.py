"""Where does K6's over-crash come from? Crash = 10+ of the last 20 (read order) below Q15."""
import pysam, numpy as np, sys
BAM = "/home/parlar_ai/dev/spike/data/giab_hg38/HG002/HG002.GRCh38.chr20.bam"
WIN = [("chr20", 38547500, 38552500), ("chr20", 38895000, 38915000)]
def crash(r):
    q = np.array(r.query_qualities); q = q[::-1] if r.is_reverse else q
    return int((q[-20:] < 15).sum() >= 10)
def tally(name, reads):
    t = {1: [0, 0], 2: [0, 0]}
    for r in reads:
        m = 1 if r.is_read1 else 2; t[m][0] += 1; t[m][1] += crash(r)
    n = t[1][0] + t[2][0]; c = t[1][1] + t[2][1]
    print(f"{name:34s} reads {n:6d} crash {100*c/max(n,1):5.2f}%  R1 {100*t[1][1]/max(t[1][0],1):5.2f}% ({t[1][0]})  R2 {100*t[2][1]/max(t[2][0],1):5.2f}% ({t[2][0]})")
names = set(l.strip() for l in open("k6/replaced_reads.txt"))
src = pysam.AlignmentFile(BAM)
def win_reads(f, keep):
    seen = set()
    for c, a, b in WIN:
        for r in f.fetch(c, a, b):
            if r.flag & 0xF0C: continue
            k = (r.query_name, r.is_read1)
            if k in seen or not keep(r): continue
            seen.add(k); yield r
tally("orig: all reads in windows", win_reads(src, lambda r: True))
tally("orig: replaced reads in windows", win_reads(src, lambda r: r.query_name in names))
tally("orig: kept reads in windows", win_reads(src, lambda r: r.query_name not in names))
m = pysam.AlignmentFile("k6/merged.bam")
tally("merged: spike in windows", win_reads(m, lambda r: r.query_name.startswith("SPIKE_")))
tally("merged: real in windows", win_reads(m, lambda r: not r.query_name.startswith("SPIKE_")))
s = pysam.AlignmentFile("k6/sim.bam")
tally("sim.bam: all spike reads", (r for r in s.fetch(until_eof=True) if not r.flag & 0xF0C))
for c, a, b in WIN:
    tally(f"merged spike {a}", (r for r in m.fetch(c, a, b) if r.query_name.startswith("SPIKE_") and not r.flag & 0xF0C))
    tally(f"orig replaced {a}", (r for r in src.fetch(c, a, b) if r.query_name in names and not r.flag & 0xF0C))
