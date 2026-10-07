"""Which run of one letter does each of K6's bin-3 spike reads actually stand over?

For every `SPIKE_` read whose template falls in run bin 3, this lists the runs of
12+ (or C 7+) in the reference inside the read's own span, and clusters the read
starts. It exists because the first write-up of the K6 correction named a single
T13 as the cause when there are five such runs
(docs/superpowers/plans/2026-10-07-quality-model-v2.md, Correction).

Usage: bin3_runs.py SIM_BAM REFERENCE_FASTA
  SIM_BAM: spike's aligned reads for the K6 events (align.sh's sim.bam).
Needs pysam. `run_bin` and `template_bin` are k6_runbins.py's, imported from it.
"""
import collections
import sys

import pysam

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from k6_runbins import run_bin, template_bin  # noqa: E402

sim, ref = sys.argv[1:3]
fa = pysam.FastaFile(ref)

def runs_in(chrom, a, b):
    """Every run of 12+ (or C 7+) of one letter overlapping [a, b) in the reference."""
    t = fa.fetch(chrom, max(0, a - 30), b + 30).upper()
    out, base, n, start = [], "", 0, 0
    for i, ch in enumerate(t + "$"):
        if ch == base:
            n += 1
        else:
            if base and run_bin(base, n) == 3:
                out.append((base, n, max(0, a - 30) + start))
            base, n, start = ch, 1, i
    return tuple(out)

cnt = collections.Counter()
starts = []
n3 = 0
for r in pysam.AlignmentFile(sim).fetch(until_eof=True):
    if not r.query_name.startswith("SPIKE_") or r.flag & 0xF0C or r.is_unmapped:
        continue
    if template_bin(fa, r) != 3:
        continue
    n3 += 1
    lo = r.reference_start if not r.is_reverse else r.reference_end - r.query_length
    cnt[runs_in(r.reference_name, lo, lo + r.query_length)] += 1
    starts.append(lo)
print(f"bin-3 spike reads: {n3}")
for k, v in cnt.most_common():
    print(f"    {v:3d} reads  runs in span: {k}")
print("read-start clusters:")
starts.sort()
grp = [starts[0]]
for s in starts[1:]:
    if s - grp[-1] > 400:
        print(f"    {grp[0]}-{grp[-1]} ({len(grp)})"); grp = []
    grp.append(s)
print(f"    {grp[0]}-{grp[-1]} ({len(grp)})")
