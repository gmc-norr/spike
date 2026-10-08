"""K2 of the duplicates fix (docs/superpowers/plans/2026-10-04-duplicates.md): only duplicates are
added, and the right ones. Usage: k2.py OLD_OUT NEW_OUT SOURCE_BAM EVENTS_VCF"""
import gzip
import subprocess
import sys

import pysam

old, new, source, events = sys.argv[1:5]
FLANK, PAD = 10_000, 1_500
ok = True


def check(label, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"{'PASS' if cond else 'FAIL'}  {label}{(': ' + detail) if detail else ''}", flush=True)


def lines(path):
    with (gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")) as f:
        return f.read()


for f in ("R1.fq.gz", "R2.fq.gz", "truth.vcf"):
    check(f"{f} identical", lines(f"{old}/{f}") == lines(f"{new}/{f}"))
views = [subprocess.run(["samtools", "view", f"{d}/sim.bam"], capture_output=True, check=True).stdout for d in (old, new)]
check("sim.bam records identical", views[0] == views[1], f"{views[0].count(b'\n')} records")

names = lambda d, f: set(open(f"{d}/{f}").read().split())  # noqa: E731
old_rep, new_rep = names(old, "replaced_reads.txt"), names(new, "replaced_reads.txt")
old_fq, new_fq = names(old, "fastq_removed_reads.txt"), names(new, "fastq_removed_reads.txt")
added = new_rep - old_rep
check("replaced_reads.txt new = old plus added", old_rep <= new_rep and len(added) >= 1,
      f"old {len(old_rep)}, new {len(new_rep)}, added {len(added)}")
check("fastq_removed_reads.txt new = old plus the same names", new_fq == old_fq | added,
      f"old {len(old_fq)}, new {len(new_fq)}")

# The spans spike reads: each event's extraction window (--flank 10000), widened by 1,500 bp.
spans = []
with pysam.VariantFile(events) as vcf:
    for rec in vcf:
        start0 = rec.pos if rec.alts[0] == "<DEL>" else rec.pos - 1
        end0 = rec.stop if rec.alts[0] == "<DEL>" else rec.pos - 1 + len(rec.ref)
        spans.append((rec.chrom, max(0, start0 - FLANK - PAD), end0 + FLANK + PAD))


def clipped(ops):
    n = 0
    for op, length in ops:
        if op not in (4, 5):
            break
        n += length
    return n


def five_prime(r):
    """Unclipped 5' end, written independently of spike's: a forward read's start minus its
    leading clips, a reverse read's end plus its trailing clips."""
    if r.is_reverse:
        return (r.reference_name, r.reference_end + clipped(reversed(r.cigartuples)), True)
    return (r.reference_name, r.reference_start - clipped(r.cigartuples), False)


mates, seen = {}, set()
with pysam.AlignmentFile(source) as bam:
    for chrom, s, e in spans:
        for r in bam.fetch(chrom, s, e):
            if r.flag & (0x4 | 0x100 | 0x200 | 0x800) or (r.query_name, r.is_read1) in seen:
                continue
            seen.add((r.query_name, r.is_read1))
            mates.setdefault(r.query_name, []).append((five_prime(r), r.is_duplicate, r.is_read1))

removed_families, copies, source_dups = set(), {}, set()
for name, ms in mates.items():
    if any(d for _, d, _ in ms):
        source_dups.add(name)
    if len(ms) != 2 or ms[0][2] == ms[1][2]:
        continue
    key = tuple(sorted(m[0] for m in ms))
    if all(d for _, d, _ in ms):
        copies[name] = key
    elif not any(d for _, d, _ in ms) and name in old_fq:
        removed_families.add(key)
expected = {n for n, k in copies.items() if k in removed_families}
check("every added name is a source duplicate", added <= source_dups, f"{len(added - source_dups)} not")
check("100% of added names are duplicates of a removed kept pair (independent key)", added <= expected,
      f"{len(added & expected)} of {len(added)}")
share = len(expected & added) / len(expected) if expected else 0.0
check(">= 99% of such duplicate pairs are added", share >= 0.99, f"{len(expected & added)} of {len(expected)} = {share:.4f}")
print(f"K2: {'PASS' if ok else 'FAIL'}", flush=True)
