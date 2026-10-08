"""Gate B mini: write recipient and spiked-view slices around background indel sites.

usage: gateb_view.py RUN RECIP SITES_TSV OUTDIR
SITES_TSV: chrom<TAB>pos (1-based) per line.
Writes OUTDIR/recip.bam and OUTDIR/view.bam (recipient minus replaced plus sim), sorted+indexed.
"""
import sys, pysam

RUN, RECIP, SITES, OUT = sys.argv[1:5]
replaced = set(l.strip() for l in open(f"{RUN}/replaced_reads.txt"))
recip = pysam.AlignmentFile(RECIP)
sim = pysam.AlignmentFile(f"{RUN}/sim.bam")
rg = recip.header.to_dict()["RG"][0]["ID"]
sites = [l.split() for l in open(SITES) if l.strip()]

ro = pysam.AlignmentFile(f"{OUT}/recip.unsorted.bam", "wb", template=recip)
vo = pysam.AlignmentFile(f"{OUT}/view.unsorted.bam", "wb", template=recip)
seen_r, seen_v = set(), set()
for chrom, pos in sites:
    p = int(pos) - 1
    for r in recip.fetch(chrom, max(0, p - 5), p + 60):
        key = (r.query_name, r.flag, r.reference_start)
        if key in seen_r:
            continue
        seen_r.add(key)
        ro.write(r)
        if r.query_name not in replaced:
            vo.write(r)
    for r in sim.fetch(chrom, max(0, p - 5), p + 60):
        key = (r.query_name, r.flag, r.reference_start)
        if key in seen_v:
            continue
        seen_v.add(key)
        a = pysam.AlignedSegment.from_dict(r.to_dict(), vo.header)
        a.set_tag("RG", rg)
        vo.write(a)
ro.close(); vo.close()
for n in ("recip", "view"):
    pysam.sort("-o", f"{OUT}/{n}.bam", f"{OUT}/{n}.unsorted.bam")
    pysam.index(f"{OUT}/{n}.bam")
print(len(sites), "sites", len(seen_r), "recip records", len(seen_v), "sim records")
