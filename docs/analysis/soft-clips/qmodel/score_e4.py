import pysam, collections, math
def wilson(k,n,z=1.96):
    c=(k+z*z/2)/(n+z*z); h=z*math.sqrt(k*(n-k)/n+z*z/4)/(n+z*z); return 100*(c-h),100*(c+h)
st=collections.defaultdict(lambda: collections.Counter())
for r in pysam.AlignmentFile("e4.bam"):
    if r.is_secondary or r.is_supplementary: continue
    m=r.query_name.split("_")[0]; s=st[m]; s["reads"]+=1
    if r.is_unmapped: s["unmapped"]+=1; continue
    clip=any(op==4 for op,_ in r.cigartuples)
    s["clipped"]+=clip; s["nm"]+=r.get_tag("NM"); s["mapped"]+=1; s["mapq0"]+=r.mapping_quality==0
    s["sa"]+=r.has_tag("SA")
print(f"{'model':8s} {'reads':>6s} {'soft-clipped % (95% CI)':>28s} {'mean NM':>8s} {'with SA %':>9s}")
for m in ("real","chain","copyQ","copyQE"):
    s=st[m]; lo,hi=wilson(s["clipped"],s["mapped"])
    print(f"{m:8s} {s['reads']:6d} {100*s['clipped']/s['mapped']:10.2f} ({lo:.2f}-{hi:.2f}) n={s['clipped']:4d} {s['nm']/s['mapped']:8.3f} {100*s['sa']/s['mapped']:9.2f}")
