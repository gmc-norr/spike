# How far from the event do replaced reads reach? (read-only)
import pysam, sys
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full/forward/normal'
REC='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam'
R=set(l.rstrip('\n') for l in open(f'{S}/run0/replaced_reads.txt'))
rec=pysam.AlignmentFile(REC)
ev=[]
for l in open(f'{S}/run0/truth.vcf'):
    if l[0]=='#': continue
    f=l.split('\t'); ev.append((f[0],int(f[1]),len(f[3])))
for c,p,L in ev[:8]:
    lo=hi=None; n=0; tot=0
    for a in rec.fetch(c,max(0,p-6000),p+6000):
        if a.is_secondary or a.is_supplementary: continue
        if a.mapping_quality<20 or a.is_duplicate: continue
        if a.query_name in R:
            n+=1
            lo=a.reference_start if lo is None else min(lo,a.reference_start)
            hi=a.reference_end if hi is None else max(hi,a.reference_end)
    # fraction replaced in bins of distance
    bins={}
    for a in rec.fetch(c,max(0,p-6000),p+6000):
        if a.is_secondary or a.is_supplementary or a.mapping_quality<20 or a.is_duplicate: continue
        d=abs((a.reference_start+a.reference_end)//2-p)//500*500
        t=bins.setdefault(d,[0,0]); t[1]+=1; t[0]+= a.query_name in R
    print(c,p,'replaced span',lo-p if lo else None,hi-p if hi else None,' frac by |dist| bin:',' '.join(f'{d}:{x[0]/x[1]:.2f}' for d,x in sorted(bins.items())))
