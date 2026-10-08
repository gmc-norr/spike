import pysam
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full/forward/normal'
sim=pysam.AlignmentFile(f'{S}/run0/sim.bam')
ev=[]
for l in open(f'{S}/run0/truth.vcf'):
    if l[0]=='#': continue
    f=l.split('\t'); ev.append((f[0],int(f[1])))
import itertools
print('example names', [a.query_name for a in itertools.islice(sim.fetch(ev[0][0],ev[0][1]-10,ev[0][1]+10),3)])
for c,p in ev[:8]:
    bins={}
    lo=hi=None
    for a in sim.fetch(c,max(0,p-6000),p+6000):
        if a.is_secondary or a.is_supplementary or a.mapping_quality<20: continue
        syn=a.query_name.startswith('ev')
        if syn:
            lo=a.reference_start if lo is None else min(lo,a.reference_start); hi=a.reference_end if hi is None else max(hi,a.reference_end)
        d=abs((a.reference_start+a.reference_end)//2-p)//500*500
        t=bins.setdefault(d,[0,0]); t[1]+=1; t[0]+=syn
    print(c,p,'synthetic span',lo-p if lo else None,hi-p if hi else None,' syn frac:',' '.join(f'{d}:{x[0]/x[1]:.2f}' for d,x in sorted(bins.items())))
