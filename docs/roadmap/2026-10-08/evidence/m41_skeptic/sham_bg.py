# Background indel read AF: before / after(sim) / after(sham), round 2b forward normal arm. Read-only.
import pysam, collections, statistics, sys
V='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/valid_lens'
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full/forward/normal'
REC='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam'
runs={}
for r in ('run0','run1','run2'):
    for l in open(f'{S}/{r}/truth.vcf'):
        if l[0]=='#': continue
        f=l.split('\t'); runs[(f[0],int(f[1]))]=r
fp=[l.rstrip('\n').split('\t') for l in open(f'{V}/fwd_fp.bed')]
ns=[l.rstrip().split('\t') for l in open(f'{V}/fwd_nonsnp.tsv')]
todo=[]
for c,s,e,tag in fp:
    g,ch,pos=tag.split(':'); pos=int(pos)
    run=runs.get((ch,pos)) or runs.get((ch,pos-1)) or runs.get((ch,pos+1))
    if not run: continue
    for c2,p,ref,alt,gt in ns:
        if c2==c and int(s)<=int(p)<=int(e) and ',' not in alt and len(ref)!=len(alt):
            todo.append((run,c,int(p),ref,alt,gt.replace('|','/'),pos))
names={}
def repl(run):
    if run not in names: names[run]=set(l.rstrip('\n') for l in open(f'{S}/{run}/replaced_reads.txt'))
    return names[run]
rec=pysam.AlignmentFile(REC)
sims={r:pysam.AlignmentFile(f'{S}/{r}/sim.bam') for r in ('run0','run1','run2')}
shams={r:pysam.AlignmentFile(f'{S}/{r}/sham.bam') for r in ('run0','run1','run2')}
def frac(reads, pos0, L, isdel):
    carry=span=0
    for a in reads:
        if a.is_secondary or a.is_supplementary or a.is_duplicate or a.is_unmapped or a.mapping_quality<20: continue
        if a.reference_start> pos0-5 or a.reference_end < pos0+L+6: continue
        span+=1
        rp=a.reference_start; hit=False
        for op,n in a.cigartuples:
            if op in (0,7,8): rp+=n
            elif op==2:
                if isdel and n==L and abs(rp-(pos0+1))<=25: hit=True
                rp+=n
            elif op==1:
                if (not isdel) and n==L and abs(rp-(pos0+1))<=25: hit=True
            elif op==3: rp+=n
        carry+=hit
    return carry,span
out=collections.defaultdict(list)
for run,c,p,ref,alt,gt,evpos in todo:
    L=abs(len(ref)-len(alt)); isdel=len(ref)>len(alt); p0=p-1
    if abs(p-evpos)<60: continue
    before=list(rec.fetch(c,p0-200,p0+200)); R=repl(run)
    kept=[a for a in before if a.query_name not in R]
    cb,sb=frac(before,p0,L,isdel)
    ca,sa=frac(kept+list(sims[run].fetch(c,p0-200,p0+200)),p0,L,isdel)
    cs,ss=frac(kept+list(shams[run].fetch(c,p0-200,p0+200)),p0,L,isdel)
    if sb>=10 and sa>=10 and ss>=10:
        out[gt].append((cb/sb, ca/sa, cs/ss, abs(p-evpos), L))
for gt in ('1/1','0/1'):
    v=out[gt]
    print(gt,'n',len(v),'median before %.3f sim %.3f sham %.3f'%tuple(statistics.median(x[i] for x in v) for i in range(3)))
    print('   |sim-before|>0.2: %d   |sham-before|>0.2: %d'%(sum(abs(x[1]-x[0])>0.2 for x in v),sum(abs(x[2]-x[0])>0.2 for x in v)))
    if gt=='1/1':
        print('   after<0.75: sim %d sham %d before %d'%(sum(x[1]<0.75 for x in v),sum(x[2]<0.75 for x in v),sum(x[0]<0.75 for x in v)))
    for lo,hi in ((60,1000),(1000,1500),(1500,2001)):
        w=[x for x in v if lo<=x[3]<hi]
        if w: print('   dist %d-%d n %d median before %.3f sim %.3f sham %.3f'%(lo,hi,len(w),*[statistics.median(x[i] for x in w) for i in range(3)]))
    for lo,hi in ((1,2),(2,5),(5,20),(20,1000)):
        w=[x for x in v if lo<=x[4]<hi]
        if w: print('   len %d-%d n %d median before %.3f sim %.3f'%(lo,hi-1,len(w),statistics.median(x[0] for x in w),statistics.median(x[1] for x in w)))
