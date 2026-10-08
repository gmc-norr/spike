import pysam, collections, statistics, sys
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full/forward/normal'
REC='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam'
runs={}
for r in ('run0','run1','run2'):
    for l in open(f'{S}/{r}/truth.vcf'):
        if l[0]=='#': continue
        f=l.split('\t'); runs[(f[0],int(f[1]))]=r
fp=[l.rstrip('\n').split('\t') for l in open('fwd_fp.bed')]
snps=collections.defaultdict(list)
for l in open('fwd_snp.tsv'):
    c,p,r,a,gt=l.rstrip().split('\t')
    if len(r)==1 and len(a)==1: snps[c].append((int(p),a,gt.replace('|','/')))
names={}
def repl(run):
    if run not in names: names[run]=set(l.rstrip('\n') for l in open(f'{S}/{run}/replaced_reads.txt'))
    return names[run]
rec=pysam.AlignmentFile(REC); sims={r:pysam.AlignmentFile(f'{S}/{r}/sim.bam') for r in runs.values()}
def af(reads,p0,alt):
    c=n=0
    for a in reads:
        if a.is_secondary or a.is_supplementary or a.is_duplicate or a.is_unmapped or a.mapping_quality<20: continue
        qp=None
        for q,rp in a.get_aligned_pairs(matches_only=True):
            if rp==p0: qp=q;break
        if qp is None: continue
        n+=1; c+= a.query_sequence[qp]==alt
    return c,n
out=collections.defaultdict(list); k=0
for c,s,e,tag in fp:
    g,ch,pos=tag.split(':'); pos=int(pos)
    run=runs.get((ch,pos)) or runs.get((ch,pos-1))
    if not run: continue
    for p,alt,gt in snps[c]:
        if not (int(s)<=p<=int(e)) or abs(p-pos)<60: continue
        k+=1
        if k>400: break
        p0=p-1
        before=list(rec.fetch(c,p0,p0+1)); R=repl(run)
        after=[x for x in before if x.query_name not in R]+list(sims[run].fetch(c,p0,p0+1))
        cb,nb=af(before,p0,alt); ca,na=af(after,p0,alt)
        if nb>=10 and na>=10: out[gt].append((cb/nb,ca/na))
for gt,v in out.items():
    print('SNP',gt,'n',len(v),'median AF before %.3f after %.3f'%(statistics.median(x[0] for x in v),statistics.median(x[1] for x in v)),'|change|>0.2:',sum(1 for x in v if abs(x[0]-x[1])>0.2))
