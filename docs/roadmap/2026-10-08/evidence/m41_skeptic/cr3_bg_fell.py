# Background non-SNP allele fraction before vs after spiking, at HG002 v4.2.1 indels within +-2 kb
# of round 2b forward (normal arm) events. Reads: recipient = HG002 hospital BAM; after = recipient
# minus run's replaced_reads.txt plus run's sim.bam. Seen-not-judged measurement for an idea list.
import pysam, collections, statistics, sys
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full/forward/normal'
REC='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam'
runs={}
for r in ('run0','run1','run2'):
    for l in open(f'{S}/{r}/truth.vcf'):
        if l[0]=='#': continue
        f=l.split('\t'); runs[(f[0],int(f[1]))]=r
fp=[l.rstrip('\n').split('\t') for l in open('../valid_lens/fwd_fp.bed')]
ns=[l.rstrip().split('\t') for l in open('../valid_lens/fwd_nonsnp.tsv')]
# assign each background record to a footprint whose event is in some run
todo=[]
for c,s,e,tag in fp:
    g,ch,pos=tag.split(':'); pos=int(pos)
    # truth POS for indels is preceding base == events.tsv pos; match within 1
    run=runs.get((ch,pos)) or runs.get((ch,pos-1)) or runs.get((ch,pos+1))
    if not run: continue
    for c2,p,ref,alt,gt in ns:
        if c2==c and int(s)<=int(p)<=int(e) and ',' not in alt and len(ref)!=len(alt):
            todo.append((run,c,int(p),ref,alt,gt.replace('|','/'),pos))
print('records',len(todo),file=sys.stderr)
names={}
def repl(run):
    if run not in names: names[run]=set(l.rstrip('\n') for l in open(f'{S}/{run}/replaced_reads.txt'))
    return names[run]
rec=pysam.AlignmentFile(REC)
sims={r:pysam.AlignmentFile(f'{S}/{r}/sim.bam') for r in ('run0','run1','run2')}
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
    if abs(p-evpos)<60: continue  # skip records overlapping the planted event itself
    before=list(rec.fetch(c,p0-200,p0+200))
    cb,sb=frac(before,p0,L,isdel)
    R=repl(run)
    after=[a for a in before if a.query_name not in R]+list(sims[run].fetch(c,p0-200,p0+200))
    ca,sa=frac(after,p0,L,isdel)
    if sb>=10 and sa>=10:
        out[gt].append((cb/sb, ca/sa, sb, sa, c, p, L, evpos))
for gt,v in out.items():
    fb=[x[0] for x in v]; fa=[x[1] for x in v]
    print(gt,'n',len(v),'median AF before %.3f after %.3f'%(statistics.median(fb),statistics.median(fa)),
          'drop>0.2:',sum(1 for x in v if x[0]-x[1]>0.2),'median depth before %d after %d'%(statistics.median(x[2] for x in v),statistics.median(x[3] for x in v)))
    if gt=='1/1':
        print('  hom with after-AF < 0.75:',sum(1 for x in v if x[1]<0.75),'of',len(v),'; before-AF<0.75:',sum(1 for x in v if x[0]<0.75))
homfp=set((x[4],x[7]) for x in out['1/1'] if x[1]<0.75)
hetfp=set((x[4],x[7]) for x in out['0/1'] if x[0]-x[1]>0.2)
nev=sum(1 for c,s,e,tag in fp if (runs.get((tag.split(':')[1],int(tag.split(':')[2]))) or runs.get((tag.split(':')[1],int(tag.split(':')[2])-1))))
print('events found in runs',nev,'footprints with a hom bg indel now <0.75:',len(homfp),'with a het bg indel dropped >0.2:',len(hetfp),'union',len(homfp|hetfp))
fell=set((x[4],x[7]) for x in out['1/1'] if x[0]>=0.75 and x[1]<0.75)
already=set((x[4],x[7]) for x in out['1/1'] if x[0]<0.75)
recs=[x for x in out['1/1'] if x[0]>=0.75 and x[1]<0.75]
print('footprints with a hom bg indel going from >=0.75 to <0.75:',len(fell),'records',len(recs),'; footprints with a hom already <0.75 before:',len(already), '; of 137, ones with no fell record:',len(homfp-fell))
import collections
print('after AF of fell records: <0.6',sum(1 for x in recs if x[1]<0.6),'of',len(recs))
