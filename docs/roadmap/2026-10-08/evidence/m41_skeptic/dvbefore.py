import subprocess, collections
V='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/valid_lens'
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full/forward/normal'
DV='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/call_snv/genome/GM24385_seracare_cancer_snv.vcf.gz'
runs={}
for r in ('run0','run1','run2'):
    for l in open(f'{S}/{r}/truth.vcf'):
        if l[0]=='#': continue
        f=l.split('\t'); runs[(f[0],int(f[1]))]=r
fp=[l.rstrip('\n').split('\t') for l in open(f'{V}/fwd_fp.bed')]
ns=[l.rstrip().split('\t') for l in open(f'{V}/fwd_nonsnp.tsv')]
recs=set()
for c,s,e,tag in fp:
    g,ch,pos=tag.split(':'); pos=int(pos)
    if not (runs.get((ch,pos)) or runs.get((ch,pos-1)) or runs.get((ch,pos+1))): continue
    for c2,p,ref,alt,gt in ns:
        if c2==c and int(s)<=int(p)<=int(e) and ',' not in alt and len(ref)!=len(alt) and abs(int(p)-pos)>=60:
            recs.add((c,int(p),ref,alt,gt.replace('|','/')))
cnt=collections.Counter()
for c,p,ref,alt,gt in sorted(recs):
    o=subprocess.run(['bcftools','query','-r',f'{c}:{p}-{p}','-f','%POS\t%REF\t%ALT\t%FILTER\t[%GT]\n',DV],capture_output=True,text=True).stdout
    m=[l.split('\t') for l in o.splitlines() if int(l.split('\t')[0])==p]
    call='none'
    for x in m:
        if x[1]==ref and alt in x[2].split(','): call=x[4].strip()+('' if x[3]=='PASS' else '_'+x[3]); break
    else:
        if m: call='other'
    cnt[(gt,call.replace('|','/'))]+=1
for k,v in sorted(cnt.items()): print(k,v)
