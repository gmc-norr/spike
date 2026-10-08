import subprocess,collections,statistics as st
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full'
R='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample'
V={'forward':R+'/D25-7403_Seq25-7598_30x_resample/raredisease_results/call_snv/genome/NA12878_snv.vcf.gz',
   'reverse':R+'/D24-14230_Seq25-7600_30x_resample/raredisease_results/call_snv/genome/GM24385_seracare_cancer_snv.vcf.gz'}
evid=collections.defaultdict(dict)
for l in open(S+'/evidence.tsv'):
    f=l.rstrip('\n').split('\t')
    if f[0]=='set': continue
    evid[(f[0],f[1],f[4],f[5],f[6],f[7])][f[2]]=f
ev=[l.rstrip('\n').split('\t') for l in open(S+'/events.tsv')][1:]
def fl(x):
    try: return float(x)
    except: return None
for s in ('forward','reverse'):
  for g in ('SNV','DEL20-49','INS5-19','INS20-49'):
    rows=[e for e in ev if e[0]==s and e[1]==g]
    with open('reg2.tsv','w') as r:
        for e in rows: r.write(f"{e[3]}\t{e[4]}\t{e[4]}\n")
    out=subprocess.run(['bcftools','query','-R','reg2.tsv','-f','%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n',V[s]],capture_output=True,text=True).stdout
    called=set()
    for l in out.splitlines():
        c,p,r_,a,gt=l.split('\t')
        if gt not in ('0/0','./.'): called.add((c,p,r_,a))
    tab=collections.defaultdict(list)
    for e in rows:
        k=(e[3],e[4],e[5],e[6]); d=evid.get((s,g)+k,{})
        real=d.get('real'); fake=d.get('fake_normal')
        if not real or not fake: continue
        tab['called' if k in called else 'uncalled'].append((fl(real[12]),fl(fake[12]),fl(real[13]),fl(fake[13]),int(real[8]),int(fake[8])))
    for kk,v in tab.items():
        def med(i): 
            xs=[x[i] for x in v if x[i] is not None]; return st.median(xs) if xs else float('nan')
        fake_ge2=sum(1 for x in v if x[5]>=2); real_ge2=sum(1 for x in v if x[4]>=2)
        print(s,g,kk,'n',len(v),'medA real %.3f fake %.3f | medE real %.3f fake %.3f | carriers>=2 real %d fake %d'%(med(0),med(1),med(2),med(3),real_ge2,fake_ge2))
