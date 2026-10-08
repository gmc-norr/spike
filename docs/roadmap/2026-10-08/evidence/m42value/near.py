import sys,subprocess,collections,statistics as st
S='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant2/full'
R='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample'
V={'forward':R+'/D25-7403_Seq25-7598_30x_resample/raredisease_results/call_snv/genome/NA12878_snv.vcf.gz',
   'reverse':R+'/D24-14230_Seq25-7600_30x_resample/raredisease_results/call_snv/genome/GM24385_seracare_cancer_snv.vcf.gz'}
groups=sys.argv[1].split(',')
ev=[l.rstrip('\n').split('\t') for l in open(S+'/events.tsv')][1:]
evid=collections.defaultdict(dict)
for l in open(S+'/evidence.tsv'):
    f=l.rstrip('\n').split('\t')
    if f[0]=='set': continue
    evid[(f[0],f[4],f[5],f[6],f[7])][f[2]]=f
for s in ('forward','reverse'):
  for g in groups:
    rows=[e for e in ev if e[0]==s and e[1]==g]
    reg=open('reg.tsv','w')
    for e in rows: reg.write(f"{e[3]}\t{max(1,int(e[4])-50)}\t{int(e[4])+50}\n")
    reg.close()
    out=subprocess.run(['bcftools','query','-R','reg.tsv','-f','%CHROM\t%POS\t%REF\t%ALT\t[%GT\t%GQ]\n',V[s]],capture_output=True,text=True).stdout
    calls=collections.defaultdict(list)
    for l in out.splitlines():
        c,p,r,a,gt,gq=l.split('\t'); calls[c].append((int(p),r,a,gt,gq))
    cnt=collections.Counter(); miss=[]
    for e in rows:
        c,p,r,a=e[3],int(e[4]),e[5],e[6]
        L=len(a)-len(r)
        cs=[x for x in calls[c] if abs(x[0]-p)<=50 and x[3] not in ('0/0','./.')]
        exact=[x for x in cs if x[0]==p and x[1]==r and x[2]==a]
        sametype=[x for x in cs if (len(x[2])-len(x[1]))*L>0 and (L==0 or abs((len(x[2])-len(x[1]))-L)<=max(2,abs(L)//3))]
        anyc=cs
        if exact: cnt['exact']+=1
        elif sametype: cnt['near_same_type']+=1
        elif anyc: cnt['other_call_only']+=1
        else:
            cnt['nothing']+=1
        if not exact:
            ev_=evid.get((s,c,str(p),r,a),{})
            real=ev_.get('real'); fake=ev_.get('fake_normal')
            miss.append((c,p,L, 'near' if sametype else ('other' if anyc else 'none'), real[12] if real else '?', fake[12] if fake else '?', real[11] if real else '?', fake[11] if fake else '?'))
    print(s,g,len(rows),dict(cnt))
    if len(sys.argv)>2:
        for m in miss: print('   miss',m)
