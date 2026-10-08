import sys,subprocess
ev=[l.rstrip('\n').split('\t') for l in open('../transplant2/full/events.tsv')][1:]
sets={'forward':sys.argv[1],'reverse':sys.argv[2]}
def load(fn):
    return {tuple(l.split('\t')[:4]) for l in open(fn)}
H={'forward':load('H1.calls.tsv'),'reverse':load('H2.calls.tsv')}
for s in ('forward','reverse'):
  for g in ('INS20-49','DUP50-299'):
    miss=[e for e in ev if e[0]==s and e[1]==g and tuple(e[3:7]) not in H[s]]
    anyc=0; ins=0; lens=[]
    for e in miss:
        c,p,r,a=e[3:7]; p=int(p); L=len(a)-len(r); w=abs(L)+100
        out=subprocess.run(['bcftools','query','-r',f'{c}:{max(1,p-w)}-{p+w}','-f','%POS\t%REF\t%ALT\t[%GT]\n',sets[s]],capture_output=True,text=True).stdout.strip().splitlines()
        alts=[o.split('\t') for o in out if o]
        alts=[o for o in alts if o[3] not in ('0/0','./.')]
        if alts: anyc+=1
        ii=[len(o[2])-len(o[1]) for o in alts if len(o[2])-len(o[1])>=10]
        if ii: ins+=1; lens.append((L,ii))
    print(s,g,'missed',len(miss),'any_alt_within_len+100',anyc,'ins>=10bp_within',ins,lens[:8])
