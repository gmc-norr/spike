import sys,subprocess
ev=[l.rstrip('\n').split('\t') for l in open('../transplant2/full/events.tsv')][1:]
sets={'forward':sys.argv[1],'reverse':sys.argv[2]}
def load(fn):
    d=set()
    for l in open(fn):
        c,p,r,a=l.split('\t')[:4]; d.add((c,p,r,a))
    return d
H={'forward':load('H1.calls.tsv'),'reverse':load('H2.calls.tsv')}
for s in ('forward','reverse'):
  for g in ('INS20-49','DEL20-49','DUP50-299','INS5-19'):
    miss=[e for e in ev if e[0]==s and e[1]==g and tuple(e[3:7]) not in H[s]]
    near=0; nearlen=0; none=0
    for e in miss:
        c,p,r,a=e[3:7]; p=int(p); L=len(a)-len(r)
        out=subprocess.run(['bcftools','query','-r',f'{c}:{max(1,p-50)}-{p+50}','-f','%POS\t%REF\t%ALT\t[%GT]\n',sets[s]],capture_output=True,text=True).stdout.strip().splitlines()
        alts=[o.split('\t') for o in out if o]
        alts=[o for o in alts if o[3] not in ('0/0','./.')]
        if not alts: none+=1; continue
        near+=1
        if any(abs((len(o[2])-len(o[1]))-L)<=max(2,abs(L)*0.2) and (len(o[2])-len(o[1]))*L>0 for o in alts): nearlen+=1
    print(s,g,'missed_exact',len(miss),'any_alt_call_within50',near,'similar_len_indel_within50',nearlen,'nothing',none)
