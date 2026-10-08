import sys, collections, bisect
T='/tmp/claude-1066/-home-parlar-ai-dev-spike/23a35c29-9600-4da2-a2c9-618425a1e0d3/scratchpad/transplant/full/events.tsv'
# Alu intervals with running max end for exact containment
alu=collections.defaultdict(list)
for l in open('alu.bed'):
    c,s,e=l.split('\t'); alu[c].append((int(s),int(e)))
idx={}
for c,v in alu.items():
    v.sort(); st=[a for a,b in v]; me=[]; m=-1
    for a,b in v: m=max(m,b); me.append(m)
    idx[c]=(st,v,me)
def inalu(c,p):
    if c not in idx: return False
    st,v,me=idx[c]; i=bisect.bisect_right(st,p)-1
    while i>=0 and me[i]>p:
        if v[i][0]<=p<v[i][1]: return True
        i-=1
    return False
def load(f):
    d=collections.defaultdict(list)
    for l in open(f):
        c,p,e,flt,src,svt,gt=l.rstrip('\n').split('\t')
        if svt!='DEL': continue
        d[c].append((int(p),int(e),flt,src,gt))
    return d
H={'HG001':load('h1_all.tsv'),'HG002':load('h2_all.tsv')}
def status(sample,c,st,en):
    L=en-st; anyc=man=False
    for p,e,flt,src,gt in H[sample][c]:
        if flt not in ('PASS',): continue
        ov=min(e,en)-max(p,st)
        if ov>0 and ov/max(e-p,L)>=0.5:
            anyc=True
            if 'manta' in src or src=='Intersection': man=True
    return anyc,man
ev=[l.rstrip('\n').split('\t') for l in open(T)][1:]
donor={'forward':'HG001','reverse':'HG002'}
rate=collections.defaultdict(collections.Counter)
strat=collections.defaultdict(collections.Counter)
for s,b,c,st,en in ev:
    st,en=int(st),int(en)
    k=inalu(c,st)+inalu(c,en)
    cls={0:'noAlu',1:'oneAlu',2:'AluAlu'}[k]
    if s in donor:
        a,m=status(donor[s],c,st,en)
        r=rate[(s,b)]; r['n']+=1; r['any']+=a; r['manta']+=m
        if en-st>=300:
            q=strat[cls]; q['n']+=1; q['any']+=a; q['manta']+=m
            q2=strat[(s,b,cls)]; q2['n']+=1; q2['manta']+=m
    else:
        a1,m1=status('HG001',c,st,en); a2,m2=status('HG002',c,st,en)
        r=rate[(s,b)]; r['n']+=1; r['disc_any']+=(a1!=a2); r['disc_manta']+=(m1!=m2)
        r['h1_any']+=a1; r['h2_any']+=a2; r['h1_manta']+=m1; r['h2_manta']+=m2
        q=strat[('shared',b,cls)]; q['n']+=1; q['disc_manta']+=(m1!=m2); q['h1_manta']+=m1; q['h2_manta']+=m2
for k in sorted(rate): print(k, dict(rate[k]))
print()
for k in sorted(strat, key=str): print(k, dict(strat[k]))
