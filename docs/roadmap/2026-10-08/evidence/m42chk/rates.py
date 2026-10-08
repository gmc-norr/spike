import sys, statistics as st
ev=[l.rstrip('\n').split('\t') for l in open('../transplant2/full/events.tsv')][1:]
def load(fn):
    d={}
    for l in open(fn):
        c,p,r,a,q,f,gt,gq,ad=l.rstrip('\n').split('\t')
        d[(c,p,r,a)]=(gt,gq,ad,q)
    return d
H={'H1':load('H1.calls.tsv'),'H2':load('H2.calls.tsv')}
def het(gt): return gt.replace('|','/') in ('0/1','1/0')
def vaf(ad):
    x=[int(v) for v in ad.split(',') if v!='.']
    return x[1]/sum(x) if len(x)>1 and sum(x)>0 else None
from collections import defaultdict
res=defaultdict(list)
for s,g,k,c,p,r,a in ev:
    res[(s,g)].append((c,p,r,a))
order=['SNV','DEL1-4','DEL5-19','DEL20-49','INS1-4','INS5-19','INS20-49','DUP50-299']
for s,donor in [('forward','H1'),('reverse','H2')]:
    for g in order:
        L=res[(s,g)]; n=len(L)
        calls=[H[donor].get(x) for x in L]
        nh=sum(1 for c in calls if c and het(c[0]))
        nany=sum(1 for c in calls if c and c[0] not in ('0/0','./.','0|0'))
        gqs=[int(c[1]) for c in calls if c and het(c[0]) and c[1]!='.']
        print(s,g,n,'het',nh,'anycall',nany,'medGQ',st.median(gqs) if gqs else None)
# shared
for g in order:
    L=res[('shared',g)]; n=len(L)
    disc=0; dgq=[]; dv=[]; both=0; h1n=h2n=0
    for x in L:
        a=H['H1'].get(x); b=H['H2'].get(x)
        ha=bool(a and het(a[0])); hb=bool(b and het(b[0]))
        h1n+=ha; h2n+=hb
        if ha!=hb: disc+=1
        if ha and hb:
            both+=1
            if a[1]!='.' and b[1]!='.': dgq.append(abs(int(a[1])-int(b[1])))
            va,vb=vaf(a[2]),vaf(b[2])
            if va is not None and vb is not None: dv.append(abs(va-vb))
    print('shared',g,n,'H1het',h1n,'H2het',h2n,'disc',disc,'medDGQ',st.median(dgq) if dgq else None,'medDVAF',round(st.median(dv),3) if dv else None)
