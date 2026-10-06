"""Controls for the read-to-read quality spread (qvar.py).
1. Shuffle: per cycle, permute real qualities across reads. Keeps every per-cycle
   distribution, destroys any read-level structure. The read-mean SD must collapse.
2. First-order chain: fit P(q_c | q_c-1 bin, cycle, mate) on the real reads (the same
   state spike's chain uses, without the base term) and sample new reads. If the
   one-step memory is the cause, its read-mean SD must land near spike's, not real's.
"""
import numpy as np
seen=set(); Q={k:{1:[],2:[]} for k in ("real","spike")}
for line in open("snv_reads.sam"):
    f=line.split("\t"); name,flag,qual=f[0],int(f[1]),f[10]
    key=(name,flag&0xC0)
    if key in seen: continue
    seen.add(key)
    kind="spike" if name.startswith("SPIKE_") else "real"
    q=np.frombuffer(qual.encode(),dtype=np.uint8).astype(np.int16)-33
    if flag&0x10: q=q[::-1]
    Q[kind][1 if flag&0x40 else 2].append(q)
for k in Q:
    for m in (1,2): Q[k][m]=np.array(Q[k][m])
rng=np.random.default_rng(7)
def rmsd(X): return X.mean(axis=1).std(ddof=1)
def all37(X): return 100*(X.min(axis=1)==37).mean()
def lt33(X): return 100*(X.mean(axis=1)<33).mean()
def show(label,X): print(f"{label:34s} read-mean SD {rmsd(X):.3f}  all-Q37 {all37(X):5.2f}%  mean<Q33 {lt33(X):5.2f}%")
real=np.vstack([Q["real"][1],Q["real"][2]]); spk=np.vstack([Q["spike"][1],Q["spike"][2]])
show("real reads",real); show("spike's reads",spk)
sh=np.vstack([np.column_stack([rng.permutation(Q["real"][m][:,c]) for c in range(151)]) for m in (1,2)])
show("control 1: real, shuffled per cycle",sh)
def pbin(q): return np.minimum(q//10,3)
vals=np.array([2,11,25,37])
sims=[]
for m in (1,2):
    X=Q["real"][m]; n=len(X)
    # cycle 0 marginal
    p0=np.array([(X[:,0]==v).mean() for v in vals])
    out=np.zeros((n,151),dtype=np.int16)
    out[:,0]=vals[rng.choice(4,n,p=p0)]
    for c in range(1,151):
        T=np.zeros((4,4))
        for b in range(4):
            sel=pbin(X[:,c-1])==b
            for j,v in enumerate(vals): T[b,j]=((X[sel,c]==v).sum()) if sel.any() else 0
        marg=np.array([(X[:,c]==v).mean() for v in vals])
        for b in range(4):
            T[b]=T[b]/T[b].sum() if T[b].sum()>=30 else marg
        prev=pbin(out[:,c-1]); u=rng.random(n); cdf=np.cumsum(T[prev],axis=1)
        out[:,c]=vals[(u[:,None]>cdf).sum(axis=1).clip(0,3)]
    sims.append(out)
show("control 2: first-order chain fit",np.vstack(sims))
