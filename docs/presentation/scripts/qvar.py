"""Spread of base quality, real vs spike's reads, from snv_reads.sam (same read set as Figure 11)."""
import numpy as np, csv, json
seen=set(); Q={k:{1:[],2:[]} for k in ("real","spike")}
for line in open("snv_reads.sam"):
    f=line.split("\t"); name,flag,qual=f[0],int(f[1]),f[10]
    key=(name,flag&0xC0)
    if key in seen: continue
    seen.add(key)
    kind="spike" if name.startswith("SPIKE_") else "real"
    mate=1 if flag&0x40 else 2
    q=np.frombuffer(qual.encode(),dtype=np.uint8).astype(np.int16)-33
    if flag&0x10: q=q[::-1]
    Q[kind][mate].append(q)
for k in Q:
    for m in (1,2): Q[k][m]=np.array(Q[k][m],dtype=float)
rng=np.random.default_rng(1)
rows=[]
for c in range(151):
    r={"cycle":c+1}
    for k in Q:
        for m in (1,2):
            col=Q[k][m][:,c]
            r[f"{k}_R{m}_sd"]=round(col.std(ddof=1),3)
            for v in (25,11,2):
                r[f"{k}_R{m}_q{v}"]=round((col==v).mean()*100,3)
    rows.append(r)
with open("qual_var.csv","w") as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def boot(fn,A,B=1000):
    n=len(A); s=[fn(A[rng.integers(0,n,n)]) for _ in range(B)]
    return np.percentile(s,[2.5,97.5])
res={}
for k in Q:
    A=np.vstack([Q[k][1],Q[k][2]])
    allb=A.ravel()
    sd_all=allb.std(ddof=1)
    rm=A.mean(axis=1)            # per-read mean
    wsd=A.std(axis=1,ddof=1)     # within-read SD
    pc_sd=lambda X: np.sqrt(np.mean(X.var(axis=0,ddof=1)))  # pooled per-cycle SD
    res[k]=dict(n_reads=len(A),
        sd_all_bases=round(sd_all,3), sd_all_ci=[round(x,3) for x in boot(lambda X:X.ravel().std(ddof=1),A,300)],
        pooled_cycle_sd=round(pc_sd(A),3),
        sd_read_mean=round(rm.std(ddof=1),3), sd_read_mean_ci=[round(x,3) for x in boot(lambda X:X.mean(axis=1).std(ddof=1),A)],
        read_mean_pct=[round(x,2) for x in np.percentile(rm,[1,5,25,50,75,95,99])],
        within_read_sd_mean=round(wsd.mean(),3), within_read_sd_ci=[round(x,3) for x in boot(lambda X:X.std(axis=1,ddof=1).mean(),A)],
        within_read_sd_pct=[round(x,2) for x in np.percentile(wsd,[50,90,99])],
        share_reads_all37=round((A.min(axis=1)==37).mean()*100,2),
        share_reads_mean_lt30=round((rm<30).mean()*100,3),
        share_reads_mean_lt33=round((rm<33).mean()*100,3),
        share_reads_mean_lt33_n=int((rm<33).sum()))
    # last 20 cycles sd
    for m in (1,2):
        res[k][f"R{m}_sd_c1_10"]=round(np.sqrt(Q[k][m][:,:10].var(axis=0,ddof=1).mean()),3)
        res[k][f"R{m}_sd_c142_151"]=round(np.sqrt(Q[k][m][:,141:].var(axis=0,ddof=1).mean()),3)
# variance decomposition: between-read share of total variance
for k in Q:
    A=np.vstack([Q[k][1],Q[k][2]])
    tot=A.var(); between=A.mean(axis=1).var()
    res[k]["between_read_share_of_var_pct"]=round(between/tot*100,2)
json.dump(res,open("qual_var.json","w"),indent=1)
print(json.dumps(res,indent=1))
