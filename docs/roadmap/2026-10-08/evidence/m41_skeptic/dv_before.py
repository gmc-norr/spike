# Hospital DV 1.6.1 genome calls (before spiking) at the hom background indels of the 121 'fell' footprints
import pysam, collections, sys
exec(open('cr3_bg_fell.py').read().split("for gt,v in out.items():")[0])
V=pysam.VariantFile('/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/call_snv/genome/GM24385_seracare_cancer_snv.vcf.gz')
recs=[x for x in out['1/1'] if x[0]>=0.75 and x[1]<0.75]
cnt=collections.Counter()
for fb,fa,sb,sa,c,p,L,ev in recs:
    best=None
    for r in V.fetch(c,max(0,p-30),p+30):
        if any(len(a)!=len(r.ref) for a in (r.alts or ())) and any(abs(len(a)-len(r.ref))==L for a in r.alts):
            s=r.samples[0]; gt=s.get('GT'); flt=';'.join(r.filter.keys())
            best=('/'.join(map(str,sorted(x for x in gt if x is not None))),flt)
    cnt[best]+=1
print(len(recs),cnt.most_common())
