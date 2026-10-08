import pysam, random, collections
B='/home/parlar_ai/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam'
bam=pysam.AlignmentFile(B); random.seed(1)
one=both=dupcov=allcov=0
for _ in range(200):
    pos=random.randint(30_000_000, 40_000_000)
    seen=collections.Counter(); 
    for a in bam.fetch('chr20',pos,pos+1):
        if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
        if a.reference_start<=pos<a.reference_end:
            allcov+=1
            if a.is_duplicate: dupcov+=1; continue
            if a.mapping_quality<20 or not a.is_proper_pair: continue
            seen[a.query_name]+=1
    one+=len(seen); both+=sum(1 for v in seen.values() if v==2)
print('pairs covering site',one,'both mates cover',both,both/one,'dup share of covering reads',dupcov/allcov)
