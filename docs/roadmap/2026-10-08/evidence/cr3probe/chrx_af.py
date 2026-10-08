"""Read-only: allele fraction and depth of HG002's own hemizygous chrX SNVs
(Q100 GT '1', non-PAR) against autosomal het SNVs, on the hospital BAM."""
import sys, pysam

bam = pysam.AlignmentFile(sys.argv[1])
vcf = pysam.VariantFile(sys.argv[2])


def af_at(chrom, pos0, alt):
    n = k = 0
    for col in bam.pileup(chrom, pos0, pos0 + 1, truncate=True, min_mapping_quality=20,
                          min_base_quality=13, stepper="nofilter"):
        for pr in col.pileups:
            r = pr.alignment
            if pr.is_del or pr.is_refskip or r.is_duplicate or r.is_secondary or r.is_supplementary:
                continue
            n += 1
            k += r.query_sequence[pr.query_position] == alt
    return n, k


def run(region, want_gt, limit):
    out = []
    for rec in vcf.fetch(*region):
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        if rec.samples[0]["GT"] != want_gt:
            continue
        n, k = af_at(rec.chrom, rec.pos - 1, rec.alts[0])
        if n >= 5:
            out.append((n, k / n))
        if len(out) >= limit:
            break
    return out

for label, region, gt in [("chrX nonPAR GT=1", ("chrX", 20000000, 30000000), (1,)),
                          ("chr20 het", ("chr20", 20000000, 30000000), (0, 1))]:
    res = run(region, gt, 200)
    afs = sorted(a for _, a in res)
    deps = sorted(n for n, _ in res)
    m = len(res)
    print(label, "n", m, "median depth", deps[m // 2], "median AF", round(afs[m // 2], 3),
          "AF p10", round(afs[m // 10], 3), "share AF<0.8", round(sum(a < 0.8 for a in afs) / m, 3))
