"""Are same-strand recurrent error sites made of distinct molecules (distinct 5' ends), or of
unflagged duplicates? Usage: python3 dupcheck.py BAM REF VCF BED regions..."""
import sys, collections, numpy as np, pysam
bam_f, ref_f, vcf_f, bed_f = sys.argv[1:5]
fa = pysam.FastaFile(ref_f); vcf = pysam.VariantFile(vcf_f); bam = pysam.AlignmentFile(bam_f)
bed = collections.defaultdict(list)
for line in open(bed_f):
    c, s, e = line.split()[:3]; bed[c].append((int(s), int(e)))
res = collections.Counter(); dist = collections.Counter()
for reg in sys.argv[5:]:
    c, se = reg.split(':'); s, e = map(int, se.split('-'))
    ref = fa.fetch(c, s, e).upper(); n = len(ref)
    mask = np.zeros(n, bool); inb = np.zeros(n, bool)
    for rec in vcf.fetch(c, s, e):
        a = rec.pos - 1 - s; b = rec.pos - 1 + len(rec.ref) - s; mask[max(a - 5, 0):max(b + 5, 0)] = True
    for bs, be in bed.get(c, ()):
        if be > s and bs < e: inb[max(bs - s, 0):min(be - s, n)] = True
    sup = collections.defaultdict(list)
    for r in bam.fetch(c, s, e):
        if r.flag & 0xF04 or r.mapping_quality < 20 or not r.is_proper_pair: continue
        if any(o not in (0, 4) for o, _ in r.cigartuples): continue
        five = r.reference_end if r.is_reverse else r.reference_start
        for q, p in r.get_aligned_pairs(matches_only=True):
            i = p - s
            if i < 0 or i >= n or mask[i] or not inb[i]: continue
            b = r.query_sequence[q]
            if b != ref[i] and b != 'N' and r.query_qualities[q] >= 10:
                sup[(p, r.is_reverse, b)].append((five, r.next_reference_start))
    for k, v in sup.items():
        if len(v) >= 3:
            res['sites>=3'] += 1
            fives = collections.Counter(x[0] for x in v)
            if max(fives.values()) == 1: res['all distinct 5p ends'] += 1
            dist[len(set(x[0] for x in v)) / len(v) >= 0.99] += 1
print(bam_f.split('/')[-1], dict(res))
