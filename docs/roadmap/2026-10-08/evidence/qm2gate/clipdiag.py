"""Clip test, why (after the verdict): by the template's longest-run bin, per set: reads with an
indel in the CIGAR, crashed share, mismatches per 100 aligned bases. Usage: clipdiag.py REF NAME=BAM ..."""
import sys, collections, pysam, numpy as np
exec(open(sys.path[0] + "/clipscore.py").read().split("def score(path):")[0].split('"""', 2)[2].replace("ref_path, vcf_path, *sets = sys.argv[1:]", "ref_path, *sets = sys.argv[1:]").replace("vcf = pysam.VariantFile(vcf_path)", ""))
for item in sets:
    name, path = item.split("=", 1)
    acc = collections.defaultdict(lambda: np.zeros(5))
    for r in pysam.AlignmentFile(path):
        if r.flag & 0x904 or not r.is_proper_pair or r.mapping_quality < 20 or abs(r.template_length) < 151 or r.query_length != 151:
            continue
        cig = r.cigartuples
        lead = cig[0][1] if cig[0][0] == 4 else 0; trail = cig[-1][1] if cig[-1][0] == 4 else 0
        span = fa.fetch(r.reference_name, max(0, r.reference_start - lead), r.reference_end + trail).upper()
        h = read_hp(span.translate(comp)[::-1] if r.is_reverse else span)
        q = np.array(r.query_qualities); qs = q[::-1] if r.is_reverse else q
        ins = sum(l for op, l in cig if op == 1); dl = sum(l for op, l in cig if op == 2); m = sum(l for op, l in cig if op in (0, 7, 8))
        a = acc[h]; a[0] += 1; a[1] += (ins + dl) > 0; a[2] += (qs[-20:] < 15).sum() >= 10; a[3] += r.get_tag("NM") - ins - dl; a[4] += m
    print(f"{name:5s} " + "  ".join(f"bin {h}: indel {100*acc[h][1]/acc[h][0]:5.2f}% crashed {100*acc[h][2]/acc[h][0]:5.2f}% mm/100 {100*acc[h][3]/acc[h][4]:.3f}" for h in range(4)))
