"""Does a read's own (called) bases show longer one-letter runs than its template because the read
crashed? Real held-out reads (real.bam): longest-run bin from the called bases against the bin from the
reference under the read, for crashed and other reads. Usage: hpsource.py REF real.bam"""
import sys, collections, pysam, numpy as np
exec(open(sys.path[0] + "/clipscore.py").read().split("def score(path):")[0].split('"""', 2)[2].replace("ref_path, vcf_path, *sets = sys.argv[1:]", "ref_path, path = sys.argv[1:]").replace("vcf = pysam.VariantFile(vcf_path)", ""))
tab = {False: collections.Counter(), True: collections.Counter()}
for r in pysam.AlignmentFile(path):
    if r.flag & 0x904 or not r.is_proper_pair or r.mapping_quality < 20 or abs(r.template_length) < 151 or r.query_length != 151:
        continue
    cig = r.cigartuples
    lead = cig[0][1] if cig[0][0] == 4 else 0; trail = cig[-1][1] if cig[-1][0] == 4 else 0
    span = fa.fetch(r.reference_name, max(0, r.reference_start - lead), r.reference_end + trail).upper()
    seq = r.query_sequence
    if r.is_reverse:
        span, seq = span.translate(comp)[::-1], seq.translate(comp)[::-1]
    q = np.array(r.query_qualities); qs = q[::-1] if r.is_reverse else q
    crashed = bool((qs[-20:] < 15).sum() >= 10)
    tab[crashed][(read_hp(span), read_hp(seq))] += 1
for crashed in (False, True):
    t = tab[crashed]; n = sum(t.values())
    higher = sum(v for (a, b), v in t.items() if b > a); lower = sum(v for (a, b), v in t.items() if b < a)
    print(f"{'crashed' if crashed else 'other  '} reads {n:7d}: called-base bin above the template's {100*higher/n:5.2f}%, below {100*lower/n:5.2f}%")
    for a in range(4):
        row = [t[(a, b)] for b in range(4)]; s = max(1, sum(row))
        print(f"    template bin {a}: called-base bin 0-3 " + " ".join(f"{100*x/s:5.1f}%" for x in row) + f"  ({sum(row)})")
