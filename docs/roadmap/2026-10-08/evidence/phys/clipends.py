import sys, re, collections, pysam
for a in sys.argv[1:]:
    name, f = a.split('=')
    c = collections.Counter(); n = 0
    for r in pysam.AlignmentFile(f):
        if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 20: continue
        n += 1
        cig = r.cigartuples
        lc = cig[0][1] if cig[0][0] == 4 else 0
        rc = cig[-1][1] if cig[-1][0] == 4 else 0
        five, three = (rc, lc) if r.is_reverse else (lc, rc)
        for end, L in (('5p', five), ('3p', three)):
            if L == 0: continue
            c[(end, '1-4' if L <= 4 else ('5-20' if L <= 20 else '21+'))] += 1
        if any(o in (1, 2) for o, _ in cig): c['indel'] += 1
    print(name, n, ' '.join('%s %s %.3f%%' % (k[0], k[1], 100 * v / n) for k, v in sorted((k, v) for k, v in c.items() if k != 'indel')), 'reads with I/D %.3f%%' % (100 * c['indel'] / n))
