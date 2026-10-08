import sys, glob, collections, pysam
comp = str.maketrans('ACGTN', 'TGCAN')
def tally(files, want_spike):
    c = [collections.Counter() for _ in range(12)]; n = 0
    for f in files:
        for r in pysam.AlignmentFile(f):
            if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.mapping_quality < 20: continue
            if r.query_name.startswith('SPIKE_') != want_spike: continue
            s = r.query_sequence
            if r.is_reverse: s = s.translate(comp)[::-1]
            n += 1
            for i in range(12): c[i][s[i]] += 1
    return n, c
S = sys.argv[1]
for label, files, sp in [('SPIKE_ reads (spiked.bam)', glob.glob(S + '/*/spiked.bam'), True), ('sample reads (base.bam)', glob.glob(S + '/*/base.bam'), False)]:
    n, c = tally(files, sp)
    print(label, 'n=%d' % n)
    for i in range(0, 11, 1):
        t = sum(c[i][b] for b in 'ACGT')
        print('  cycle %2d ' % (i + 1) + ' '.join('%s %.3f' % (b, c[i][b] / t) for b in 'ACGT'))
