import sys, gzip, subprocess
fa, chrom, pos1, refb, altb = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
out = subprocess.run(['samtools', 'faidx', fa, f'{chrom}:{pos1-12}-{pos1+12}'], capture_output=True, text=True, check=True).stdout.split('\n', 1)[1].replace('\n', '').upper()
assert out[12] == refb, out
ref_k = out; alt_k = out[:12] + altb + out[13:]
rc = lambda s: s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
counts = {}
for d in sys.argv[6:]:
    c = {'spike_ref': 0, 'spike_alt': 0, 'kept_ref': 0, 'kept_alt': 0}
    for fq in ('R1.fq.gz', 'R2.fq.gz'):
        with gzip.open(f'{d}/{fq}', 'rt') as f:
            while True:
                h = f.readline()
                if not h: break
                s = f.readline().strip(); f.readline(); f.readline()
                kind = 'spike' if h.startswith('@SPIKE_') else 'kept'
                if ref_k in s or rc(ref_k) in s: c[kind + '_ref'] += 1
                if alt_k in s or rc(alt_k) in s: c[kind + '_alt'] += 1
    print(d.split('/')[-2], c)
