import sys, re, bisect, collections
calls = collections.defaultdict(list)
for line in open(sys.argv[2]):
    c, p, t, i = line.rstrip('\n').split('\t')
    calls[c].append((int(p), t, i))
for c in calls: calls[c].sort()
keys = {c: [x[0] for x in v] for c, v in calls.items()}
ALU5 = 'GGCCGGGCGCGGTGGCTCA'
def rc(s):
    return s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
stats = collections.defaultdict(lambda: collections.Counter())
for line in open(sys.argv[1]):
    c, p, l, seq, gt = line.rstrip('\n').split('\t')
    p = int(p); seq = seq.upper()
    alu = (ALU5 in seq) or (rc(ALU5) in seq)
    polya = bool(re.search('A{10,}', seq[-60:])) or bool(re.search('T{10,}', seq[:60]))
    cls = 'alu_head' if alu else ('polyA_only' if polya else 'other')
    g = gt.split(':')[0]
    zyg = 'hom' if g in ('1|1', '1/1') else 'het'
    hits = []
    if c in keys:
        lo = bisect.bisect_left(keys[c], p - 300); hi = bisect.bisect_right(keys[c], p + 300)
        hits = calls[c][lo:hi]
    anyc = len(hits) > 0
    insc = any(h[1] == 'INS' for h in hits)
    manta_ins = any(h[1] == 'INS' and 'Manta' in h[2] for h in hits)
    s = stats[cls]
    s['n'] += 1; s['any_sv_call'] += anyc; s['ins_call'] += insc; s['manta_ins'] += manta_ins
    s[zyg] += 1
for k, s in stats.items():
    print(k, dict(s))
