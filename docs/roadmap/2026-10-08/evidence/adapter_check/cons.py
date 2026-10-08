import sys, collections, re
comp = str.maketrans("ACGTN","TGCAN")
cnt = {1: collections.defaultdict(collections.Counter), 2: collections.defaultdict(collections.Counter)}
n = {1:0, 2:0}
for line in sys.stdin:
    f = line.split('\t')
    flag = int(f[1]); tlen = abs(int(f[8])); cig = f[5]; seq = f[9]
    if not (0 < tlen < 150): continue
    mate = 1 if flag & 64 else 2
    rev = flag & 16
    if rev:
        seq = seq.translate(comp)[::-1]
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cig)
        lead = int(ops[-1][0]) if ops[-1][1] == 'S' else 0
    else:
        ops = re.findall(r'(\d+)([MIDNSHP=X])', cig)
        lead = int(ops[0][0]) if ops[0][1] == 'S' else 0
    if lead: continue
    n[mate] += 1
    for i, b in enumerate(seq[tlen:tlen+45]):
        cnt[mate][i][b] += 1
for m in (1, 2):
    s = ''; ag = []
    for i in range(45):
        c = cnt[m][i]
        if not c: break
        b, k = c.most_common(1)[0]; s += b; ag.append(k / sum(c.values()))
    print(m, n[m], s)
    print('  min agreement first 33:', round(min(ag[:33]), 3), 'first:', round(ag[0], 3))
