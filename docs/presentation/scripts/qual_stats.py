import re, json, collections
seen = set(); reads = []
for line in open("snv_reads.sam"):
    f = line.rstrip("\n").split("\t")
    name, flag = f[0], int(f[1])
    key = (name, flag & 0xC0)
    if key in seen: continue
    seen.add(key)
    md = next((t[5:] for t in f[11:] if t.startswith("MD:Z:")), None)
    reads.append(dict(kind="spike" if name.startswith("SPIKE_") else "real", flag=flag, pos=int(f[3]),
                      cigar=f[5], seq=f[9], qual=[ord(c)-33 for c in f[10]], md=md, mapq=int(f[4])))
print("reads", collections.Counter(r["kind"] for r in reads))
print("Q values", sorted(collections.Counter(q for r in reads for q in r["qual"]).items()))

def bin_(q): return 0 if q < 10 else 1 if q < 20 else 2 if q < 30 else 3
# 1. transitions in sequencing order, cycles 2..end
trans = {k: collections.Counter() for k in ("real", "spike")}
for r in reads:
    q = r["qual"][::-1] if r["flag"] & 0x10 else r["qual"]
    for a, b in zip(q, q[1:]):
        trans[r["kind"]][(bin_(a), b)] += 1
# 2. low-Q count in last 20 cycles
lowtail = {k: collections.Counter() for k in ("real", "spike")}
for r in reads:
    q = r["qual"][::-1] if r["flag"] & 0x10 else r["qual"]
    if len(q) >= 20:
        lowtail[r["kind"]][sum(1 for v in q[-20:] if v < 15)] += 1
# 3. mismatch rate by Q on aligned bases, excluding variant positions
def aligned(r):
    # yields (refpos, readidx, is_mismatch) using CIGAR + MD
    if r["md"] is None: return []
    md_tokens = re.findall(r"(\d+)|(\^[A-Z]+)|([A-Z])", r["md"])
    # expand MD into per-aligned-base match/mismatch list (excluding deletions)
    md_seq = []
    for num, dele, mm in md_tokens:
        if num: md_seq += [False]*int(num)
        elif mm: md_seq.append(True)
    out = []; rp = r["pos"]; qi = 0; mi = 0
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", r["cigar"]):
        n = int(n)
        if op in "M=X":
            for j in range(n):
                out.append((rp + j, qi + j, md_seq[mi + j] if mi + j < len(md_seq) else False))
            rp += n; qi += n; mi += n
        elif op in "IS": qi += n
        elif op in "DN": rp += n
    return out
pos_tot = collections.Counter(); pos_mm = collections.Counter(); per_read = []
for r in reads:
    if r["mapq"] < 20: per_read.append([]); continue
    al = aligned(r); per_read.append(al)
    for rp, qi, mm in al:
        pos_tot[rp] += 1; pos_mm[rp] += mm
variant = {p for p in pos_tot if pos_tot[p] >= 10 and pos_mm[p] / pos_tot[p] >= 0.1}
print("variant positions excluded", len(variant))
mm = {k: collections.Counter() for k in ("real", "spike")}; tot = {k: collections.Counter() for k in ("real", "spike")}
for r, al in zip(reads, per_read):
    for rp, qi, m in al:
        if rp in variant: continue
        q = r["qual"][qi]
        tot[r["kind"]][q] += 1; mm[r["kind"]][q] += m
out = {"trans": {k: {f"{a}>{b}": v for (a, b), v in c.items()} for k, c in trans.items()},
       "lowtail": {k: dict(c) for k, c in lowtail.items()},
       "mm": {k: {q: [mm[k][q], tot[k][q]] for q in tot[k]} for k in tot}}
json.dump(out, open("qual_stats.json", "w"), indent=1)
for k in ("real", "spike"):
    print(k, "mismatch by Q:", {q: (mm[k][q], tot[k][q], round(mm[k][q]/tot[k][q], 5)) for q in sorted(tot[k])})
    print(k, "lowtail:", sorted(lowtail[k].items()))
    rows = collections.defaultdict(dict)
    for (a, b), v in trans[k].items(): rows[a][b] = v
    for a in sorted(rows):
        s = sum(rows[a].values()); print(k, "prev bin", a, "n", s, {b: round(v/s, 4) for b, v in sorted(rows[a].items())})
