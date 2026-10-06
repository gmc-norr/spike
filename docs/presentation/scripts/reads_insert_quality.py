import csv, collections
seen = set()
tl = {"spike": [], "real": []}
qsum = {k: {1: [0]*151, 2: [0]*151} for k in tl}
qn = {k: {1: [0]*151, 2: [0]*151} for k in tl}
crash = {k: [0, 0] for k in tl}
for line in open("snv_reads.sam"):
    f = line.split("\t")
    name, flag, tlen, qual = f[0], int(f[1]), int(f[8]), f[10]
    key = (name, flag & 0xC0)
    if key in seen: continue
    seen.add(key)
    kind = "spike" if name.startswith("SPIKE_") else "real"
    mate = 1 if flag & 0x40 else 2
    q = [ord(c) - 33 for c in qual]
    if flag & 0x10: q = q[::-1]          # sequencing order
    for i, v in enumerate(q[:151]):
        qsum[kind][mate][i] += v; qn[kind][mate][i] += 1
    if len(q) >= 20:
        crash[kind][0] += sum(1 for v in q[-20:] if v < 15) >= 10
        crash[kind][1] += 1
    if mate == 1 and flag & 0x2 and 0 < abs(tlen) <= 1500:
        tl[kind].append(abs(tlen))
for k in tl:
    v = sorted(tl[k]); n = len(v)
    print(k, "pairs", n, "median", v[n//2], "mean", round(sum(v)/n, 1), "crash", crash[k][0], "/", crash[k][1], round(100*crash[k][0]/crash[k][1], 2), "%")
with open("insert.csv", "w") as f:
    w = csv.writer(f); w.writerow(["kind", "tlen"])
    for k in tl:
        for x in tl[k]: w.writerow([k, x])
with open("qual_cycle.csv", "w") as f:
    w = csv.writer(f); w.writerow(["cycle", "real_R1", "spike_R1", "real_R2", "spike_R2"])
    for i in range(151):
        w.writerow([i+1] + [round(qsum[k][m][i]/qn[k][m][i], 2) if qn[k][m][i] else "" for m in (1, 2) for k in ("real", "spike")][0:0] +
                   [round(qsum["real"][1][i]/qn["real"][1][i],2), round(qsum["spike"][1][i]/qn["spike"][1][i],2),
                    round(qsum["real"][2][i]/qn["real"][2][i],2), round(qsum["spike"][2][i]/qn["spike"][2][i],2)])
