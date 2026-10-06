import csv
def load(p):
    return [(int(l.split()[1]), int(l.split()[2])) for l in open(p)]
o, h, m = load("depth_orig.txt"), load("depth_het.txt"), load("depth_hom.txt")
assert [x[0] for x in o] == [x[0] for x in h] == [x[0] for x in m]
B = 250
with open("depth_bins.csv", "w") as f:
    w = csv.writer(f); w.writerow(["start", "original", "het", "hom"])
    for i in range(0, len(o), B):
        s = o[i][0]
        mean = lambda a: sum(x[1] for x in a[i:i+B]) / len(a[i:i+B])
        w.writerow([s, round(mean(o),2), round(mean(h),2), round(mean(m),2)])
# summary inside vs outside
ins = lambda a: [d for p,d in a if 38900000 < p <= 38910000]
out = lambda a: [d for p,d in a if p <= 38898000 or p > 38912000]
avg = lambda x: sum(x)/len(x)
for name, a in [("original", o), ("het", h), ("hom", m)]:
    print(name, "inside", round(avg(ins(a)),1), "outside", round(avg(out(a)),1), "ratio", round(avg(ins(a))/avg(out(a)),3))
