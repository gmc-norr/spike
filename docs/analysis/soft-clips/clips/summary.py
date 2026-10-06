"""One cause per soft-clipped read, from classify2.py, foreign.py and foreign_sw.py."""
import pickle, collections
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; total = D["total"]
cls = pickle.load(open("cls2.pkl", "rb")); F = dict(pickle.load(open("foreign_cls.pkl", "rb"))); SW = pickle.load(open("foreign_sw.pkl", "rb"))
group = {"badend": "bad end", "hiQerr": "bad end", "polyG": "bad end", "short": "bad end",
         "adapter": "adapter", "adapter2": "adapter",
         "site": "sample's own variant", "samplevar": "sample's own variant",
         "chimera": "chimeric fragment", "chimericpair": "chimeric fragment",
         "lowcomp": "simple repeat", "repeatmap": "simple repeat"}
rest = {"local": "bad end", "local-inverted": "chimeric fragment", "not local": "not in the reference"}
def cause(i, k):
    if k == "foreign":
        f = F[i]
        return rest[SW[i]] if f == "rest" else group[f]
    return group[k]
order = ["bad end", "adapter", "sample's own variant", "chimeric fragment", "simple repeat", "not in the reference"]
prio = {k: j for j, k in enumerate(order)}
per_read = {}
for i, (c, k) in enumerate(zip(C, cls)):
    g = cause(i, k); key = (c["name"], c["r1"])
    if key not in per_read or prio[g] > prio[per_read[key]]: per_read[key] = g   # a read with two kinds counts as the less ordinary one
cnt = collections.Counter(per_read.values()); n = len(per_read)
print(f"reads {total}; soft-clipped {n} ({100*n/total:.2f}%)")
print(f"{'cause':22s} {'reads':>6s} {'of clipped':>11s} {'of all reads':>13s}")
for k in order: print(f"{k:22s} {cnt[k]:6d} {100*cnt[k]/n:10.1f}% {100*cnt[k]/total:12.3f}%")
pickle.dump(per_read, open("cause_per_read.pkl", "wb"))
