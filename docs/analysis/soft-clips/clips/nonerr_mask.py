import pickle, numpy as np
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]; cls = pickle.load(open("cls2.pkl", "rb"))
nonerr = {"adapter", "site", "foreign", "chimera"}
mask = {}
for c, k in zip(C, cls):
    if k not in nonerr: continue
    key = (c["name"], c["r1"]); m = mask.setdefault(key, np.zeros(151, bool))
    if c["three"]: m[151 - c["L"]:] = True
    else: m[:c["L"]] = True
pickle.dump(mask, open("nonerr_mask.pkl", "wb")); print("reads with non-error clipped bases:", len(mask))
