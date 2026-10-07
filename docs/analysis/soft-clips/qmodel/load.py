"""Real reads (not SPIKE_) in the 25 windows of the SNV run, with qualities in
sequencing order and, per base, whether it mismatches the reference -- clipped
bases included, placed where they would have aligned (case file 2026-10-06)."""
import os
import pysam, numpy as np, collections, pickle
BAM = os.environ["SNV_BAM"]
REF = os.environ["REF"]
fa = pysam.FastaFile(REF)
wins = [(p - 2500, p + 2500) for p in (38550000 + i * 60000 for i in range(25))]
reads = []; seen = set()
mm_pos = collections.Counter(); cov_pos = collections.Counter()
bam = pysam.AlignmentFile(BAM)
for a, b in wins:
    ref = fa.fetch("chr20", a - 1000, b + 1000).upper(); off = a - 1000
    for r in bam.fetch("chr20", a, b):
        if r.flag & 0xF0C or r.query_name.startswith("SPIKE_"): continue
        key = (r.query_name, r.is_read1)
        if key in seen: continue
        seen.add(key)
        seq = r.query_sequence; q = np.array(r.query_qualities, dtype=np.int16)
        n = len(seq); refpos = [None] * n
        for qp, rp in r.get_aligned_pairs():
            if qp is not None and rp is not None: refpos[qp] = rp
        cig = r.cigartuples
        if cig[0][0] == 4:   # left clip: place before the first aligned base
            L = cig[0][1]
            for i in range(L): refpos[i] = r.reference_start - (L - i)
        if cig[-1][0] == 4:
            R = cig[-1][1]
            for i in range(R): refpos[n - R + i] = r.reference_end + i
        clipped = np.zeros(n, bool)
        if cig[0][0] == 4: clipped[:cig[0][1]] = True
        if cig[-1][0] == 4: clipped[n - cig[-1][1]:] = True
        mm = np.zeros(n, bool); ok = np.zeros(n, bool); pos = np.full(n, -1)
        for i, rp in enumerate(refpos):
            if rp is None or not (0 <= rp - off < len(ref)): continue
            ok[i] = True; pos[i] = rp; mm[i] = seq[i] != ref[rp - off] and seq[i] != "N"
            if not clipped[i]:
                cov_pos[rp] += 1; mm_pos[rp] += mm[i]
        bases = np.frombuffer(seq.encode(), np.uint8)
        if r.is_reverse:  # sequencing order
            q, mm, ok, pos, clipped = q[::-1], mm[::-1], ok[::-1], pos[::-1], clipped[::-1]
            bases = np.frombuffer(seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1].encode(), np.uint8)
        reads.append(dict(name=r.query_name, mate=1 if r.is_read1 else 2, q=q, mm=mm, ok=ok, pos=pos, clip=clipped, base=bases))
var = {p for p, c in cov_pos.items() if c >= 5 and mm_pos[p] / c >= 0.10}
for d in reads:
    d["ok"] &= ~np.isin(d["pos"], list(var))
pickle.dump(reads, open("reads.pkl", "wb"))
print("reads", len(reads), "variant positions masked", len(var), "lengths", collections.Counter(len(d["q"]) for d in reads).most_common(3))
