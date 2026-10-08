"""What happens in the bad tails that make real reads soft-clip? Step 1: collect per-base data.

Every primary, mapped, non-duplicate, MAPQ >= 20 read of the HG002 35x chr20 slice whose fragment
is at least a read long (no adapter), in SEQUENCING orientation (reverse reads flipped and
complemented). Each base gets: quality, called base, the reference base where it sits (soft-clipped
bases placed where they would have aligned), the reference bases one before and one after in
sequencing order (for phasing), whether it is soft-clipped, its placed reference position, and whether it is counted (placed on
an A/C/G/T reference base, called base not N, not inserted, not at a site of the sample's own
variant: >= 5 reads and >= 10% non-reference).

Groups: badclip = the read's clip cause is 'bad end' (clips/summary.py) and it has a 3' clip;
crash = >= 10 of the last 20 qualities below Q15; good reads are kept 1 in 10.
With a 4th argument "all", 1 in 8 of all those reads instead (an unbiased sample).
Usage: collect.py SLICE_BAM REF OUT.pkl [all]"""
import sys, pickle, zlib, numpy as np, pysam

bam_path, ref_path, out = sys.argv[1:4]
SAMPLE_ALL = len(sys.argv) > 4 and sys.argv[4] == "all"  # keep 1 in 8 of every read instead, whatever its group
CODE = np.full(256, 4, np.uint8)
for i, b in enumerate(b"ACGT"):
    CODE[b] = i
COMP = np.array([3, 2, 1, 0, 4], np.uint8)  # A<->T, C<->G, N
cause = pickle.load(open("clips/cause_per_read.pkl", "rb"))

bam = pysam.AlignmentFile(bam_path)
reads = [r for r in bam.fetch("chr20") if not r.flag & 0xF0C]
lo = min(r.reference_start for r in reads) - 500
hi = max(r.reference_end for r in reads) + 500
ref = CODE[np.frombuffer(pysam.FastaFile(ref_path).fetch("chr20", lo, hi).upper().encode(), np.uint8)]

# the sample's own variants: >= 5 reads and >= 10% of them non-reference
cov = np.array(bam.count_coverage("chr20", lo, hi, quality_threshold=0), np.int64)  # 4 x len
tot = cov.sum(0)
refc = np.where(ref < 4, cov[np.minimum(ref, 3), np.arange(len(ref))], 0)
masked = (tot >= 5) & ((tot - refc) >= 0.10 * tot)
print(f"reference {lo}-{hi}; masked variant positions {masked.sum()}")

recs = []
n_seen = 0
for r in reads:
    if r.mapping_quality < 20 or abs(r.template_length) < r.query_length or r.query_length != 151:
        continue
    n_seen += 1
    q = np.array(r.query_qualities, np.uint8)
    L = len(q)
    if r.is_reverse:
        q_seq = q[::-1]
    else:
        q_seq = q
    crashed = (q_seq[-20:] < 15).sum() >= 10
    key = (r.query_name, r.is_read1)
    cig = r.cigartuples
    clip3 = (cig[0][1] if cig[0][0] == 4 else 0) if r.is_reverse else (cig[-1][1] if cig[-1][0] == 4 else 0)
    clip5 = (cig[-1][1] if cig[-1][0] == 4 else 0) if r.is_reverse else (cig[0][1] if cig[0][0] == 4 else 0)
    badclip = cause.get(key) == "bad end" and clip3 > 0
    anyclip = key in cause
    good = not anyclip and not crashed
    if SAMPLE_ALL:
        if zlib.crc32(r.query_name.encode()) % 8:
            continue
    elif good and zlib.crc32(r.query_name.encode()) % 10:
        continue
    elif not (badclip or crashed or good):
        continue
    placed = np.full(L, -1, np.int64)
    clipped = np.zeros(L, bool)
    rp, qp = r.reference_start, 0
    for i, (op, l) in enumerate(cig):
        if op == 4:
            placed[qp:qp + l] = np.arange(rp - l, rp) if i == 0 else np.arange(rp, rp + l)
            clipped[qp:qp + l] = True
            qp += l
        elif op in (0, 7, 8):
            placed[qp:qp + l] = np.arange(rp, rp + l)
            qp += l
            rp += l
        elif op == 1:
            qp += l
        elif op in (2, 3):
            rp += l
    called = CODE[np.frombuffer(r.query_sequence.encode(), np.uint8)]
    idx = placed - lo
    ok = (placed >= 0) & (idx >= 1) & (idx < len(ref) - 1)
    safe = np.where(ok, idx, 1)
    rb = np.where(ok, ref[safe], 4)
    rprev, rnext = np.where(ok, ref[safe - 1], 4), np.where(ok, ref[safe + 1], 4)  # reference order
    counted = ok & (rb < 4) & (called < 4) & ~np.where(ok, masked[safe], True)
    if r.is_reverse:  # to sequencing order: flip and complement; the base "after" in sequencing is the one before on the reference
        called, rb = COMP[called[::-1]], COMP[rb[::-1]]
        rprev, rnext = COMP[rnext[::-1]], COMP[rprev[::-1]]
        clipped, counted = clipped[::-1], counted[::-1]
        placed = placed[::-1]
    recs.append(dict(name=r.query_name, r1=r.is_read1, q=q_seq.copy(), called=called, ref=rb, rprev=rprev, rnext=rnext,
                     clipped=clipped.copy(), counted=counted.copy(), pos=placed.copy(), rev=r.is_reverse, clip3=clip3, clip5=clip5, crashed=bool(crashed),
                     badclip=bool(badclip), anyclip=anyclip, good=bool(good), cause=cause.get(key)))
pickle.dump(recs, open(out, "wb"))
print(f"reads looked at {n_seen}; kept {len(recs)}: badclip {sum(x['badclip'] for x in recs)}, "
      f"crashed {sum(x['crashed'] for x in recs)}, crashed+badclip {sum(x['crashed'] and x['badclip'] for x in recs)}, "
      f"good (1 in 10) {sum(x['good'] for x in recs)}")
