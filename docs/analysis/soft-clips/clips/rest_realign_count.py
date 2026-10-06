import pysam, pickle, collections
D = pickle.load(open("clips.pkl", "rb")); C = D["clips"]
st = collections.Counter(); far = collections.Counter()
for r in pysam.AlignmentFile("rest_ge20.bam"):
    if r.is_secondary or r.is_supplementary: continue
    if r.is_unmapped: st["unmapped"] += 1; continue
    c = C[int(r.query_name)]; L = len(r.query_sequence); al = r.query_alignment_length; nm = r.get_tag("NM")
    ident = 1 - nm / max(1, al)
    if al < 0.6 * L: st["maps only partly"] += 1; continue
    near = r.reference_name == "chr20" and abs(r.reference_start - c["boundary"]) < 1000
    key = ("near the clip" if near else ("chr20, far" if r.reference_name == "chr20" else "other chromosome")) + (", MAPQ>=20" if r.mapping_quality >= 20 else ", MAPQ<20")
    st[key] += 1
    if not near: far[r.reference_name] += 1
n = sum(st.values())
print(f"relaxed bwa-mem2 (-k 13 -T 20 -B 2), {n} clips:")
for k, v in st.most_common(): print(f"  {k:32s} {v:4d} ({100*v/n:.0f}%)")
print("  far hits by chromosome:", far.most_common(6))
