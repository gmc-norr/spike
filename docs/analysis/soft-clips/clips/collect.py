"""Every soft clip in the real chr20 slice (HG002 35x, bwa-mem2 2.2.1), with what is needed to
say why it is there: side (3'/5' in sequencing order), length, qualities, identity to the
reference where the bases would have aligned, adapter prefix, poly-G, SA tag, fragment length,
and how many other reads clip at the same place."""
import os
import pysam, pickle, collections, numpy as np
REF = os.environ["REF"]
fa = pysam.FastaFile(REF); chrom_seq = None
comp = str.maketrans("ACGTN", "TGCAN")
bam = pysam.AlignmentFile(os.environ["SLICE"])
lo = min(r.reference_start for r in bam.fetch("chr20", 38500000, 38501000)) - 2000
ref = fa.fetch("chr20", lo, 40210000).upper()
clips = []; total = 0
for r in pysam.AlignmentFile(os.environ["SLICE"]):
    if r.flag & 0xF0C: continue
    total += 1
    cig = r.cigartuples
    for side, (op, L) in (("left", cig[0]), ("right", cig[-1])):
        if op != 4: continue
        seq = r.query_sequence; q = r.query_qualities
        if side == "left":
            s, qq = seq[:L], list(q[:L]); rpos = list(range(r.reference_start - L, r.reference_start)); boundary = r.reference_start
        else:
            s, qq = seq[-L:], list(q[-L:]); rpos = list(range(r.reference_end, r.reference_end + L)); boundary = r.reference_end
        refs = "".join(ref[p - lo] if 0 <= p - lo < len(ref) else "N" for p in rpos)
        ident = sum(a == b for a, b in zip(s, refs)) / L
        three = (side == "right") != r.is_reverse          # 3' end in sequencing order
        sq = s.translate(comp)[::-1] if r.is_reverse else s  # sequencing orientation
        qs = qq[::-1] if r.is_reverse else qq
        if three: first = sq            # bases run 5'->3' away from the insert: adapter starts at the clip
        else: first = sq[::-1]
        clips.append(dict(name=r.query_name, r1=r.is_read1, rev=r.is_reverse, side=side, three=three, L=L,
                          seq_seqorient=sq, q_seqorient=qs, meanq=float(np.mean(qq)), lowfrac=float(np.mean(np.array(qq) < 15)),
                          ident=ident, boundary=boundary, tlen=abs(r.template_length), proper=r.is_proper_pair,
                          sa=r.has_tag("SA"), mapq=r.mapping_quality, clip_ref_orient=s))
pickle.dump(dict(clips=clips, total=total), open("clips.pkl", "wb"))
reads = {(c["name"], c["r1"]) for c in clips}
print(f"primary mapped reads {total}; reads with a soft clip {len(reads)} ({100*len(reads)/total:.2f}%); clips {len(clips)}")
with open("clips_ge20.fa", "w") as f:
    for i, c in enumerate(clips):
        if c["L"] >= 20: f.write(f">{i}\n{c['clip_ref_orient']}\n")
