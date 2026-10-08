"""Read-only: in-repeat reads at RFC1 (chr4:39348424, AAAAG) in a BAM.
A read is in-repeat when >= 90% of its 5-mers are rotations of the motif or
its reverse complement. Reports how many such reads sit at the locus, their
MAPQ, and whether their mate is in-repeat too (an IRR pair, which spike's
reference-overlapping tiling cannot make)."""
import sys, pysam

MOTIF = "AAAAG"
def rots(m):
    return {m[i:] + m[:i] for i in range(len(m))}
RC = MOTIF[::-1].translate(str.maketrans("ACGT", "TGCA"))
KM = rots(MOTIF) | rots(RC)

def purity(seq):
    if seq is None or len(seq) < 50:
        return 0.0
    k = [seq[i:i + 5] for i in range(len(seq) - 4)]
    return sum(x in KM for x in k) / len(k)

bam = pysam.AlignmentFile(sys.argv[1])
irr = []
for r in bam.fetch("chr4", 39347900, 39349000):
    if r.is_secondary or r.is_supplementary or r.is_duplicate:
        continue
    if purity(r.query_sequence) >= 0.9:
        irr.append(r)
names = {r.query_name for r in irr}
mate_pure = 0
mate_elsewhere = 0
mq = [r.mapping_quality for r in irr]
for r in irr:
    if r.mate_is_unmapped:
        continue
    if r.next_reference_name != "chr4" or abs(r.next_reference_start - 39348424) > 2000:
        mate_elsewhere += 1
    # mate in-repeat? look it up if local
    if r.next_reference_name == "chr4" and abs(r.next_reference_start - 39348424) < 2000:
        try:
            m = bam.mate(r)
            if purity(m.query_sequence) >= 0.9:
                mate_pure += 1
        except ValueError:
            pass
print("in-repeat reads at locus", len(irr), "MAPQ0", sum(q == 0 for q in mq),
      "mate in-repeat (local)", mate_pure, "mate placed >2kb away/other chrom", mate_elsewhere)
