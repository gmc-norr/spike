"""Read-only CR3 probe on an EXISTING transplant run (round 2, forward, run0).

For planted events, find HG002 Q100 hom-alt non-SNV records 300-2000 bp away,
and compare the indel's allele fraction in the recipient BAM against the
spiked view (recipient minus replaced_reads plus sim.bam).
A read carries the indel if it has an I/D op starting within 3 bp of POS+1
(left-normalised Q100 records; bwa may shift in repeats, so +-3).
"""
import sys, subprocess, pysam

RUN = sys.argv[1]
RECIP = sys.argv[2]
Q100 = sys.argv[3]
MAXEV = int(sys.argv[4])

replaced = set()
with open(f"{RUN}/replaced_reads.txt") as fh:
    for line in fh:
        replaced.add(line.strip())

events = []
with open(f"{RUN}/truth.vcf") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        events.append((f[0], int(f[1]), len(f[3]), len(f[4])))
events = events[:MAXEV]

vcf = pysam.VariantFile(Q100)
recip = pysam.AlignmentFile(RECIP)
sim = pysam.AlignmentFile(f"{RUN}/sim.bam")


def carries(read, pos0, is_del, ln):
    if read.cigartuples is None:
        return None
    rpos = read.reference_start
    covered = False
    hit = False
    for op, n in read.cigartuples:
        if op in (0, 7, 8):
            if rpos <= pos0 < rpos + n:
                covered = True
            rpos += n
        elif op == 2:
            if abs(rpos - (pos0 + 1)) <= 3 and is_del:
                hit = True
            if rpos <= pos0 < rpos + n:
                covered = True
            rpos += n
        elif op == 1:
            if abs(rpos - (pos0 + 1)) <= 3 and not is_del:
                hit = True
        elif op == 3:
            rpos += n
    # require the read to reach 10 bp past the site on both sides
    if read.reference_start > pos0 - 10 or read.reference_end < pos0 + ln + 10:
        return None
    return hit


def count(bam, chrom, pos0, is_del, ln, skip=None, only_hap=False):
    n = k = 0
    for r in bam.fetch(chrom, pos0, pos0 + 1):
        if r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_unmapped:
            continue
        if skip is not None and r.query_name in skip:
            continue
        if only_hap and "_hap_" not in r.query_name:
            continue
        c = carries(r, pos0, is_del, ln)
        if c is None:
            continue
        n += 1
        k += int(c)
    return n, k

print("event\tbg_indel\tdist\trecip_n\trecip_af\tspiked_n\tspiked_af\thap_n\thap_carry")
for chrom, epos, rl, al in events:
    for rec in vcf.fetch(chrom, max(0, epos - 2000), epos + 2000):
        gt = rec.samples[0]["GT"]
        if gt is None or len(gt) != 2 or gt != (1, 1):
            continue
        r, a = rec.ref, rec.alts[0]
        if len(r) == len(a):
            continue
        d = abs(rec.pos - epos)
        if d < 300:
            continue
        is_del = len(r) > len(a)
        ln = abs(len(r) - len(a))
        pos0 = rec.pos - 1
        rn, rk = count(recip, chrom, pos0, is_del, ln)
        kn, kk = count(recip, chrom, pos0, is_del, ln, skip=replaced)
        sn, sk = count(sim, chrom, pos0, is_del, ln)
        hn, hk = count(sim, chrom, pos0, is_del, ln, only_hap=True)
        if rn < 10:
            continue
        tn, tk = kn + sn, kk + sk
        print(f"{chrom}:{epos}\t{rec.pos}{r[:6]}>{a[:6]}\t{d}\t{rn}\t{rk/rn:.2f}\t{tn}\t{(tk/tn if tn else float('nan')):.2f}\t{hn}\t{hk}")
