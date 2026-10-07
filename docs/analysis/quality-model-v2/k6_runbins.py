"""K6 by template run bin: crashed-read share of spike's remade reads against the donors in K6's windows.

A read is crashed when 10+ of its last 20 qualities (sequencing order) are below Q15, as in K6.
The run bin is the longest one-letter run in the reference under the read, in sequencing
order, binned as spike's run_bin (A/G/T 7-8 -> 1, 9-11 -> 2, 12+ -> 3; C 5-6 -> 2, 7+ -> 3).
'expected' gives spike's reads the donors' crash share in each run bin.
Usage: k6_runbins.py SIM_BAM ORIGINAL_BAM REFERENCE_FASTA
  SIM_BAM: spike's aligned reads (align.sh's sim.bam) for the K6 events; SPIKE_ reads are counted.
  ORIGINAL_BAM: the 31-value HG002 chr20 BAM spike ran on; every read in the windows is counted."""
import sys
import numpy as np
import pysam

WINDOWS = [("chr20", 38547500, 38552500), ("chr20", 38895000, 38915000)]
COMP = str.maketrans("ACGTN", "TGCAN")


def run_bin(base, n):
    if base == "C":
        return 3 if n >= 7 else 2 if n >= 5 else 0
    if base in "AGT":
        return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
    return 0


def template_bin(fa, r):
    n = r.query_length
    if r.is_reverse:
        end = r.reference_end
        t = fa.fetch(r.reference_name, max(0, end - n), end).upper().translate(COMP)[::-1]
    else:
        t = fa.fetch(r.reference_name, r.reference_start, r.reference_start + n).upper()
    best, base, run = 0, "N", 0
    for b in t:
        if b == base and b != "N":
            run += 1
        else:
            base, run = b, int(b != "N")
        best = max(best, run_bin(base, run))
    return best


def crashed(r):
    q = np.array(r.query_qualities)
    q = q[::-1] if r.is_reverse else q
    return int((q[-20:] < 15).sum() >= 10)


def tally(fa, reads):
    t = np.zeros((2, 4, 2), int)
    for r in reads:
        m, h = 0 if r.is_read1 else 1, template_bin(fa, r)
        t[m, h, 0] += 1
        t[m, h, 1] += crashed(r)
    return t


def in_windows(bam, keep):
    seen = set()
    for chrom, a, b in WINDOWS:
        for r in bam.fetch(chrom, a, b):
            if r.flag & 0xF0C or r.is_unmapped or not keep(r):
                continue
            key = (r.query_name, r.is_read1)
            if key not in seen:
                seen.add(key)
                yield r


def main() -> None:
    sim, orig, ref = sys.argv[1:4]
    fa = pysam.FastaFile(ref)
    spike = tally(fa, (r for r in pysam.AlignmentFile(sim).fetch(until_eof=True)
                       if r.query_name.startswith("SPIKE_") and not r.flag & 0xF0C and not r.is_unmapped))
    donors = tally(fa, in_windows(pysam.AlignmentFile(orig), lambda r: True))
    print("crashed share by run bin 0 / 1 / 2 / 3 (reads in brackets)")
    for name, t in (("donors in windows", donors), ("spike remade reads", spike)):
        for m in range(2):
            cells = "  ".join(f"{100 * t[m, h, 1] / t[m, h, 0]:5.1f}% ({t[m, h, 0]:4d})" if t[m, h, 0] else "    -        "
                              for h in range(4))
            print(f"{name:20s} R{m + 1}: {cells}")
    rate = donors[:, :, 1] / np.maximum(donors[:, :, 0], 1)
    n_spike = spike[:, :, 0].sum()
    for name, t in (("donors in windows", donors), ("spike remade reads", spike)):
        print(f"{name:20s} overall {100 * t[:, :, 1].sum() / t[:, :, 0].sum():5.2f}% of {t[:, :, 0].sum()}; "
              f"run bin 3 holds {100 * t[:, 3, 0].sum() / t[:, :, 0].sum():.1f}%")
    print(f"expected for spike's reads at the donors' per-bin rates: {100 * (spike[:, :, 0] * rate).sum() / n_spike:.2f}%")


if __name__ == "__main__":
    main()
