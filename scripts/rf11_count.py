#!/usr/bin/env python3
"""RF11 plan: `ins_reads`'s read count, today's rule and the proposed one, from a BAM.

A copy of `check_ins_reads` in `src/validate.rs`, read through `samtools view`:

- the reads are those `reads_with_inserted_sequence` queries, overlapping
  `[pos - 100, pos + 100]`, that `usable_alignment` admits: not unmapped,
  secondary, supplementary, duplicate or QC-fail, and MAPQ >= 20 (`--min-mapq`'s
  default; 255, "unavailable", counts as 0 as noodles reads it);
- `cigar_shows_insertion_near` walks the CIGAR from the 0-based alignment start:
  an `I` of at least `min_len` counts where it sits within 100 of `pos`, and so
  does a soft clip, but only once `min_len` reaches the clip floor -- 50 today
  (`INS_MAX_EVIDENCE_LEN`), `F` under the proposed fix;
- `pos` is the truth record's VCF POS, compared with 0-based reference positions
  as the Rust does; the count is of distinct read names.

Per read name it keeps the longest `I` and the longest soft clip within the pad,
so every rule below is one comparison:

  count(L, floor) = #names with  max_I >= L  or  (L >= floor and max_S >= L)

Usage:
  rf11_count.py sites BAM SITES_TSV > null_counts.tsv   (null-* rows of rf11_sites.py)
  rf11_count.py one BAM CHROM POS MIN_LEN FLOOR          (prints one count)
"""
import re
import subprocess
import sys

PAD = 100
MIN_MAPQ = 20
LENGTHS = (20, 30, 40, 45, 49, 50)
UNUSABLE = 0x4 | 0x100 | 0x800 | 0x400 | 0x200
CIGAR_OP = re.compile(r"(\d+)([MIDNSHP=X])")


def per_read(bam, chrom, pos, pad=PAD, min_mapq=MIN_MAPQ):
    """{read name: (longest I, longest soft clip)} within `pad` of `pos`."""
    region = f"{chrom}:{max(0, pos - pad) + 1}-{pos + pad}"
    out = subprocess.run(
        ["samtools", "view", bam, region], capture_output=True, text=True, check=True
    ).stdout
    reads = {}
    for line in out.splitlines():
        f = line.split("\t", 6)
        flag, mapq = int(f[1]), int(f[4])
        if flag & UNUSABLE or f[5] == "*":
            continue
        if (0 if mapq == 255 else mapq) < min_mapq:
            continue
        ref_pos = int(f[3]) - 1
        best_i = best_s = 0
        for n, op in CIGAR_OP.findall(f[5]):
            n = int(n)
            if op == "I":
                if abs(ref_pos - pos) <= pad:
                    best_i = max(best_i, n)
            elif op == "S":
                if abs(ref_pos - pos) <= pad:
                    best_s = max(best_s, n)
            elif op in "MDN=X":
                ref_pos += n
        old = reads.get(f[0], (0, 0))
        reads[f[0]] = (max(old[0], best_i), max(old[1], best_s))
    return reads


def count(reads, min_len, floor):
    return sum(
        1
        for best_i, best_s in reads.values()
        if best_i >= min_len or (min_len >= floor and best_s >= min_len)
    )


def main(argv):
    if argv[0] == "one":
        bam, chrom, pos, min_len, floor = argv[1:6]
        print(count(per_read(bam, chrom, int(pos)), int(min_len), int(floor)))
        return
    bam, sites = argv[1], argv[2]
    header = ["set", "chrom", "pos"]
    for length in LENGTHS:
        header += [f"I{length}", f"IS{length}"]
    print("\t".join(header))
    with open(sites) as f:
        for line in f:
            name, chrom, pos = line.split()
            if not name.startswith("null"):
                continue
            reads = per_read(bam, chrom, int(pos))
            row = [name, chrom, pos]
            for length in LENGTHS:
                # I only (floor above any length), and I or clip (floor 0).
                row += [str(count(reads, length, 10**9)), str(count(reads, length, 0))]
            print("\t".join(row))


if __name__ == "__main__":
    main(sys.argv[1:])
