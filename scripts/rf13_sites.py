#!/usr/bin/env python3
"""RF13 plan: the kill test's sites, on fresh seeds (RF11's sites have been seen).

Prints `<set>\t<chrom>\t<pos>\t<seq or ->` lines:

  rand   6 random positions in the HG002 stvar benchmark, drawn as rf11_sites.py
         draws `k1` (isolated, >= 250 kb apart) but with RF13's own seed
  hg     12 real HG002 insertions of 50-300 bp: single ALT, het or hom-alt,
         pure (1-base REF, ALT starting with it), isolated (no other variant of
         >= 10 bp within 1 kb), >= 1 kb inside a benchmark interval. `seq` is the
         inserted sequence (ALT past its anchor); `pos` is the VCF POS, which is
         where `ins:chr20:POS:SEQ` puts it.
  null   200 random positions drawn as rf11_sites.py draws `null-rand`, own seed

Usage: rf13_sites.py BENCHMARK_BED HG002_VCF_GZ RMSK_BED > rf13_sites.tsv
"""
import bisect
import gzip
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rf11_sites as base  # noqa: E402

CHROM = "chr20"
SEEDS = {"rand": 20261001, "hg": 20261002, "null": 20261003}


def real_insertions(vcf_path, bench, excl):
    bench_starts = [s for s, _ in bench]
    starts = [s for s, _ in excl]
    out = []
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != CHROM or "," in f[4]:
                continue
            pos, ref, alt = int(f[1]), f[3], f[4]
            gt = f[9].split(":")[0].replace("|", "/")
            if gt not in ("0/1", "1/0", "1/1"):
                continue
            if not (len(ref) == 1 and alt[0] == ref and 50 <= len(alt) - 1 <= 300):
                continue
            if set(alt[1:].upper()) - set("ACGT"):
                continue
            b = bisect.bisect_right(bench_starts, pos) - 1
            if b < 0 or not (bench[b][0] + base.EDGE <= pos < bench[b][1] - base.EDGE):
                continue
            # Isolated: the only exclusion interval over pos is this record's own.
            i = bisect.bisect_right(starts, pos) - 1
            own = (pos - base.VARIANT_GAP, pos + 1 + base.VARIANT_GAP)
            if i >= 0 and excl[i] != own and pos < excl[i][1]:
                continue
            out.append((pos, alt[1:].upper()))
    return out


def main(bench_path, vcf_path, rmsk_path):
    bench = base.read_bed(bench_path, CHROM)
    spans = base.indel_spans(vcf_path, CHROM)
    excl = base.merged_exclusion(spans)
    starts = [s for s, _ in excl]
    base.SEEDS.update({"rand": SEEDS["rand"], "null": SEEDS["null"]})
    rand = base.draw("rand", bench, bench, 6, base.K1_GAP, excl, starts)
    null = base.draw("null", bench, bench, 200, base.NULL_GAP, excl, starts)
    hg = sorted(random.Random(SEEDS["hg"]).sample(real_insertions(vcf_path, bench, excl), 12))
    for pos in rand:
        print(f"rand\t{CHROM}\t{pos}\t-")
    for pos, seq in hg:
        print(f"hg\t{CHROM}\t{pos}\t{seq}")
    for pos in null:
        print(f"null\t{CHROM}\t{pos}\t-")


if __name__ == "__main__":
    main(*sys.argv[1:])
