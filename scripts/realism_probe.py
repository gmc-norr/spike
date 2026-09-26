#!/usr/bin/env python3
"""Realism probe: today's `spike validate` checks on variants HG002 really carries.

spike's checks hold a spike-in to fixed rules about the aligner's output. This
asks what those rules say about *real* variants, in HG002's own BAM, with GIAB's
own truth: if they fail real variants too, a FAIL on a spike-in says nothing
about spike.

Selects from the T2T-Q100 chr20 VCF, inside the stvar benchmark (>= 1 kb from its
edges), single-ALT, genotype het or hom-alt, and isolated -- no other variant whose
REF and ALT differ by >= 10 bp within 1 kb of the event:
  ins  pure insertions (1-base REF, ALT starting with it) of 20-49 bp
  del  pure deletions (1-base ALT, REF starting with it) of >= 500 bp

Writes <out>/real_ins.vcf and <out>/real_del.vcf as spike truth VCFs (SIM_VAF 0.5
for het, 1.0 for hom), plus <out>/events.tsv with each event's class, genotype and
whether it sits in a RepeatMasker Simple_repeat or Low_complexity interval.

Usage: realism_probe.py OUT_DIR BENCHMARK_BED HG002_VCF_GZ RMSK_BED REFERENCE_FAI
"""
import bisect
import gzip
import os
import sys

CHROM = "chr20"
EDGE = 1_000
GAP = 1_000
MIN_INDEL = 10


def intervals(path, keep=None):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if f[0] == CHROM and (keep is None or keep(f)):
            out.append((int(f[1]), int(f[2])))
    return sorted(out)


def inside(pos, ivs, starts, margin=0):
    i = bisect.bisect_right(starts, pos) - 1
    return i >= 0 and ivs[i][0] + margin <= pos < ivs[i][1] - margin


def main(out, bench_path, vcf_path, rmsk_path, fai_path):
    os.makedirs(out, exist_ok=True)
    bench = intervals(bench_path)
    bstarts = [s for s, _ in bench]
    rep = intervals(rmsk_path, lambda f: f[3].split("/")[1] in ("Simple_repeat", "Low_complexity"))
    # Repeat intervals overlap, so test membership by scanning the few before pos.
    rstarts = [s for s, _ in rep]

    def in_repeat(lo, hi):
        i = bisect.bisect_right(rstarts, hi)
        return any(e > lo for s, e in rep[max(0, i - 50):i])

    records, indels = [], []
    with gzip.open(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != CHROM:
                continue
            pos, ref, alts = int(f[1]), f[3], f[4].split(",")
            if any(a.startswith("<") or abs(len(a) - len(ref)) >= MIN_INDEL for a in alts):
                indels.append((pos, pos + len(ref)))
            if len(alts) != 1:
                continue
            gt = f[9].split(":")[0].replace("|", "/")
            zyg = {"0/1": "het", "1/0": "het", "1/1": "hom"}.get(gt)
            if zyg is None:
                continue
            alt = alts[0]
            if len(ref) == 1 and alt[0] == ref and 20 <= len(alt) - 1 <= 49:
                records.append(("ins", pos, ref, alt, zyg))
            elif len(alt) == 1 and ref[0] == alt and len(ref) - 1 >= 500:
                records.append(("del", pos, ref, alt, zyg))
    indels.sort()
    istarts = [s for s, _ in indels]

    def isolated(lo, hi):
        i = bisect.bisect_right(istarts, hi + GAP)
        near = [(s, e) for s, e in indels[max(0, i - 200):i] if e > lo - GAP]
        return len(near) == 1  # only the event itself

    length = next(int(l.split("\t")[1]) for l in open(fai_path) if l.split("\t")[0] == CHROM)
    header = (
        "##fileformat=VCFv4.3\n"
        f"##contig=<ID={CHROM},length={length}>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
        '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="length">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n'
        '##INFO=<ID=SIM_VAF,Number=1,Type=Float,Description="fraction">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    vcfs = {k: open(os.path.join(out, f"real_{k}.vcf"), "w") for k in ("ins", "del")}
    for v in vcfs.values():
        v.write(header)
    table = open(os.path.join(out, "events.tsv"), "w")
    table.write("class\tpos\tlen\tzyg\trepeat\n")
    for kind, pos, ref, alt, zyg in records:
        span_end = pos + len(ref)
        if not (inside(pos, bench, bstarts, EDGE) and inside(span_end, bench, bstarts, EDGE)):
            continue
        if not isolated(pos, span_end):
            continue
        vaf = "0.5" if zyg == "het" else "1.0"
        gt = "0/1" if zyg == "het" else "1/1"
        if kind == "ins":
            n = len(alt) - 1
            vcfs["ins"].write(
                f"{CHROM}\t{pos}\treal_ins_{pos}\t{ref}\t{alt}\t999\tPASS\t"
                f"SVTYPE=INS;SVLEN={n};SIM_VAF={vaf}\tGT\t{gt}\n"
            )
        else:
            n = len(ref) - 1
            vcfs["del"].write(
                f"{CHROM}\t{pos}\treal_del_{pos}\t{ref[0]}\t<DEL>\t999\tPASS\t"
                f"SVTYPE=DEL;END={pos + n};SVLEN=-{n};SIM_VAF={vaf}\tGT\t{gt}\n"
            )
        rep_flag = "rep" if in_repeat(pos - 50, span_end + 50) else "norep"
        table.write(f"{kind}\t{pos}\t{n}\t{zyg}\t{rep_flag}\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
