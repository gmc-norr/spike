"""Draw the duplicates check's events (docs/superpowers/plans/2026-10-04-duplicates.md) and write
them as a VCF spike reads: 20 het SNVs, 20 hom SNVs, 3 hom 300 bp deletions.

Usage: draw_events.py --bam SOURCE_BAM --reference FASTA --calls VCF --region chr:start-end --seed 7 --out VCF
"""
import argparse
import bisect
import random
import sys

import pysam

import measure

TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
SPACING, CALL_FLANK, N_FLANK, DEL_LEN, DEL_STEP = 20_000, 1_000, 1_000, 300, 50


def draw(bam, fasta, calls, chrom, start, end, seed, n_het=20, n_hom=20, n_del=3, max_candidates=20_000):
    """Events in 1-based [start, end]: dicts with kind, chrom, start, end (1-based, inclusive), and
    ref/alt for SNVs. `calls` are the 1-based positions raredisease's VCF holds."""
    rng = random.Random(seed)
    calls = sorted(calls)
    length = fasta.get_reference_length(chrom)
    events, tried = [], 0

    def spaced(lo, hi):
        return all(max(lo, e["start"]) - min(hi, e["end"]) >= SPACING for e in events)

    def no_call(lo, hi):
        i = bisect.bisect_left(calls, lo - CALL_FLANK)
        return i == len(calls) or calls[i] > hi + CALL_FLANK

    def clean(lo, hi):
        a, b = max(0, lo - 1 - N_FLANK), min(length, hi + N_FLANK)
        return "N" not in fasta.fetch(chrom, a, b).upper()

    def base_ok(p):
        return fasta.fetch(chrom, p - 1, p).upper() in TRANSITION and measure.coverage_ok(bam, chrom, p)

    def candidate(width):
        nonlocal tried
        while tried < max_candidates:
            tried += 1
            lo = rng.randint(start, end - width + 1)
            hi = lo + width - 1
            if not (spaced(lo, hi) and no_call(lo, hi) and clean(lo, hi)):
                continue
            if all(base_ok(p) for p in range(lo, hi + 1, DEL_STEP if width > 1 else 1)):
                return lo, hi
        sys.exit(f"draw_events: {len(events)} events accepted from {max_candidates} candidates, "
                 f"{n_het + n_hom + n_del} wanted")

    for kind, n, width in (("het", n_het, 1), ("hom", n_hom, 1), ("del", n_del, DEL_LEN)):
        for _ in range(n):
            lo, hi = candidate(width)
            e = {"kind": kind, "chrom": chrom, "start": lo, "end": hi}
            if width == 1:
                e["ref"] = fasta.fetch(chrom, lo - 1, lo).upper()
                e["alt"] = TRANSITION[e["ref"]]
            events.append(e)
    return events


def write_vcf(events, fasta, path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##INFO=<ID=SIM_VAF,Number=1,Type=Float,Description="Allele fraction to simulate">\n')
        f.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
        f.write('##ALT=<ID=DEL,Description="Deletion">\n')
        for chrom in sorted({e["chrom"] for e in events}):
            f.write(f"##contig=<ID={chrom},length={fasta.get_reference_length(chrom)}>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for e in sorted(events, key=lambda e: (e["chrom"], e["start"])):
            if e["kind"] == "del":
                anchor = e["start"] - 1
                ref = fasta.fetch(e["chrom"], anchor - 1, anchor).upper()
                info = f"SVTYPE=DEL;END={e['end']};SIM_VAF=1.0"
                f.write(f"{e['chrom']}\t{anchor}\t.\t{ref}\t<DEL>\t.\tPASS\t{info}\n")
            else:
                vaf = "0.5" if e["kind"] == "het" else "1.0"
                f.write(f"{e['chrom']}\t{e['start']}\t.\t{e['ref']}\t{e['alt']}\t.\tPASS\tSIM_VAF={vaf}\n")


def call_positions(path, chrom, start, end):
    """Every 1-based reference base a record of the calls VCF covers near the region."""
    out = set()
    with pysam.VariantFile(path) as vcf:
        for rec in vcf.fetch(chrom, max(0, start - 1 - CALL_FLANK), end + CALL_FLANK):
            out.update(range(rec.pos, rec.pos + len(rec.ref)))
    return out


def main():
    ap = argparse.ArgumentParser()
    for a in ("bam", "reference", "calls", "region", "seed", "out"):
        ap.add_argument(f"--{a}", required=True)
    a = ap.parse_args()
    chrom, start, end = measure.region(a.region)
    fasta, bam = pysam.FastaFile(a.reference), pysam.AlignmentFile(a.bam)
    events = draw(bam, fasta, call_positions(a.calls, chrom, start, end), chrom, start, end, int(a.seed))
    write_vcf(events, fasta, a.out)
    for e in events:
        print(f"{e['kind']}\t{e['chrom']}:{e['start']}-{e['end']}\t{e.get('ref', '')}>{e.get('alt', '')}", flush=True)


if __name__ == "__main__":
    main()
