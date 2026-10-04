"""The duplicates check (docs/superpowers/plans/2026-10-04-duplicates.md): after raredisease's own
duplicate marking, do spiked hom SNVs read the old allele from leftover source duplicates?

Usage: measure.py --spiked BAM --baseline BAM --source BAM --sim BAM --reference FASTA --calls VCF
                  --events VCF --replaced TXT --source-region chr:start-end --region chr:start-end
                  --spiked-metrics TXT --baseline-metrics TXT --c1 TXT
"""
import argparse
import sys

import pysam

COUNTED_MAPQ, COUNTED_BQ = 5, 10           # DeepVariant 1.6.1's make_examples defaults
UNCOUNTED = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
NOT_PRIMARY = 0x4 | 0x100 | 0x800
MARGIN, LEAK = 0.02, 0.5                   # the locked D rule


def counted(read):
    return not (read.flag & UNCOUNTED) and read.mapping_quality >= COUNTED_MAPQ


def bases_at(bam, chrom, pos1):
    """(name, base) of each counted read with a base of quality >= 10 at 1-based pos1."""
    pos0 = pos1 - 1
    out = []
    for r in bam.fetch(chrom, pos0, pos0 + 1):
        if not counted(r):
            continue
        for q, p in r.get_aligned_pairs(matches_only=True):
            if p == pos0:
                if r.query_qualities[q] >= COUNTED_BQ:
                    out.append((r.query_name, r.query_sequence[q].upper()))
                break
            if p > pos0:
                break
    return out


def allele_counts(bam, fasta, chrom, pos1):
    """(names of the reads with the reference base, all reads with a base) at 1-based pos1."""
    ref = fasta.fetch(chrom, pos1 - 1, pos1).upper()
    bases = bases_at(bam, chrom, pos1)
    return [n for n, b in bases if b == ref], len(bases)


def pooled(sites):
    """Sum of (ref names, total) over sites."""
    return sum(len(r) for r, _ in sites), sum(t for _, t in sites)


def share(num, den):
    return num / den if den else float("nan")


def verdict(r_spk, r_real, leak):
    if r_spk - r_real <= MARGIN:
        return "does not matter"
    return "matters" if leak >= LEAK else "something else"


def split_ref_reads(names, source_dups, replaced):
    out = {"source duplicate": 0, "not taken": 0, "replaced": 0, "SPIKE_": 0}
    for n in names:
        if n.startswith("SPIKE_"):
            out["SPIKE_"] += 1
        elif n in source_dups:
            out["source duplicate"] += 1
        elif n in replaced:
            out["replaced"] += 1
        else:
            out["not taken"] += 1
    return out


def coverage_ok(bam, chrom, pos1, lo=20, hi=45, good_share=0.95):
    """20-45 primary, mapped, non-duplicate, non-QC-fail records over pos1, >= 95% MAPQ >= 20 and proper."""
    n = good = 0
    for r in bam.fetch(chrom, pos1 - 1, pos1):
        if r.flag & UNCOUNTED:
            continue
        n += 1
        if r.mapping_quality >= 20 and r.is_proper_pair:
            good += 1
    return lo <= n <= hi and good >= good_share * n


def inside(bam, chrom, start1, end1):
    """Names of the counted reads lying wholly inside 1-based [start1, end1]."""
    return [r.query_name for r in bam.fetch(chrom, start1 - 1, end1)
            if counted(r) and r.reference_start >= start1 - 1 and r.reference_end <= end1]


def dup_counts(bam, chrom, pos1, flank=300):
    """(duplicate-flagged, all) primary mapped records overlapping 1-based [pos1 - flank, pos1 + flank]."""
    dups = total = 0
    for r in bam.fetch(chrom, pos1 - 1 - flank, pos1 + flank):
        if r.flag & NOT_PRIMARY:
            continue
        total += 1
        dups += bool(r.flag & 0x400)
    return dups, total


# --- the run ---------------------------------------------------------------------------

def region(text):
    chrom, span = text.split(":")
    start, end = (int(x.replace(",", "")) for x in span.split("-"))
    return chrom, start, end


def read_events(path):
    """[(kind, chrom, start, end)] from draw_events.py's VCF."""
    out = []
    with pysam.VariantFile(path) as vcf:
        for rec in vcf:
            vaf = float(rec.info["SIM_VAF"])
            if rec.alts[0] == "<DEL>":
                out.append(("del", rec.chrom, rec.pos + 1, rec.stop))
            else:
                out.append(("hom" if vaf == 1.0 else "het", rec.chrom, rec.pos, rec.pos))
    return out


def real_sites(calls, source, chrom, start, end, genotype):
    """raredisease's 1-bp-for-1-bp records with this sorted GT that share no position and pass the
    coverage test in the source BAM."""
    seen, recs = {}, []
    with pysam.VariantFile(calls) as vcf:
        sample = list(vcf.header.samples)[0]
        for rec in vcf.fetch(chrom, start - 1, end):
            seen[rec.pos] = seen.get(rec.pos, 0) + 1
            gt = rec.samples[sample]["GT"]
            if (len(rec.ref) == 1 and rec.alts and len(rec.alts) == 1 and len(rec.alts[0]) == 1
                    and None not in gt and tuple(sorted(gt)) == genotype):
                recs.append(rec.pos)
    return [p for p in recs if seen[p] == 1 and coverage_ok(source, chrom, p)], set(seen)


def names_with(bam, chrom, start, end, flag):
    return {r.query_name for r in bam.fetch(chrom, start - 1, end)
            if not (r.flag & NOT_PRIMARY) and r.flag & flag}


def percent_duplication(path):
    lines = open(path).read().splitlines()
    i = next(i for i, l in enumerate(lines) if l.startswith("LIBRARY"))
    head, row = lines[i].split("\t"), lines[i + 1].split("\t")
    return float(row[head.index("PERCENT_DUPLICATION")])


def main():
    ap = argparse.ArgumentParser()
    for a in ("spiked", "baseline", "source", "sim", "reference", "calls", "events", "replaced",
              "source-region", "region", "spiked-metrics", "baseline-metrics", "c1"):
        ap.add_argument(f"--{a}", required=True)
    a = ap.parse_args()
    out = lambda s="": print(s, flush=True)  # noqa: E731

    fasta = pysam.FastaFile(a.reference)
    spiked, baseline, source, sim = (pysam.AlignmentFile(p) for p in (a.spiked, a.baseline, a.source, a.sim))
    chrom, start, end = region(a.region)
    s_chrom, s_start, s_end = region(a.source_region)
    events = read_events(a.events)
    replaced = set(open(a.replaced).read().split())
    source_dups = names_with(source, s_chrom, s_start, s_end, 0x400)
    out(f"events: {sum(e[0] == 'het' for e in events)} het, {sum(e[0] == 'hom' for e in events)} hom, "
        f"{sum(e[0] == 'del' for e in events)} del; source duplicate names in {a.source_region}: {len(source_dups)}")

    hom = [(c, s) for k, c, s, _ in events if k == "hom"]
    het = [(c, s) for k, c, s, _ in events if k == "het"]
    real_hom, call_pos = real_sites(a.calls, source, chrom, start, end, (1, 1))
    real_het, _ = real_sites(a.calls, source, chrom, start, end, (0, 1))
    out(f"real sites passing the coverage test: {len(real_hom)} hom, {len(real_het)} het")

    spk = [allele_counts(spiked, fasta, c, p) for c, p in hom]
    real = [allele_counts(baseline, fasta, chrom, p) for p in real_hom]
    shifted = [allele_counts(baseline, fasta, chrom, p + 1) for p in real_hom if p + 1 not in call_pos]
    ref_spk, n_spk = pooled(spk)
    ref_real, n_real = pooled(real)
    ref_shift, n_shift = pooled(shifted)
    r_spk, r_real, r_shift = share(ref_spk, n_spk), share(ref_real, n_real), share(ref_shift, n_shift)
    ref_names = [n for names, _ in spk for n in names]
    split = split_ref_reads(ref_names, source_dups, replaced)
    leak = share(split["source duplicate"], len(ref_names)) if ref_names else 0.0

    out()
    out("per spiked hom SNV: position, old-allele reads / counted reads, of which source duplicates")
    for (c, p), (names, total) in zip(hom, spk):
        s = split_ref_reads(names, source_dups, replaced)
        out(f"  {c}:{p}  {len(names)}/{total}  {s['source duplicate']}")
    out()
    out(f"R_spk  = {ref_spk}/{n_spk} = {r_spk:.4f}  (20 spiked hom SNVs, spiked BAM)")
    out(f"R_real = {ref_real}/{n_real} = {r_real:.4f}  ({len(real_hom)} real hom SNVs, baseline BAM)")
    out(f"R_spk - R_real = {r_spk - r_real:.4f}  (locked margin {MARGIN})")
    out(f"L = {split['source duplicate']}/{len(ref_names)} = {leak:.4f}  (locked floor {LEAK})")
    out(f"old-allele reads at spiked hom SNVs, split: {split}")

    c1 = open(a.c1).read().strip()
    c1_pass = c1.endswith("C1: PASS")
    c2_pass = r_real <= 0.10 and r_shift >= 0.90
    c3_spk = sum(t >= 15 for _, t in spk)
    c3_real = sum(t >= 15 for _, t in real)
    c3_pass = c3_spk >= 15 and c3_real >= 100
    out()
    out(f"C1: {c1}")
    out(f"C2: R_real {r_real:.4f} <= 0.10 and ref share at +1 bp {ref_shift}/{n_shift} = {r_shift:.4f} "
        f"({len(shifted)} sites) >= 0.90: {'PASS' if c2_pass else 'FAIL'}")
    out(f"C3: spiked hom SNVs with >= 15 reads {c3_spk}/20 (>= 15), real hom SNVs {c3_real} (>= 100): "
        f"{'PASS' if c3_pass else 'FAIL'}")

    out()
    out("reported, not judged:")
    bam_route = []
    for c, p in hom:
        names, _ = allele_counts(source, fasta, c, p)
        kept = [n for n, _ in bases_at(source, c, p) if n not in replaced]
        kept_ref = [n for n in names if n not in replaced]
        sim_ref, sim_total = allele_counts(sim, fasta, c, p)
        bam_route.append((kept_ref + sim_ref, len(kept) + sim_total))
    r, n = pooled(bam_route)
    out(f"  BAM route (source minus replaced_reads.txt, plus sim.bam), hom SNVs: old allele {r}/{n} = {share(r, n):.4f}")
    het_spk = [allele_counts(spiked, fasta, c, p) for c, p in het]
    het_real = [allele_counts(baseline, fasta, chrom, p) for p in real_het]
    r, n = pooled(het_spk)
    out(f"  het SNVs, alt share: spiked {n - r}/{n} = {share(n - r, n):.4f}", )
    r, n = pooled(het_real)
    out(f"  het SNVs, alt share: real {n - r}/{n} = {share(n - r, n):.4f} ({len(real_het)} sites)")
    for tag, bam in (("spiked", spiked), ("baseline", baseline)):
        d = [dup_counts(bam, c, p) for c, p in hom + het]
        dd, tt = sum(x for x, _ in d), sum(y for _, y in d)
        out(f"  duplicate share within 300 bp of the 40 spiked SNVs, {tag} BAM: {dd}/{tt} = {share(dd, tt):.4f}")
    for k, c, s, e in events:
        if k != "del":
            continue
        in_spk = inside(spiked, c, s, e)
        in_base = inside(baseline, c, s, e)
        out(f"  hom DEL {c}:{s}-{e}: counted reads wholly inside: spiked {len(in_spk)} "
            f"(source duplicates {sum(n in source_dups for n in in_spk)}), baseline {len(in_base)}")
    out(f"  Picard PERCENT_DUPLICATION: spiked {percent_duplication(a.spiked_metrics):.6f}, "
        f"baseline {percent_duplication(a.baseline_metrics):.6f}")

    out()
    if not (c1_pass and c2_pass and c3_pass):
        out("D: inconclusive (a control failed)")
    else:
        out(f"D: {verdict(r_spk, r_real, leak)}")


if __name__ == "__main__":
    sys.exit(main())
