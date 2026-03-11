#!/usr/bin/env python3
"""Quantify realism metrics for spiked tandem duplications.

Inputs:
- A spiked BAM file.
- A truth VCF containing DUP events (Spike-style with SVTYPE=DUP and SIM_VAF).

Outputs one row per DUP in TSV:
- Event interval and expected VAF
- Coverage metrics (dup region, flanks, ratio)
- Split-read / soft-clip indicators near the junction window
- Optional BAF checks from user-specified SNP sites

The truth VCF from spike writes DUP start/end in 0-based coordinates in INFO/END.
Set --truth-1based if your VCF uses standard 1-based coordinates.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from dataclasses import dataclass
from typing import List, Tuple

try:
    import pysam
except Exception as exc:  # pragma: no cover - dependency guard
    print("ERROR: pysam is required. Install with `pip install pysam`.", file=sys.stderr)
    raise SystemExit(1) from exc


@dataclass
class DupEvent:
    idx: int
    chrom: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    vaf: float


def parse_truth_vcfs(path: str, use_zero_based: bool) -> List[DupEvent]:
    events: List[DupEvent] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom = f[0]
            try:
                pos = int(f[1])
            except ValueError:
                continue

            info = f[7]
            svtype = None
            end = None
            vaf = 0.5
            for entry in info.split(";"):
                if not entry:
                    continue
                if entry.startswith("SVTYPE="):
                    svtype = entry.split("=", 1)[1]
                elif entry.startswith("END="):
                    end = int(entry.split("=", 1)[1])
                elif entry.startswith("SIM_VAF="):
                    try:
                        vaf = float(entry.split("=", 1)[1])
                    except ValueError:
                        pass

            if svtype != "DUP":
                continue
            if end is None:
                continue

            start = pos if use_zero_based else max(pos - 1, 0)
            events.append(DupEvent(len(events), chrom, start, end, vaf))
    return events


def parse_snp_list(path: str) -> List[Tuple[str, int, str, str]]:
    sites = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt = parts[:4]
            try:
                p = int(pos)
            except ValueError:
                continue
            # accept either 0-based or 1-based from caller by auto-normalizing
            # if looks like common 1-based coordinate, keep as-is and adjust later if needed
            sites.append((chrom, p, ref.upper(), alt.upper()))
    return sites


def mean_depth(bam: pysam.AlignmentFile, chrom: str, start: int, end: int, min_mapq: int, min_bq: int) -> float:
    if end <= start:
        return 0.0
    total = 0
    for col in bam.pileup(
        chrom,
        start,
        end,
        min_base_quality=min_bq,
        min_mapping_quality=min_mapq,
        truncate=True,
        max_depth=100000,
    ):
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            total += 1
    return total / float(end - start)


def median(vals: List[float]) -> str:
    if not vals:
        return "NA"
    return f"{statistics.median(vals):.4f}"


def count_split_softreads(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    min_mapq: int,
) -> Tuple[int, int, int]:
    if end <= start:
        return 0, 0, 0

    sa_reads = set()
    soft_reads = set()
    total_reads = 0
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.mapping_quality < min_mapq:
            continue
        total_reads += 1
        name = read.query_name
        if read.has_tag("SA"):
            sa_reads.add(name)

        cigartuples = read.cigartuples or ()
        if cigartuples and ((cigartuples[0][0] == 4) or (cigartuples[-1][0] == 4)):
            soft_reads.add(name)

    return len(sa_reads), len(soft_reads), total_reads


def ref_alt_counts_at_pos(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos0: int,
    ref: str,
    alt: str,
    min_mapq: int,
    min_bq: int,
) -> Tuple[int, int, int]:
    # pos0 is 0-based genomic coordinate
    refc = 0
    altc = 0
    depth = 0

    for col in bam.pileup(
        chrom,
        pos0,
        pos0 + 1,
        min_base_quality=min_bq,
        min_mapping_quality=min_mapq,
        truncate=True,
        max_depth=100000,
    ):
        if col.reference_pos != pos0:
            continue
        for pr in col.pileups:
            if pr.is_del or pr.is_refskip or pr.query_position is None:
                continue
            read = pr.alignment
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if read.mapping_quality < min_mapq:
                continue

            b = read.query_sequence[pr.query_position].upper()
            if b not in ("A", "C", "G", "T"):
                continue
            q = read.query_qualities[pr.query_position]
            if q < min_bq:
                continue
            depth += 1
            if b == ref:
                refc += 1
            elif b == alt:
                altc += 1
    return refc, altc, depth


def baf_from_sites(
    bam: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
    sites: List[Tuple[str, int, str, str]],
    min_mapq: int,
    min_bq: int,
    min_depth: int,
    sites_are_one_based: bool,
) -> Tuple[int, float, float, float, int]:
    """Return (used, mean_af, median_af, max_abs_dev_from_half, tested).

    max_abs_dev_from_half is the largest |AF - 0.5| across usable SNP sites.
    For a het DUP this shifts toward ~0.33 or ~0.67, giving a value near 0.17.
    For a diploid region it stays near 0.0.
    """
    af_values: List[float] = []
    tested = 0
    used = 0

    for s_chrom, s_pos, ref, alt in sites:
        if s_chrom != chrom:
            continue
        pos0 = (s_pos - 1) if sites_are_one_based else s_pos
        if pos0 < start or pos0 >= end:
            continue
        tested += 1
        refc, altc, depth = ref_alt_counts_at_pos(
            bam,
            chrom,
            pos0,
            ref,
            alt,
            min_mapq=min_mapq,
            min_bq=min_bq,
        )
        usable = refc + altc
        if usable < min_depth:
            continue
        used += 1
        af_values.append(altc / usable)

    if not af_values:
        return 0, float("nan"), float("nan"), float("nan"), tested
    max_dev = max(abs(af - 0.5) for af in af_values)
    return (
        used,
        float(statistics.mean(af_values)),
        statistics.median(af_values),
        max_dev,
        tested,
    )


def main() -> int:
    p = argparse.ArgumentParser(description="Run DUP realism QA checks against spike truth.")
    p.add_argument("--bam", required=True, help="Spiked BAM/CRAM")
    p.add_argument("--truth", required=True, help="Truth VCF containing DUP entries")
    p.add_argument("--out", default="dup_realism_report.tsv", help="Output TSV path")
    p.add_argument("--ref", default=None, help="Optional: reference FASTA for CRAM input")
    p.add_argument("--junction-window", type=int, default=500, help="Window around event for split-read scan")
    p.add_argument("--coverage-flank", type=int, default=5000, help="Flank bp for depth comparison")
    p.add_argument("--min-mapq", type=int, default=20, help="Minimum mapping quality")
    p.add_argument("--min-baseq", type=int, default=20, help="Minimum base quality for depth/AF")
    p.add_argument("--min-baf-depth", type=int, default=10, help="Minimum usable depth per SNP for AF")
    p.add_argument("--snp-list", default=None, help="Optional SNP list: chrom pos ref alt")
    p.add_argument("--snp-list-1based", action="store_true", help="Interpret SNP positions as 1-based")
    p.add_argument("--truth-1based", action="store_true", help="Truth VCF uses 1-based positions")
    p.add_argument("--max-events", type=int, default=0, help="Process only first N events")
    p.add_argument("--event", default=None, help="Optional filter: chr:start-end")

    args = p.parse_args()

    use_zero_based_truth = not args.truth_1based
    events = parse_truth_vcfs(args.truth, use_zero_based=use_zero_based_truth)
    if not events:
        print("No DUP events found in truth VCF", file=sys.stderr)
        return 1

    if args.event:
        try:
            e_chr, e_rng = args.event.split(":", 1)
            e_start, e_end = e_rng.split("-", 1)
            e_start_i = int(e_start)
            e_end_i = int(e_end)
        except Exception:
            print("Invalid --event format. Expected chr:start-end", file=sys.stderr)
            return 1
        events = [ev for ev in events if ev.chrom == e_chr and ev.start >= e_start_i and ev.end <= e_end_i]
        if not events:
            print("No events match --event", file=sys.stderr)
            return 1

    if args.max_events > 0:
        events = events[: args.max_events]

    snps: List[Tuple[str, int, str, str]] = []
    if args.snp_list:
        snps = parse_snp_list(args.snp_list)

    with pysam.AlignmentFile(args.bam, "rb", reference_filename=args.ref) as bam:
        with open(args.out, "w", newline="", encoding="utf-8") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow([
                "event_id",
                "chrom",
                "start",
                "end",
                "len",
                "expected_vaf",
                "event_depth",
                "left_depth",
                "right_depth",
                "flank_depth",
                "coverage_ratio",
                "junction_sa_reads",
                "junction_softclip_reads",
                "window_reads",
                "baf_sites",
                "baf_used",
                "baf_mean",
                "baf_median",
                "baf_max_dev",
            ])

            for ev in events:
                start = max(ev.start, 0)
                end = max(ev.end, start + 1)

                evt_depth = mean_depth(
                    bam,
                    ev.chrom,
                    start,
                    end,
                    min_mapq=args.min_mapq,
                    min_bq=args.min_baseq,
                )

                left_start = max(0, start - args.coverage_flank)
                left_end = start
                right_start = end
                right_end = end + args.coverage_flank

                left_depth = mean_depth(
                    bam,
                    ev.chrom,
                    left_start,
                    left_end,
                    min_mapq=args.min_mapq,
                    min_bq=args.min_baseq,
                )
                right_depth = mean_depth(
                    bam,
                    ev.chrom,
                    right_start,
                    right_end,
                    min_mapq=args.min_mapq,
                    min_bq=args.min_baseq,
                )
                flank_depth = (left_depth + right_depth) / 2.0 if (left_depth + right_depth) > 0 else 0.0
                ratio = evt_depth / flank_depth if flank_depth > 0 else float("nan")

                q_start = max(0, start - args.junction_window)
                q_end = end + args.junction_window
                sa_cnt, soft_cnt, win_reads = count_split_softreads(
                    bam,
                    ev.chrom,
                    q_start,
                    q_end,
                    min_mapq=args.min_mapq,
                )

                baf_sites = "NA"
                baf_used = "0"
                baf_mean = "NA"
                baf_median = "NA"
                baf_max_dev = "NA"

                if snps:
                    used, mean_af, median_af, max_dev, tested = baf_from_sites(
                        bam,
                        ev.chrom,
                        start,
                        end,
                        snps,
                        min_mapq=args.min_mapq,
                        min_bq=args.min_baseq,
                        min_depth=args.min_baf_depth,
                        sites_are_one_based=args.snp_list_1based,
                    )
                    baf_sites = str(tested)
                    baf_used = str(used)
                    if not math.isnan(mean_af):
                        baf_mean = f"{mean_af:.4f}"
                        baf_median = f"{median_af:.4f}"
                        baf_max_dev = f"{max_dev:.4f}"

                writer.writerow(
                    [
                        f"dup_{ev.idx + 1}",
                        ev.chrom,
                        start,
                        end,
                        end - start,
                        f"{ev.vaf:.4f}",
                        f"{evt_depth:.4f}",
                        f"{left_depth:.4f}",
                        f"{right_depth:.4f}",
                        f"{flank_depth:.4f}",
                        f"{ratio:.4f}" if ratio == ratio else "NA",
                        sa_cnt,
                        soft_cnt,
                        win_reads,
                        baf_sites,
                        baf_used,
                        baf_mean,
                        baf_median,
                        baf_max_dev,
                    ]
                )

    print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
