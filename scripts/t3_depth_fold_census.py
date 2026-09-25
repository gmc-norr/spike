#!/usr/bin/env python3
"""T3: what CR2's depth-fold warnings are -- mappability or real depth.

For every event in EVENTS_FILE, measures two windows -- the worst bin spike named
in its own warning (only present when it warned), and the 2 kb anchor window the
tiling was scaled at -- four ways:

  pool_frag_depth  mean fragment depth over proper pairs whose both mates pass
                   MAPQ >= 20 and the standard flag filters (the donor pool's rule,
                   and the estimator SIM_DEPTH_FOLD itself uses)
  any_read_depth   mean read depth over mapped, primary, non-duplicate,
                   non-QC-fail records at ANY MAPQ
  lowmapq_share    the share of those records with MAPQ < 20
  gc               the GC fraction of the reference over the window

The anchor's position is not in the log, so it is identified: of the event's two
breakpoints, the one whose pool_frag_depth is within 10% of the `scaled_by` spike
printed. Neither within 10% -> UNIDENTIFIED, and the event is left out of the
classification.

Classifies each warning event by `fold_any = max(r, 1/r)`, `r` = the any-MAPQ
read-depth ratio bin/anchor, against 1.5 -- `census::DEPTH_FOLD_WARN_ABOVE`, the
threshold spike already warns at.

The control (the plan's amendment): for EVERY event, `fold_any` over every 1 kb bin
of a uniform grid across [span_start - 2000, span_end + 2000) -- the replacement
footprint -- against that event's anchor, reported as the per-event maximum.

Judges nothing: the plan's verdict rule is applied by whoever reads this output.

Usage: t3_depth_fold_census.py RUNS_DIR EVENTS_FILE BAM REFERENCE
RUNS_DIR holds <n>/log per event, numbered from 1 in EVENTS_FILE's order.
"""
import pathlib
import re
import subprocess
import sys

WARN_ABOVE = 1.5           # census::DEPTH_FOLD_WARN_ABOVE
ANCHOR_WINDOW = 2000       # estimate_coverage_at(..., 2000)
ANCHOR_TOLERANCE = 0.10    # the plan's 10%
HAP_FLANK = 2000           # the footprint the grid spans
GRID_BIN = 1000            # depth_fold's own BIN
# UNMAP | SECONDARY | QCFAIL | DUP | SUPPLEMENTARY, for `samtools view -F`.
EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
# `samtools depth` has no `--ff`; its `-G` ADDS to a default filter-out list that
# already holds UNMAP, SECONDARY, QCFAIL and DUP, so SUPPLEMENTARY is all that is
# missing. Its default min-MQ is 0, which is the any-MAPQ this measures. `-J`
# counts a position a read's CIGAR deletes as covered, so the depth is the read's
# reference span -- the same semantics as validate's `count_depth_in_region`.
DEPTH_FLAGS = "-a -Q 0 -q 0 -J -G SUPPLEMENTARY"

WARNING = re.compile(
    r"the donor's depth over (\S+):(\d+)-(\d+) is ([\d.]+)x, but every fragment this "
    r"event tiles is scaled by the ([\d.]+)x .*?\(([\d.]+)-fold\)"
)


def sh(command):
    """Run a shell pipeline and return stdout, refusing to return silence.

    An unrecognised samtools option makes samtools print usage to stderr, exit
    non-zero and write nothing to stdout -- which a caller summing stdout reads as
    a depth of 0.0. That happened here (`--ff` is `samtools view`'s spelling, not
    `samtools depth`'s) and every any-MAPQ depth came back 0.0 while the pool
    depths beside them read 55x. Fail loudly instead.
    """
    done = subprocess.run(command, shell=True, capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(f"command failed ({done.returncode}): {command}\n{done.stderr}")
    return done.stdout


def any_read_depth(bam, chrom, start, end):
    """Mean read depth over [start, end), 0-based half-open, at any MAPQ."""
    region = f"{chrom}:{start + 1}-{end}"
    out = sh(f"samtools depth {DEPTH_FLAGS} -r '{region}' '{bam}'")
    total = sum(int(line.split("\t")[2]) for line in out.splitlines() if line)
    return total / (end - start)


def mapq_counts(bam, chrom, start, end):
    """(all, mapq>=20) primary non-duplicate records overlapping the window."""
    region = f"{chrom}:{start + 1}-{end}"
    allc = sh(f"samtools view -c -F {EXCLUDE} '{bam}' '{region}'").strip()
    hi = sh(f"samtools view -c -F {EXCLUDE} -q 20 '{bam}' '{region}'").strip()
    return int(allc or 0), int(hi or 0)


def pool_fragments(bam, chrom, lo, hi):
    """Pool-eligible fragment spans overlapping [lo, hi), read once.

    One fragment per pair, taken from the leftmost mate (TLEN > 0) so no pair is
    counted twice; its span is [POS-1, POS-1+TLEN). Requires proper pair (0x2),
    mate mapped and this read at MAPQ >= 20; samtools has no mate-MAPQ filter, so
    the mate is checked through MQ:i where the aligner wrote one. The query window
    is widened so a pair starting before `lo` still counts.
    """
    pad = 2000
    region = f"{chrom}:{max(0, lo - pad) + 1}-{hi}"
    out = sh(
        f"samtools view -F {EXCLUDE} -f 0x2 -q 20 '{bam}' '{region}' "
        "| awk -F'\\t' '$9 > 0 {mq=\"\"; for (i=12; i<=NF; i++) if ($i ~ /^MQ:i:/) "
        "{split($i,a,\":\"); mq=a[3]} if (mq == \"\" || mq+0 >= 20) print $4, $9}'"
    )
    spans = []
    for line in out.splitlines():
        if not line:
            continue
        pos, tlen = line.split()
        start = int(pos) - 1
        spans.append((start, start + int(tlen)))
    return spans


def mean_over(spans, start, end):
    covered = sum(max(0, min(e, end) - max(s, start)) for s, e in spans)
    return covered / (end - start)


def gc(reference, chrom, start, end):
    out = sh(f"samtools faidx '{reference}' '{chrom}:{start + 1}-{end}'")
    seq = "".join(out.splitlines()[1:]).upper()
    acgt = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / acgt if acgt else float("nan")


def window(bam, reference, chrom, start, end, frags=None):
    spans = frags if frags is not None else pool_fragments(bam, chrom, start, end)
    allc, hi = mapq_counts(bam, chrom, start, end)
    return {
        "region": f"{chrom}:{start}-{end}",
        "start": start, "end": end,
        "pool_frag_depth": mean_over(spans, start, end),
        "any_read_depth": any_read_depth(bam, chrom, start, end),
        "lowmapq_share": (allc - hi) / allc if allc else float("nan"),
        "gc": gc(reference, chrom, start, end),
    }


def fold(ratio):
    if not ratio or ratio != ratio:
        return float("nan")
    return max(ratio, 1 / ratio)


def main(runs_dir, events_file, bam, reference):
    runs = pathlib.Path(runs_dir)
    events = pathlib.Path(events_file).read_text().splitlines()
    print("n\tevent\twarned\tspike_fold\tbin\tbin_pool\tbin_any\tbin_lowmq\tbin_gc"
          "\tanchor\tanc_pool\tanc_any\tanc_lowmq\tanc_gc\tfold_any\tclass\tgrid_max_fold_any")
    for n, spec in enumerate(events, start=1):
        log = runs / str(n) / "log"
        text = log.read_text() if log.exists() else ""
        found = WARNING.search(text)
        chrom = spec.split(":")[1]
        span_start, span_end = (int(x) for x in spec.split(":")[2].split("-"))

        # One pool read covering the whole footprint, reused by every window.
        lo, hi = max(0, span_start - HAP_FLANK - ANCHOR_WINDOW), span_end + HAP_FLANK + ANCHOR_WINDOW
        frags = pool_fragments(bam, chrom, lo, hi)

        # The anchor: whichever breakpoint's pool depth matches `scaled_by`.
        scaled_by = float(found.group(5)) if found else None
        candidates = []
        anchor = None
        for which, pos in (("start", span_start), ("end", span_end)):
            a = window(bam, reference, chrom,
                       max(0, pos - ANCHOR_WINDOW // 2), pos + ANCHOR_WINDOW // 2, frags)
            a["which"], a["pos"] = which, pos
            candidates.append(a)
            if scaled_by is not None and anchor is None and \
                    abs(a["pool_frag_depth"] - scaled_by) <= ANCHOR_TOLERANCE * scaled_by:
                anchor = a
        # A non-warning event has no `scaled_by` to match, so the control uses the
        # candidate the pool is deeper at -- stated, not silently chosen: the
        # control statistic is a ratio, and the deeper anchor is the conservative
        # one (it makes a thin bin's fold larger, never smaller).
        control_anchor = anchor or max(candidates, key=lambda c: c["pool_frag_depth"])

        # The control grid over the footprint.
        grid_max = 0.0
        g = max(0, span_start - HAP_FLANK)
        while g < span_end + HAP_FLANK:
            e = min(g + GRID_BIN, span_end + HAP_FLANK)
            if e - g >= GRID_BIN // 2:
                d = any_read_depth(bam, chrom, g, e)
                r = d / control_anchor["any_read_depth"] if control_anchor["any_read_depth"] else 0
                grid_max = max(grid_max, fold(r))
            g = e

        if not found:
            print(f'{n}\t{spec}\tno\t-\t-\t-\t-\t-\t-\t'
                  f'{control_anchor["which"]}@{control_anchor["pos"]}\t'
                  f'{control_anchor["pool_frag_depth"]:.1f}\t{control_anchor["any_read_depth"]:.1f}\t'
                  f'{control_anchor["lowmapq_share"]:.3f}\t{control_anchor["gc"]:.3f}\t-\t-\t'
                  f'{grid_max:.2f}')
            continue

        bchrom, bs, be = found.group(1), int(found.group(2)), int(found.group(3))
        b = window(bam, reference, bchrom, bs, be, frags)
        if anchor is None:
            cands = " ".join(f'{c["which"]}={c["pool_frag_depth"]:.1f}' for c in candidates)
            print(f'{n}\t{spec}\tyes\t{found.group(6)}\t{b["region"]}\t'
                  f'{b["pool_frag_depth"]:.1f}\t{b["any_read_depth"]:.1f}\t'
                  f'{b["lowmapq_share"]:.3f}\t{b["gc"]:.3f}\t'
                  f'UNIDENTIFIED(scaled_by={scaled_by:.1f}; {cands})\t-\t-\t-\t-\t-\t'
                  f'anchor_unidentified\t{grid_max:.2f}')
            continue
        r = b["any_read_depth"] / anchor["any_read_depth"] if anchor["any_read_depth"] else 0
        f_any = fold(r)
        klass = "mappability" if f_any <= WARN_ABOVE else "real_depth"
        print(f'{n}\t{spec}\tyes\t{found.group(6)}\t{b["region"]}\t'
              f'{b["pool_frag_depth"]:.1f}\t{b["any_read_depth"]:.1f}\t'
              f'{b["lowmapq_share"]:.3f}\t{b["gc"]:.3f}\t'
              f'{anchor["which"]}@{anchor["pos"]}\t{anchor["pool_frag_depth"]:.1f}\t'
              f'{anchor["any_read_depth"]:.1f}\t{anchor["lowmapq_share"]:.3f}\t'
              f'{anchor["gc"]:.3f}\t{f_any:.2f}\t{klass}\t{grid_max:.2f}')


if __name__ == "__main__":
    main(*sys.argv[1:5])
