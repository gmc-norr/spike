#!/usr/bin/env python3
"""T7: the sample's own non-SNP variants inside each event's footprint (CR3).

For every event, counts the records in a sample VCF that
  1. overlap the footprint `[start - HAP_FLANK, end + HAP_FLANK)`,
  2. are non-SNP -- not every ALT a single base against a single-base REF; a
     symbolic `<...>` ALT counts as non-SNP,
  3. and the sample carries (GT not 0/0, 0|0, ./. or .|.).
Each counted record goes in one size class by |len(longest ALT) - len(REF)|, or by
SVLEN when the ALT is symbolic: 1, 2-5, 6-20, 21-50, >50 bp.

Also reports, per footprint, the count WITHOUT the non-SNP test -- the plan's first
control: if the two totals are close, the filter is excluding nothing.

Overlap is computed here from POS and len(REF), not delegated to a bcftools flag:
bcftools 1.9 has no `--regions-overlap` and the default moved between versions.

Usage: t7_footprint_variants.py EVENTS_FILE VCF [FLANK]
EVENTS_FILE holds one `<type>:<chrom>:<start>-<end>` spec per line.
"""
import subprocess
import sys

HAP_FLANK = 2000
CLASSES = [(1, 1, "1bp"), (2, 5, "2-5bp"), (6, 20, "6-20bp"),
           (21, 50, "21-50bp"), (51, None, ">50bp")]
NO_CALL = {"0/0", "0|0", "./.", ".|."}


def sh(command):
    """Run a command, refusing to return silence on a failure.

    T3's census returned 0.0 everywhere because a samtools option did not exist:
    the tool printed usage to stderr, exited non-zero and wrote nothing, and the
    caller summed nothing to zero. Fail loudly.
    """
    done = subprocess.run(command, capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(f"command failed ({done.returncode}): {' '.join(command)}\n"
                           f"{done.stderr[:2000]}")
    return done.stdout


def size_class(delta):
    for lo, hi, name in CLASSES:
        if delta >= lo and (hi is None or delta <= hi):
            return name
    return None


def info_int(info, key):
    for field in info.split(";"):
        if field.startswith(key + "="):
            try:
                return abs(int(field.split("=", 1)[1].split(",")[0]))
            except ValueError:
                return None
    return None


def classify(ref, alts, info):
    """(is_non_snp, size class) for one record, or (False, None) for a pure SNP."""
    symbolic = any(a.startswith("<") for a in alts)
    if symbolic:
        svlen = info_int(info, "SVLEN")
        return True, size_class(svlen if svlen else 51)
    if len(ref) == 1 and all(len(a) == 1 for a in alts):
        return False, None
    delta = max(abs(len(a) - len(ref)) for a in alts)
    # An MNV (equal lengths, more than one base) changes no length but is not a
    # SNP either; SampleCopies cannot represent it any better than an indel.
    return True, size_class(delta if delta > 0 else 1)


def scan(chrom, lo, hi, vcf):
    """Every record whose reference span overlaps [lo, hi), 0-based half-open."""
    # Query a window widened by the longest REF this VCF holds so a record
    # starting before `lo` still comes back; 1-based inclusive for bcftools.
    pad = 60000
    region = f"{chrom}:{max(1, lo - pad + 1)}-{hi}"
    out = sh(["bcftools", "view", "-H", "-r", region, vcf])
    kept = []
    for line in out.splitlines():
        if not line:
            continue
        f = line.split("\t")
        pos = int(f[1])
        ref, alt, info = f[3], f[4], f[7]
        start = pos - 1
        end = start + len(ref)
        if end <= lo or start >= hi:
            continue
        gt = f[9].split(":")[0] if len(f) > 9 else "./."
        if gt in NO_CALL:
            continue
        alts = alt.split(",")
        non_snp, klass = classify(ref, alts, info)
        kept.append((non_snp, klass))
    return kept


def main(events_file, vcf, flank=str(HAP_FLANK)):
    flank = int(flank)
    events = [l for l in open(events_file).read().splitlines() if l.strip()]
    header = ["n", "event", "footprint", "all_carried", "non_snp"] + \
             [name for _, _, name in CLASSES]
    print("\t".join(header))
    rows = []
    for n, spec in enumerate(events, start=1):
        parts = spec.split(":")
        chrom = parts[1]
        start, end = (int(x) for x in parts[2].split("-"))
        lo, hi = max(0, start - flank), end + flank
        kept = scan(chrom, lo, hi, vcf)
        per_class = {name: 0 for _, _, name in CLASSES}
        non_snp = 0
        for is_non_snp, klass in kept:
            if is_non_snp:
                non_snp += 1
                if klass:
                    per_class[klass] += 1
        rows.append((n, spec, len(kept), non_snp, per_class))
        print("\t".join([str(n), spec, f"{chrom}:{lo}-{hi}", str(len(kept)), str(non_snp)] +
                        [str(per_class[name]) for _, _, name in CLASSES]))

    n_events = len(rows)
    fired = [r for r in rows if r[3] > 0]
    print()
    print(f"events={n_events} footprints with at least one carried non-SNP record="
          f"{len(fired)} ({100.0 * len(fired) / n_events:.1f}%)")
    all_carried = sum(r[2] for r in rows)
    all_non_snp = sum(r[3] for r in rows)
    print(f"control 1, is the filter doing anything: carried records in all footprints="
          f"{all_carried}, of which non-SNP={all_non_snp} "
          f"(SNP share {100.0 * (all_carried - all_non_snp) / all_carried:.1f}%)"
          if all_carried else "control 1: no carried records at all")
    counts = sorted(r[3] for r in rows)
    print(f"non-SNP records per footprint: min {counts[0]} median "
          f"{counts[len(counts) // 2]} max {counts[-1]}")
    for _, _, name in CLASSES:
        hit = sum(1 for r in rows if r[4][name] > 0)
        tot = sum(r[4][name] for r in rows)
        print(f"  class {name:>8}: fires on {hit} of {n_events} footprints "
              f"({100.0 * hit / n_events:.1f}%), {tot} records in all")


if __name__ == "__main__":
    main(*sys.argv[1:4])
