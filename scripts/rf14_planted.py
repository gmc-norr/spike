#!/usr/bin/env python3
"""RF14: `del_planted` -- are spike's own reads for a deletion there, carrying its join?

The replica of the proposed row, run before any Rust. The same idea as RF13's
`ins_planted`, for DEL:

- **Whose reads.** The truth ID is `sim_del_N`; spike names that event's reads
  `evNNNN_hap_*` (`simulate.rs`, `truth.rs`: both number `i + 1` over one event
  list). No other read counts.
- **Which records.** Those over START +/- 500 or END +/- 500 (the windows
  `split_reads` reads) that are not secondary or supplementary. MAPQ, duplicate,
  QC-fail and unmapped flags do not exclude a record.
- **Carrying.** The probe is the event haplotype's 31 bases across the join,
  `ref[START-15, START) + ref[END, END+16)`, in either orientation. A read carries
  it if some 31-base window of the read is within 2 substitutions of it: spike
  writes the sample's own SNPs onto the event copy and errors onto every read
  (case file, RF6). Every probe base is a reference base, so the tolerance covers
  all 31.
- **Verdict.** PASS at 1 or more carrying reads.
- **Not evaluable (FAIL).** The ID is not `sim_del_N`, or START < 15.

`START` is the truth record's POS read as `load_truth_events` reads a DEL (the
0-based first deleted base, which is the 1-based anchor), `END` its INFO END.

Negatives (for K):
  N1   the unspiked slice.bam
  N2a  END + 50, a truth whose deletion runs 50 bp too far
  N2b  START - 50, one that starts 50 bp too early
  N3   START + 1000 and END + 1000, the event moved
  N4   the ID renumbered to the next event

A negative is **indistinguishable** when its own probe, in either orientation, is
within 2 substitutions of some 31-base window of what spike planted, read from the
reference alone: `ref[START-2000, START) + ref[END, END+2000)`, the haplotype spike
tiles, to its 2 kb flanks. Then the planted reads spell that truth's join too, and
no rule reading 31 bases can tell the two truths apart. It is decided from the
reference before any read is counted.

Usage:
  rf14_planted.py k K_DIR REFERENCE           score every run under K_DIR
  rf14_planted.py null BAM SITES REFERENCE    N5 on the unspiked BAM
  rf14_planted.py one BAM REFERENCE TRUTH_VCF
"""
import json
import os
import re
import subprocess
import sys

PAD = 500
FLANK = 15
K = 31
MAX_MISMATCH = 2
HAP_FLANK = 2000
ID_RE = re.compile(r"^sim_del_(\d+)$")
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def run(cmd):
    """stdout of `cmd`, raising on a non-zero exit (case file: T3's census run)."""
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


def fetch(ref, chrom, start0, end0):
    start0 = max(0, start0)
    if end0 <= start0:
        return ""
    out = run(["samtools", "faidx", ref, f"{chrom}:{start0 + 1}-{end0}"])
    return "".join(out.split("\n")[1:]).upper()


def probe(ref, chrom, start, end):
    j = fetch(ref, chrom, start - FLANK, start) + fetch(ref, chrom, end, end + K - FLANK)
    return [j] if revcomp(j) == j else [j, revcomp(j)]


def within(window, p, limit=None):
    limit = MAX_MISMATCH if limit is None else limit
    mism = 0
    for a, b in zip(window, p):
        if a != b:
            mism += 1
            if mism > limit:
                return False
    return True


def carries(seq, ps):
    seq = seq.upper()
    return any(within(seq[i:i + K], p) for p in ps for i in range(len(seq) - K + 1))


def carriers(bam, ref, chrom, start, end, rec_id, any_name=False):
    m = ID_RE.match(rec_id)
    if not m:
        return f"id {rec_id} is not sim_del_N"
    if start < FLANK:
        return "start < 15"
    prefix = f"ev{int(m.group(1)):04}_hap_"
    ps = probe(ref, chrom, start, end)
    regions = [f"{chrom}:{max(0, p - PAD) + 1}-{p + PAD}" for p in (start, end)]
    names = set()
    for region in regions:
        for line in run(["samtools", "view", bam, region]).splitlines():
            f = line.split("\t", 10)
            if (not any_name and not f[0].startswith(prefix)) or int(f[1]) & (0x100 | 0x800):
                continue
            if f[0] in names:
                continue
            if carries(f[9], ps):
                names.add(f[0])
    return names


def indistinguishable(ref, chrom, start, end, neg_start, neg_end):
    """The negative's probe is spelled, within 2, by what spike planted for (start, end)."""
    hap = fetch(ref, chrom, start - HAP_FLANK, start) + fetch(ref, chrom, end, end + HAP_FLANK)
    ps = probe(ref, chrom, neg_start, neg_end)
    return any(within(hap[i:i + K], p) for p in ps for i in range(len(hap) - K + 1))


def del_records(truth):
    for line in open(truth):
        if not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            if "SVTYPE=DEL" in f[7]:
                end = int(re.search(r"(?:^|;)END=(\d+)", f[7]).group(1))
                yield f[0], int(f[1]), end, f[2], f


def count(c):
    return 0 if isinstance(c, str) else len(c)


def split_reads_verdict(run_dir):
    """Master's `split_reads` row on the run's own truth, from validate.json."""
    path = os.path.join(run_dir, "validate.json")
    try:
        rows = json.load(open(path))["checks"]
    except (OSError, ValueError, KeyError):
        return "no validate.json"
    r = [c for c in rows if c["check"] == "split_reads"]
    return ("PASS" if r[0]["pass"] else "FAIL") + f" {r[0]['observed']}" if r else "no row"


def main(argv):
    mode = argv[0]
    if mode == "one":
        bam, ref, truth = argv[1:4]
        for chrom, start, end, rec_id, _ in del_records(truth):
            c = carriers(bam, ref, chrom, start, end, rec_id)
            print(f"{chrom}:{start}-{end} {rec_id} -> {c if isinstance(c, str) else len(c)}")
        return
    if mode == "sites":
        # From the reference alone, before any read is counted: how close each
        # K positive's join is to the reference near its breakpoints (RF6's guard
        # distance), and which of its negatives are indistinguishable.
        sites, ref = argv[1:3]
        for line in open(sites):
            name, chrom, pos, end = line.rstrip("\n").split("\t")
            if name == "null":
                continue
            pos = int(pos)
            spans = [(pos, int(end))] if name == "hg" else [(pos, pos + n) for n in (50, 300, 1000, 10000)]
            for s, e in spans:
                j = probe(ref, chrom, s, e)
                near_ref = fetch(ref, chrom, s - 1000, s + 1000) + "N" + fetch(ref, chrom, e - 1000, e + 1000)
                dist = min(
                    sum(a != b for a, b in zip(near_ref[i:i + K], p))
                    for p in j for i in range(len(near_ref) - K + 1)
                )
                blind = [
                    label for label, (ns, ne) in (
                        ("N2a", (s, e + 50)), ("N2b", (s - 50, e)), ("N3", (s + 1000, e + 1000)),
                    ) if indistinguishable(ref, chrom, s, e, ns, ne)
                ]
                print(f"{name}\t{s}\t{e}\t{e - s}\tguard_distance={dist}\tindistinguishable={','.join(blind) or '-'}")
        return
    if mode == "null":
        bam, sites, ref = argv[1:4]
        passed = n = 0
        for line in open(sites):
            name, chrom, pos, _ = line.rstrip("\n").split("\t")
            if name != "null":
                continue
            c = carriers(bam, ref, chrom, int(pos), int(pos) + 1000, "sim_del_1")
            n += 1
            passed += count(c) >= 1
        print(f"N5 null sites passing: {passed} of {n}")
        return

    kdir, ref = argv[1:3]
    pos_rows, judged, blind, unfiltered = [], [], [], []
    n_judged = 0
    for name in sorted(os.listdir(kdir), key=lambda s: (len(s), s)):
        d = os.path.join(kdir, name)
        if not os.path.isdir(d):
            continue
        truth = os.path.join(d, "run", "truth.vcf")
        merged = os.path.join(d, "run", "merged.bam")
        if not (os.path.exists(truth) and os.path.exists(merged)):
            pos_rows.append((name, "did not reach validate", ""))
            continue
        chrom, start, end, rec_id, _ = next(del_records(truth))
        c = carriers(merged, ref, chrom, start, end, rec_id)
        pos_rows.append((name, c if isinstance(c, str) else len(c), split_reads_verdict(d)))
        n_id = int(ID_RE.match(rec_id).group(1)) if ID_RE.match(rec_id) else 0
        negs = {
            "N1 empty": (os.path.join(d, "slice.bam"), start, end, rec_id, False),
            "N2a END+50": (merged, start, end + 50, rec_id, True),
            "N2b START-50": (merged, start - 50, end, rec_id, True),
            "N3 moved +1000": (merged, start + 1000, end + 1000, rec_id, True),
            "N4 wrong id": (merged, start, end, f"sim_del_{n_id + 1}", False),
        }
        for label, (bam, s, e, i, can_blind) in negs.items():
            n = count(carriers(bam, ref, chrom, s, e, i))
            if can_blind and indistinguishable(ref, chrom, start, end, s, e):
                blind.append((name, label, n))
                continue
            n_judged += 1
            if n >= 1:
                judged.append((name, label, n))
        if name.startswith("h_"):
            u = carriers(os.path.join(d, "slice.bam"), ref, chrom, start, end, rec_id, any_name=True)
            unfiltered.append((name, count(u)))
    for name, c, sr in pos_rows:
        print(f"  {name}: {c}   (master split_reads: {sr})")
    ran = [c for _, c, _ in pos_rows if isinstance(c, int)]
    ok = sum(1 for c in ran if c >= 1)
    print(f"K+ runs reaching validate: {len(ran)} of {len(pos_rows)}")
    print(f"K+ carriers >= 1: {ok} of {len(ran)}  (min {min(ran) if ran else '-'}, "
          f"median {sorted(ran)[len(ran) // 2] if ran else '-'})")
    print(f"K- judged negatives with a carrier: {len(judged)} of {n_judged} {judged}")
    print(f"indistinguishable negatives (not judged): {len(blind)} {blind}")
    print(f"check of the check, h_ sites' unspiked slice with the name filter off: {unfiltered}")


if __name__ == "__main__":
    main(sys.argv[1:])
