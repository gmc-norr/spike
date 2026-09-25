#!/usr/bin/env python3
"""N18: does validate under-count small-indel carriers because of read ends?

Re-implements `spike validate`'s allele-fraction counting (the N13 pad rule
for an indel, the base at POS for an SNV, one vote per fragment and none for
mates that disagree, MAPQ >= 20, N14's binomial grade at an expected 0.5),
then recounts keeping only reads that reach F bases past the site's repeat
region on each side. REVIEW.md's N18 plan says what is computed and what
counts as supported or refuted.

The site lists are plain truth VCFs in spike's format (SIM_VAF=0.5):

    # het indels on chr20 (N15's 6,663 sites)
    bcftools view -H -r chr20 -f PASS -g het -m2 -M2 -v indels \\
        HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz   # REF, ALT <= 11 bp, shared anchor
    # SNVs in chr20:38-40 Mb (N16's set)

Usage:
    n18_indel_flank.py --bam B --ref R --indels I.vcf --snvs S.vcf \\
        [--check-log spike-debug.err --check-indel-json I.json --check-snv-json S.json]
"""

import argparse
import json
import math
import re
import statistics

import pysam

MIN_MAPQ = 20
INDEL_POS_PAD = 10
MIN_PILEUP_DEPTH = 5
ALT_ERROR_RATE = 0.001
REF_ERROR_RATE = 0.01
AF_TAIL = 0.005
AF_POWER = 0.99
MIN_ALT_READS = 3
FLANKS = [None, 0, 5, 10, 20, 30]  # None: validate's own rule, no flank filter


# --- N14's grade, as validate.rs has it ------------------------------------

def binomial_pmf(n, p):
    step = math.log(p) - math.log(1 - p)
    ln = n * math.log(1 - p)
    out = []
    for k in range(n + 1):
        out.append(math.exp(ln))
        if k < n:
            ln += math.log(n - k) - math.log(k + 1) + step
    return out


def alt_error_floor(n):
    at_least = 1.0
    for k, prob in enumerate(binomial_pmf(n, ALT_ERROR_RATE)):
        if at_least < AF_TAIL:
            return max(k, MIN_ALT_READS)
        at_least -= prob
    return n + 1


def grade(alt, total, vaf=0.5):
    """'nocov', 'low', 'shallow' (not evaluable), or 'pass' / 'fail'."""
    if total == 0:
        return "nocov"
    if total < MIN_PILEUP_DEPTH:
        return "low"
    p = vaf * (1 - REF_ERROR_RATE) + (1 - vaf) * ALT_ERROR_RATE
    floor = alt_error_floor(total)
    pmf = binomial_pmf(total, p)
    if sum(pmf[floor:]) < AF_POWER:
        return "shallow"
    x = min(alt, total)
    at_most, at_least = sum(pmf[: x + 1]), sum(pmf[x:])
    in_range = at_most >= AF_TAIL and (vaf >= 1.0 or at_least >= AF_TAIL)
    return "pass" if alt >= floor and in_range else "fail"


# --- reads -----------------------------------------------------------------

def usable(r):
    return not (
        r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail
    ) and r.mapping_quality >= MIN_MAPQ and r.query_name


def indel_vote(r, pos, ref_len, kind, length):
    """validate.rs cigar_indel_vote (pad rule): 'C', 'S' or None."""
    ref_pos = r.reference_start
    carries = cov_anchor = cov_far = False
    far = pos + ref_len
    for op, n in r.cigartuples:
        if op in (0, 7, 8):
            cov_anchor |= ref_pos <= pos < ref_pos + n
            cov_far |= ref_pos <= far < ref_pos + n
            ref_pos += n
        elif (op == 2 and kind == "D") or (op == 1 and kind == "I"):
            if n == length and abs(ref_pos - (pos + 1)) <= INDEL_POS_PAD:
                carries = True
            if op == 2:
                ref_pos += n
        elif op in (2, 3):
            ref_pos += n
    if carries:
        return "C"
    if cov_anchor and cov_far:
        return "S"
    return None


def fragment_votes(reads):
    """One vote per read name; mates that disagree give none (N16)."""
    by_name = {}
    for name, vote in reads:
        by_name.setdefault(name, []).append(vote)
    out = []
    for votes in by_name.values():
        if all(v == votes[0] for v in votes):
            out.append(votes[0])
    return out


# --- sites -----------------------------------------------------------------

def load_sites(path):
    sites = []
    for line in open(path):
        if line.startswith("#"):
            continue
        c = line.split("\t")
        sites.append((c[0], int(c[1]), c[3].upper(), c[4].upper()))
    return sites


def repeat_region(fasta, chrom, pos, ref, alt):
    """[s, e): the bases an indel at 0-based anchor `pos` can slide along."""
    lo = max(0, pos - 2000)
    seq = fasta.fetch(chrom, lo, pos + 2000).upper()
    at = lambda i: seq[i - lo] if 0 <= i - lo < len(seq) else None
    if len(ref) > len(alt):  # deletion of ref[1:]
        length = len(ref) - len(alt)
        s, e = pos + 1, pos + 1 + length
        while at(s - 1) is not None and at(s - 1) == at(s - 1 + length):
            s -= 1
        while at(e) is not None and at(e) == at(e - length):
            e += 1
    else:  # insertion of alt[1:] before pos + 1
        ins = alt[1:]
        length = len(ins)
        s = e = pos + 1
        k = 0
        while at(e) is not None and at(e) == ins[k % length]:
            e += 1
            k += 1
        k = 0
        while at(s - 1) is not None and at(s - 1) == ins[(length - 1 - k) % length]:
            s -= 1
            k += 1
    return s, e


def qualifies(start, end, s, e, flank, indel):
    """Aligned span reaches `flank` bases past [s, e) on each side, anchors
    included for an indel. `flank` None: no filter."""
    if flank is None:
        return True
    if indel:
        return start <= s - 1 - flank and end >= e + 1 + flank
    return start <= s - flank and end >= e + flank


def measure_indels(bam, fasta, sites):
    """Per site: its region and every usable read's (name, vote, start, end)."""
    out = []
    for chrom, vcf_pos, ref, alt in sites:
        pos = vcf_pos - 1
        kind = "D" if len(ref) > len(alt) else "I"
        length = abs(len(ref) - len(alt))
        s, e = repeat_region(fasta, chrom, pos, ref, alt)
        reads = []
        for r in bam.fetch(chrom, pos, pos + len(ref) + 1):
            if not usable(r):
                continue
            v = indel_vote(r, pos, len(ref), kind, length)
            if v is not None:
                reads.append((r.query_name, v, r.reference_start, r.reference_end))
        out.append((vcf_pos, s, e, reads))
    return out


def measure_snvs(bam, sites):
    out = []
    for chrom, vcf_pos, ref, alt in sites:
        pos = vcf_pos - 1
        reads = []
        for r in bam.fetch(chrom, pos, pos + 1):
            if not usable(r):
                continue
            for qp, rp in r.get_aligned_pairs(matches_only=True):
                if rp == pos:
                    base = r.query_sequence[qp].upper()
                    if base in "ACGT":
                        reads.append((r.query_name, base, r.reference_start, r.reference_end))
                    break
        out.append((vcf_pos, pos, pos + 1, alt, reads))
    return out


def count_indel(site, flank):
    _, s, e, reads = site
    votes = fragment_votes(
        (n, v) for n, v, a, b in reads if qualifies(a, b, s, e, flank, indel=True)
    )
    return votes.count("C"), votes.count("S")


def count_snv(site, flank):
    _, s, e, alt, reads = site
    votes = fragment_votes(
        (n, v) for n, v, a, b in reads if qualifies(a, b, s, e, flank, indel=False)
    )
    return votes.count(alt), len(votes)


def summarise(label, counts):
    """counts: list of (alt, total). Returns the metric row."""
    grades = [grade(a, t) for a, t in counts]
    evaluable = [g for g in grades if g in ("pass", "fail")]
    fails = sum(g == "fail" for g in evaluable)
    fracs = [a / t for a, t in counts if t >= 10]
    rate = 100.0 * fails / len(evaluable) if evaluable else float("nan")
    mean = statistics.mean(fracs) if fracs else float("nan")
    kept = sum(t for _, t in counts)
    print(
        f"{label}\tsites={len(counts)}\tevaluable={len(evaluable)}\tout_of_range={fails} "
        f"({rate:.2f}%)\tmean_fraction={mean:.4f} (n>=10: {len(fracs)})\tfragments={kept}"
    )
    return rate, mean


# --- the ruler check against validate itself --------------------------------

def check_ruler(indels, snvs, log, indel_json, snv_json):
    ok = True
    if log:
        spike = {}
        for line in open(log):
            m = re.search(r"allele_freq SNP \S+:(\d+)-\d+ \(unknown\): (\d+) carry, (\d+) span", line)
            if m:
                spike[int(m.group(1)) + 1] = (int(m.group(2)), int(m.group(3)))
        mine = {site[0]: count_indel(site, None) for site in indels}
        same = sum(mine[p] == spike.get(p) for p in mine)
        share = same / len(mine)
        print(f"ruler\tindel carries/spans match validate at {same}/{len(mine)} sites ({100*share:.2f}%)")
        ok &= share >= 0.99
    if indel_json:
        spike = verdicts(indel_json)
        mine = {site[0]: grade(*to_alt_total(count_indel(site, None))) == "pass" for site in indels}
        same = sum(mine[p] == spike.get(p) for p in mine)
        share = same / len(mine)
        print(f"ruler\tindel verdicts match validate at {same}/{len(mine)} sites ({100*share:.2f}%)")
        ok &= share >= 0.99
    if snv_json:
        d = json.load(open(snv_json))
        spike = {
            int(re.search(r":(\d+)-", c["event"]).group(1)) + 1: c["observed"]
            for c in d["checks"] if c["check"] == "allele_freq"
        }
        same = n = 0
        for site in snvs:
            alt, total = count_snv(site, None)
            want = spike.get(site[0])
            if want is None:
                continue
            n += 1
            g = grade(alt, total)
            if g == "nocov":
                got = "no coverage"
            elif g == "low":
                got = f"low depth ({total})"
            elif g == "shallow":
                got = f"too shallow ({alt / total:.2f} at {total} reads"
            else:
                got = f"{alt / total:.2f}"
            same += want == got or (g == "shallow" and want.startswith(got))
        share = same / n
        print(f"ruler\tSNV observed fractions match validate at {same}/{n} sites ({100*share:.2f}%)")
        ok &= share >= 0.99
    return ok


def to_alt_total(ct):
    c, s = ct
    return c, c + s


def verdicts(path):
    d = json.load(open(path))
    return {
        int(re.search(r":(\d+)-", c["event"]).group(1)) + 1: c["pass"]
        for c in d["checks"] if c["check"] == "allele_freq"
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--indels", required=True)
    ap.add_argument("--snvs", required=True)
    ap.add_argument("--check-log")
    ap.add_argument("--check-indel-json")
    ap.add_argument("--check-snv-json")
    a = ap.parse_args()

    bam = pysam.AlignmentFile(a.bam)
    fasta = pysam.FastaFile(a.ref)
    indels = measure_indels(bam, fasta, load_sites(a.indels))
    snv_sites = [x for x in load_sites(a.snvs) if len(x[2]) == 1 and len(x[3]) == 1]
    snvs = measure_snvs(bam, snv_sites)

    if not check_ruler(indels, snvs, a.check_log, a.check_indel_json, a.check_snv_json):
        print("RULER CHECK FAILED: nothing below is to be read")
        return 1

    for flank in FLANKS:
        tag = "unfiltered" if flank is None else f"F={flank}"
        summarise(f"indel\t{tag}", [to_alt_total(count_indel(s, flank)) for s in indels])
    for flank in FLANKS:
        tag = "unfiltered" if flank is None else f"F={flank}"
        summarise(f"snv\t{tag}", [count_snv(s, flank) for s in snvs])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
