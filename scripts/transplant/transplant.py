#!/usr/bin/env python3
"""Transplant test, round 1: spiked deletions against real ones, HG001 <-> HG002.

The rule, the events and the controls are locked in
docs/superpowers/plans/2026-09-28-transplant-round1.md; this script only
carries them out.

  forward   HG001-only het deletions, spiked into HG002; real = HG001
  reverse   HG002-only het deletions, spiked into HG001; real = HG002
  shared    deletions het in both: HG001's reads against HG002's

d = fake minus real, per transplanted event. c, per shared event, is taken
the same way round as d -- recipient minus donor -- so forward compares with
HG002 minus HG001 (the plan's wording) and reverse with its mirror,
HG001 minus HG002. Settled before any evidence was read.

Usage: transplant.py pilot|full --out DIR --spike SPIKE --bench DIR [--threads N]
"""
import argparse
import gzip
import os
import random
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evidence  # noqa: E402

CV = os.path.expanduser("~/dev/cnv_validation/from_rv/raredisease_seracare_na12878_resample")
BAMS = {
    "HG001": f"{CV}/D25-7403_Seq25-7598_30x_resample/raredisease_results/alignment/D25-7403_Seq25-7598_30x_sorted_md.bam",
    "HG002": f"{CV}/D24-14230_Seq25-7600_30x_resample/raredisease_results/alignment/D24-14230_Seq25-7600_30x_sorted_md.bam",
}
REF = os.path.expanduser("~/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta")
SV_TRUTH = {
    "HG001": os.path.expanduser("~/dev/cnv_validation/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.svs.vcf.gz"),
    "HG002": os.path.expanduser("~/dev/spike/data/giab_hg38/HG002/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"),
}
# The real BAMs' own command (their @PG), so spike's reads are aligned as theirs were.
ALIGNER = "bwa-mem2 mem -M -K 100000000 -t {threads} -R '@RG\\tID:sim\\tPL:ILLUMINA\\tSM:sim'"
SETS = {"forward": ("HG001", "HG002"), "reverse": ("HG002", "HG001")}  # set: (donor, recipient)
BINS = ["50-299", "300-999", "1k-10k"]
MIN_DEL = 50
NEIGHBOUR = 1_000
RUN_GAP = 100_000
SEED = 1
MAX_RECIPIENT_EVIDENCE = 1
B1_VAF = 0.25
B2_SHIFT = 200
REFUSED = re.compile(r"^\s*DEL\s+(\S+):(\d+)-(\d+) \(\d+bp\): \d+ of \d+ reads over it", re.M)


# --- the logic the tests pin down ---------------------------------------------

def deletion_of(chrom, pos, ref, alt):
    """A pure deletion of >= 50 bp as (chrom, start, end), 0-based half-open; else None."""
    if len(alt) != 1 or alt != ref[0] or len(ref) - 1 < MIN_DEL:
        return None
    return (chrom, pos, pos + len(ref) - 1)


def size_bin(length):
    return "50-299" if length < 300 else "300-999" if length < 1000 else "1k-10k" if length < 10000 else "10k+"


def is_het(gt):
    return gt.replace("|", "/") in ("0/1", "1/0")


def isolated(events, others, margin=NEIGHBOUR):
    """The events with no other SV (not themselves) within `margin` bp."""
    by_chrom = {}
    for o in others:
        by_chrom.setdefault(o[0], []).append(o)
    keep = []
    for e in events:
        near = [o for o in by_chrom.get(e[0], []) if o != e and o[1] < e[2] + margin and o[2] > e[1] - margin]
        if not near:
            keep.append(e)
    return keep


def draw(cands, n, seed, keep):
    """Up to n of `cands`, in a seeded shuffle of their sorted order, taking those `keep` accepts."""
    order = sorted(cands)
    random.Random(seed).shuffle(order)
    out = []
    for e in order:
        if len(out) == n:
            break
        if keep(e):
            out.append(e)
    return out


def batches(events, gap=RUN_GAP):
    """Events split into spike runs, each event >= `gap` bp from every other in its run."""
    runs = []
    for e in events:
        for run in runs:
            if all(o[0] != e[0] or o[1] >= e[2] + gap or e[1] >= o[2] + gap for o in run):
                run.append(e)
                break
        else:
            runs.append([e])
    return runs


def parse_refused(stderr):
    """The events spike's RF8 refusal lists, as (chrom, start, end), 0-based half-open."""
    return [(m[0], int(m[1]) - 1, int(m[2])) for m in REFUSED.findall(stderr)]


def percentile(values, q):
    """Linear interpolation between order statistics (numpy's default)."""
    v = sorted(values)
    pos = (len(v) - 1) * q / 100
    lo = int(pos)
    hi = min(lo + 1, len(v) - 1)
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


def judge(d, c):
    """The locked rule: median(d) within c's 25th-75th, and d's 10th-90th width <= 1.5x c's."""
    med = percentile(d, 50)
    c25, c75 = percentile(c, 25), percentile(c, 75)
    wd = percentile(d, 90) - percentile(d, 10)
    wc = percentile(c, 90) - percentile(c, 10)
    bias_ok = c25 <= med <= c75
    spread_ok = wd <= 1.5 * wc
    return {"n_d": len(d), "n_c": len(c), "median_d": med, "c25": c25, "c75": c75,
            "width_d": wd, "width_c": wc, "bias_ok": bias_ok, "spread_ok": spread_ok,
            "pass": bias_ok and spread_ok}


# --- inputs -------------------------------------------------------------------

def proper_bound(bam_path):
    """bwa's proper-pair upper bound for this library, by bwa's own rule (mem_pestat), from chr1:50-60 Mb."""
    tl = []
    with pysam.AlignmentFile(bam_path) as bam:
        for r in bam.fetch("chr1", 50_000_000, 60_000_000):
            if r.flag & 0xF0C or not (r.is_paired and r.is_read1) or r.mapping_quality < 20:
                continue
            if r.next_reference_id != r.reference_id:
                continue
            t = r.template_length
            if (not r.is_reverse and r.mate_is_reverse and t > 0) or (r.is_reverse and not r.mate_is_reverse and t < 0):
                tl.append(abs(t))
    tl.sort()
    n = len(tl)
    p25, p75 = tl[int(.25 * n + .499)], tl[int(.75 * n + .499)]
    iqr = p75 - p25
    lo, hi = max(1, int(p25 - 2 * iqr + .499)), int(p75 + 2 * iqr + .499)
    kept = [x for x in tl if lo <= x <= hi]
    avg = sum(kept) / len(kept)
    std = (sum(x * x for x in kept) / len(kept) - avg * avg) ** 0.5
    high = int(p75 + 3 * iqr + .499)
    if high < avg + 4 * std:
        high = int(avg + 4 * std + .499)
    return high


def sv_intervals(sample):
    """Every truth SV (>= 50 bp either way) of `sample`, as (chrom, start, end)."""
    out = subprocess.run(["bcftools", "query", "-i", "abs(strlen(REF)-strlen(ALT))>=50",
                          "-f", "%CHROM\t%POS\t%REF\t%ALT\n", SV_TRUTH[sample]],
                         check=True, capture_output=True, text=True).stdout
    ivs = []
    for line in out.splitlines():
        chrom, pos, ref, alt = line.split("\t")
        pos = int(pos)
        ivs.append((chrom, pos, pos + max(1, len(ref) - 1)))
    return ivs


def candidates(bench):
    """forward, reverse and shared deletions, pure and het where the plan needs het."""
    def read(name):
        with pysam.VariantFile(os.path.join(bench, name)) as vf:
            for rec in vf:
                gt = "/".join("." if a is None else str(a) for a in rec.samples[0]["GT"])
                yield rec, gt
    sets = {"forward": [], "reverse": [], "shared": []}
    for key, name in (("forward", "fn.vcf.gz"), ("reverse", "fp.vcf.gz")):
        for rec, gt in read(name):
            e = deletion_of(rec.chrom, rec.pos, rec.ref, rec.alts[0]) if len(rec.alts) == 1 else None
            if e and is_het(gt):
                sets[key].append(e)
    comp = {rec.info["MatchId"]: (rec, gt) for rec, gt in read("tp-comp.vcf.gz")}
    for rec, gt in read("tp-base.vcf.gz"):
        other = comp.get(rec.info["MatchId"])
        e = deletion_of(rec.chrom, rec.pos, rec.ref, rec.alts[0]) if len(rec.alts) == 1 else None
        if e and other and is_het(gt) and is_het(other[1]):
            sets["shared"].append(e)
    return sets


# --- running spike, the aligner and the sham --------------------------------------

def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, **kw)


def event_spec(e, arm):
    chrom, start, end = e
    if arm == "B2":
        start, end = start + B2_SHIFT, end + B2_SHIFT
    vaf = B1_VAF if arm == "B1" else 0.5
    return f"del:{chrom}:{start}-{end};af={vaf}"


def spike_run(spike, recipient, events, arm, out_dir, threads, log):
    """spike on `events` until it accepts them; returns (accepted, refused)."""
    todo, refused = list(events), []
    while todo:
        if os.path.exists(out_dir):
            sh(["rm", "-rf", out_dir])
        cmd = [spike, "--bam", BAMS[recipient], "--reference", REF, "--seed", str(SEED),
               "--threads", str(threads), "--aligner", ALIGNER.format(threads=threads), "--align", "-o", out_dir]
        for e in todo:
            cmd += ["--event", event_spec(e, arm)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(cmd)}\nexit {p.returncode}\n{p.stderr}\n")
        if p.returncode == 0:
            return todo, refused
        listed = parse_refused(p.stderr)
        if not listed:
            raise RuntimeError(f"spike failed in {out_dir} for a reason other than RF8; see {log}")
        # A B2 event is refused at its shifted place; map it back to its truth event.
        shift = B2_SHIFT if arm == "B2" else 0
        back = {(c, s - shift, en - shift) for c, s, en in listed}
        refused += [e for e in todo if e in back]
        todo = [e for e in todo if e not in back]
    return [], refused


def sham(recipient, events, run_dir, threads):
    """The originals spike replaced, realigned unchanged with the same command: sham.bam."""
    names = os.path.join(run_dir, "replaced_reads.txt")
    regions = [f"{c}:{max(1, s - 20_000)}-{en + 20_000}" for c, s, en in events]
    raw = os.path.join(run_dir, "sham_orig.bam")
    sh(["samtools", "view", "-b", "-F", "0x900", "-N", names, "-o", raw, BAMS[recipient]] + regions)
    r1, r2 = os.path.join(run_dir, "sham_R1.fq.gz"), os.path.join(run_dir, "sham_R2.fq.gz")
    sh(f"samtools collate -u -O {raw} | samtools fastq -n -1 {r1} -2 {r2} -0 /dev/null -s /dev/null -",
       shell=True, stderr=subprocess.DEVNULL)
    out = os.path.join(run_dir, "sham.bam")
    sh(f"{ALIGNER.format(threads=threads)} {REF} {r1} {r2} 2> {run_dir}/sham_align.log"
       f" | samtools sort -@ {threads} -o {out} - && samtools index {out}", shell=True, executable="/bin/bash")
    return out


def skip_set(run_dir):
    with open(os.path.join(run_dir, "replaced_reads.txt")) as fh:
        return {line.strip() for line in fh if line.strip()}


# --- the whole thing ----------------------------------------------------------

def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v)) for v in row) + "\n")


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["pilot", "full"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--spike", required=True)
    ap.add_argument("--bench", required=True)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    log = os.path.join(a.out, "spike.log")
    per_bin, bins = (10, ["300-999"]) if a.stage == "pilot" else (100, BINS)
    arms = {"forward": ["normal", "B1"] if a.stage == "pilot" else ["normal", "B1", "B2"],
            "reverse": ["normal"]}

    bound = {s: proper_bound(b) for s, b in BAMS.items()}
    write_tsv(os.path.join(a.out, "bounds.tsv"), ["sample", "proper_pair_bound"], sorted(bound.items()))

    others = sv_intervals("HG001") + sv_intervals("HG002")
    cands = {k: isolated(v, others) for k, v in candidates(a.bench).items()}

    ev_rows, chosen = [], {}
    for key in ("forward", "reverse", "shared"):
        for b in bins:
            pool = [e for e in cands[key] if size_bin(e[2] - e[1]) == b]
            if key == "shared":
                keep = lambda e: True  # noqa: E731
            else:
                donor, recipient = SETS[key]

                def keep(e, recipient=recipient):
                    row = evidence.measure(BAMS[recipient], *e, bound[recipient])
                    ev_rows.append([key, b, "recipient_before", *e, *[row[c] for c in evidence.COLUMNS]])
                    return row["n_any"] <= MAX_RECIPIENT_EVIDENCE
            chosen[(key, b)] = draw(pool, per_bin, SEED, keep)
    write_tsv(os.path.join(a.out, "events.tsv"), ["set", "bin", "chrom", "start", "end"],
              [[k, b, *e] for (k, b), es in chosen.items() for e in es])

    refused_rows = []
    for key in ("forward", "reverse"):
        donor, recipient = SETS[key]
        events = [e for b in bins for e in chosen[(key, b)]]
        for e in events:
            row = evidence.measure(BAMS[donor], *e, bound[donor])
            ev_rows.append([key, size_bin(e[2] - e[1]), "real", *e, *[row[c] for c in evidence.COLUMNS]])
        for arm in arms[key]:
            for i, run in enumerate(batches(events)):
                run_dir = os.path.join(a.out, key, arm, f"run{i}")
                accepted, refused = spike_run(a.spike, recipient, run, arm, run_dir, a.threads, log)
                refused_rows += [[key, arm, size_bin(e[2] - e[1]), *e] for e in refused]
                if not accepted:
                    continue
                skip = skip_set(run_dir)
                sim = os.path.join(run_dir, "sim.bam")
                sham_bam = sham(recipient, accepted, run_dir, a.threads) if arm == "normal" else None
                for e in accepted:
                    row = evidence.measure(BAMS[recipient], *e, bound[recipient], skip, [sim])
                    ev_rows.append([key, size_bin(e[2] - e[1]), f"fake_{arm}", *e, *[row[c] for c in evidence.COLUMNS]])
                    if sham_bam:
                        row = evidence.measure(BAMS[recipient], *e, bound[recipient], skip, [sham_bam])
                        ev_rows.append([key, size_bin(e[2] - e[1]), "sham", *e, *[row[c] for c in evidence.COLUMNS]])
    for b in bins:
        for e in chosen[("shared", b)]:
            for s in ("HG001", "HG002"):
                row = evidence.measure(BAMS[s], *e, bound[s])
                ev_rows.append(["shared", b, s, *e, *[row[c] for c in evidence.COLUMNS]])
    write_tsv(os.path.join(a.out, "evidence.tsv"), ["set", "bin", "which", "chrom", "start", "end"] + evidence.COLUMNS, ev_rows)
    write_tsv(os.path.join(a.out, "refused.tsv"), ["set", "arm", "bin", "chrom", "start", "end"], refused_rows)
    judge_all(a.out, ev_rows, bins, arms)


def judge_all(out, ev_rows, bins, arms):
    cols = ["set", "bin", "which", "chrom", "start", "end"] + evidence.COLUMNS
    rows = [dict(zip(cols, r)) for r in ev_rows]
    val = {}
    for r in rows:
        val[(r["set"], r["which"], r["chrom"], r["start"], r["end"])] = r
    results = []
    for key in ("forward", "reverse"):
        donor, recipient = SETS[key]
        for b in bins:
            shared = [r for r in rows if r["set"] == "shared" and r["bin"] == b and r["which"] == "HG001"]
            for metric in ("J", "E1"):
                c = []
                for r in shared:
                    k = ("shared", recipient, r["chrom"], r["start"], r["end"])
                    kd = ("shared", donor, r["chrom"], r["start"], r["end"])
                    if val[k][metric] is not None and val[kd][metric] is not None:
                        c.append(val[k][metric] - val[kd][metric])
                for arm in arms[key]:
                    d = []
                    for r in rows:
                        if r["set"] == key and r["bin"] == b and r["which"] == f"fake_{arm}":
                            real = val.get((key, "real", r["chrom"], r["start"], r["end"]))
                            if real and r[metric] is not None and real[metric] is not None:
                                d.append(r[metric] - real[metric])
                    if d and c:
                        j = judge(d, c)
                        results.append([key, arm, b, metric] + [j[k] for k in
                                       ("n_d", "n_c", "median_d", "c25", "c75", "width_d", "width_c", "bias_ok", "spread_ok", "pass")])
    write_tsv(os.path.join(out, "judge.tsv"), ["set", "arm", "bin", "metric", "n_d", "n_c", "median_d", "c25", "c75",
                                               "width_d", "width_c", "bias_ok", "spread_ok", "pass"], results)


if __name__ == "__main__":
    main(sys.argv[1:])
