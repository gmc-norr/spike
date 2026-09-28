#!/usr/bin/env python3
"""Transplant test, round 2: spiked SNVs, 1-49 bp indels and 50-299 bp duplications against real ones.

The rule, the events and the controls are locked in
docs/superpowers/plans/2026-09-28-transplant-round2.md; this script only
carries them out. Samples, BAMs, aligner, spiking and sham are round 1's
(transplant.py).

  forward   HG001-only het variants, spiked into HG002; real = HG001
  reverse   HG002-only het variants, spiked into HG001; real = HG002
  shared    variants het in both: c = recipient minus donor, as round 1

Usage: round2.py pilot|full --out DIR --spike SPIKE --prep DIR [--threads N] [--exclude EVENTS.tsv]
  --prep is count2.sh's output directory.
"""
import argparse
import bisect
import os
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import small_evidence as se  # noqa: E402
import transplant as t  # noqa: E402

GROUPS = ["SNV", "DEL1-4", "DEL5-19", "DEL20-49", "INS1-4", "INS5-19", "INS20-49", "DUP50-299"]
PILOT_GROUPS = ["SNV", "INS20-49", "DUP50-299"]
SMALL_MARGIN = 150
LONG_SPAN = 1_000  # records longer than this are checked for isolation apart
DRAW_SEED = {"pilot": 1, "full": 2}
UNIQUE_MAX = 1.10  # round 2b: a duplication's array over its copy, below this
# Which broken control must fail which metric (round 2's plan).
MUST_FAIL = {"A": ("B1",), "E": ("B1",), "J": ("B1", "B2")}
REFUSED_SNV = re.compile(r"^\s*SNV\s+(\S+):(\d+) (\S+)>(\S+): \d+ of \d+ reads over it", re.M)
REFUSED_DUP = re.compile(r"^\s*DUP\s+(\S+):(\d+)-(\d+) \(\d+bp\): \d+ of \d+ reads over it", re.M)
COLS = ["set", "group", "which", "kind", "chrom", "pos", "ref", "alt"]


# --- the logic the tests pin down ---------------------------------------------

def kind_of(ref, alt):
    """("SNV", 1), ("DEL", n) or ("INS", n) for a 1-49 bp pure indel with its anchor; else None."""
    if len(ref) == 1 and len(alt) == 1:
        return ("SNV", 1)
    if len(alt) == 1 and ref[0] == alt and 1 <= len(ref) - 1 <= 49:
        return ("DEL", len(ref) - 1)
    if len(ref) == 1 and alt[0] == ref and 1 <= len(alt) - 1 <= 49:
        return ("INS", len(alt) - 1)
    return None


def group_of(kind, length):
    if kind == "SNV":
        return "SNV"
    if kind == "DUP":
        return "DUP50-299" if 50 <= length <= 299 else None
    return kind + ("1-4" if length <= 4 else "5-19" if length <= 19 else "20-49")


def group_of_event(e):
    kind, _, _, ref, alt = e
    return group_of(kind, 1 if kind == "SNV" else abs(len(alt) - len(ref)))


def metrics_of(group):
    return ["A"] if group == "SNV" else ["J"] if group.startswith("DUP") else ["A", "E"]


def _inside(region, chrom, pos, ref):
    ivs = region.get(chrom, [])
    lo, hi = pos - 1, pos - 1 + len(ref)
    i = bisect.bisect_right(ivs, (lo, float("inf"))) - 1
    return i >= 0 and ivs[i][0] <= lo and hi <= ivs[i][1]


class Neighbours:
    """Every record of both samples, and the SVs, to ask what lies within a margin of a span."""

    def __init__(self, spans):  # (chrom, start, end, key); 1-based closed
        short, long_ = {}, {}
        for c, s, e, k in spans:
            (long_ if e - s > LONG_SPAN else short).setdefault(c, []).append((s, e, k))
        # By span only: a record's key and an SV's None cannot be ordered.
        self.short = {c: sorted(v, key=lambda x: x[:2]) for c, v in short.items()}
        self.short_starts = {c: [s for s, _, _ in v] for c, v in self.short.items()}
        self.long_starts, self.long_max_end = {}, {}
        for c, v in long_.items():
            v.sort(key=lambda x: x[:2])
            self.long_starts[c] = [s for s, _, _ in v]
            ends, m = [], -1
            for _, e, _ in v:
                m = max(m, e)
                ends.append(m)
            self.long_max_end[c] = ends

    def any_near(self, chrom, s, e, margin, own):
        lst, starts = self.short.get(chrom, []), self.short_starts.get(chrom, [])
        i = bisect.bisect_left(starts, s - margin - LONG_SPAN)
        while i < len(lst) and lst[i][0] <= e + margin:
            if lst[i][1] >= s - margin and lst[i][2] != own:
                return True
            i += 1
        j = bisect.bisect_right(self.long_starts.get(chrom, []), e + margin)
        return j > 0 and self.long_max_end[chrom][j - 1] >= s - margin


def classify(h1, h2, region, margin, svs=()):
    """forward / reverse / shared small variants: het in one sample and absent from the other, or het in both;
    inside one region interval; no other record of either sample, nor an SV, within `margin` bp."""
    g1 = {(c, p, r, a): g for c, p, r, a, g in h1}
    g2 = {(c, p, r, a): g for c, p, r, a, g in h2}
    spans = [(c, p, p + len(r) - 1, (c, p, r, a)) for c, p, r, a, _ in list(h1) + list(h2)]
    spans += [(c, s, e, None) for c, s, e in svs]
    near = Neighbours(spans)
    sets = {"forward": [], "reverse": [], "shared": []}
    for key in sorted(set(g1) | set(g2)):
        c, p, r, a = key
        k = kind_of(r, a)
        if k is None or not _inside(region, c, p, r):
            continue
        a1, a2 = g1.get(key), g2.get(key)
        if t.is_het(a1 or "") and a2 is None:
            which = "forward"
        elif t.is_het(a2 or "") and a1 is None:
            which = "reverse"
        elif t.is_het(a1 or "") and t.is_het(a2 or ""):
            which = "shared"
        else:
            continue
        if near.any_near(c, p, p + len(r) - 1, margin, key):
            continue
        sets[which].append((k[0], c, p, r, a))
    return sets


def isolated_dups(dups, svs, margin, segment):
    """The duplications with no SV of either sample, other than their own records, within `margin` bp."""
    by_chrom = {}
    for o in svs:
        by_chrom.setdefault(o[0], []).append(o)
    keep = []
    for event, own in dups:
        s, e = segment(event)
        if not any(o not in own and o[1] < e + margin and o[2] > s - margin for o in by_chrom.get(event[1], [])):
            keep.append(event)
    return keep


def unique_dups(dups, fetch, segment):
    """Round 2b: the duplications whose copy is not part of a longer tandem array (array < 1.10x the copy)."""
    keep = []
    for e in dups:
        s, en = segment(e)
        rs, re_ = se.repeat_region(fetch, e[1], s, en, e[4][len(e[3]):].upper())
        if (re_ - rs) / (en - s) < UNIQUE_MAX:
            keep.append(e)
    return keep


def event_spec(event, arm, segment=None):
    kind, chrom, pos, ref, alt = event
    vaf = t.B1_VAF if arm == "B1" else 0.5
    if kind == "DUP":
        s, e = segment
        if arm == "B2":
            s, e = s + t.B2_SHIFT, e + t.B2_SHIFT
        return f"dup:{chrom}:{s}-{e};af={vaf}"
    if arm == "B2":
        raise ValueError("B2 is run on duplications only")
    return f"snp:{chrom}:{pos}:{ref}:{alt};af={vaf}"


def parse_refused(stderr):
    """The events spike's RF8 refusal lists, as spike labels them."""
    return ([("SNV", c, int(p), r, a) for c, p, r, a in REFUSED_SNV.findall(stderr)]
            + [("DUP", c, int(s) - 1, int(e)) for c, s, e in REFUSED_DUP.findall(stderr)])


def refused_of(events, listed, arm, segment):
    """The events of a run that `listed` names; a listed event that is none of them stops the run."""
    def label(e):
        if e[0] == "DUP":
            s, en = segment(e)
            shift = t.B2_SHIFT if arm == "B2" else 0
            return ("DUP", e[1], s + shift, en + shift)
        return ("SNV", e[1], e[2], e[3], e[4])
    by_label = {label(e): e for e in events}
    unknown = [x for x in listed if x not in by_label]
    if unknown:
        raise RuntimeError(f"spike refused events that are not this run's: {unknown}")
    return [by_label[x] for x in listed]


def verdict(result, controls, metric):
    """pass / fail, or inconclusive when a control that must fail `metric` passed or is missing."""
    if any(controls.get(arm) is None or controls[arm]["pass"] for arm in MUST_FAIL[metric]):
        return "inconclusive"
    return "pass" if result["pass"] else "fail"


def pilot_stops(rows):
    """The pilot's stop reasons: B1 passing any metric. rows: set, arm, group, metric, ..., pass."""
    return [f"B1 passes {r[3]} in {r[2]}" for r in rows if r[1] == "B1" and r[-1]]


# --- inputs -------------------------------------------------------------------

def read_bed(path):
    region = {}
    with open(path) as fh:
        for line in fh:
            c, s, e = line.split("\t")[:3]
            region.setdefault(c, []).append((int(s), int(e)))
    return {c: sorted(v) for c, v in region.items()}


def read_table(path):
    with open(path) as fh:
        return [(c, int(p), r, a, g) for c, p, r, a, g in (line.rstrip("\n").split("\t") for line in fh)]


def dup_candidates(bench, fetch):
    """forward / reverse / shared 50-299 bp duplications: insertions that copy the reference beside them,
    each with its own truth records (for isolation)."""
    def read(name):
        with pysam.VariantFile(os.path.join(bench, name)) as vf:
            for rec in vf:
                gt = "/".join("." if a is None else str(a) for a in rec.samples[0]["GT"])
                yield rec, gt

    def dup(rec):
        if len(rec.alts) != 1 or len(rec.ref) != 1 or rec.alts[0][0] != rec.ref:
            return None
        e = ("DUP", rec.chrom, rec.pos, rec.ref, rec.alts[0])
        if group_of_event(e) is None or se.segment_of(fetch, rec.chrom, rec.pos, rec.ref, rec.alts[0]) is None:
            return None
        return e

    def own(rec):
        return (rec.chrom, rec.pos, rec.pos + max(1, len(rec.ref) - 1))
    sets = {"forward": [], "reverse": [], "shared": []}
    for key, name in (("forward", "fn.vcf.gz"), ("reverse", "fp.vcf.gz")):
        for rec, gt in read(name):
            e = dup(rec)
            if e and t.is_het(gt):
                sets[key].append((e, [own(rec)]))
    comp = {rec.info["MatchId"]: (rec, gt) for rec, gt in read("tp-comp.vcf.gz")}
    for rec, gt in read("tp-base.vcf.gz"):
        other = comp.get(rec.info["MatchId"])
        e = dup(rec)
        if e and other and t.is_het(gt) and t.is_het(other[1]):
            sets["shared"].append((e, [own(rec), own(other[0])]))
    return sets


# --- running spike ------------------------------------------------------------

def spike_run(spike, recipient, events, arm, out_dir, threads, log, segment):
    """spike on `events` until it accepts them; returns (accepted, refused)."""
    todo, refused = list(events), []
    while todo:
        if os.path.exists(out_dir):
            t.sh(["rm", "-rf", out_dir])
        cmd = [spike, "--bam", t.BAMS[recipient], "--reference", t.REF, "--seed", str(t.SEED),
               "--threads", str(threads), "--aligner", t.ALIGNER.format(threads=threads), "--align", "-o", out_dir]
        for e in todo:
            cmd += ["--event", event_spec(e, arm, segment(e) if e[0] == "DUP" else None)]
        p = subprocess.run(cmd, capture_output=True, text=True)
        with open(log, "a") as fh:
            fh.write(f"$ {' '.join(cmd)}\nexit {p.returncode}\n{p.stderr}\n")
        if p.returncode == 0:
            return todo, refused
        listed = parse_refused(p.stderr)
        if not listed:
            raise RuntimeError(f"spike failed in {out_dir} for a reason other than RF8; see {log}")
        gone = refused_of(todo, listed, arm, segment)
        refused += gone
        todo = [e for e in todo if e not in gone]
    return [], refused


# --- the whole thing ----------------------------------------------------------

def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["pilot", "full"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--spike", required=True)
    ap.add_argument("--prep", required=True, help="count2.sh's output directory")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--exclude", help="an events.tsv whose events may not be drawn (the pilot's)")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    log = os.path.join(a.out, "spike.log")
    per_group, groups = (10, PILOT_GROUPS) if a.stage == "pilot" else (100, GROUPS)
    fasta = pysam.FastaFile(t.REF)

    def fetch(chrom, start, end):
        return fasta.fetch(chrom, max(0, start), max(0, end)).upper()

    seg_cache = {}

    def segment(e):
        if e not in seg_cache:
            seg_cache[e] = se.segment_of(fetch, e[1], e[2], e[3], e[4])
        return seg_cache[e]

    def span(e):
        if e[0] == "DUP":
            return (e[1], *segment(e))
        return (e[1], e[2] - 1, e[2] - 1 + len(e[3]))

    excluded = set()
    if a.exclude:
        with open(a.exclude) as fh:
            next(fh)
            excluded = {(k, c, int(p), r, al) for _, _, k, c, p, r, al in (line.rstrip("\n").split("\t") for line in fh)}
    svs = t.sv_intervals("HG001") + t.sv_intervals("HG002")
    small = classify(read_table(os.path.join(a.prep, "h1.all.tsv")), read_table(os.path.join(a.prep, "h2.all.tsv")),
                     read_bed(os.path.join(a.prep, "small.bed")), SMALL_MARGIN, svs)
    dups = {k: unique_dups(isolated_dups(v, svs, t.NEIGHBOUR, segment), fetch, segment)
            for k, v in dup_candidates(os.path.join(a.prep, "insbench"), fetch).items()}
    cands = {k: small[k] + dups[k] for k in small}
    t.write_tsv(os.path.join(a.out, "pools.tsv"), ["set", "group", "n"],
                [[k, g, sum(1 for e in cands[k] if group_of_event(e) == g)] for k in cands for g in GROUPS])

    ev_rows, chosen = [], {}

    def row_of(key, which, e, row):
        return [key, group_of_event(e), which, *e, *[row[c] for c in se.COLUMNS]]

    for key in ("forward", "reverse", "shared"):
        for g in groups:
            pool = t.without([e for e in cands[key] if group_of_event(e) == g], excluded)
            if key == "shared":
                keep = lambda e: True  # noqa: E731
            else:
                recipient = t.SETS[key][1]

                def keep(e, key=key, recipient=recipient):
                    row = se.measure(t.BAMS[recipient], fetch, e)
                    ev_rows.append(row_of(key, "recipient_before", e, row))
                    return se.recipient_count(row, e[0]) <= t.MAX_RECIPIENT_EVIDENCE
            chosen[(key, g)] = t.draw(pool, per_group, DRAW_SEED[a.stage], keep)
    t.write_tsv(os.path.join(a.out, "events.tsv"), ["set", "group", "kind", "chrom", "pos", "ref", "alt"],
                [[k, g, *e] for (k, g), es in chosen.items() for e in es])

    arms = {"forward": ["normal", "B1", "B2"], "reverse": ["normal"]}
    refused_rows = []
    for key in ("forward", "reverse"):
        donor, recipient = t.SETS[key]
        events = [e for g in groups for e in chosen[(key, g)]]
        for e in events:
            ev_rows.append(row_of(key, "real", e, se.measure(t.BAMS[donor], fetch, e)))
        for arm in arms[key]:
            arm_events = [e for e in events if arm != "B2" or e[0] == "DUP"]
            spans = [(*span(e), i) for i, e in enumerate(arm_events)]
            for n, run in enumerate(t.batches(spans)):
                run_events = [arm_events[x[3]] for x in run]
                run_dir = os.path.join(a.out, key, arm, f"run{n}")
                accepted, refused = spike_run(a.spike, recipient, run_events, arm, run_dir, a.threads, log, segment)
                refused_rows += [[key, arm, group_of_event(e), *e] for e in refused]
                if not accepted:
                    continue
                skip = t.skip_set(run_dir)
                sim = os.path.join(run_dir, "sim.bam")
                sham_bam = t.sham(recipient, [span(e) for e in accepted], run_dir, a.threads) if arm == "normal" else None
                for e in accepted:
                    ev_rows.append(row_of(key, f"fake_{arm}", e, se.measure(t.BAMS[recipient], fetch, e, skip, [sim])))
                    if sham_bam:
                        ev_rows.append(row_of(key, "sham", e, se.measure(t.BAMS[recipient], fetch, e, skip, [sham_bam])))
    for g in groups:
        for e in chosen[("shared", g)]:
            for s in ("HG001", "HG002"):
                ev_rows.append(row_of("shared", s, e, se.measure(t.BAMS[s], fetch, e)))
    t.write_tsv(os.path.join(a.out, "evidence.tsv"), COLS + se.COLUMNS, ev_rows)
    t.write_tsv(os.path.join(a.out, "refused.tsv"), ["set", "arm", "group", "kind", "chrom", "pos", "ref", "alt"],
                refused_rows)
    results = judge_all(a.out, ev_rows, groups)
    if a.stage == "pilot":
        stops = pilot_stops(results)
        with open(os.path.join(a.out, "pilot_stops.txt"), "w") as fh:
            fh.write("\n".join(stops) + "\n" if stops else "none\n")


def judge_all(out, ev_rows, groups):
    rows = [dict(zip(COLS + se.COLUMNS, r)) for r in ev_rows]
    val = {(r["set"], r["which"], r["kind"], r["chrom"], r["pos"], r["ref"], r["alt"]): r for r in rows}
    fields = ("n_d", "n_c", "median_d", "c25", "c75", "width_d", "width_c", "bias_ok", "spread_ok", "pass")
    results = []
    for key in ("forward", "reverse"):
        donor, recipient = t.SETS[key]
        arms = ["normal", "B1", "B2"] if key == "forward" else ["normal"]
        for g in groups:
            shared = [r for r in rows if r["set"] == "shared" and r["group"] == g and r["which"] == "HG001"]
            ev = lambda r: (r["kind"], r["chrom"], r["pos"], r["ref"], r["alt"])  # noqa: E731
            for metric in metrics_of(g):
                c = []
                for r in shared:
                    hi, lo = val[("shared", recipient, *ev(r))][metric], val[("shared", donor, *ev(r))][metric]
                    if hi is not None and lo is not None:
                        c.append(hi - lo)
                for arm in arms:
                    d = []
                    for r in rows:
                        if r["set"] == key and r["group"] == g and r["which"] == f"fake_{arm}":
                            real = val.get((key, "real", *ev(r)))
                            if real and r[metric] is not None and real[metric] is not None:
                                d.append(r[metric] - real[metric])
                    if d and c:
                        j = t.judge(d, c)
                        results.append([key, arm, g, metric] + [j[f] for f in fields])
    t.write_tsv(os.path.join(out, "judge.tsv"), ["set", "arm", "group", "metric", *fields], results)
    judged = {(r[0], r[1], r[2], r[3]): {"pass": r[-1]} for r in results}
    verdicts = []
    for key in ("forward", "reverse"):
        for g in groups:
            statuses = []
            for metric in metrics_of(g):
                res = judged.get((key, "normal", g, metric))
                if res is None:
                    continue
                controls = {arm: judged.get(("forward", arm, g, metric)) for arm in ("B1", "B2")}
                statuses.append(verdict(res, controls, metric))
                verdicts.append([key, g, metric, statuses[-1]])
            verdicts.append([key, g, "overall", t.overall(statuses)])
    t.write_tsv(os.path.join(out, "verdicts.tsv"), ["set", "group", "metric", "verdict"], verdicts)
    return results


if __name__ == "__main__":
    main(sys.argv[1:])
