#!/usr/bin/env python3
"""Transplant test, round 3: fake-vs-real mismatch against split-half counting noise at each site.

The rule and its stop checks are locked in
docs/superpowers/plans/2026-09-28-transplant-round3.md; this script only
carries them out. It runs no spike: it re-reads round 2b's reads (round2.py's
--out directory) with small_evidence.py's metrics, on all reads and on random
halves of the read pairs.

Usage: round3.py --full DIR --out DIR [--processes N] [--salts K] [--draws B]
"""
import argparse
import hashlib
import multiprocessing
import os
import random
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import round2 as r2  # noqa: E402
import small_evidence as se  # noqa: E402
import transplant as t  # noqa: E402

LIMIT = 2.25  # 1.5^2: rounds 1b and 2's spread allowance, on a variance scale
CALIBRATION = (0.8, 1.25)
BOOT_SEED = 3
SALTS = 10
DRAWS = 2_000
ARMS = {"forward": ["normal", "B1", "B2"], "reverse": ["normal"]}


# --- the logic the tests pin down ---------------------------------------------

def half_of(name, salt):
    """The half (0 or 1) a read pair goes to under `salt`; mates share a name."""
    return hashlib.blake2b(f"{salt}:{name}".encode()).digest()[0] & 1


class OtherHalf:
    """Names to skip: those in `skip`, and every pair outside half `h` of `salt`."""

    def __init__(self, salt, h, skip=frozenset()):
        self.salt, self.h, self.skip = salt, h, skip

    def __contains__(self, name):
        return name in self.skip or half_of(name, self.salt) != self.h


def paths_for(bam, skip=frozenset(), extras=(), half=None):
    """(BAM, names to skip) pairs as small_evidence reads them; `half` = (salt, h) keeps that half only."""
    if half is None:
        return [(bam, frozenset(skip))] + [(p, frozenset()) for p in extras]
    salt, h = half
    return [(bam, OtherHalf(salt, h, skip))] + [(p, OtherHalf(salt, h)) for p in extras]


def measure_on(paths, fetch, event):
    """small_evidence.measure, on the given (BAM, names to skip) pairs."""
    kind, chrom, pos, ref, alt = event
    row = dict.fromkeys(se.COLUMNS)
    if kind == "SNV":
        row.update(se.measure_snv(paths, chrom, pos, alt))
    elif kind in ("DEL", "INS"):
        row.update(se.measure_indel(paths, fetch, kind, chrom, pos, ref, alt))
    elif kind == "DUP":
        row.update(se.measure_dup(paths, fetch, chrom, pos, ref, alt))
    else:
        raise ValueError(f"unknown kind {kind!r}")
    return row


def split_noise(pairs):
    """Mean of (x0 - x1)^2 / 4 over the salts where both halves are defined; None if none is."""
    gaps = [(x0 - x1) ** 2 / 4 for x0, x1 in pairs if x0 is not None and x1 is not None]
    return sum(gaps) / len(gaps) if gaps else None


def mismatch_ratio(items):
    """Sum of squared differences over summed noise; items are (difference, noise)."""
    noise = sum(v for _, v in items)
    return sum(x * x for x, _ in items) / noise if noise > 0 else None


def items_of(fake, real, metric):
    """(fake - real, v_fake + v_real) per event both have, with value and noise defined.
    fake, real: {event: {metric: (value, v)}}."""
    out = []
    for e, row in fake.items():
        if e not in real:
            continue
        (x, vx), (y, vy) = row[metric], real[e][metric]
        if None not in (x, vx, y, vy):
            out.append((x - y, vx + vy))
    return out


def rho_interval(f_items, c_items, draws, seed):
    """5th and 95th percentiles of R_f / R_c over bootstrap draws, each side resampled apart."""
    rng = random.Random(seed)
    rhos = []
    for _ in range(draws):
        rf = mismatch_ratio([rng.choice(f_items) for _ in f_items])
        rc = mismatch_ratio([rng.choice(c_items) for _ in c_items])
        if rf is not None and rc:
            rhos.append(rf / rc)
    return interval(rhos)


def interval(values):
    """The 5th and 95th percentiles, or (None, None) with no values."""
    if not values:
        return (None, None)
    return (t.percentile(values, 5), t.percentile(values, 95))


def judge3(lo, hi):
    if lo is None or hi is None:
        return "unsure"
    if hi <= LIMIT:
        return "pass"
    if lo > LIMIT:
        return "fail"
    return "unsure"


def verdict3(status, controls, metric):
    """The metric's status, or inconclusive unless every control that must break it failed."""
    if any(controls.get(arm) != "fail" for arm in r2.MUST_FAIL[metric]):
        return "inconclusive"
    return status


def calibration(rows):
    """Split-half noise over binomial noise A(1-A)/N, summed over SNV rows {A, N, v}."""
    usable = [r for r in rows if r["A"] is not None and r["v"] is not None and r["N"] > 0]
    return sum(r["v"] for r in usable) / sum(r["A"] * (1 - r["A"]) / r["N"] for r in usable)


def calibrated(ratio):
    return CALIBRATION[0] <= ratio <= CALIBRATION[1]


def as_written(v):
    """A value as round 2b's evidence.tsv writes it (transplant.write_tsv)."""
    return "" if v is None else (f"{v:.4f}" if isinstance(v, float) else str(v))


def same_as_written(v, cell):
    return as_written(v) == cell


def runs_of(events, span):
    """event -> the run number round2.py put it in (transplant.batches, in the given order)."""
    spans = [(*span(e), i) for i, e in enumerate(events)]
    return {events[x[3]]: n for n, run in enumerate(t.batches(spans)) for x in run}


def in_run_bed(event, bed):
    """Whether the event's POS lies inside one of its run's events.bed intervals (chrom, start, end)."""
    _, chrom, pos, _, _ = event
    return any(c == chrom and s <= pos - 1 < e for c, s, e in bed)


# --- measuring (one worker process) -------------------------------------------

_fasta = None
_skip = {}


def _fetch(chrom, start, end):
    global _fasta
    if _fasta is None:
        _fasta = pysam.FastaFile(t.REF)
    return _fasta.fetch(chrom, max(0, start), max(0, end)).upper()


def _fresh():
    """Pool initializer: a worker opens its own FASTA, never the parent's handle."""
    global _fasta
    _fasta = None
    _skip.clear()


def _skip_of(run_dir):
    if run_dir not in _skip:
        _skip.clear()  # one run's names at a time
        _skip[run_dir] = frozenset(t.skip_set(run_dir))
    return _skip[run_dir]


def measure_task(task):
    """(set, which, event, bam, run_dir or None, salts) -> (set, which, event, full row, {metric: (value, v, salts used)})."""
    key, which, event, bam, run_dir, salts = task
    skip, extras = (_skip_of(run_dir), [os.path.join(run_dir, "sim.bam")]) if run_dir else (frozenset(), [])
    full = measure_on(paths_for(bam, skip, extras), _fetch, event)
    halves = [[measure_on(paths_for(bam, skip, extras, half=(k, h)), _fetch, event) for h in (0, 1)]
              for k in range(salts)]
    noise = {}
    for m in r2.metrics_of(r2.group_of_event(event)):
        pairs = [(a[m], b[m]) for a, b in halves]
        used = sum(1 for x0, x1 in pairs if x0 is not None and x1 is not None)
        noise[m] = (full[m], split_noise(pairs), used)
    return key, which, event, full, noise


# --- the whole thing ----------------------------------------------------------

def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def event_of(row):
    return (row["kind"], row["chrom"], int(row["pos"]), row["ref"], row["alt"])


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True, help="round2.py's --out directory of round 2b's full run")
    ap.add_argument("--out", required=True)
    ap.add_argument("--processes", type=int, default=16)
    ap.add_argument("--salts", type=int, default=SALTS)
    ap.add_argument("--draws", type=int, default=DRAWS)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    chosen = {}
    for r in read_tsv(os.path.join(a.full, "events.tsv")):
        chosen.setdefault((r["set"], r["group"]), []).append(event_of(r))
    written = {(r["set"], r["which"], event_of(r)): r for r in read_tsv(os.path.join(a.full, "evidence.tsv"))}
    refused = {(r["set"], r["arm"], event_of(r)) for r in read_tsv(os.path.join(a.full, "refused.tsv"))}

    def span(e):
        if e[0] == "DUP":
            return (e[1], *se.segment_of(_fetch, e[1], e[2], e[3], e[4]))
        return (e[1], e[2] - 1, e[2] - 1 + len(e[3]))

    # Check 3, while the tasks are built: each event lies in the run it is said to be in.
    tasks, misplaced = [], []
    for key in ("forward", "reverse"):
        donor, recipient = t.SETS[key]
        events = [e for g in r2.GROUPS for e in chosen.get((key, g), [])]
        tasks += [(key, "real", e, t.BAMS[donor], None, a.salts) for e in events]
        for arm in ARMS[key]:
            arm_events = [e for e in events if arm != "B2" or e[0] == "DUP"]
            runs = runs_of(arm_events, span)
            beds = {}
            for e in arm_events:
                if (key, arm, e) in refused:
                    continue
                run_dir = os.path.join(a.full, key, arm, f"run{runs[e]}")
                if run_dir not in beds:
                    with open(os.path.join(run_dir, "events.bed")) as fh:
                        beds[run_dir] = [(c, int(s), int(en)) for c, s, en in (ln.split("\t")[:3] for ln in fh)]
                if not in_run_bed(e, beds[run_dir]):
                    misplaced.append((key, arm, e, run_dir))
                tasks.append((key, f"fake_{arm}", e, t.BAMS[recipient], run_dir, a.salts))
    for g in r2.GROUPS:
        for e in chosen.get(("shared", g), []):
            tasks += [("shared", s, e, t.BAMS[s], None, a.salts) for s in ("HG001", "HG002")]
    if misplaced:
        stop(a.out, f"check 3 failed: {len(misplaced)} events not inside their run's events.bed, first {misplaced[0]}")
    print(f"round3: {len(tasks)} measurements, {a.salts} salts, {a.processes} processes", flush=True)

    tasks.sort(key=lambda x: (x[4] or "", x[0], x[1]))  # a worker keeps one run's names at a time
    with multiprocessing.Pool(a.processes, initializer=_fresh) as pool:
        results = pool.map(measure_task, tasks, chunksize=8)

    # Check 1: every all-reads value is round 2b's, as written.
    mismatches = []
    for key, which, e, full, _ in results:
        row = written.get((key, which, e))
        if row is None:
            mismatches.append((key, which, e, "no row in evidence.tsv"))
            continue
        for c in se.COLUMNS:
            if not same_as_written(full[c], row[c]):
                mismatches.append((key, which, e, c, as_written(full[c]), row[c]))
    rows = [[key, r2.group_of_event(e), which, *e, m, value, v, used]
            for key, which, e, _, noise in results for m, (value, v, used) in noise.items()]
    t.write_tsv(os.path.join(a.out, "evidence3.tsv"),
                ["set", "group", "which", "kind", "chrom", "pos", "ref", "alt", "metric", "value", "v", "salts_used"], rows)
    if mismatches:
        stop(a.out, f"check 1 failed: {len(mismatches)} values differ from round 2b, first {mismatches[0]}")

    # Check 2: split-half noise against binomial noise on real SNVs.
    snv = [{"A": noise["A"][0], "N": full["n_carry"] + full["n_ref"], "v": noise["A"][1]}
           for key, which, e, full, noise in results if e[0] == "SNV" and not which.startswith("fake")]
    cal = calibration(snv)
    with open(os.path.join(a.out, "checks.txt"), "w") as fh:
        fh.write(f"check 1: {len(results)} measurements, 0 differ from round 2b\n")
        fh.write(f"check 2: split-half / binomial noise on {len(snv)} real SNV measurements = {cal:.4f} "
                 f"(must lie within {CALIBRATION[0]}-{CALIBRATION[1]})\n")
        fh.write("check 3: 0 events outside their run's events.bed\n")
    if not calibrated(cal):
        stop(a.out, f"check 2 failed: split-half / binomial noise on SNVs is {cal:.4f}")

    judge_all(a.out, results, a.draws)
    print("round3: done", flush=True)


def stop(out, why):
    with open(os.path.join(out, "STOPPED.txt"), "w") as fh:
        fh.write(why + "\n")
    print(f"round3: STOPPED: {why}", flush=True)
    sys.exit(2)


def judge_all(out, results, draws):
    by = {}  # (set, which, group) -> {event: {metric: (value, v)}}
    for key, which, e, _, noise in results:
        by.setdefault((key, which, r2.group_of_event(e)), {})[e] = {m: (x, v) for m, (x, v, _) in noise.items()}
    judged, table = {}, []
    for key in ("forward", "reverse"):
        donor, recipient = t.SETS[key]
        for g in r2.GROUPS:
            shared_rec, shared_don = by.get(("shared", recipient, g), {}), by.get(("shared", donor, g), {})
            real = by.get((key, "real", g), {})
            for m in r2.metrics_of(g):
                c = items_of(shared_rec, shared_don, m)
                for arm in ARMS[key]:
                    fake = by.get((key, f"fake_{arm}", g), {})
                    if not fake:
                        continue
                    f = items_of(fake, real, m)
                    rf, rc = mismatch_ratio(f), mismatch_ratio(c)
                    lo, hi = rho_interval(f, c, draws, BOOT_SEED) if f and c else (None, None)
                    status = judge3(lo, hi)
                    judged[(key, arm, g, m)] = status
                    table.append([key, arm, g, m, len(f), len(fake) - len(f), len(c), len(shared_rec) - len(c),
                                  rf, rc, rf / rc if rf is not None and rc else None, lo, hi, status])
    t.write_tsv(os.path.join(out, "judge3.tsv"),
                ["set", "arm", "group", "metric", "n_f", "left_out_f", "n_c", "left_out_c",
                 "R_f", "R_c", "rho", "rho5", "rho95", "status"], table)
    verdicts = []
    for key in ("forward", "reverse"):
        for g in r2.GROUPS:
            statuses = []
            for m in r2.metrics_of(g):
                status = judged.get((key, "normal", g, m))
                if status is None:
                    continue
                controls = {arm: judged.get(("forward", arm, g, m)) for arm in ("B1", "B2")}
                statuses.append(verdict3(status, controls, m))
                verdicts.append([key, g, m, statuses[-1]])
            verdicts.append([key, g, "overall", t.overall(statuses)])
    t.write_tsv(os.path.join(out, "verdicts3.tsv"), ["set", "group", "metric", "verdict"], verdicts)


if __name__ == "__main__":
    main(sys.argv[1:])
