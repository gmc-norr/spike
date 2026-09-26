#!/usr/bin/env python3
"""One TSV line per event from a real_events.sh output directory.

Scores `spike validate`'s `--json` report, not its text table: the table is
space-padded, so a parser has to guess where the event label ends.

Usage: real_events_score.py OUT_DIR [CHECK ...]
With no CHECK, every check is reported. With one or more, the summary at the end
counts how many events had each named check FAIL.
"""
import json
import pathlib
import statistics
import sys


def harness_failure(run_dir, spike_exit):
    """Why this event was not scored, if the reason is the harness's own.

    `real_events.sh`'s slice step writes its exit status to `slice.exit` and,
    when it fails, no `spike.exit` at all -- a wrong BAM path or a contig name
    the BAM does not carry is this harness failing, not spike declining the
    event. Older run directories wrote the string `slice failed` into
    `spike.exit` itself, so a non-numeric `spike.exit` is read the same way.
    Returns a reason string, or None when nothing says the harness failed.
    """
    slice_exit = run_dir / "slice.exit"
    if slice_exit.exists():
        text = slice_exit.read_text().strip()
        if text != "0":
            return f"slice failed (slice.exit={text})"
    if spike_exit == "NA":
        return "no spike.exit: spike never ran"
    if not spike_exit.lstrip("-").isdigit():
        return f"spike.exit is not a number: {spike_exit!r}"
    return None


def main(out_dir, checks):
    out = pathlib.Path(out_dir)
    events = (out / "events.txt").read_text().splitlines()
    # A blank line here would number every later event one lower than
    # real_events.sh numbered its directories, and nothing would look wrong.
    # Refuse the file instead; real_events.sh refuses it too.
    blank = [i for i, line in enumerate(events, start=1) if not line.strip()]
    if blank:
        raise SystemExit(
            f"{out / 'events.txt'}: blank line(s) at {blank}. Every later event "
            f"would be attributed to the wrong run directory. Remove them and "
            f"re-run real_events.sh."
        )
    rows = []
    print("n\tevent\tspike_exit\tvalidate_exit\tchecks")
    for n, event in enumerate(events, start=1):
        d = out / str(n)
        spike_exit = (d / "spike.exit").read_text().strip() if (d / "spike.exit").exists() else "NA"
        val_exit = (d / "validate.exit").read_text().strip() if (d / "validate.exit").exists() else "NA"
        broken = harness_failure(d, spike_exit)
        report = d / "validate.json"
        if broken is not None:
            print(f"{n}\t{event}\t{spike_exit}\t{val_exit}\tHARNESS_FAILED: {broken}")
            rows.append((n, event, spike_exit, val_exit, None, broken))
            continue
        if not report.exists():
            print(f"{n}\t{event}\t{spike_exit}\t{val_exit}\tNO_REPORT")
            rows.append((n, event, spike_exit, val_exit, None, None))
            continue
        try:
            data = json.loads(report.read_text())
        except Exception as error:               # a truncated report is not a pass
            print(f"{n}\t{event}\t{spike_exit}\t{val_exit}\tUNPARSEABLE: {error}")
            rows.append((n, event, spike_exit, val_exit, None, None))
            continue
        fields = []
        by_check = {}
        for check in data.get("checks", []):
            name = check["check"]
            status = "PASS" if check["pass"] else "FAIL"
            advisory = ":ADV" if check.get("advisory") else ""
            fields.append(f'{name}={check["observed"]}:{status}{advisory}')
            by_check.setdefault(name, []).append(check)
        print(f"{n}\t{event}\t{spike_exit}\t{val_exit}\t{';'.join(fields)}")
        rows.append((n, event, spike_exit, val_exit, by_check, None))

    # Three reasons an event is not scored, kept apart: the harness failed
    # before spike ran, spike refused the event, or spike ran and left no
    # readable report. Filing the first under the second reported a broken
    # harness as "spike declined this event".
    broken = [r for r in rows if r[5] is not None]
    rest = [r for r in rows if r[5] is None]
    scored = [r for r in rest if r[4] is not None and r[2] == "0"]
    refused = [r for r in rest if r[2] != "0"]
    unscored = [r for r in rest if r[4] is None and r[2] == "0"]
    print()
    print(f"events={len(rows)} scored={len(scored)} refused_by_spike={len(refused)} "
          f"no_report={len(unscored)} harness_failed={len(broken)}")
    for r in broken:
        print(f"  HARNESS FAILED: {r[0]} {r[1]} {r[5]}")
    for r in refused:
        print(f"  refused: {r[0]} {r[1]} spike_exit={r[2]}")
    for r in unscored:
        print(f"  no report: {r[0]} {r[1]}")
    for check in checks:
        present = [r for r in scored if check in r[4]]
        failed = [r for r in present if any(not c["pass"] for c in r[4][check])]
        print(f"{check}: present on {len(present)} of {len(scored)} scored, "
              f"FAIL on {len(failed)}")
        for r in failed:
            obs = ",".join(c["observed"] for c in r[4][check])
            print(f"  FAIL: {r[0]} {r[1]} observed={obs}")
        values = []
        for r in present:
            for c in r[4][check]:
                try:
                    values.append(float(c["observed"]))
                except ValueError:
                    pass
        if values:
            # A true median: `values[len(values) // 2]` is the upper of the two
            # middle values on an even n, which is not the median and was
            # printed under that name.
            values.sort()
            mid = statistics.median(values)
            print(f"  {check} observed: min {values[0]:.2f} median {mid:.2f} "
                  f"max {values[-1]:.2f} (n={len(values)})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
