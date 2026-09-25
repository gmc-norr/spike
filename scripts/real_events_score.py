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
import sys


def main(out_dir, checks):
    out = pathlib.Path(out_dir)
    events = (out / "events.txt").read_text().splitlines()
    rows = []
    print("n\tevent\tspike_exit\tvalidate_exit\tchecks")
    for n, event in enumerate(events, start=1):
        d = out / str(n)
        spike_exit = (d / "spike.exit").read_text().strip() if (d / "spike.exit").exists() else "NA"
        val_exit = (d / "validate.exit").read_text().strip() if (d / "validate.exit").exists() else "NA"
        report = d / "validate.json"
        if not report.exists():
            print(f"{n}\t{event}\t{spike_exit}\t{val_exit}\tNO_REPORT")
            rows.append((n, event, spike_exit, val_exit, None))
            continue
        try:
            data = json.loads(report.read_text())
        except Exception as error:               # a truncated report is not a pass
            print(f"{n}\t{event}\t{spike_exit}\t{val_exit}\tUNPARSEABLE: {error}")
            rows.append((n, event, spike_exit, val_exit, None))
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
        rows.append((n, event, spike_exit, val_exit, by_check))

    scored = [r for r in rows if r[4] is not None and r[2] == "0"]
    refused = [r for r in rows if r[2] != "0"]
    unscored = [r for r in rows if r[4] is None and r[2] == "0"]
    print()
    print(f"events={len(rows)} scored={len(scored)} refused_by_spike={len(refused)} "
          f"no_report={len(unscored)}")
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
            values.sort()
            mid = values[len(values) // 2]
            print(f"  {check} observed: min {values[0]:.2f} median {mid:.2f} "
                  f"max {values[-1]:.2f} (n={len(values)})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
