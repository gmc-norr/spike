#!/usr/bin/env python3
"""RF6 criteria C1 and C3: master's and the new validate on the same merged BAMs.

Reads a real_events.sh output directory run with BEFORE_SPIKE=<master's binary>:
per event, validate.before.exit / validate.before.json (master) and
validate.json.exit / validate.json (new). Prints one line per event and the counts
the plan's criteria are judged on.

Usage: rf6_score.py OUT_DIR
"""
import json
import pathlib
import sys

OUT = pathlib.Path(sys.argv[1])
events = (OUT / "events.txt").read_text().splitlines()


def read(path):
    return path.read_text().strip() if path.exists() else "NA"


def row(report, check):
    for r in report.get("checks", []):
        if r["check"] == check:
            flag = "PASS" if r["pass"] else "FAIL"
            return f"{r['observed']} {flag}{' (advisory)' if r['advisory'] else ''}"
    return "absent"


scored, refused, now_fail, now_pass = 0, [], [], []
for n, spec in enumerate(events, start=1):
    d = OUT / str(n)
    spike_exit = read(d / "spike.exit")
    if spike_exit != "0":
        refused.append((n, spec, spike_exit))
        print(f"{n:3} {spec:30} spike exit {spike_exit}: not scored")
        continue
    before, after = read(d / "validate.before.exit"), read(d / "validate.json.exit")
    rb = json.loads((d / "validate.before.json").read_text())
    ra = json.loads((d / "validate.json").read_text())
    scored += 1
    if before == "0" and after != "0":
        now_fail.append((n, spec))
    if before != "0" and after == "0":
        now_pass.append((n, spec))
    print(
        f"{n:3} {spec:30} exit master={before} new={after}  "
        f"split_reads master=[{row(rb, 'split_reads')}] new=[{row(ra, 'split_reads')}]  "
        f"junction_sequence=[{row(ra, 'junction_sequence')}]"
    )

print(f"\nscored {scored}, not scored {len(refused)}")
print(f"master exit 0 -> new exit != 0 (must be none): {len(now_fail)} {now_fail}")
print(f"master exit != 0 -> new exit 0: {len(now_pass)} {now_pass}")
