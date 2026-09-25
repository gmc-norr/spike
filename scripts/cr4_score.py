#!/usr/bin/env python3
"""CR4 criterion C4: how often the warning fires on the 40 seeded chr20 DELs.

Usage: cr4_score.py OUT_DIR  (after cr4_run.sh; reads OUT_DIR/cr4-c4-events.txt
and OUT_DIR/cr4-c4/<n>/{exit,log,run/truth.vcf})
"""
import pathlib
import re
import sys

S = pathlib.Path(sys.argv[1])
events = (S / "cr4-c4-events.txt").read_text().split()
rows, refused = [], []
for i, spec in enumerate(events, start=1):
    d = S / "cr4-c4" / str(i)
    code = int((d / "exit").read_text())
    log = (d / "log").read_text()
    if code != 0:
        last = [l for l in log.splitlines() if l.strip()][-1:]
        refused.append((i, spec, last[0][:160] if last else ""))
        continue
    truth = (d / "run" / "truth.vcf").read_text()
    resist = float(re.search(r"SIM_RESIST=([0-9.]+)", truth).group(1))
    census = re.search(r"cannot edit: (\d+) of (\d+)", log)
    warned = "WARN" in log and "spike cannot edit" in log
    rows.append((i, spec, resist, int(census.group(1)), int(census.group(2)), warned))

for i, spec, r, res, cnt, w in rows:
    print(f"{i:3} {spec:28} R={r:.3f} ({res}/{cnt}) {'WARN' if w else ''}")
for i, spec, why in refused:
    print(f"{i:3} {spec:28} REFUSED: {why}")
n = len(rows)
warned = sum(1 for row in rows if row[5])
values = sorted(row[2] for row in rows)
print(f"\nran {n}, refused {len(refused)}, warned {warned} of {n}")
if values:
    mid = values[len(values) // 2]
    print(f"R: min {values[0]:.3f}, median {mid:.3f}, max {values[-1]:.3f}")
verdict = (
    "INCONCLUSIVE (more than 4 refusals)" if len(refused) > 4
    else "PASS" if warned <= 8 else "FAIL"
)
print(f"C4 (warned <= 8): {verdict}")
