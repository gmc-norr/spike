#!/usr/bin/env python3
"""RF8 kill test K: how many of the 40 seeded deletions would the refusal stop.

Counts the runs whose SIM_RESIST is above 0.5, the threshold RF8's plan locked.
Usage: rf8_score.py OUT_DIR  (after `cr4_run.sh OUT_DIR ... chr1` with master's
binary; reads OUT_DIR/cr4-c4-events.txt and OUT_DIR/cr4-c4/<n>/{exit,log,run/truth.vcf})
"""
import pathlib
import re
import sys

REFUSE_ABOVE = 0.5

S = pathlib.Path(sys.argv[1])
events = (S / "cr4-c4-events.txt").read_text().split()
rows, failed = [], []
for i, spec in enumerate(events, start=1):
    d = S / "cr4-c4" / str(i)
    code = int((d / "exit").read_text())
    log = (d / "log").read_text()
    if code != 0:
        last = [l for l in log.splitlines() if l.strip()][-1:]
        failed.append((i, spec, last[0][:160] if last else ""))
        continue
    truth = (d / "run" / "truth.vcf").read_text()
    resist = float(re.search(r"SIM_RESIST=([0-9.]+)", truth).group(1))
    census = re.search(r"cannot edit: (\d+) of (\d+)", log)
    rows.append((i, spec, resist, int(census.group(1)), int(census.group(2))))

for i, spec, r, res, cnt in rows:
    print(f"{i:3} {spec:28} R={r:.3f} ({res}/{cnt}) {'ABOVE' if r > REFUSE_ABOVE else ''}")
for i, spec, why in failed:
    print(f"{i:3} {spec:28} EXIT!=0: {why}")
above = sum(1 for row in rows if row[2] > REFUSE_ABOVE)
values = sorted(row[2] for row in rows)
print(f"\nran {len(rows)}, non-zero exit {len(failed)}, above {REFUSE_ABOVE}: {above} of {len(rows)}")
if values:
    print(f"R: min {values[0]:.3f}, median {values[len(values) // 2]:.3f}, max {values[-1]:.3f}")
verdict = (
    "INCONCLUSIVE (more than 4 non-zero exits)" if len(failed) > 4
    else "PASS" if above <= 1 else "FAIL"
)
print(f"K (above {REFUSE_ABOVE} on at most 1): {verdict}")
