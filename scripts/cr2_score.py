#!/usr/bin/env python3
"""CR2 criterion C4: how often the depth-fold warning fires on the 40 chr20 DUPs.

Usage: cr2_score.py OUT_DIR  (after cr2_run.sh; reads OUT_DIR/cr2-c4-events.txt
and OUT_DIR/cr2-c4/<n>/{exit,log,run/truth.vcf})
"""
import pathlib
import re
import sys

S = pathlib.Path(sys.argv[1])
events = (S / "cr2-c4-events.txt").read_text().split()
rows, refused = [], []
for i, spec in enumerate(events, start=1):
    d = S / "cr2-c4" / str(i)
    code = int((d / "exit").read_text())
    log = (d / "log").read_text()
    if code != 0:
        last = [l for l in log.splitlines() if l.strip()][-1:]
        refused.append((i, spec, last[0][:160] if last else ""))
        continue
    truth = (d / "run" / "truth.vcf").read_text()
    fold = float(re.search(r"SIM_DEPTH_FOLD=([0-9.]+)", truth).group(1))
    warn = [l for l in log.splitlines() if "WARN" in l and "SIM_DEPTH_FOLD" in l]
    where = re.search(r"depth over (\S+) is ([0-9.]+x), but .* by the ([0-9.]+x)", warn[0]) if warn else None
    rows.append((i, spec, fold, bool(warn), where.groups() if where else None))

for i, spec, fold, warned, where in rows:
    extra = f"  WARN {where[0]} {where[1]} vs {where[2]}" if warned and where else ""
    print(f"{i:3} {spec:28} fold={fold:.2f}{extra}")
for i, spec, why in refused:
    print(f"{i:3} {spec:28} REFUSED: {why}")
n = len(rows)
warned = sum(1 for row in rows if row[3])
values = sorted(row[2] for row in rows)
print(f"\nran {n}, refused {len(refused)}, warned {warned} of {n}")
if values:
    print(f"fold: min {values[0]:.2f}, median {values[len(values) // 2]:.2f}, max {values[-1]:.2f}")
verdict = (
    "INCONCLUSIVE (more than 4 refusals)" if len(refused) > 4
    else "PASS" if warned <= 8 else "FAIL"
)
print(f"C4 (warned <= 8): {verdict}")
