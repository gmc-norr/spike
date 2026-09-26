#!/usr/bin/env python3
"""RF13 second attempt, C2 and C3: the Rust row against its replica, and the exit status.

For every run under K_DIR that reached validate, runs the NEW spike's
`validate --json` on the run's merged.bam three times -- with the run's own truth,
with N2's wrong-letters truth, and with N3's POS + 1000 truth -- and compares the
`ins_planted` observed count with rf13b_planted.py's on the same input.

C3: the new validate's exit status on the run's own truth, beside the one master's
validate wrote to validate.json.exit when the run was made.

Usage: rf13b_c2.py K_DIR NEW_SPIKE REFERENCE
"""
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rf13b_planted import carriers, ins_records, wrong_letters  # noqa: E402


def rust_count(spike, bam, ref, truth_path):
    out = subprocess.run(
        [spike, "validate", "--bam", bam, "--truth", truth_path, "--reference", ref, "--json"],
        capture_output=True, text=True,
    )
    rows = [c for c in json.loads(out.stdout)["checks"] if c["check"] == "ins_planted"]
    return rows[0]["observed"], out.returncode


def main(kdir, spike, ref):
    mismatches, n, exits = [], 0, []
    tmp = tempfile.mkdtemp()
    for run in sorted(os.listdir(kdir)):
        d = os.path.join(kdir, run)
        truth = os.path.join(d, "run", "truth.vcf")
        merged = os.path.join(d, "run", "merged.bam")
        if not (os.path.exists(truth) and os.path.exists(merged)):
            continue
        header = [l for l in open(truth) if l.startswith("#")]
        f = next(ins_records(truth))
        chrom, pos, rec_id, alt = f[0], int(f[1]), f[2], f[4]
        cases = {
            "own": (pos, alt),
            "N2": (pos, wrong_letters(alt, pos)),
            "N3": (pos + 1000, alt),
        }
        for label, (p, a) in cases.items():
            g = f.copy()
            g[1], g[4] = str(p), a
            path = os.path.join(tmp, f"{run}_{label}.vcf")
            with open(path, "w") as out:
                out.writelines(header)
                out.write("\t".join(g) + "\n")
            observed, rc = rust_count(spike, merged, ref, path)
            c = carriers(merged, ref, chrom, p, a, rec_id)
            replica = c if isinstance(c, str) else str(len(c))
            n += 1
            if observed != replica:
                mismatches.append((run, label, observed, replica))
            if label == "own":
                master_rc = open(os.path.join(d, "validate.json.exit")).read().strip()
                exits.append((run, master_rc, str(rc)))
    print(f"C2 Rust ins_planted == replica: {n - len(mismatches)} of {n}  {mismatches}")
    moved = [e for e in exits if e[1] != e[2]]
    print(f"C3 exit status, master -> new, on {len(exits)} runs: "
          f"{sum(1 for e in exits if e[1] == '0')} -> {sum(1 for e in exits if e[2] == '0')} exit 0")
    for run, old, new in moved:
        print(f"  {run}: {old} -> {new}")


if __name__ == "__main__":
    main(*sys.argv[1:])
