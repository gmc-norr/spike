#!/usr/bin/env python3
"""RF14, C2 and C3: the Rust row against its replica, and the exit status.

For every run directory under each RUNS_DIR that reached validate (run/truth.vcf and
run/merged.bam present):

C2  runs the NEW spike's `validate --json` on the run's merged.bam four times -- the
    run's own truth, N2a (END + 50), N2b (START - 50) and N3 (both + 1000) -- and
    compares the `del_planted` observed count with rf14_planted.py's on the same input.
C3  runs MASTER's and the NEW spike's `validate` on the run's own truth and compares
    their exit statuses. Pass: no run where master exits 0 and the new one does not.

Usage: rf14_c2.py NEW_SPIKE MASTER_SPIKE REFERENCE RUNS_DIR...
"""
import concurrent.futures
import json
import os
import re
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rf14_planted as row  # noqa: E402


def validate(spike, bam, ref, truth):
    out = subprocess.run(
        [spike, "validate", "--bam", bam, "--truth", truth, "--reference", ref, "--json"],
        capture_output=True, text=True,
    )
    try:
        checks = json.loads(out.stdout)["checks"]
    except ValueError:
        raise SystemExit(f"no JSON from {spike} on {bam} {truth}: {out.stderr[-500:]}")
    return checks, out.returncode


def shifted(truth_lines, start, end):
    """The truth VCF with its one DEL record moved to (start, end)."""
    out = []
    for line in truth_lines:
        if line.startswith("#"):
            out.append(line)
            continue
        f = line.rstrip("\n").split("\t")
        f[1] = str(start)
        f[7] = re.sub(r"(^|;)END=\d+", rf"\g<1>END={end}", f[7])
        f[7] = re.sub(r"(^|;)SVLEN=-?\d+", rf"\g<1>SVLEN=-{end - start}", f[7])
        out.append("\t".join(f) + "\n")
    return out


def one_run(args):
    d, new, master, ref, tmp = args
    truth = os.path.join(d, "run", "truth.vcf")
    merged = os.path.join(d, "run", "merged.bam")
    lines = open(truth).readlines()
    chrom, start, end, rec_id, _ = next(row.del_records(truth))
    cases = {
        "own": (start, end),
        "N2a": (start, end + 50),
        "N2b": (start - 50, end),
        "N3": (start + 1000, end + 1000),
    }
    tag = d.strip("/").replace("/", "_")
    c2, c3 = [], None
    for label, (s, e) in cases.items():
        path = os.path.join(tmp, f"{tag}_{label}.vcf")
        with open(path, "w") as fh:
            fh.writelines(shifted(lines, s, e))
        checks, rc = validate(new, merged, ref, path)
        rust = [c["observed"] for c in checks if c["check"] == "del_planted"]
        replica = str(row.count(row.carriers(merged, ref, chrom, s, e, rec_id)))
        c2.append((label, rust[0] if rust else "no row", replica))
        if label == "own":
            _, master_rc = validate(master, merged, ref, path)
            c3 = (master_rc, rc)
    return d, c2, c3


def main(new, master, ref, *dirs):
    runs = []
    for top in dirs:
        for name in sorted(os.listdir(top)):
            d = os.path.join(top, name)
            if os.path.exists(os.path.join(d, "run", "truth.vcf")) and \
                    os.path.exists(os.path.join(d, "run", "merged.bam")):
                runs.append(d)
    tmp = tempfile.mkdtemp()
    n2, mismatches, exits = 0, [], []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        for d, c2, c3 in pool.map(one_run, [(d, new, master, ref, tmp) for d in runs]):
            for label, rust, replica in c2:
                n2 += 1
                if rust != replica:
                    mismatches.append((d, label, rust, replica))
            exits.append((d, c3[0], c3[1]))
    print(f"C2 Rust del_planted == replica: {n2 - len(mismatches)} of {n2}  {mismatches}")
    broke = [(d, m, n) for d, m, n in exits if m == 0 and n != 0]
    fixed = [(d, m, n) for d, m, n in exits if m != 0 and n == 0]
    print(f"C3 on {len(exits)} runs: master exit 0 on {sum(m == 0 for _, m, _ in exits)}, "
          f"new exit 0 on {sum(n == 0 for _, _, n in exits)}")
    print(f"C3 master 0 -> new non-zero: {len(broke)} {broke} -> {'PASS' if not broke else 'FAIL'}")
    print(f"C3 master non-zero -> new 0: {len(fixed)}")
    for top in dirs:
        sub = [(m, n) for d, m, n in exits if d.startswith(top)]
        print(f"  {top}: {len(sub)} runs, master 0 on {sum(m == 0 for m, _ in sub)}, "
              f"new 0 on {sum(n == 0 for _, n in sub)}")


if __name__ == "__main__":
    main(*sys.argv[1:])
