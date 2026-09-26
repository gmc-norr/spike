#!/usr/bin/env python3
"""RF8 criteria C1b and C3's probe half: the lowmap and uniform probes.

Runs CR4's probe event on each BAM three ways: master's binary, the new binary,
and the new binary with --allow-resistant. Prints each run's exit, SIM_RESIST
and output md5s, so the refusal and the byte-identity can be read off.

Usage: rf8_probes.py OUT_DIR MASTER_SPIKE NEW_SPIKE  (needs samtools on PATH)
"""
import hashlib
import pathlib
import re
import shutil
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from review_sv_model import Probe  # noqa: E402

S = pathlib.Path(sys.argv[1])
MASTER, NEW = pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
EVENT = "del:chrT:10000-14000;af=1"
RUNS = (("master", MASTER, []), ("new", NEW, []), ("new+allow", NEW, ["--allow-resistant"]))


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest() if path.exists() else "absent"


shutil.rmtree(S, ignore_errors=True)
S.mkdir(parents=True)
probe = Probe(MASTER, S)
for bam_label, kwargs in (("uniform", {}), ("lowmap", {"lowmap": True})):
    bam, _ = probe.make_bam(bam_label, **kwargs)
    print(f"== {bam_label}")
    for tag, binary, extra in RUNS:
        dest = S / f"{bam_label}-{tag}"
        command = [
            str(binary.resolve()), "--bam", str(bam), "--reference", str(probe.fasta),
            "--flank", "2000", "--seed", "17", "-o", str(dest), "--event", EVENT, *extra,
        ]
        done = subprocess.run(command, capture_output=True, text=True)
        (S / f"{bam_label}-{tag}.log").write_text(done.stdout + done.stderr)
        truth = dest / "truth.vcf"
        resist = None
        if truth.exists():
            found = re.search(r"SIM_RESIST=([^;\t]+)", truth.read_text())
            resist = found.group(1) if found else None
        # ##fileDate is the day's date, not the run's content.
        body = (
            hashlib.md5("".join(
                l for l in truth.read_text().splitlines(keepends=True)
                if not l.startswith("##fileDate")
            ).encode()).hexdigest() if truth.exists() else "absent"
        )
        print(
            f"  {tag:10} exit={done.returncode} SIM_RESIST={resist} "
            f"R1={md5(dest / 'R1.fq.gz')} R2={md5(dest / 'R2.fq.gz')} "
            f"replaced={md5(dest / 'replaced_reads.txt')} truth={body}"
        )
        if done.returncode != 0:
            last = [l for l in done.stderr.splitlines() if l.strip()][-1:]
            print(f"    stderr: {last[0][:200] if last else ''}")
