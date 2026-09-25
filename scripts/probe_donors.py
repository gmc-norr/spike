#!/usr/bin/env python3
"""Build the Codex review's three donor BAMs and a bwa-mem2-indexed reference.

`probe_loop.sh` needs a merged BAM, so the probe reference has to be indexed for
bwa-mem2 -- `review_sv_model.Probe` only writes and faidxes it. This makes both.

Usage: probe_donors.py OUT_DIR SPIKE
Prints `ref=<path>` and one `<label>=<bam path>` line per donor.
"""
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from review_sv_model import Probe  # noqa: E402


def main(out_dir, spike):
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    probe = Probe(pathlib.Path(spike), out)
    subprocess.run(["bwa-mem2", "index", str(probe.fasta)],
                   check=True, capture_output=True)
    print(f"ref={probe.fasta}")
    for label, kwargs in (("uniform", {}), ("lowmap", {"lowmap": True}),
                          ("variable", {"variable": True}),
                          ("homdel", {"homdel": True})):
        bam, _ = probe.make_bam(label, **kwargs)
        print(f"{label}={bam}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
