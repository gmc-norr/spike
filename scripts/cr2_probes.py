#!/usr/bin/env python3
"""CR2 criteria C1-C3: the variable and uniform DUP probes, before vs depth-fold binary.

Usage: cr2_probes.py OUT_DIR BEFORE_SPIKE NEW_SPIKE  (needs samtools on PATH)

BEFORE_SPIKE is the CR4 census binary (2b5b193), NEW_SPIKE the depth-fold one.
"""
import hashlib
import pathlib
import re
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from review_sv_model import Probe  # noqa: E402

S = pathlib.Path(sys.argv[1])
BINARIES = {"before": pathlib.Path(sys.argv[2]), "new": pathlib.Path(sys.argv[3])}
EVENT = "dup:chrT:10000-28000;af=0.5"


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


results = {}
for tag, binary in BINARIES.items():
    out = S / f"cr2-probe-{tag}"
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    probe = Probe(binary, out)
    for bam_label, kwargs in (("uniform", {}), ("variable", {"variable": True})):
        bam, _ = probe.make_bam(bam_label, **kwargs)
        dest = probe.run(f"dup_{bam_label}", bam, [EVENT])
        log = (out / f"dup_{bam_label}.log").read_text()
        truth = (dest / "truth.vcf").read_text()
        record = [l for l in truth.splitlines() if not l.startswith("#")][0]
        fold = re.search(r"SIM_DEPTH_FOLD=([^;\t]+)", record)
        warning = [l for l in log.splitlines() if "WARN" in l and "SIM_DEPTH_FOLD" in l]
        results[(tag, bam_label)] = {
            "R1": md5(dest / "R1.fq.gz"),
            "R2": md5(dest / "R2.fq.gz"),
            "replaced": md5(dest / "replaced_reads.txt"),
            "fold": fold.group(1) if fold else None,
            "warning": warning[0] if warning else None,
            "truth_lines": truth.splitlines(),
        }

for bam_label in ("variable", "uniform"):
    b, n = results[("before", bam_label)], results[("new", bam_label)]
    print(f"== {bam_label}")
    print(f"  SIM_DEPTH_FOLD before={b['fold']} new={n['fold']}")
    print(f"  warning: {n['warning']}")
    for key in ("R1", "R2", "replaced"):
        print(f"  {key}: before {b[key]}  new {n[key]}  same={b[key] == n[key]}")
    only_b = [l for l in b["truth_lines"] if l not in n["truth_lines"]]
    only_n = [l for l in n["truth_lines"] if l not in b["truth_lines"]]
    print(f"  truth lines only before: {only_b}")
    print(f"  truth lines only new:    {only_n}")
