#!/usr/bin/env python3
"""CR4 criteria C1-C3: the lowmap and uniform probes, master vs census binary.

Usage: cr4_probes.py OUT_DIR MASTER_SPIKE NEW_SPIKE  (needs samtools on PATH)
"""
import hashlib
import pathlib
import re
import shutil
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from review_sv_model import Probe  # noqa: E402

S = pathlib.Path(sys.argv[1])
BINARIES = {"master": pathlib.Path(sys.argv[2]), "new": pathlib.Path(sys.argv[3])}
EVENT = "del:chrT:10000-14000;af=1"


def md5(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


results = {}
for tag, binary in BINARIES.items():
    out = S / f"cr4-probe-{tag}"
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    probe = Probe(binary, out)
    for bam_label, kwargs in (("uniform", {}), ("lowmap", {"lowmap": True})):
        bam, _ = probe.make_bam(bam_label, **kwargs)
        dest = probe.run(f"del_{bam_label}", bam, [EVENT])
        log = (out / f"del_{bam_label}.log").read_text()
        truth = (dest / "truth.vcf").read_text()
        record = [l for l in truth.splitlines() if not l.startswith("#")]
        resist = re.search(r"SIM_RESIST=([^;\t]+)", record[0]) if record else None
        results[(tag, bam_label)] = {
            "R1": md5(dest / "R1.fq.gz"),
            "R2": md5(dest / "R2.fq.gz"),
            "replaced": md5(dest / "replaced_reads.txt"),
            "resist": resist.group(1) if resist else None,
            "warned": "spike cannot edit" in log and "WARN" in log,
            "truth_lines": truth.splitlines(),
        }

for bam_label in ("lowmap", "uniform"):
    m, n = results[("master", bam_label)], results[("new", bam_label)]
    print(f"== {bam_label}")
    print(f"  SIM_RESIST master={m['resist']} new={n['resist']}  warned new={n['warned']}")
    for key in ("R1", "R2", "replaced"):
        print(f"  {key}: master {m[key]}  new {n[key]}  same={m[key] == n[key]}")
    only_m = [l for l in m["truth_lines"] if l not in n["truth_lines"]]
    only_n = [l for l in n["truth_lines"] if l not in m["truth_lines"]]
    print(f"  truth lines only in master: {only_m}")
    print(f"  truth lines only in new:    {only_n}")
