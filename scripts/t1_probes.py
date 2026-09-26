#!/usr/bin/env python3
"""T1 criteria C1-C4: the advisory census rows on the Codex review's probes.

Builds the review script's three donor BAMs (`uniform`, `lowmap`, `variable`),
runs CR4's C1 deletion and CR2's C1 duplication on each, and runs `spike validate`
on the donor BAM against each run's truth VCF. The census rows read the truth VCF,
so the donor BAM is a fine `--bam` here: C5 and C6 use a merged real-data BAM.

C4 also strips `SIM_RESIST` and `SIM_DEPTH_FOLD` (and their header lines) from one
truth VCF and validates that copy under both binaries.

Usage: t1_probes.py OUT_DIR SPIKE_NEW SPIKE_BASE
Prints one `key=value` line per measurement; scores nothing.
"""
import pathlib
import re
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from review_sv_model import Probe  # noqa: E402

DEL = "del:chrT:10000-14000;af=1"
DUP = "dup:chrT:10000-28000;af=0.5"


def validate(spike, bam, truth, fasta, out, extra=()):
    command = [str(spike), "validate", "--bam", str(bam), "--truth", str(truth),
               "--reference", str(fasta), *extra]
    done = subprocess.run(command, capture_output=True, text=True)
    # Explicit suffixes: `Path.with_suffix` would turn both `x.validate` and
    # `x.json` into the same `x.txt` and the second run would clobber the first.
    pathlib.Path(str(out) + ".out").write_text(done.stdout)
    pathlib.Path(str(out) + ".err").write_text(done.stderr)
    return done.returncode, done.stdout


def row(table, check):
    r"""The report line naming `check`, or None.

    The name must start the line or follow whitespace and must end at
    whitespace or a `:`. A leading `\s` alone could never match `Advisory`,
    which starts at column 0, so this script's own "no advisory line" check was
    a silent no-op. The trailing lookahead keeps `split_reads` from matching
    `split_reads_each_end`.
    """
    for line in table.splitlines():
        if re.search(rf"(?:^|\s){re.escape(check)}(?=[\s:]|$)", line):
            return line.rstrip()
    return None


def main(out_dir, spike_new, spike_base):
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    probe = Probe(pathlib.Path(spike_new), out)
    bams = {
        "uniform": probe.make_bam("uniform")[0],
        "lowmap": probe.make_bam("lowmap", lowmap=True)[0],
        "variable": probe.make_bam("variable", variable=True)[0],
    }
    for name, bam in bams.items():
        for tag, event in (("del", DEL), ("dup", DUP)):
            label = f"{name}_{tag}"
            # About half the lowmap pairs are uneditable, which spike refuses
            # since RF8; T1 measures the census on the run it is told to make.
            extra = ("--allow-resistant",) if name == "lowmap" else ()
            run = probe.run(label, bam, events=(event,), extra=extra)
            truth = run / "truth.vcf"
            info = [l for l in truth.read_text().splitlines() if not l.startswith("#")]
            for field in ("SIM_RESIST", "SIM_DEPTH_FOLD"):
                found = re.search(rf"{field}=([^;\t]*)", info[0])
                print(f"{label}.{field}={found.group(1) if found else 'ABSENT'}")
            code, table = validate(spike_new, bam, truth, probe.fasta,
                                   out / f"{label}.validate")
            print(f"{label}.exit={code}")
            for check in ("resistant", "depth_fold"):
                print(f"{label}.row.{check}={row(table, check)}")
            code, _ = validate(spike_new, bam, truth, probe.fasta,
                               out / f"{label}.json", extra=("--json",))
            print(f"{label}.json.exit={code}")

    # C4: the same truth VCF with the two census fields and header lines removed.
    src = out / "lowmap_del" / "truth.vcf"
    stripped = out / "lowmap_del_nocensus.vcf"
    kept = []
    for line in src.read_text().splitlines():
        if line.startswith("##INFO=<ID=SIM_RESIST") or line.startswith("##INFO=<ID=SIM_DEPTH_FOLD"):
            continue
        kept.append(re.sub(r"SIM_RESIST=[^;]*;SIM_DEPTH_FOLD=[^;]*;", "", line))
    stripped.write_text("\n".join(kept) + "\n")
    for tag, spike in (("new", spike_new), ("base", spike_base)):
        code, table = validate(spike, bams["lowmap"], stripped, probe.fasta,
                               out / f"nocensus_{tag}.validate")
        print(f"nocensus.{tag}.exit={code}")
        print(f"nocensus.{tag}.rows={len([l for l in table.splitlines() if 'PASS' in l or 'FAIL' in l])}")
        for check in ("resistant", "depth_fold", "Advisory"):
            print(f"nocensus.{tag}.row.{check}={row(table, check)}")
        (out / f"nocensus_{tag}.table").write_text(table)


if __name__ == "__main__":
    main(*sys.argv[1:4])
