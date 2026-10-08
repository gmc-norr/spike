"""Mutation checks for the duplicates fix: the unmutated suite must be green, then each mutant red."""
import os, subprocess, sys
D, T = sys.argv[1], sys.argv[2]
M = [
 ("--fastq-prefix alone accepted", "src/main.rs", '#[arg(long, value_name = "NAME", requires = "into_fastq")]', '#[arg(long, value_name = "NAME")]'),
 ("a / allowed", "src/main.rs", 'p == ".." || p.contains(\'/\')', 'p == ".."'),
 ("'..' allowed", "src/main.rs", 'p == "." || p == ".." ||', 'p == "." ||'),
 ("default renamed", "src/main.rs", 'const DEFAULT_FASTQ_PREFIX: &str = "spiked";', 'const DEFAULT_FASTQ_PREFIX: &str = "spike";'),
 ("prefix after the mate", "src/main.rs", 'format!("{}_R{}.fastq.gz", prefix, mate)', 'format!("R{1}_{0}.fastq.gz", prefix, mate)'),
 ("run_fastq ignores the prefix", "src/main.rs", "    let [r1, r2] = whole_sample_paths(output_dir, prefix);\n    let status", "    let [r1, r2] = whole_sample_paths(output_dir, DEFAULT_FASTQ_PREFIX);\n    let status"),
 ("old name kept as alias", "src/main.rs", '    #[arg(long, num_args = 2, value_names = ["RAW_R1", "RAW_R2"])]\n    into_fastq', '    #[arg(long, alias = "raw-fastq", num_args = 2, value_names = ["RAW_R1", "RAW_R2"])]\n    into_fastq'),
 ("readme rows swapped", "src/main.rs", "    match whole_sample {\n        Some(prefix) =>", "    match whole_sample.filter(|_| false) {\n        Some(prefix) =>"),
 ("next steps branches swapped", "src/main.rs", "    if into_fastq {\n        format!(\n            \"Next steps: the whole", "    if !into_fastq {\n        format!(\n            \"Next steps: the whole"),
 ("R1/R2 line drops the warning", "src/main.rs", 'around the events only (not the whole sample): {} and {}', 'around the events: {} and {}'),
 ("check message keeps the old name", "src/main.rs", 'bail!("--into-fastq {} is not a file', 'bail!("--raw-fastq {} is not a file'),
]
def test():
    r = subprocess.run(["cargo", "test", "-j", "16", "--release"], cwd=D, capture_output=True, text=True,
                       env={**os.environ, "CARGO_TARGET_DIR": T})
    out = r.stdout + r.stderr
    failed = [l.split()[1] for l in out.splitlines() if l.startswith("test ") and l.endswith("FAILED")]
    summary = [l for l in out.splitlines() if l.startswith("test result")]
    return r.returncode, failed, summary, out
rc, failed, summary, out = test()
print(f"unmutated: exit {rc}, {summary}", flush=True)
if rc != 0:
    print(out[-3000:]); sys.exit("the unmutated suite is not green; no score")
caught = 0
for label, f, old, new in M:
    p = os.path.join(D, f)
    src = open(p).read()
    assert src.count(old) == 1, label
    open(p, "w").write(src.replace(old, new))
    try:
        rc, failed, summary, out = test()
    finally:
        open(p, "w").write(src)
    built = any(l.startswith("test result") for l in out.splitlines())
    red = rc != 0 and built and failed
    caught += bool(red)
    print(f"{'CAUGHT' if red else 'SURVIVED'}  {label}: {failed if built else 'BUILD FAILED'}", flush=True)
print(f"{caught} of {len(M)} caught", flush=True)
