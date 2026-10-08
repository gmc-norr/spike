"""Mutation checks for the duplicates fix: the unmutated suite must be green, then each mutant red."""
import os, subprocess, sys
D, T = sys.argv[1], sys.argv[2]
M = [
 ("strand left out of the key", "src/origin.rs", "        let family = Fragment { name, mates }.family();\n",
  "        let family: Vec<FivePrime> = Fragment { name, mates }.family().into_iter().map(|f| FivePrime { reverse: false, ..f }).collect();\n"),
 ("any family returned", "src/origin.rs", ".filter(|(family, _)| removed_families.contains(family))", ".filter(|(family, _)| removed_families.contains(family) || true)"),
 ("a one-mate fragment used", "src/origin.rs", "let pair = mates.len() == 2 && mates.iter().any(|m| m.first) && mates.iter().any(|m| !m.first);", "let pair = !mates.is_empty();"),
 ("flag not required on both mates", "src/origin.rs", "        if all_flagged {\n            copies.push", "        if !none_flagged {\n            copies.push"),
 ("an unflagged pair returned", "src/origin.rs", "        if all_flagged {\n            copies.push", "        if all_flagged || !removed.contains(name) {\n            copies.push"),
 ("QC-failed records kept", "src/origin.rs", "for r in records.iter().filter(|r| !r.qc_fail) {", "for r in records.iter() {"),
 ("no dedup across spans", "src/origin.rs", "&& seen.insert((r.name.clone(), r.first))", "&& { seen.insert((r.name.clone(), r.first)); true }"),
 ("no widening", "src/main.rs", "    let pad = crate::stats::MAX_FRAGMENT_LEN as u64;\n    let mut sides", "    let pad = 0u64;\n    let mut sides"),
 ("a fusion's second side dropped", "src/main.rs", "            sides.push((chrom_b, *bp_b, *bp_b));\n", ""),
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
