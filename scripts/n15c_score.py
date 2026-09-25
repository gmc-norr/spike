#!/usr/bin/env python3
"""N15, third attempt: score the locked criteria from `spike validate` runs.

"Out of range" is an allele_freq FAIL among the listed sites whose observed
value is a fraction (sites below depth 5 are left out, as in N18 and N19).
Only sites in --sites are scored; the context records are not.

Usage:
    n15c_score.py --sites S.vcf --isolated I.txt --base B.json --second X.json --new N.json \\
        [--spellings a.err b.err c.err]
"""
import argparse
import json
import re

NUM = re.compile(r"^\d+\.\d+$")
COUNT = re.compile(r"allele_freq SNP (\S+):(\d+)-\d+ \(unknown\): (\d+) carry, (\d+) span")


def verdicts(path, keys):
    out = {}
    for c in json.load(open(path))["checks"]:
        if c["check"] == "allele_freq":
            m = re.search(r" (\S+):(\d+)-", c["event"])
            k = f"{m.group(1)}:{int(m.group(2)) + 1}"  # 0-based label -> VCF POS
            if k in keys:
                assert k not in out, f"two allele_freq rows for {k}"
                out[k] = c
    assert set(out) == keys, f"{path}: {len(keys - set(out))} sites have no allele_freq row"
    return out


def rate(v, ks):
    graded = [k for k in ks if NUM.match(v[k]["observed"])]
    fail = sum(not v[k]["pass"] for k in graded)
    return fail, len(graded), 100 * fail / len(graded)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    for f in ("sites", "isolated", "base", "second", "new"):
        ap.add_argument(f"--{f}", required=True)
    ap.add_argument("--spellings", nargs=3)
    a = ap.parse_args()
    keys = {f"{c[0]}:{c[1]}" for c in (l.split("\t") for l in open(a.sites) if not l.startswith("#"))}
    iso = set(open(a.isolated).read().split())
    assert iso <= keys
    classes = {"isolated": iso, "not isolated": keys - iso, "overall": keys}
    print(f"sites {len(keys)}: isolated {len(iso)}, not isolated {len(keys - iso)}")
    r = {}
    for run in ("base", "second", "new"):
        v = verdicts(getattr(a, run), keys)
        for cls, ks in classes.items():
            f, g, x = rate(v, ks)
            r[(run, cls)] = x
            print(f"  {run:6s} {cls:13s} out of range {f}/{g} = {x:.2f}%")
    c = {
        "C1 not isolated falls >= 5 pts below master": r[("base", "not isolated")] - r[("new", "not isolated")] >= 5,
        "C2 isolated rises <= 0.5 pts above master": r[("new", "isolated")] - r[("base", "isolated")] <= 0.5,
        "C3 overall below master": r[("new", "overall")] < r[("base", "overall")],
        "C4 not isolated below the second attempt": r[("new", "not isolated")] < r[("second", "not isolated")],
    }
    if a.spellings:
        got = []
        for p in a.spellings:
            got.append([(int(m.group(3)), int(m.group(4))) for m in map(COUNT.search, open(p)) if m])
        print(f"  spellings: {got}")
        c["C5 the three spellings agree"] = all(len(g) == 1 for g in got) and got[0] == got[1] == got[2]
    for name, ok in c.items():
        print(f"{name}: {'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
