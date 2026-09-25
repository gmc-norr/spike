#!/usr/bin/env python3
"""N15, third attempt: build site lists the way N15 and N18 did.

sites:     GIAB PASS, biallelic, het indels, REF and ALT <= 11 bp, in spike's
           truth format at SIM_VAF=0.5 (the lists N15 and N18 used).
isolated:  sites with no other GIAB record's POS in [POS - 25, POS + 25 + len(REF)].
context:   (optional) every other GIAB PASS record within --context bp of a
           site, one line per ALT, SIM_VAF 1.0 if hom-alt else 0.5.

Usage: n15c_sites.py GIAB.vcf.gz HEADER_FROM.vcf OUT_PREFIX [--context N] chr...

HEADER_FROM is any truth VCF spike wrote; its header is copied.
"""
import argparse
import bisect
import subprocess


def records(giab, chroms, extra=()):
    cmd = ["bcftools", "view", "-H", "-f", "PASS", *extra, giab, *chroms]
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        c = line.split("\t")
        yield c[0], int(c[1]), c[3], c[4], c[9].split(":")[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("giab"); ap.add_argument("header_from"); ap.add_argument("out")
    ap.add_argument("--context", type=int, default=0)
    ap.add_argument("chroms", nargs="+")
    a = ap.parse_args()
    header = [l for l in open(a.header_from) if l.startswith("#")]
    sites = [(ch, p, r, al) for ch, p, r, al, _ in records(a.giab, a.chroms, ["-g", "het", "-m2", "-M2", "-v", "indels"])
             if len(r) <= 11 and len(al) <= 11]
    allrec = list(records(a.giab, a.chroms))
    by_chrom = {}
    for ch, p, r, al, gt in allrec:
        by_chrom.setdefault(ch, []).append((p, r, al, gt))
    with open(a.out + "-sites.vcf", "w") as f:
        f.writelines(header)
        for i, (ch, p, r, al) in enumerate(sites, 1):
            f.write(f"{ch}\t{p}\tv{i}\t{r}\t{al}\t999\tPASS\tSIM_VAF=0.500;SIM_GENE=unknown\tGT\t0/1\n")
    site_keys = {(ch, p, r, al) for ch, p, r, al in sites}
    # As N15 built it: no other record's POS in [POS - 25, POS + 25 + len(REF)].
    allpos = {ch: sorted(q for q, _, _, _ in v) for ch, v in by_chrom.items()}
    with open(a.out + "-isolated.txt", "w") as f:
        for ch, p, r, al in sites:
            lo = bisect.bisect_left(allpos[ch], p - 25)
            hi = bisect.bisect_right(allpos[ch], p + 25 + len(r))
            if hi - lo == 1:
                f.write(f"{ch}:{p}\n")
    if a.context:
        pos_by = {ch: sorted(p for c2, p, _, _ in sites if c2 == ch) for ch in a.chroms}
        out, n = [], len(sites)
        for ch, p, r, al, gt in allrec:
            if (ch, p, r, al) in site_keys:
                continue
            ps = pos_by.get(ch, [])
            i = bisect.bisect_left(ps, p - a.context)
            if i < len(ps) and ps[i] <= p + a.context:
                for alt in al.split(","):
                    n += 1
                    hom = gt.replace("|", "/") in ("1/1", "2/2")
                    vaf = "1.000" if hom else "0.500"
                    out.append(f"{ch}\t{p}\tv{n}\t{r}\t{alt}\t999\tPASS\tSIM_VAF={vaf};SIM_GENE=unknown\tGT\t{'1/1' if hom else '0/1'}\n")
        with open(a.out + "-context.vcf", "w") as f:
            f.writelines(header)
            lines = sorted(
                [l for l in open(a.out + "-sites.vcf") if not l.startswith("#")] + out,
                key=lambda l: (a.chroms.index(l.split("\t")[0]), int(l.split("\t")[1])))
            f.writelines(lines)


if __name__ == "__main__":
    main()
