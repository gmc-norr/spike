"""Which Platinum haplotype (GT first = hap1) do spike's synthetic reads carry, per event?

Counts 25-mers centred on each Platinum het SNP in the footprint, in SPIKE_ reads only.
"""
import gzip, subprocess, sys

FA = "/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
PV = "/home/parlar_ai/dev/cnv_validation/platinum_pedigree_truthset_v1.2/NA12878_hq_v1.2.1.latest-smallvar.vcf.gz"
EVENTS = [11100300, 11113400]
K = 12


def ref(start, end):
    out = subprocess.run(["samtools", "faidx", FA, f"chr19:{start}-{end}"], capture_output=True, text=True).stdout
    return "".join(out.split("\n")[1:]).upper()


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


probes = {}
for ev in EVENTS:
    out = subprocess.run(["bcftools", "view", "-H", "-v", "snps", "-g", "het", "-r",
                          f"chr19:{ev-2000}-{ev+2000}", PV], capture_output=True, text=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        pos, r, a, gt = int(f[1]), f[3], f[4], f[9].split(":")[0]
        h1, h2 = (a, r) if gt == "1|0" else (r, a)
        left, right = ref(pos - K, pos - 1), ref(pos + 1, pos + K)
        probes[(ev, pos)] = {"hap1": left + h1 + right, "hap2": left + h2 + right}

for d in sys.argv[1:]:
    tally = {ev: {"hap1": 0, "hap2": 0} for ev in EVENTS}
    for fq in ("R1.fq.gz", "R2.fq.gz"):
        with gzip.open(f"{d}/{fq}", "rt") as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline().strip()
                f.readline(); f.readline()
                if not h[1:].startswith("SPIKE_"):
                    continue
                for (ev, pos), pr in probes.items():
                    for lab, kmer in pr.items():
                        if kmer in s or rc(kmer) in s:
                            tally[ev][lab] += 1
    calls = []
    for ev in EVENTS:
        t = tally[ev]
        calls.append("hap1" if t["hap1"] > t["hap2"] else "hap2")
    rel = "cis" if calls[0] == calls[1] else "trans"
    print(d, {ev: tally[ev] for ev in EVENTS}, calls, rel)
