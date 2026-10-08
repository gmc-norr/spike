import gzip, subprocess, sys

FA = "/home/parlar_ai/dev/spike/data/giab_hg38/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
# (1-based pos, allele on Platinum hap1 (GT first), allele on hap2)
SNPS = [(11103857, "T", "G"), (11105885, "C", "A"), (11105941, "A", "G")]
EVENT = (11105500, "G", "A")
K = 12


def ref(start, end):
    out = subprocess.run(["samtools", "faidx", FA, f"chr19:{start}-{end}"], capture_output=True, text=True).stdout
    return "".join(out.split("\n")[1:]).upper()


def rc(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


probes = {}
for pos, h1, h2 in SNPS + [EVENT]:
    left = ref(pos - K, pos - 1)
    right = ref(pos + 1, pos + K)
    probes[pos] = {"hap1" if (pos, h1, h2) in SNPS else "REF": left + h1 + right,
                   "hap2" if (pos, h1, h2) in SNPS else "ALT": left + h2 + right}

for d in sys.argv[1:]:
    counts = {}
    alt_reads = []
    names_seqs = []
    for fq in ("R1.fq.gz", "R2.fq.gz"):
        with gzip.open(f"{d}/{fq}", "rt") as f:
            while True:
                h = f.readline()
                if not h:
                    break
                s = f.readline().strip()
                f.readline(); f.readline()
                names_seqs.append((h.split()[0][1:].split(":")[0].split("/")[0], s))
    # fragment = name; pool alleles by fragment
    frag = {}
    for name, s in names_seqs:
        spike = name.startswith("SPIKE_")
        for pos, pr in probes.items():
            for lab, kmer in pr.items():
                if kmer in s or rc(kmer) in s:
                    key = ("spike" if spike else "orig", pos, lab)
                    counts[key] = counts.get(key, 0) + 1
                    frag.setdefault(name, set()).add((pos, lab))
    print("==", d)
    for k in sorted(counts):
        print("  ", k, counts[k])
    # among fragments carrying the event ALT, which SNP haplotype do they carry
    link = {}
    for name, al in frag.items():
        if (EVENT[0], "ALT") in al:
            for pos, lab in al:
                if pos != EVENT[0]:
                    link[(pos, lab)] = link.get((pos, lab), 0) + 1
    print("   ALT-fragment linked SNP alleles:", dict(sorted(link.items())))
