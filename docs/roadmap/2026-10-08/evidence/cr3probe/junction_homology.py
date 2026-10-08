"""Read-only reference arithmetic: junction homology of deletions.
homology = left + right, where left = k with ref[S-k:S] == ref[E-k:E] and
right = k with ref[S:S+k] == ref[E:E+k] (0-based half-open deletion [S,E))."""
import sys, pysam

fa = pysam.FastaFile(sys.argv[1])


def hom(chrom, s, e, cap=200):
    left = 0
    while left < cap and fa.fetch(chrom, s - left - 1, s - left) == fa.fetch(chrom, e - left - 1, e - left):
        left += 1
    right = 0
    while right < cap and fa.fetch(chrom, s + right, s + right + 1) == fa.fetch(chrom, e + right, e + right + 1):
        right += 1
    return left + right


def summary(label, dels):
    hs = sorted(hom(*d) for d in dels)
    n = len(hs)
    print(label, "n", n, "median", hs[n // 2], ">=2bp", sum(h >= 2 for h in hs), "blunt(0)", sum(h == 0 for h in hs), "max", hs[-1])

# demo VCF (repo data): POS = base before; END = last deleted (1-based)
demo = []
for line in open(sys.argv[2]):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if "SVTYPE=DEL" not in f[7]:
        continue
    end = int([x for x in f[7].split(";") if x.startswith("END=")][0][4:])
    demo.append((f[0], int(f[1]), end))
summary("demo LDLR DELs (repo VCF)", demo)

# ClinVar pathogenic >=50 bp deletions (BED, 0-based) on chr19 LDLR +- 50 kb
cv = []
for line in open(sys.argv[3]):
    f = line.split("\t")
    if f[0] == "19" and 11040000 < int(f[1]) and int(f[2]) < 11190000 and f[3].startswith("Deletion"):
        if int(f[2]) - int(f[1]) <= 200000:
            cv.append(("chr19", int(f[1]), int(f[2])))
summary("ClinVar LDLR-region DELs", cv)
