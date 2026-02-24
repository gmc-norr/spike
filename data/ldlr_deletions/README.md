# LDLR Known Pathogenic Deletions (GRCh38/hg38)

Curated set of 10 known large structural deletions in the LDLR gene that cause
Familial Hypercholesterolemia (FH). These are real mutations documented in the
literature and clinical databases, formatted for use with `sv_exon_sim`.

## Gene: LDLR (chr19p13.2)

- **Transcript:** ENST00000252444 (Ensembl canonical) / NM_000527.5 (RefSeq)
- **Gene span (hg38):** chr19:11,090,578-11,133,816 (0-based half-open)
- **Exons:** 18
- **OMIM:** 606945
- **Alu content:** ~55% of gene sequence, making it highly susceptible to
  Alu-mediated non-allelic homologous recombination (NAHR)

Large rearrangements account for 5-10% of all pathogenic LDLR variants
depending on population.

## LDLR Exon Map (Ensembl ENST00000252444, GRCh38, 0-based half-open)

```
Exon  1: 11,090,578 - 11,090,919  (341 bp)   5' UTR + signal peptide
          --- intron 1: 9,303 bp ---
Exon  2: 11,100,222 - 11,100,345  (123 bp)   Ligand-binding repeat 1
          --- intron 2: 2,318 bp ---
Exon  3: 11,102,663 - 11,102,786  (123 bp)   Ligand-binding repeat 2
          --- intron 3: 2,433 bp ---
Exon  4: 11,105,219 - 11,105,600  (381 bp)   Ligand-binding repeats 3-4
          --- intron 4: 964 bp ---
Exon  5: 11,106,564 - 11,106,687  (123 bp)   Ligand-binding repeat 5
          --- intron 5: 704 bp ---
Exon  6: 11,107,391 - 11,107,514  (123 bp)   Ligand-binding repeat 6-7
          --- intron 6: 3,137 bp ---
Exon  7: 11,110,651 - 11,110,771  (120 bp)   EGF precursor homology (EGF-A)
          --- intron 7: 742 bp ---
Exon  8: 11,111,513 - 11,111,639  (126 bp)   EGF precursor homology (EGF-B)
          --- intron 8: 1,638 bp ---
Exon  9: 11,113,277 - 11,113,449  (172 bp)   Beta-propeller (YWTD 1)
          --- intron 9: 85 bp ---
Exon 10: 11,113,534 - 11,113,762  (228 bp)   Beta-propeller (YWTD 2)
          --- intron 10: 2,331 bp ---
Exon 11: 11,116,093 - 11,116,212  (119 bp)   Beta-propeller (YWTD 3)
          --- intron 11: 646 bp ---
Exon 12: 11,116,858 - 11,116,998  (140 bp)   Beta-propeller (YWTD 4-5)
          --- intron 12: 3,093 bp ---
Exon 13: 11,120,091 - 11,120,233  (142 bp)   Beta-propeller (YWTD 6)
          --- intron 13: 136 bp ---
Exon 14: 11,120,369 - 11,120,522  (153 bp)   EGF precursor homology (EGF-C)
          --- intron 14: 2,651 bp ---
Exon 15: 11,123,173 - 11,123,344  (171 bp)   O-linked sugar domain
          --- intron 15: 4,663 bp ---
Exon 16: 11,128,007 - 11,128,085  (78 bp)    Membrane-spanning domain
          --- intron 16: 1,427 bp ---
Exon 17: 11,129,512 - 11,129,670  (158 bp)   Cytoplasmic domain (NPXY)
          --- intron 17: 1,610 bp ---
Exon 18: 11,131,280 - 11,133,816  (2536 bp)  3' UTR
```

## Variants

### 1. Promoter + Exon 1 (French Canadian founder)

- **ID:** `LDLR_del_promoter_ex1_FrenchCanadian`
- **Size:** 15,888 bp
- **POS:** 11,077,742 (upstream of gene)
- **END:** 11,093,630 (intron 1: 11,090,919 - 11,100,222)
- **Exons deleted:** 1
- **Mechanism:** Alu NAHR (AluY elements)
- **Prevalence:** ~60% of French Canadian FH cases (Quebec founder, traceable to Kamouraska)
- **Effect:** Null allele — complete loss of transcription (promoter + exon 1 removed)
- **ClinVar:** [RCV000003899](https://www.ncbi.nlm.nih.gov/clinvar/RCV000003899/)
  (NC_000019.10 coordinates — native GRCh38)
- **Reference:** Boileau et al.; Hobbs et al. 1992; Bhm et al. 2004 (PMID: 14756670)

### 2. Exon 5 (Danish)

- **ID:** `LDLR_del_ex5_Danish`
- **Size:** 1,042 bp
- **POS:** 11,106,050 (intron 4: 11,105,600 - 11,106,564)
- **END:** 11,107,092 (intron 5: 11,106,687 - 11,107,391)
- **Exons deleted:** 5
- **Mechanism:** Alu NAHR (AluSx / AluJo)
- **Effect:** In-frame deletion removing ligand-binding repeat 5
- **Reference:** Nyegaard et al. BMC Med Genet 2006;7:55 (PMID: 16725053)

### 3. Exons 4-6 (ClinVar)

- **ID:** `LDLR_del_ex4-6_ClinVar`
- **Size:** 3,800 bp
- **POS:** 11,104,400 (intron 3: 11,102,786 - 11,105,219)
- **END:** 11,108,200 (intron 6: 11,107,514 - 11,110,651)
- **Exons deleted:** 4, 5, 6
- **Mechanism:** Alu NAHR
- **Effect:** In-frame deletion of ligand-binding repeats 3-7; disrupts LDL binding
- **ClinVar:** [VCV000526767](https://www.ncbi.nlm.nih.gov/clinvar/variation/526767/)

### 4. Exons 7-8 (Danish)

- **ID:** `LDLR_del_ex7-8_Danish`
- **Size:** 3,012 bp
- **POS:** 11,109,600 (intron 6: 11,107,514 - 11,110,651)
- **END:** 11,112,612 (intron 8: 11,111,639 - 11,113,277)
- **Exons deleted:** 7, 8
- **Mechanism:** Alu NAHR (AluSg/x on both ends)
- **Effect:** Disrupts EGF precursor homology domain (EGF-A and EGF-B)
- **Reference:** Nyegaard et al. BMC Med Genet 2006;7:55 (PMID: 16725053)

### 5. Exons 5-10 (Czech, NHEJ)

- **ID:** `LDLR_del_ex5-10_Czech`
- **Size:** 7,636 bp
- **POS:** 11,106,300 (intron 4: 11,105,600 - 11,106,564)
- **END:** 11,113,936 (intron 10: 11,113,762 - 11,116,093)
- **Exons deleted:** 5, 6, 7, 8, 9, 10
- **Mechanism:** Non-homologous end joining (NHEJ) — unlike most LDLR deletions
- **Prevalence:** 4 unrelated probands in Czech cohort (1,441 FH patients)
- **Effect:** Frameshift; removes ligand-binding repeats 5-7 and EGF-A/B domains
- **Reference:** Goldmann et al. BMC Med Genet 2010;11:115 (PMID: 20727145)

### 6. Exons 3-12 (Czech)

- **ID:** `LDLR_del_ex3-12_Czech`
- **Size:** 17,604 bp
- **POS:** 11,101,200 (intron 2: 11,100,345 - 11,102,663)
- **END:** 11,118,804 (intron 12: 11,116,998 - 11,120,091)
- **Exons deleted:** 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
- **Mechanism:** Alu NAHR (FLAM_A / AluY)
- **Prevalence:** 1 proband
- **Effect:** Loss of most ligand-binding and EGF precursor homology domains
- **Reference:** Goldmann et al. BMC Med Genet 2010;11:115 (PMID: 20727145)

### 7. Exons 9-14 (Czech, most common)

- **ID:** `LDLR_del_ex9-14_Czech`
- **Size:** 9,713 bp
- **POS:** 11,112,200 (intron 8: 11,111,639 - 11,113,277)
- **END:** 11,121,913 (intron 14: 11,120,522 - 11,123,173)
- **Exons deleted:** 9, 10, 11, 12, 13, 14
- **Mechanism:** Alu NAHR (AluSq on both ends)
- **Prevalence:** 10 unrelated probands — most common large deletion in Czech population
- **Effect:** Disrupts EGF precursor homology domain; likely class 2 (transport-defective)
- **Reference:** Goldmann et al. BMC Med Genet 2010;11:115 (PMID: 20727145)

### 8. Exons 9-15 (Czech)

- **ID:** `LDLR_del_ex9-15_Czech`
- **Size:** 12,400 bp
- **POS:** 11,112,100 (intron 8: 11,111,639 - 11,113,277)
- **END:** 11,124,500 (intron 15: 11,123,344 - 11,128,007)
- **Exons deleted:** 9, 10, 11, 12, 13, 14, 15
- **Mechanism:** Alu NAHR (AluJb / AluSx1)
- **Prevalence:** 8 unrelated probands
- **Effect:** Disrupts EGF precursor homology and O-linked sugar domains
- **Reference:** Goldmann et al. BMC Med Genet 2010;11:115 (PMID: 20727145)

### 9. Exons 13-15 (Danish)

- **ID:** `LDLR_del_ex13-15_Danish`
- **Size:** 6,298 bp
- **POS:** 11,118,800 (intron 12: 11,116,998 - 11,120,091)
- **END:** 11,125,098 (intron 15: 11,123,344 - 11,128,007)
- **Exons deleted:** 13, 14, 15
- **Mechanism:** Alu NAHR (AluSg/x on both ends)
- **Note:** Original breakpoint characterization showed 15 bp insertion at junction
- **Effect:** Disrupts EGF precursor homology domain (EGF-C) and O-linked sugar domain
- **Reference:** Nyegaard et al. BMC Med Genet 2006;7:55 (PMID: 16725053)

### 10. Exons 4-18 (Italian, near-complete gene loss)

- **ID:** `LDLR_del_ex4-18_Italian`
- **Size:** 30,400 bp
- **POS:** 11,104,000 (intron 3: 11,102,786 - 11,105,219)
- **END:** 11,134,400 (past gene end at 11,133,816)
- **Exons deleted:** 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18
- **Mechanism:** Alu NAHR (AluSx in intron 3 / AluJb in exon 18)
- **Effect:** Loss of nearly all functional domains (ligand-binding, EGF precursor,
  O-linked glycan, transmembrane, cytosolic). Null allele.
- **Reference:** Genes 2023;14(6):1275 (MDPI, Italian family)

## Usage

Simulate all deletions (each as heterozygous):

```bash
sv_exon_sim \
  --bam <BAM> \
  --reference <FASTA> \
  --vcf crates/sv_exon_sim/data/ldlr_deletions/ldlr_known_deletions_hg38.vcf \
  -o /tmp/ldlr_sim
```

## Coordinate Conventions

- **VCF POS**: 1-based preceding base (numerically equal to 0-based start)
- **VCF END**: 1-based inclusive (numerically equal to 0-based exclusive end)
- **Internal (sv_exon_sim)**: 0-based half-open `[start, end)`
- **SIM_VAF=0.5**: heterozygous (standard for autosomal dominant FH)

## Breakpoint Verification

All breakpoints were verified to fall within the correct introns based on
Ensembl ENST00000252444 exon coordinates on GRCh38. The French Canadian
deletion coordinates (variant 1) are native GRCh38 from ClinVar
(NC_000019.10). All other breakpoints are placed in intronic regions
consistent with the Alu-mediated NAHR mechanism reported in the literature.

The exact base-pair positions of intronic breakpoints are approximate — real
Alu-Alu recombination breakpoints occur within the homologous Alu sequence
and vary by tens of bp between patients carrying the "same" deletion.

## References

1. Nyegaard et al. "Genomic characterization of five deletions in the LDLR gene
   in Danish Familial Hypercholesterolemia patients." BMC Med Genet 2006;7:55.
   PMID: 16725053. DOI: 10.1186/1471-2350-7-55

2. Goldmann et al. "Genomic characterization of large rearrangements of the LDLR
   gene in Czech patients with familial hypercholesterolemia." BMC Med Genet
   2010;11:115. PMID: 20727145. DOI: 10.1186/1471-2350-11-115

3. Bhm et al. "LDLR gene rearrangement analysis in French-Canadian familial
   hypercholesterolemia subjects." J Lipid Res 2004;45:1456-1461.
   PMID: 14756670.

4. ClinVar RCV000003899 — French Canadian >15kb deletion (NC_000019.10,
   native GRCh38 coordinates).

5. ClinVar VCV000526767 — Exon 4-6 deletion.

6. LDLR exon 4-18 deletion in Italian family. Genes 2023;14(6):1275.
   DOI: 10.3390/genes14061275
