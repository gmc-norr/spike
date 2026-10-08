# Other tools: facts and sources behind Table 8 (slides 23-24)

Gathered 2026-10-08 by a research agent from each tool's paper, documentation, repository code and release pages. Three claims were re-checked by hand: BAMSurgeon's `doc/Manual.md` (SV reads simulated with wgsim over a reference template), wgsim's `wgsim.c` (one quality for every base) and NEAT's `neat/read_simulator/utils/vcf_func.py` (SNV, insertion and deletion only; other types become `UnknownVariant`). The rest is as the agent reported it. These tools change; re-check before reusing the slides.

Every cell below has a primary source. Where a source did not say something, the cell reads "not found". Two cautions for the slide:

- **BAMSurgeon has changed since its last release.** Release 1.4.1 (2022) builds SV reads by Velvet local assembly of the sample's reads plus wgsim. The current master, rewritten in August 2026 and not released, builds SV reads from reference slices plus wgsim. `doc/Manual.pdf` no longer exists; it is now `doc/Manual.md`.
- **NEAT's README overstates its SV support.** The README lists inversions, translocations and duplications. The code and the 2026 JOSS paper support only SNVs, insertions and deletions.

## Fact table

| Tool | 1 Family | 2 Types | 3 How variant reads are made | 4 Quality of new reads | 5 Keeps sample's own haplotype | 6 Focus / AF | 7 Truth | 8 Built-in check | 9 Output | 10 Language; release; last commit |
|---|---|---|---|---|---|---|---|---|---|---|
| BAMSurgeon | Existing BAM/CRAM | SNV, indel, DEL, DUP, INV, INS, TRN/BND; `<CNV>` skipped (copy-number file only adjusts VAF) | SNV/indel: bases edited in place, then realigned. SV in 1.4.1: Velvet assembly, edit contig, wgsim reads. SV in master: reference-slice template, wgsim reads | wgsim; `--simerr` default 0 gives one fixed quality, 'I' (Q40). Optional real "donor" reads for DUP interiors | SNV/indel: skips sites near existing het alleles (`--snvfrac`). Master's SV template "carries no small variants" | Somatic; per-variant VAF | Truth VCF (achieved VAF, SOMATIC); optional BS tag on altered reads | Depth and coverage-ratio checks drop failing sites; master adds `validate_spikein.py` | BAM only (no FASTQ output found) | Python 3; 1.4.1, 2022-08-24 (GitHub and Bioconda); 2026-08-05 |
| SomatoSim | Existing BAM | SNV only | Base changed in place, no realignment | N/A (keeps original BQ/MQ) | Not found (keeps the strand ratio) | Somatic; VAF range or per site | `simulation_output.txt` (no VCF) | Yes: reports output VAF and coverage | BAM | Python 3; v1.0.0 (setup.py; no tags, not on Bioconda); 2021-03-29 |
| VarBen | Existing BAM (Illumina, Ion Torrent) | SNV, ins, del, complex delins; SV: del, inv, tandem dup, 3 translocation kinds; CNV ("working on testing") | Original reads edited, then remapped. For SVs, code splices reference sequence into the read past the breakpoint | N/A (code keeps original qualities) | `--snpfrac` avoids spiking on top of existing het alleles; `--haplosize` | Somatic (oncology); per-variant AF | `success_list.txt`, `invalid_mutation.txt` | Depth, minimum-mutated-reads and coverage-difference filters; achieved-VAF report not found | BAM | Python 2; 1.0 (setup.py, no releases); 2021-03-16 |
| SVEngine | New reads from haplotype FASTA (BAM input option commented out in code) | DEL, INS (domestic/foreign), INV, DUP, inverted DUP, TRA | Contigs edited, reads made by xwgsim (modified wgsim), aligned with bwa mem | Uniform error rate, one fixed quality per run (xwgsim code) | Variant placed on a chosen haplotype (VARHAP); one FASTA per haplotype; using the sample's variants not found | Germline, somatic, clonal tree; per-variant fraction | Altered FASTA + VAR/BED truth; `var2vcf.py` in repo | Not found | FASTA, FASTQ (default) or BAM | Python 2 + C; 1.0.0 (setup.py, no tags); 2019-11-14 (Bitbucket) |
| NEAT | New reads from reference + random and/or VCF variants | Code: SNV, ins, del; other VCF types become a placeholder | Own simulator samples reads from the mutated reference | Default model (human Illumina) or learned offline from a real FASTQ with a separate utility | `include_vcf` (README example uses NA12878.vcf); phasing not found | Mainly germline; genotype by ploidy; README contradicts itself on tumour/normal | Golden VCF (coverage, allele balance), golden BAM | Not found (`compare-vcfs` scores a caller's VCF) | FASTQ (whole genome possible), BAM, VCF | Python ≥3.10; 4.7.1, 2026-09-15 (Bioconda 4.7.1); 2026-09-15 |
| ART | New reads from any FASTA | None (read simulator only) | N/A | Built-in empirical profiles, or user profiles via `art_profiler_illumina` | N/A | N/A | ALN, SAM (incl. error-free SAM) | N/A | FASTQ | C++; MountRainier 2016-06-05 ("the latest version"); not found |
| VarSim | New reads from a simulated diploid genome | SNV, ins, del, MNP, complex; SV ins, del, dup, inv; translocations from v0.7.0 | Modified vcf2diploid builds the genome; reads by ART (default) or DWGSIM | ART profiles (`--profile_1/2`) | Diploid; `--vc_prop_het`; phasing not found | Germline + somatic (COSMIC); normal contamination by mixing | Truth VCF; true locations in read names; map file | Downstream validation only | FASTQ | Java + Python; v0.8.6, 2020-09-26; 2022-01-28 |
| VISOR | New reads from HACk-built FASTA haplotypes | DEL, INV, tandem/inverted DUP, INS, STR expansion/contraction, new tandem repeats, 3 translocation kinds, SNP, MNP | Haplotype FASTA; short reads by pywgsim (wgsim core), long reads by Badread, 10x reads (beta); minimap2 alignment | Short: one fixed quality (Q17 at default `--error 0.02`). Long: Badread presets or user models | Haplotype-specific; "optionally with nearby single-nucleotide variants"; HP tag | Purity, clone fractions | BAM with optional HP/CL tags; VCF not found | Not found | Sorted BAM; FASTQ with `--fastq` | Python 3; v1.1.2.1, 2023-11-27 (Bioconda 1.1.2.1); 2024-10-21 |

## Sources and quotes

**BAMSurgeon**
- Paper: Ewing AD, Houlahan KE, Hu Y, et al. Nat Methods 2015;12(7):623–630, doi:10.1038/nmeth.3407 (europepmc PMC4856034).
- SV paper: Lee AY, Ewing AD, et al. Genome Biol 2018;19:188, doi:10.1186/s13059-018-1539-5.
- Repo docs: github.com/adamewing/bamsurgeon `doc/Manual.md` and `doc/design-notes.md` (master); `doc/Manual.tex` and `bin/addsv.py` at tag 1.4.1.
- wgsim quality code: github.com/lh3/wgsim `wgsim.c` line 248, `Q = (ERR_RATE == 0.0)? 'I' : …`.
- 3, paper: "modifying reads covering the selected sites, realigning a requisite number…"
- 3, Manual.tex 1.4.1: "regions… are assembled using Velvet… Read coverage is simulated over the contig using wgsim… -e 0 -r 0 -R 0".
- 3, Manual.md (master): "builds the mutated haplotype directly from reference sequence, simulates reads over it with `wgsim`".
- 4: addsv.py 1.4.1 passes `'-e', str(err_rate)`, with `--simerr` default 0.0. Manual.md: with `--donorbam` "the interior comes from real reads with a real error profile."
- 5, Manual.tex: "we suggest phasing BAMs beforehand and simulating mutations on a per-haplotype basis."
- 9, Manual.md: "writes the mutated BAM and a truth VCF".

**SomatoSim**
- Paper: Hawari MA, Hong CS, Biesecker LG. BMC Bioinformatics 2021;22:109, doi:10.1186/s12859-021-04024-8 (PMC7936459). Repo README: github.com/BieseckerLab/SomatoSim.
- 3/4: "The reads are not realigned… the read and the variant allele retain their original read MQ and BQ".
- 9: "SomatoSim outputs a BAM file containing simulated variants."

**VarBen**
- Paper: Li Z, Fang S, Zhang R, et al. J Mol Diagn 2021;23(3):285–299, doi:10.1016/j.jmoldx.2020.11.010. Full text blocked (403), so I used the abstract, README, `docs/VarBen_manual.pdf` and code (`varben/deal_sv/dealReadsType.py`).
- 3, abstract: "directly editing the original sequencing reads". Manual: "remap edited reads to reference genome".
- 4, code: `qual = read.query_qualities; read.query_sequence = new_seq; read.query_qualities = qual`.
- 9, manual: "output directory name for edited bam file".

**SVEngine**
- Paper: Xia LC, Ai D, Lee H, Andor N, Li C, Zhang NR, Ji HP. GigaScience 2018, giy081, doi:10.1093/gigascience/giy081 (PMC6057526).
- Repo: bitbucket.org/charade/svengine (`README.rst`, `doc/Manual.rst`, `svengine/mf/mutforge.py`, `xwgsim/xwgsim.c`).
- 3, paper: "spikes structural variants into the contigs, samples short reads from the altered contigs… performs the alignment."
- 4, xwgsim.c: `Q = (ERR_RATE == 0.0)? 'I' : …` and `qstr[i] = Q`.
- 9, Manual: `-x {fasta,fastq,bam} … (default: fastq)`.

**NEAT**
- Papers: Stephens ZD, Hudson ME, Mainzer LS, Taschuk M, Weber MR, Iyer RK. PLoS One 2016;11(11):e0167047, doi:10.1371/journal.pone.0167047. Allen JM, Gandhi KR, Alhazmy R, Wasnik Y, Fliege CE. JOSS 2026;11(121):9056, doi:10.21105/joss.09056. The PDF title is "Enhancing short-read sequencing simulation: Updates to NEAT"; the README cites it as "next-generation".
- Repo: github.com/ncsa/NEAT (README, `neat/read_simulator/utils/vcf_func.py`, `neat/variants/unknown_variant.py`).
- 2: `UnknownVariant` is "a placeholder type of variant… used for input variants of unsupported type". JOSS Table 1: "Two variant types supported" → "Framework to expand variant types".
- 3, README: "FASTQ files with reads sampled from a provided reference genome".
- 4, README: "Leave empty to use default model (default model based on human, sequenced by Illumina)"; `neat model-qual-score -i input_reads.fastq`.
- 9, README: "FASTQ, BAM, and VCF".

**ART**
- Paper: Huang W, Li L, Myers JR, Marth GT. Bioinformatics 2012;28(4):593–594, doi:10.1093/bioinformatics/btr708. Official page: niehs.nih.gov/research/resources/software/biostatistics/art.
- 4, abstract: "built-in, technology-specific read error models and base quality value profiles parameterized empirically… customized… quality profiles."
- 9, page: "ART outputs reads in the FASTQ format, and alignments in the ALN format… SAM".

**VarSim**
- Paper: Mu JC, Mohiyuddin M, Li J, Bani Asadi N, Gerstein MB, Abyzov A, Wong WH, Lam HYK. Bioinformatics 2015;31(9):1469–1471, doi:10.1093/bioinformatics/btu828 (PMC4410653).
- Docs: bioinform.github.io/varsim. Release notes: GitHub v0.7.0.
- 3/4, paper: "supports DWGSIM and ART… uses ART as the default since ART learns an error profile based on real sequencing reads."
- 9, docs: "reads… out/lane*.fq.gz and out/simu.truth.vcf".

**VISOR**
- Paper: Bolognini D, Sanders A, Korbel JO, Magi A, Benes V, Rausch T. Bioinformatics 2020;36(4):1267–1269, doi:10.1093/bioinformatics/btz719.
- Docs: davidebolo1993.github.io/visordoc. Code: `VISOR/VISOR.py`, `VISOR/SHORtS/SHORtS.py`; pywgsim `src/wgsim_mod.c` line 288 (same fixed-quality rule as wgsim).
- 3, abstract: "SVs are implanted into FASTA haplotypes… reads are drawn at random from these haplotypes using standard error profiles."
- 9, VISOR.py: `--fastq` "store synthetic read pairs in FASTQ format"; docs: "SHORtS and LASeR store in the output folder a sorted BAM."

## Published comparisons

- **Duncavage et al., J Mol Diagn 2023;25(1):3–16 (AMP/API/CAP).** The full text was blocked (403/Cloudflare), so its statements about specific tools are not found. From the AMP press release (amp.org, Oct 18 2022): in silico files "cannot supplant the use of physical samples". Labs should know the limits "for… variant types susceptible to systematic sequencing and mapping errors".
- **Milhaven & Pfeifer, Heredity 2023 (PMC9905089).** Independent test of ART, DWGSIM, InSilicoSeq, Mason, NEAT and wgsim against real data: "quality scores of reads simulated with DWGSIM and wgsim were a poor representation of the empirical data". This matters for the slide because the SV reads of BAMSurgeon and SVEngine, and VISOR's short reads, all come from wgsim.
- **Lee et al. 2018 (BAMSurgeon's own authors, so not independent).** "the simulated reads do not necessarily reflect the non-uniform coverage" of real samples.
- **Claims by competing tools' authors (not independent).** The SVEngine paper says BAMSurgeon's input BAM "has to be large (typically >30x) in order to successfully assemble local contigs". The SomatoSim paper says reference-only simulation "fails to capture the instrument- and experiment-specific sequencing error".