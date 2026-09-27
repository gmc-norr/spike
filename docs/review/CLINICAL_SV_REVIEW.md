**Assessment: Spike has a useful haplotype representation, but the current simulator is not ready to supply trusted germline WGS SV truth sets without substantial restrictions and additional validation.** Several problems affect emitted molecules and dosage, even on a unique, error-free reference. I would use it for exploratory caller tests today, and develop it into a supplementary clinical validation tool after addressing the findings below.

Reviewed 2026-09-25, source commit `4efa0f4`. Intended use: paired-end short-read WGS, germline structural variants. This is a new assessment of the current code, including the fixes recorded in `REVIEW.md`; its historical findings are not automatically current findings. No simulator implementation was changed during this review.

The central distinction is between **representing a rearranged sequence correctly** and **sampling a realistic sequencing library from that sequence in a real diploid background**. The segment model largely solves the first problem for simple events. The remaining issues are concentrated in the second problem, event composition, and the truth/validation contract.

**Evidence and limits.** I inspected the simulation, haplotype, synthesis, extraction, phasing, input/output, validation, and validation-harness code; built the current binary; ran its tests and Clippy; and constructed independent synthetic BAM probes. The probes use a 40 kb unique random reference, 150 bp paired reads, 400 bp fragments, Q60, a baseline of 75× base depth, and Spike seed 17. They intentionally remove biological and mapping ambiguity so that the model can be tested directly. Four cases were also aligned with `bwa-mem2` and merged using Spike's own generated scripts; `samtools depth` reproduced the sequence-based measurements. These are simulator correctness experiments, not measurements of clinical caller sensitivity. I did not run a new human whole-genome or multi-caller clinical benchmark.

Reproduce the measurements with [scripts/review_sv_model.py](../../scripts/review_sv_model.py):

```bash
cargo build --offline
python3 scripts/review_sv_model.py --output /tmp/spike-review-new
# Also check the emitted sequences through alignment and merge:
python3 scripts/review_sv_model.py --with-alignment --output /tmp/spike-review-aligned-new
```

The output directory must be new. Requirements: Python 3.9+, `samtools`, and optionally `bwa-mem2`, on PATH. The script writes the toy inputs, run logs, truth, FASTQs, `summary.json`, and optionally `alignment-summary.json`. Its exit status reports execution success, not whether these known model defects are present. Original review artifacts are under `/tmp/spike-clinical-review/`; the script is the durable reproduction.

**Current build health.** `cargo build --offline` succeeds, with one unused-method warning. The suite has **421 passing tests and one ignored test** when `bcftools` is on PATH. Without it, one gVCF warning test fails: 420 pass, one fails, one is ignored. Clippy completes with 13 binary warnings and 14 test-target warnings, mostly style. The missing test dependency should be declared explicitly. A passing unit suite currently does not establish the properties measured below.

| ID | Priority | Finding | Evidence |
| --- | --- | --- | --- |
| R1 | High | Nearby, non-overlapping events restore each other's deleted sequence | Measured in FASTQ and aligned/merged BAM |
| R2 | High | One depth estimate flattens donor coverage and distorts dosage | Measured in FASTQ and aligned/merged BAM |
| R3 | High | Synthetic haplotypes erase background indels | Measured with a homozygous background deletion |
| R4 | High for difficult loci | Filtered donor molecules remain resistant to the event | Measured with mixed-MAPQ input and aligned/merged BAM |
| R5 | High for long INS | Exhausted placement retries admit novel-only fragments into a reference-overlap budget | Measured over four insertion sizes |
| R6 | High for translocations | Additive fusion evidence does not represent a balanced germline rearrangement | Code and molecule-count analysis |
| R7 | High for truth integrity | Genotypes, ploidy, and inserted sequence are not faithfully represented in truth | VCF genotype round trip measured; other paths inspected |
| R8 | Medium | Mate recovery discards unmatched R1 before the recovery pass | Code plus BAM boundary-pair reproduction |
| R9 | High for interpreting a benchmark | Current QC and harness results cannot establish SV sequence/genotype correctness or clinical precision | Code plus measured QC blind spots |

**R1. Nearby non-overlapping events can cancel each other's biological effect.** Relevant code: [event overlap check](../../src/main.rs#L626), [per-event replacement](../../src/simulate.rs#L166), and [output combination](../../src/simulate.rs#L301).

The overlap check considers event intervals, but each event replaces reads across a larger haplotype footprint with 2 kb flanks. Events are simulated independently against the original donor. Combining outputs removes originals suppressed by either event, but unconditionally retains the synthetic reads from both. Synthetic flanks from event B can therefore contain reference sequence that event A deleted. The previous H1 fix removes a different problem—retained *originals* undoing suppression—and does not solve this one.

Reproduction: `del:chrT:10000-11000` and `del:chrT:12000-13000`, without `--allow-overlap`. These events are accepted by default. Depth was measured inside each deleted interval, away from the junctions.

| Requested events | First deletion interior | Second deletion interior |
| --- | ---: | ---: |
| First deletion only, AF=1 | 0.00× | Not deleted |
| Both deletions, AF=1 | **70.21×** | **74.34×** |
| First deletion only, AF=0.5 | 34.59× | Not deleted |
| Both deletions, AF=0.5 | **52.27×** | **53.33×** |

Baseline depth is 75×. A homozygous deletion should not recover essentially normal depth because another deletion was added nearby. For unclassified copies, independent suppression plus overlapping synthetic coverage can also inflate otherwise unchanged sequence: two AF=0.5 replacements have approximately `(1-0.5)^2 + 0.5 + 0.5 = 1.25` times the original depth in a shared unchanged interior.

**Action:** immediately reject intersecting replacement footprints, including a fragment margin, unless a shared simulation handles them. The complete solution is to group connected events, assign them to explicit haplotypes, and generate/suppress molecules once per group. Simply dropping synthetic reads by name cannot compose their sequence or phase correctly.

**R2. Uniform tiling changes the donor's spatial coverage profile.** Relevant code: [first covered breakpoint selection](../../src/simulate.rs#L480), [fragment count](../../src/simulate.rs#L548), and [uniform placement](../../src/simulate.rs#L725).

Spike measures fragment depth in a 2 kb window at the first covered breakpoint side, suppresses local originals, and uses that one estimate across the entire variant haplotype. This replaces real variation in coverage with uniform synthetic coverage. It matters for WGS too: sequence composition, library preparation, and mapping produce spatial variation even without capture.

Reproduction: a heterozygous tandem duplication of `chrT:10000-28000`. The donor has 75× depth near the first breakpoint and an interior section at 18.75×. That low section becomes **81.09×**, a **4.32×** increase. A locally proportional CN2→CN3 change should produce about **28.13×**, a 1.5× increase. The high-depth section becomes 110.92×, close to the expected 112.5×. With a uniformly covered donor, both sections instead behave approximately as expected.

The approximate rule implemented inside a duplication is `D_out(x) = (1-v) D_donor(x) + 2v D_anchor`, rather than `(1+v) D_donor(x)`. Sampling noise explains small deviations from these expectations, not the large low-depth distortion. Inversions and replacement flanks can also have their profile changed. Zero-depth sides are sometimes filled using depth from the other side; the existing warning documents that behavior but does not make it realistic.

**Action:** learn a spatial fragment-start intensity, preferably with library and sequence-context terms, and project that process onto the altered haplotypes. Avoid simply treating low aligned depth as a molecular sampling bias: some of it is mappability and should arise naturally when simulated molecules are aligned. Establish dosage preservation in well-mappable control windows and characterize difficult regions separately. Add local, binned donor/output comparisons instead of relying on one event-average ratio.

**R3. Background SNPs are preserved, but background indels and existing rearrangements are not.** Relevant code: [sample copies](../../src/loh.rs#L34), [SNP-only gVCF parsing](../../src/loh.rs#L513), and [application to reference-derived haplotypes](../../src/simulate.rs#L156).

`SampleCopies` stores a base per reference position. The gVCF path explicitly excludes non-SNP alleles; its bookkeeping for spanning deletions does not install those deletions into the synthetic sequence. Thus the simulator preserves many SNP alleles while replacing the sample's indel/SV sequence with reference sequence. The pileup path cannot repair that representation limitation.

Reproduction: the donor has a homozygous 2 bp deletion at reference positions `[18000,18002)`. Spike adds a heterozygous tandem duplication spanning it. All three resulting genomic copies should still carry that background deletion. Instead, sequence-specific read counting finds **32 deletion-supporting reads and 57 reference-supporting reads: AF=0.360**, down from 1.0. Retained donor molecules carry the deletion; the two synthetic copies resurrect the reference bases.

This can alter assembly, repeat context, phasing, and apparent allele balance near an SV. It can also erase or combine an already-present background event unpredictably. It is especially relevant when moving events between samples: absence from a caller's output is not proof that the donor lacks the event.

**Action:** support sample-specific haplotypes containing SNPs, indels, and existing SVs, ideally from phased calls or a suitable assembly. Until then, restrict trusted tests to independently characterized backgrounds and explicitly reject or label footprints containing unsupported variation. Fail closed in a strict mode when the requested sample-variant input cannot be read; [the current error handler](../../src/simulate.rs#L93) logs a warning and proceeds with empty copies.

**R4. The donor-training filter also determines which molecules can be edited.** Relevant code: [proper-pair requirement](../../src/extract.rs#L104), [MAPQ/flag filters](../../src/extract.rs#L777), [replacement names](../../src/main.rs#L553), and [merge operation](../../src/main.rs#L1456).

Only qualifying proper pairs enter the donor pool. Low-MAPQ pairs, discordant pairs, unmapped mates, and excluded flagged records are generally absent from the replacement-name list and remain in the original BAM after merging. Those molecules cannot undergo the event. Filtering reads to learn a reliable library model is reasonable; using the same selection as the complete editable population is a separate and consequential assumption.

Reproduction: half the donor pairs have MAPQ 60 and half MAPQ 0. A requested AF=1 deletion leaves **37.5× inside the deletion**, entirely from the MAPQ-0 half of the original 75× input. This is not a whole-sample homozygous deletion. `spike validate` nevertheless passes the event's coverage check at observed ratio 0.00 because its default MAPQ filter also hides those retained molecules.

Some clinical callers deliberately ignore these reads, which limits the impact for those specific tests. Others use low-MAPQ, unmapped, or discordant evidence in assembly or SV discovery. Lowering `--min-mapq` addresses only part of the selection.

**Action:** separate model-training eligibility from replacement eligibility. Track all affected molecules, their mates, and supplementary/secondary representations; decide consistently how the altered genome changes them. Retain an explicit exclusion census and quantify event-resistant depth. Avoid blindly deleting unrelated multimappers: assignment uncertainty should be modeled or the locus declared unsupported.

**R5. Long insertion placement breaks its own reference-overlap constraint.** Relevant code: [novel-only start exclusion from the count](../../src/simulate.rs#L594), [bounded redraw loop](../../src/simulate.rs#L729), and [nearest-reference coordinate fallback](../../src/synth.rs#L813).

The count formula excludes starts wholly inside inserted sequence. The placement loop redraws at most ten times, then accepts its last start even if it still lies wholly inside the insertion. The generator accepts those starts because both ends can use a nearest-reference fallback. This consumes a budget intended for reference-overlapping fragments with novel-only reads.

| Random insertion length | Synthetic pairs emitted | Pairs with a detectable reference 31-mer in either mate |
| --- | ---: | ---: |
| 1 kb | 500 | 498 |
| 10 kb | 500 | 485 |
| 100 kb | 500 | 171 |
| 1 Mb | 500 | 17 |

The 31-mer measure is a conservative anchor diagnostic, not an exact count of all reference-overlapping fragments: a fragment with only a very short reference overlap can lack such a seed. The severe size-dependent loss matches the explicit retry bug. The same run still suppresses reference molecules and writes the requested AF to truth.

**Action:** sample valid start intervals directly, with probabilities proportional to interval lengths for the sampled fragment size. Do not accept an invalid placement on retry exhaustion. Decide separately whether to simulate novel-only molecules: they are legitimate molecules in real WGS and can matter to assembly, but require their own count and truth provenance. Excluding them deliberately is itself a limitation of an insertion benchmark.

**R6. Fusion mode is a junction-evidence spike, not a balanced germline translocation model.** Relevant code: [additive choice](../../src/simulate.rs#L138), [additive count](../../src/simulate.rs#L571), and [fusion haplotype](../../src/haplotype.rs#L343).

A fusion keeps every original pair and adds only fragments crossing one new adjacency. A balanced heterozygous reciprocal translocation should alter one copy at each participating chromosome, retain the other copies, and produce both derivative adjacencies while conserving copy number. Two BND mate records describe the two ends of **one adjacency**; they are not the reciprocal derivative chromosome.

With equally covered partners and AF=0.5, the additive formula produces `C` junction fragments on top of the `C` unmodified reference fragments, whereas an ordinary heterozygous replacement would allocate roughly half the molecules to each allele. This can make local breakpoint evidence much stronger than in a copy-neutral germline sample. Even providing the reciprocal adjacency as a second event leaves the additive retention issue.

**Action:** provide explicit derivative-chromosome paths, dosage, and phase, with suppression/replacement across both partners. Label the current mode as additive junction evidence. It remains useful for testing whether a caller recognizes a specified join, but does not establish sensitivity or genotyping performance for balanced translocations. The legacy junction DUP mode shares the additive dosage problem; use the full tandem model for simple DUP work after fixing R1–R4.

**R7. Truth lacks information required for germline genotype and sequence validation.** Relevant code: [VCF input fields retained](../../src/vcf_input.rs#L261), [genotype inferred from AF](../../src/truth.rs#L321), [INS output](../../src/truth.rs#L252), [random INS construction](../../src/main.rs#L1244), and [`af=het`](../../src/main.rs#L423).

- **Input GT is not the simulated genotype.** A measured input DEL with `GT=1/1` and no simulation-specific AF produces `SIM_VAF=0.500; GT=0/1`. FILTER and homozygous-reference records are also not selection gates, though those cases have warning counters. A specification-only VCF mode is legitimate, but must not be confused with reproducing a sample's truth genotypes.
- **No ploidy or absolute copy-number model.** Both the simulation and emitted GT assume two starting copies. A male non-PAR X/Y event cannot have its haploid truth represented properly; baseline CNVs and CN>4 gains also lack an explicit model. Setting AF=1 can obtain complete deletion of eligible reads, but still writes diploid `1/1`.
- **`af=het` changes the underlying event fraction.** Beta(40,40) is applied to copy suppression and generation, not just observation noise. Values above 0.5 cause some contribution from the second copy. This creates a dosage mixture instead of a fixed diploid heterozygote with stochastic sampling. For ordinary germline probes, use an explicit 0.5 until genotype and sampling bias are separated.
- **Insertion sequence is lost from truth.** INS records contain `<INS>` and SVLEN, without the supplied/generated allele sequence. Random sequence is created as a local value during haplotype building and is not retained in the event. Feeding that truth VCF back cannot recover the original insertion sequence. Position/length alone cannot validate an assembled insertion or discriminate unrelated equal-length alleles.
- **Requested AF is presented as truth despite caps/floors.** Additive events cap AF at 0.95 while truth uses the uncapped request; a minimum of two synthetic fragments can also raise support at low depth. Warnings do not update the machine-readable truth. This is less central at normal germline depth but makes titration and robustness studies misleading.

**Action:** store input sample, GT, ploidy, phase, baseline/output CN, actual insertion sequence, derivative adjacencies, and requested molecule fraction separately. Record emitted molecule counts and post-alignment measurements separately from the intended genome. Export sequence-resolved alleles or a referenced allele FASTA plus hashes. Make specification mode versus genotype-reproduction mode explicit.

**R8. Unmatched R1 records are removed before mate recovery.** Relevant code: [BAM pass-one pairing](../../src/extract.rs#L135) and the corresponding CRAM loop around line 384.

The expression `(read1_map.remove(&name), read2_map.remove(&name))` removes R1 even when R2 is absent. That unmatched R1 is consequently unavailable to the second pass. Unmatched R2 can remain, producing an asymmetric recovery behavior.

Measured boundary pair `p003200` has R1 at 12900 and R2 at 13150, both proper and MAPQ 60. Extraction of `[8000,13000)` sees R1, and the widened second query covers R2, but the pair never appears in the extracted/replaced set. It remains in the original after merge, so this example is not a whole-BAM read-loss claim. It does contradict the documented mate recovery, truncates standalone FASTQ output at some boundaries, and can bias small donor pools.

**Action:** check that both maps contain the name before removing either entry. Test both mate orders on BAM and CRAM, including a mate outside the initial query.

**R9. QC passes are weaker than the truth claims being made.** Relevant code: [coverage tolerance](../../src/validate.rs#L587), [split-read checks](../../src/validate.rs#L636), [insertion evidence](../../src/validate.rs#L698), [SA parsing](../../src/validate.rs#L2235), [harness comparison](../../scripts/validate_pipeline.sh#L292), and [benchmark invocation](../../scripts/validate_pipeline.sh#L849).

The current checks have useful diagnostic value, but several are evidence-presence checks rather than event validation:

- A DUP/DEL event-average depth ratio can conceal large local errors. R2's 4.32× interior amplification passes the overall DUP ratio check at 1.32 versus expected 1.50. R4's resistant low-MAPQ depth is invisible to the default check.
- Split evidence counts SA alignments near the other locus, without verifying strand, the aligned breakpoint implied by CIGAR, sequence, or allele fraction. For an inversion, this is not proof of both correct junctions. Small deletions represented only by a CIGAR `D` may lack SA tags and fail despite valid sequence; this is already acknowledged in the earlier review.
- INS evidence accepts nearby insertion/soft-clipping operations above a length threshold; it does not establish the identity of the inserted sequence. The relevant unresolved N15 limitation remains present.
- The synthetic probes' global validator exits are nonzero because fragments have SD=0 and no duplicate-marking step was run. The report above claims specific event checks passed, **not** that the entire validator returned success. Fixed thresholds such as an insert-size SD range are library heuristics, not a comparison to the measured donor.
- The end-to-end harness is largely a deletion/Delly test with a gain-over-background gate. It cannot validate INS/DUP/INV/BND performance by inference. A TP-count increase does not identify which newly recovered events are true additions versus changes caused by realignment, and does not establish genotype accuracy.
- Spike's truth contains added events only. Calls matching real background variants are not necessarily false positives. Whole-callset precision against this partial truth is not a clinical precision estimate. Conversely, a negative control with no call does not prove a locus is variant-free.

**Action:** distinguish genome truth, molecular truth, alignment evidence, and caller output. Validate both junctions and orientations where applicable, sequence identity, genotype/CN, local depth and B-allele balance, and preservation of known background variants. Match outcomes by event identity. Measure clinical precision only against sufficiently complete truth within an explicit confident/evaluable region set, with separate accounting for background and newly introduced calls. The current harness should be described as an integration/regression test.

**How realistic is the underlying model?** The ordered-segment representation should be retained. It naturally produces appropriate sequence for a simple deletion, tandem duplication, inversion, or specified insertion, and actual alignment can produce split reads and discordant pairs without hard-coding their CIGARs. Empirical fragment lengths, separately learned R1/R2 quality profiles, random mate orientation, and incorporation of phased SNPs are useful foundations. The full tandem duplication layout is preferable to independent depth/junction additions.

Its current realism is uneven:

| Aspect | Assessment for germline WGS |
| --- | --- |
| Isolated, simple SV sequence on a known background | Sound representation; sampling still requires the fixes above |
| Local CN and depth | Reasonable on uniform high-quality donors; unreliable on heterogeneous coverage |
| Background alleles | SNP-aware, not a complete sample haplotype |
| Simple tandem duplication | Represented; dispersed/inverted duplications and higher CN states need explicit paths |
| Balanced translocation | Not represented by the additive fusion model |
| Inserted sequence | Explicit sequence accepted, but random A/C/G/T is not a biological model for mobile elements/repeat insertions; truth loses sequence |
| Breakpoint complexity | No general event specification for microhomology, junction insertions, templated sequence, or multiple linked rearrangements |
| Sequencing errors | Learned quality distributions, but errors assume perfectly calibrated Phred scores and largely independent, uniformly chosen substitutions; indels are optional and context-poor |
| Library structure | Fixed mean read length, one pooled model; no explicit PCR/optical duplicate families, lane/library heterogeneity, or adapter/read-through process |
| Hard loci | Real alignment helps, but filtering, reconstructed reference haplotypes, and uniform sampling change precisely the evidence that makes those loci difficult |

Quality-score realism is not equivalent to error realism. Matching mean Q or per-cycle Q does not establish a matched substitution spectrum, repeat-dependent indel errors, pair-level quality correlation, clipping, duplicate family sizes, or empirical MAPQ distribution. I would calibrate these against held-out real libraries and real SV carriers before introducing more complicated quality models.

There are also useful engineering improvements: bound memory and stream output for large event sets; use quality histograms instead of retaining multiple copies of observations; avoid rescanning an entire plain gVCF for every event; preserve per-library/read-group context through FASTQ reconstruction; and publish a reproducible environment with declared external test dependencies. Fragment normalization should use the actual sampled distribution: [statistics](../../src/stats.rs#L27) accept lengths below 10 kb, whereas [generation](../../src/simulate.rs#L697) conditions them on the read length and a 1500 bp maximum. That mismatch is code-confirmed; I did not quantify its impact on an ordinary WGS library here.

**Recommended development order.** First make the simulator's genome and dosage claims correct; improvements to superficial read quality should follow those invariants.

1. **Prevent known invalid truth sets:** reject overlapping replacement footprints; fix mate recovery and insertion placement; preserve insertion sequence; require explicit germline GT/ploidy; reject or explicitly exclude unsupported backgrounds and failed sample-haplotype loading.
2. **Introduce a shared sample-genome model:** represent the starting haplotypes and their CN, compose variants onto selected copies, build all derivative paths, and apply each affected molecule once. Keep a distinction between the high-quality training subset and the complete molecule population.
3. **Make sampling donor-aware:** model local start intensity and library strata, use the actual fragment/read-length distributions, preserve or deliberately model duplicate families, and allow unbiased sampling that can produce zero supporting reads. Record why any placement was rejected.
4. **Add independent invariant tests:** expected CN2→CN1 and CN2→CN3 ratios; homozygous deletion with no residual assigned-copy molecules; copy-neutral inversions/translocations; unchanged hom-alt SNPs and indels; phase consistency; composed nearby events; reference/novel read budgets across insertion lengths; BAM/CRAM parity; low-MAPQ/discordant molecules; sex-chromosome ploidy; input/output allele and genotype round trips. Include alignment-based held-out cases alongside unit tests.
5. **Create durable provenance:** source commit/build identifier, seed and resolved parameters, reference and input hashes, sample/library identity, allele sequences, edits and exclusions, pre/post depth, sequence support, and truth completeness/evaluable regions. The current run README is useful but insufficient as a complete validation manifest. A truth record should not imply successful molecule generation merely because its event specification parsed.

**How I would use it in a clinical validation study.** Establish real-sample performance first, then use corrected Spike to add targeted challenges and regression coverage. This is consistent with the AMP/API/CAP recommendations: in-silico material supplements physical samples, and systematic sequencing/mapping limitations need explicit consideration. [AMP's primary summary of the recommendations](https://www.amp.org/AMP/assets/File/pressreleases/2022/AMP_In_Silico_NGS_Pipeline_Validation_Press_Release101722.pdf), [joint report, 2023](https://doi.org/10.1016/j.jmoldx.2022.09.007).

Use the lab's actual library preparation, reference, alignment, duplicate handling, caller ensemble, filtering, genotyping, and reporting workflow. Run real positive controls with orthogonal support appropriate to the variant class, including DUP/INV/translocations that a deletion-heavy benchmark cannot establish. Use multiple libraries/samples as well as multiple simulation seeds: seeds alone do not test systematic assay effects.

Pin the exact GIAB release, reference build, and its matching confident/evaluation regions. NIST now lists **HG002 v5.0q** for small variants and SVs, based on the T2T-HG002 v1.1 assembly, and marks several older HG002 benchmarks deprecated. Follow the selected release's README and exclusions; the older DEL/INS benchmark should not be treated as comprehensive truth for all SV classes. I verified NIST's current listing, but did not audit or run the v5.0q release files in this review. [NIST GIAB release guidance](https://www.nist.gov/programs-projects/genome-bottle). The older benchmark's isolated DEL/INS scope is documented in its [original NIST publication](https://www.nist.gov/publications/robust-benchmark-detection-germline-large-deletions-and-insertions).

Define a challenge matrix before tuning: SV class; size relative to read and fragment length; tandem versus dispersed/inverted duplication; repeat and segmental-duplication context; medically relevant loci; GC and coverage strata; het/hom/haploid status; nearby background variation; balanced versus unbalanced rearrangements; and batch/library effects. Sequence-resolved real alleles and phased backgrounds are much more informative than uniformly random insertion strings.

For each case, retain an unmodified control and, where feasible, a control subjected to equivalent extraction/realignment without an edit. Process both through the same production steps. Report detection sensitivity, genotype/CN accuracy, breakpoint/sequence accuracy, and false-positive burden separately, stratified by the matrix and with uncertainty intervals. Establish acceptance criteria prospectively; neither a convenient number of simulations nor one overall recall threshold establishes performance in every reportable class. More simulated positives cannot by themselves establish specificity.

For immediate experiments with current Spike, limit interpretation to individually inspected, isolated simple autosomal events in well-characterized high-quality backgrounds; use explicit AF=0.5 for diploid heterozygotes and the full DUP model. Inspect the merged BAM, local dosage, background alleles, both junctions where relevant, and evidence after the production duplicate-handling step. The confirmed model defects above preclude treating an unqualified Spike truth VCF or aggregate caller recall as clinical validation evidence.
