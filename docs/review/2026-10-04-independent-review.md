**Independent review of spike — 4 October 2026**

Reviewed commit `2764dbc`. Production code was not changed. The pre-existing `.gitignore` modification and untracked `.ignore` were left alone.

The segment-based haplotype representation is a good foundation. For isolated events in reasonably uniform, diploid short-read data, the overall approach is understandable and useful. However, the current implementation has reproducible correctness bugs, and its output should not be treated as a generally faithful sample-genome simulation. In particular, an emitted truth record is not proof that the requested allele fraction reached the reads.

This review distinguishes newly reproduced failures from limitations already recorded in `REVIEW.md` and `CLINICAL_SV_DESIGN_NOTES.md`. Historical measurements below are explicitly identified; they were not rerun on the original GIAB datasets.

**Evidence and scope**

Reviewed the Rust simulation, extraction, sample-copy, reference, truth, validation, output-script and input-parsing paths, alongside the existing design/review records and Python validation helpers. Ran:

- `cargo build --locked --offline`: passed.
- `cargo test --locked --offline`: **672 passed, 2 failed, 1 ignored**. Both failures require `bcftools`, which is absent from PATH; this dependency is already documented. These failures do not establish a product regression.
- `python3 -B -m pytest -q scripts/duplicates scripts/transplant`: **75 passed**.
- `cargo clippy --locked --offline --all-targets`: completed, with 12 production warnings and 2 additional test warnings. Mostly style/API-size issues; not the main risk.
- Independent executable probes using a random 12 kb reference, indexed BAMs constructed with samtools, explicit genotypes, exact sequence probes, and generated FASTQ scripts.

Reproduce the independent probes with:

```bash
cargo build --locked --offline
python3 scripts/review_20261004.py
```

The script requires samtools and creates a fresh temporary directory. It prints that directory and writes its inputs, outputs and JSON results there. Its destructive-output test acts only on a disposable generated FASTQ. The original review run's artifacts are under `/tmp/spike-review-20261004`.

No real aligner/caller benchmark was run in this review: bwa-mem2, minimap2 and bcftools were unavailable on PATH. The deletion-validation probe uses exact, constructed alignments, so it tests the validator independently of an aligner's representation choices. Checks used the debug build; release-mode tests were not run.

**Newly reproduced findings, ordered by priority**

**1. P1 — The full-FASTQ route can destroy a raw input.**

Locations: `src/main.rs:1369`, `src/main.rs:2054`, `src/main.rs:2115`.

`validate_into_fastq` checks only whether raw paths exist. The generated `fastq.sh` streams directly into the requested output path. If an output aliases either raw input, shell redirection truncates that input before decompression. The failure trap then deletes both output paths, including the aliased raw input. It can also delete a pre-existing second output that this run never reached.

Measured: passing a disposable `alias_R1.fq.gz` as both RAW_R1 and OUT_R1 produced exit 1, a decompression error, and **the raw file no longer existed**. This also applies to CLI runs when the raw path is `<output>/<prefix>_R1.fastq.gz`. Different strings can alias through symlinks or hard links; string comparison alone is insufficient.

Fix: reject collisions between every input and output, and between outputs; write both mates into temporary files beside their destinations, validate them, and publish only on success. Cleanup must remove only temporary files owned by this run.

**2. P1 — `origin` removal and replacement use different probabilities.**

Locations: `src/origin.rs:231`, `src/origin.rs:246`, `src/origin.rs:392`, `src/origin.rs:459`.

Removal uses the most confident mate's probability for the entire fragment. Replacement depth sums each mate's independent placement probability. For a uniquely anchored pair with one MAPQ-60 mate and one MAPQ-0 mate listing a single alternative, removal is approximately 1, while the average weight used for replacement is approximately `(1 + 0.5) / 2 = 0.75`.

Measured at AF 1 for a SNP: **149 original pairs removed, 111 synthetic pairs added**, versus **149 removed, 148 added** with confident placements. This creates loss of coverage across an otherwise copy-neutral replacement footprint. The committed probe compares `clean` and `origin` on the same mixed-confidence BAM at `--min-mapq 0`. An additional check using `origin` with the default MAPQ and flank settings, retaining high-quality donors outside the footprint, reproduced the same 148-to-111 difference.

This is an internal conservation problem, separate from whether the MAPQ/XA heuristic is well calibrated. Fix: derive both removal and replacement from the same joint fragment-placement weights. Test a unique mate paired with an ambiguous mate, both mates ambiguous, alternative placements on different contigs, and duplicate-family members with different confidence.

**3. P1 — An allele already present in the donor can contradict the requested truth AF.**

Locations: `src/simulate.rs:203`, `src/simulate.rs:249`, `src/main.rs:1530`.

The REF check validates against the reference FASTA, not against the donor's genotype. Suppression retains some original reads and synthesis adds the requested ALT. When the donor already carries ALT, retained originals also support it. There is no baseline-allele check or adjustment.

Measured: a donor homozygous for `chr1:5001 T>A`, requested again at `af=0.5`, exited successfully and wrote `SIM_VAF=0.500;GT=0/1`. Exact 25-base probes in the output found **25 ALT-supporting reads and zero REF-supporting reads**. The result remained homozygous.

This matters particularly when selecting events from real variant catalogs. Fix: detect existing event support and either reject the site with an explicit reason, or define and implement a target-final-genotype/AF operation. Arbitrarily adding the same allele cannot establish the requested final fraction.

**4. P2 — A gVCF containing only hom-alt calls loses them during fallback.**

Location: `src/loh.rs:182`.

`sample_copies` discards the entire parsed `RegionSnps` whenever `het` is empty, including valid `hom_alt` entries. Pileup then replaces the discarded calls. Its minimum depth of 10 means a known hom-alt allele can disappear in a shallow region; the same issue arises when the relevant reads fail the pileup's MAPQ filter.

Measured: at an 8-read-depth hom-alt SNP in the replacement flank, providing a hom-alt-only gVCF produced **5 reference-bearing synthetic reads and zero ALT-bearing synthetic reads**. Adding an unrelated het record to the same gVCF, solely to exercise the other branch, produced **zero reference-bearing and 7 ALT-bearing synthetic reads** at that hom-alt site.

Fix: preserve usable gVCF hom-alt calls independently of whether het calls exist. If fallback fills missing information, merge it under an explicit precedence rule. Current parser-level tests that hom-alt calls are retained do not cover their subsequent loss in `sample_copies`.

**5. P2 — The coverage validator treats deleted bases as covered.**

Locations: `src/validate.rs:3016`, `src/validate.rs:3025`, `src/validate.rs:3037`.

`ref_span` includes CIGAR `D` and `N`, then the depth counter treats the entire resulting interval as covered. This is reference span, not sequenced-base depth. A read spanning a true deletion therefore contributes apparent depth to the very bases it deletes.

Measured: an exactly constructed, homozygous 50 bp deletion had **zero aligned-base coverage at every deleted position**, verified with `samtools depth -aa`. Nevertheless, `coverage_ratio` reported **1.02**, against expected **0.00**, and failed. `del_planted` independently found 28 supporting reads and passed. The fixture also lacks duplicate marking, so its global duplicate row failed separately; the finding concerns the incorrect coverage row itself.

Fix: accumulate aligned `M`, `=` and `X` blocks for a base-depth check. If a deletion-inclusive span metric is wanted, give it its own name and a compatible expectation. Samtools explicitly distinguishes the two semantics through the `depth -J` option: [samtools depth documentation](https://www.htslib.org/doc/samtools-depth.html).

**6. P2 — Different R1/R2 cycle counts are silently collapsed.**

Locations: `src/bam_stats.rs:224`, `src/synth.rs:489`.

The statistics scan merges both mates into one length histogram; ties select the longer length. `SimConfig` and synthesis then carry one cycle count for both mates.

Measured: an ordinary paired BAM with **100 bp R1 and 150 bp R2** produced synthetic **150 bp R1 and 150 bp R2**. For the added R1 cycles, there are no training observations; quality sampling reaches its Q20 fallback. Thus both sequence length and the quality tail are invented, while original R1 reads remain 100 bp.

Fix: model cycles and trimming separately for each mate, preferably per library/read group. If this input shape is outside the supported scope, reject it explicitly instead of silently changing it.

**7. P2 — Full-FASTQ success does not establish that the mates pair up.**

Locations: `src/main.rs:2069`, `src/main.rs:2112`, `src/main.rs:2128`.

The script processes raw mates independently. It checks total removed records and the number of added synthetic records, but never compares retained raw read names or total raw record counts. It also assumes four-line FASTQ records without checking sequence/quality lengths or incomplete records.

Measured: RAW_R1 names `[a,b]` and RAW_R2 names `[b,a]`, with an empty removal list and one synthetic pair, returned **exit 0** and preserved the mismatched order. An empty removal list is a normal possibility for additive events. Downstream paired alignment receives incorrectly paired reads despite the success message.

Fix: validate raw mates together as a stream, checking normalized names, record structure and lengths. Validate the specific removed-name set, not just its cardinality. This is separate from the destructive-output bug in finding 1.

**8. P2 — Invalid indel-error probabilities are accepted.**

Locations: `src/main.rs:178`, `src/main.rs:463`, `src/synth.rs:565`.

`--indel-error-rate` is described as a fraction of sequencing errors, but it has no finite/range validation. Values **NaN, 2.0 and -0.5** all exited successfully and wrote truth/output files. NaN and negative values disable the branch; values above 1 make every error selected by that branch an indel.

Fix: require a finite value in `[0,1]`, and clarify the help's conflicting “per base” versus “fraction of total error” wording. This is a configuration correctness issue, lower priority than the preceding model/output failures.

**Model assessment**

The strongest part is `VariantHaplotype`: reference segments, orientation and novel sequence give a consistent way to produce SNP/indel and SV reads. Keep it. The full tandem-duplication model is also preferable to independently adding depth and junction evidence. Empirical fragment lengths, separate quality distributions for R1/R2, randomized mate orientation and phased SNP inheritance are sensible components.

Several assumptions restrict what those components can establish:

| Area | Assessment and consequence |
| --- | --- |
| Background genome | `SampleCopies` is two maps of single-base substitutions. It cannot preserve background indels or SVs. The earlier CR3 measurement reported a hom-alt 2 bp deletion diluted to AF 0.360 under a het DUP. The present types and non-SNP filtering still have this limitation; that historical measurement was not rerun. |
| Local depth | One breakpoint depth controls uniform tiling across the event. A local donor depth `D(x)` becomes approximately `(1-v)D(x)+vC` for a copy-neutral edit, with `C` taken elsewhere. Uneven coverage is flattened and uncovered regions can gain synthetic coverage. `SIM_DEPTH_FOLD` measures and warns; it does not correct this. The previous decision to retain a warning-only model is documented. |
| Meaning of AF | AF acts as an edited-copy/replacement fraction for full DEL/DUP/INV, but as an added-junction fraction for additive modes. The same number is not one universal observed read fraction. Make edit fraction, expected copy number, and measured read evidence distinct fields. |
| Ploidy and genotype | The sample is always represented as two copies. Input GT is not used as the event model; truth GT is `1/1` at AF >=0.9 and `0/1` otherwise. That does not encode haploidy, pre-existing CN changes, or cellular mixtures. A value of 0.9 is not itself proof of a homogeneous homozygous genotype. |
| Phasing | Read-backed blocks and gVCF PS links make sense. Unconnected blocks select event copies independently. That is a defensible uncertainty assumption, but it cannot establish chromosome-wide phase or linked events. Base quality is also absent from the pileup/assignment interfaces, so low-quality bases can influence copy selection. |
| Quality/error process | The Markov quality model is useful. Sampling errors from `10^(-Q/10)` assumes calibration; it does not learn the empirical mismatch process, error-context biases or PCR errors. Matching quality marginals does not demonstrate matching variant-caller difficulty. |
| Library trimming | Inferring trimming from depletion of terminal A and then applying one TruSeq adapter prefix is a dataset-specific heuristic. It was measured on particular HG001/HG002 libraries; it should be configurable and verified on other preprocessing protocols. In the raw-FASTQ route, synthesized reads already follow the BAM's processing assumptions before entering preprocessing again. The downstream effect of that second pass was not measured here. |
| Complex/repetitive events | `origin` is appropriately experimental. MAPQ/XA weights, omitted alternatives and duplicate-family identity are approximations, in addition to finding 2. Fusion mode adds one junction and preserves originals; it does not simulate a balanced translocation. Novel-insertion interiors are deliberately not sampled when a fragment would lack a reference anchor, limiting assembly/unmapped-read realism. |

The low-count floor needs a narrower claim. In another probe, requesting SNP AF 0.001 emitted two synthetic pairs and recorded `SIM_VAF=0.014`; **neither pair carried the SNP**. Stochastic absence at low depth is legitimate, but two pairs somewhere in a roughly 4 kb footprint do not guarantee that an event was planted. `SIM_VAF` here is a model-derived expectation, not measured realization. Record informative fragment counts and distinguish expected from observed AF; do not condition on positive evidence silently, since that would bias sensitivity estimates.

Validation should have separate contracts for structural correctness, evidence detectability, statistical agreement and realism against independent data. The current `del_planted`/`ins_planted` checks are useful provenance checks, but a matching `SPIKE_` read and one junction probe cannot prove the requested AF or that a caller faces realistic evidence. `--strict` promotes advisory rows; it does not remove their measurement limitations. Existing transplantation and negative-control work is valuable and should become a repeatable acceptance suite, rather than relying only on a high unit-test count.

**Architecture and maintainability**

The project can improve incrementally; a wholesale rewrite is unnecessary. The main architectural problem is that one run does not have a single explicit model of the sample, molecules, edits and resulting evidence. Probabilities, coordinate conventions and metadata are instead reconciled across several passes. Finding 2 is a concrete consequence.

1. **Separate input parsing, planning, simulation and publication.** Turn raw CLI/VCF specifications into validated `EventPlan`s, and return an `EventResult` containing requested parameters, modeled expectations, generated/suppressed counts, evidence counts and warnings. Truth, README and logs should consume the same result. Currently parallel vectors (`adjusted_afs`, `resistant`, `depth_folds`, `event_stats`, outputs) must remain index-aligned.
2. **Represent fragments and library identity explicitly.** `ReadPair` loses alignment blocks, read groups and molecule/barcode information, while `OriginRecord` rebuilds a different fragment view. Preserve enough identity to reason about mates, duplicate families, per-library synthesis and read-group-aware output. Do not silently pool multiple samples and label replacements with the first SM.
3. **Create a shared alignment-reading layer.** Extraction, LOH, origin, census and validation repeat BAM/CRAM scans and eligibility decisions. `census` currently calls helpers owned by `validate`, reversing the useful dependency direction. Share coordinate/CIGAR mechanics while retaining named policies for training, editing and measurement; those policies should not become one indiscriminate filter.
4. **Use stronger domain types.** Replace string model selectors with enums, give AF validated construction, and distinguish reference positions, cuts and haplotype offsets. `TruthEvent` has a string tag and many independently optional fields; an enum can prevent impossible combinations. Make haplotype segment/layout mutation controlled: public segments plus separately stored concatenated sequence and total length create a synchronization obligation.
5. **Remove avoidable memory duplication.** `SharedReference::load` loads entire event chromosomes while `ReferenceReader` also caches clones of those same sequences. During loading, reference storage is approximately twice the requested sequence size, plus transient buffers. Events and kept output pairs are also retained together until final writing. Prefer an indexed window cache/shared immutable sequence storage, then stream or spool finalized event outputs. These are source-derived resource risks; peak RSS was not benchmarked here.
6. **Split the large orchestration/validation modules along responsibilities.** `main.rs` has about 2,600 lines before its tests and `validate.rs` about 3,800. Shell generation, CLI parsing, measurements, verdict policies and presentation should be independently testable. The embedded shell should become a tested output component or a paired-FASTQ Rust implementation.
7. **Make runs auditable and atomic.** Write a manifest with version/commit, seed, normalized events, input identities, model settings and completion state. Publish a completed run as a coherent set of outputs. Reusing an output directory currently permits old alignments and newly generated truth/FASTQ files to coexist without a run-identity check.
8. **Make the evidence easier to maintain.** No CI workflow or pinned Rust toolchain is present in the tracked tree; Cargo.lock is present. Add reproducible test dependencies and small end-to-end fixtures for the findings above. Keep research narratives in dated reports, and maintain a short current limitations/status document. The long README and historical review files intermingle old observations, partial fixes and current contracts, making it unnecessarily difficult to determine what remains true.

**Recommended order of work**

First fix destructive FASTQ publication, the origin conservation mismatch, existing-allele handling and gVCF fallback. Then correct the validator's base-depth semantics and per-mate library statistics. Keep nearby-event rejection and `origin`'s experimental status while those fixes are validated.

Next introduce explicit event results and one fragment-placement model used by both removal and synthesis. Those changes support the larger work: background-indel preservation, coverage-aware sampling and genotype/CN semantics. Require independent conservation and sequence checks, plus held-out real-data comparisons, before describing any of these models as broadly realistic.
