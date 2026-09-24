# spike — code review

Date: 2026-09-24 · Reviewed: working tree on `master` at `f5428ce` plus uncommitted edits to `extract.rs`, `loh.rs`, `simulate.rs`, `synth.rs`.

**Verdict:** the code is tidy and hard to crash, but in several common cases it produces quietly wrong output: the realized allele fraction is not what was asked for, or the variant lands in the wrong place. For a tool whose job is to produce truth sets, these matter most, because nothing downstream will flag them. Eight high-severity issues below; fix those before trusting any benchmark made with spike.

## How this was checked

- The code was read in four parts: CLI and input parsing; haplotype and simulation core; read synthesis and I/O; LOH, truth VCF and `validate`.
- Every finding was reproduced in a scratch copy: a unit test, a harness calling the real functions, or a real run on chr20 data (HG002 NovaSeq 35x, NA18488, GRCh38). The numbers below are measured, not estimated.
- Findings marked *plausible* are believed but not fully proven.
- Apart from adding this file, nothing in the repo was changed.

## Health check

| Check | Result |
| --- | --- |
| `cargo build --release` | OK, 1 warning (unused `primary_chrom`, `is_within_single_segment` in `haplotype.rs`) |
| `cargo test` | 128 passed, 0 failed |
| `cargo clippy --all-targets` | Style only: 6× `is_multiple_of`, 4× too many arguments, 2× use `?`, 1× no-effect op, 1× range loop, 1× manual `contains` |
| `scripts/validate_pipeline.sh` | Broken (see M17) |

The tests pass, but most would still pass with the high-severity bugs below. See [Test gaps](#test-gaps).

## Summary table

| ID | Severity | Problem | Where |
| --- | --- | --- | --- |
| H1 | High | **Fixed.** Nearby events undo each other's read suppression | `main.rs:431-437`, `main.rs:517-530` |
| H2 | High | **Fixed.** Exon numbers ignore strand | `exon.rs:99-103` |
| H3 | High | **Fixed.** BND orientation misread from VCF input | `vcf_input.rs:325-327` |
| H4 | High | **Fixed.** Inverted fusion simulated ~flank bp from truth position | `haplotype.rs:337-341`, `truth.rs:127-131` |
| H5 | High | **Fixed.** Fusion / junction-DUP get 2× junction reads | `simulate.rs:337-341` |
| H6 | High | LOH allele chosen at random per SNP | `loh.rs:112-121` |
| H7 | High | Insertions ≥ ~500 bp yield almost no insert-carrying reads | `synth.rs:745-746` |
| H8 | High | Synthetic reads lack the sample's own SNPs | `simulate.rs:77-80` |
| M1 | Medium | Allele fraction drifts at haplotype edges | `simulate.rs:123-130, 352, 437` |
| M2 | Medium | Short-insert libraries under-tiled | `simulate.rs:379` |
| M3 | Medium | `--flank` < 2000 leaves extra reads | `main.rs:379` |
| M4 | Medium | Chromosome-end segments overcount length | `haplotype.rs:81-104` |
| M5 | Medium | `merge.sh` loses/duplicates reads; adds sample `SIM` | `main.rs:1130`, `main.rs:667-671` |
| M6 | Medium | Truth VCF unsorted, no `##contig`, `REF=N` | `truth.rs:24-72, 94`; `main.rs:451` |
| M7 | Medium | **Fixed.** Same `--seed` gives different output | `extract.rs:376` |
| M8 | Medium | `--region` merged with distant events | `main.rs:172-180` |
| M9 | Medium | Fusion read pools double-counted | `main.rs:906-911`; `simulate.rs:471-503` |
| M10 | Medium | `validate` coverage check always passes on WGS | `validate.rs:753-756` |
| M11 | Medium | `validate` exits 0 when every check errors | `validate.rs:76-80, 136-158` |
| M12 | Medium | `validate` split-read check passes with no simulation | `validate.rs:503-517` |
| M13 | Medium | All synthetic pairs are F1R2 | `synth.rs:497-504, 736-742` |
| M14 | Medium | Missing base qualities → invalid FASTQ | `extract.rs:453-460` |
| M15 | Medium | CRAM extraction ~300× slower than BAM | `extract.rs:240-243, 314-317` |
| M16 | Medium | LOH pileup memory ~1 GB per Mb | `loh.rs:580-583` |
| M17 | Medium | `validate_pipeline.sh` no longer runs | `scripts/validate_pipeline.sh` |
| L1–L19 | Low | Parsing edge cases, robustness, minor I/O | see [Low](#low-severity) |

## High severity

### H1 · Nearby events undo each other's read suppression

`main.rs:431-437` collects every event's `kept_originals`, then `dedup_by_name` (`main.rs:517-530`) keeps the last copy of each read name. Each event returns *all* unsuppressed reads in its extraction window (event ± `--flank`, default 10 kb), not just its own footprint. A read suppressed by event 1 comes back through event 2's list.

- DEL chr1:[50000,52000) + SNP at 58000 (no overlap, so allowed): depth in the deleted region **1.000×** (expected 0.5×). The SNP comes out at AF 0.37 instead of ~0.5.
- Two SNPs 5 kb apart on chr20: the first drops from AF **0.556 to 0.355**. 78 of the 83 reads it suppressed reappear.
- With `--region`, every event on that chromosome interacts. Two SNPs 40 kb apart: AF 0.355 and 0.429.

**Fix:** gather suppressed read names from all events into one set and remove them from the final output, instead of keep-last dedup. Or run events in sequence over one shared pool.

### H2 · Exon numbers ignore strand

`exon.rs:99-103` numbers exons 1..n by genomic start. The BED has no strand column, and the exon number in the name (`TP53_exon1`) is ignored. Duplicate transcript lines also shift the numbering. Fusion breakpoints (`exon.rs:392-397`) always take gene A's genomic-left side.

- With a TP53 BED (minus strand), `del:TP53:exon1-exon1` resolves to `TP53_exon11` (7669608-7669690).
- README examples on minus-strand genes all hit the wrong exons: `dup:BRCA1:exon2-exon5`, `inv:TP53:exon3-exon6`, EML4-ALK (ALK).

**Fix:** read a strand column; number exons in transcript order, or parse the number from the name. Mirror the breakpoint side for minus-strand genes.

### H3 · BND orientation misread from VCF input

`vcf_input.rs:325-327` sets `inverted = open == ']' && close == ']'`. Per the VCF spec, `t[p[` and `]p]t` are same-strand joins; `t]p]` and `[p[t` are inverted. So `]p]t` is parsed as inverted and `[p[t` as forward. Both are wrong. The `Fusion` model also can't express the partner-first forms (`]p]t`, `[p[t`).

- One BCR-ABL1 mate pair, as written by `truth.rs`: with the chr22 record first → `chr22→chr9, inverted: false`. With the chr9 record first (normal sort order) → `chr9→chr22, inverted: true`, i.e. ABL1-left + reverse-complemented BCR.
- The unit tests assert the wrong behaviour, and contradict the function's own doc comment.

**Fix:** decode all four bracket forms. For partner-first forms, swap A and B, or build the event from the canonical mate only.

### H4 · Inverted fusion placed ~flank bp from the truth position

`haplotype.rs:337-341` fetches `ref_B[bp_b .. bp_b + flank]` and then reverse-complements it, so the B-side junction lands at `bp_b + flank - 1`.

- VCF input `N]chr20:34000000]`: B-side split reads at chr20:34001957-34001999. The truth VCF says 34000000.
- The mate record is written as `[chrA:p[N` (`truth.rs:127-131`). The mate of `t]p]` should also be `t]p]` form.
- `validate` only checks chrom A, so it still passes.

**Fix:** use `revcomp(ref_B[bp_b + 1 - flank .. bp_b + 1])` and write the mate as `N]chrA:pos]`.

### H5 · Fusion and junction-model DUP get 2× junction reads

`simulate.rs:339` uses `zone_per_bp = 2.0 * mean_frag`. Only fragments starting in (bp − f, bp) can cross a breakpoint, so the zone is f, not 2f. The supporting fraction becomes 2v / (1 + 2v).

| Requested VAF | Measured supporting fraction |
| --- | --- |
| 0.05 | 0.091 |
| 0.10 | 0.166 |
| 0.20 | 0.285 |
| 0.50 | 0.498 |

It is only right at 0.5, which is probably why earlier validation passed. The README's `af=0.05` fusion example comes out at ~9%.

**Fix:** n = cov · v / (1 − v) per breakpoint (capped as v → 1), or suppress originals at the breakpoint at rate v.

### H6 · LOH allele chosen at random per SNP

`loh.rs:112-121` picks the target allele by coin flip at each het SNP. A fragment spanning two SNPs whose picks disagree scores a tie. Ties go into `classified_set` but not `loh_set`, so at VAF 0.5 they are never suppressed (`simulate.rs:143-158`). The README says ties are excluded from both sets. Phase from the gVCF (`0|1`) is discarded.

- Het DEL chr20:38412500-38422500, where HG002 truth has every SNP phased `1|0`: the surviving allele is REF at 38412870 but ALT at 38419488, 38420312, 38421231. Some SNPs stay het (38414348 at 0.28 ref, 38414563 at 0.73).
- DUP: REF duplicated at 38424935, ALT at 38426138 and 38426372 — a haplotype that doesn't exist in the sample.
- Unit test with two in-phase SNPs spanned by all 40 fragments: empty target set in 93 of 200 seeds.
- Real logs: 8–21% of classified fragments are ties.

**Fix:** phase SNPs (gVCF GT phase / PS, else read-backed phasing) and pick one haplotype per phase block. Treat ties as unclassified (random suppression at the VAF rate).

### H7 · Long insertions yield almost no insert-carrying reads

`synth.rs:745-746` maps R1 start and R2 end back to the reference with `hap_to_ref(...)?`. If either falls in novel sequence the pair is dropped (accepted by `simulate.rs:401-403`). Only fragments spanning the whole insertion survive.

| Insertion length | Pairs carrying ≥ 10 inserted bases (30 runs) |
| --- | --- |
| 20 bp | 453 |
| 100 bp | 434 |
| 250 bp | 243 |
| 500 bp | 8 |
| 1000 bp | 0 |

The README's own 500 bp insertion example produces essentially no evidence. Long sequence-resolved VCF insertions are hit the same way.

**Fix:** keep these pairs. For `ReadPair` metadata, fall back to the nearest reference-mapped base or the breakpoint position.

### H8 · Synthetic reads lack the sample's own SNPs

The haplotype is built from the reference. Only DUP-interior het SNPs are substituted (`simulate.rs:77-80`; the variant map is empty for DEL, `simulate.rs:279`). Reads across the 2 kb flanks are replaced at the VAF rate with reference-only reads, so real SNPs drift toward REF.

| Site | Original ref/alt | After simulation |
| --- | --- | --- |
| DEL left flank 38411227 | 24/17 | 39/8 |
| DEL left flank 38411969 | 18/15 | 28/6 |
| DUP flank 38427232 | 18/30 | 25/11 |

Hom-alt SNPs will look het. INV, INS and SNP events share the same code path.

**Fix:** apply the sample's alleles (from pileup, on the correct haplotype) across the whole haplotype footprint. Or replace only reads that actually need to change.

## Medium severity

### M1 · Allele fraction drifts at haplotype edges
Tiling starts uniformly on `[0, H − f]` (`simulate.rs:437`), suppression uses overlap-fraction scaling (`simulate.rs:123-130`), and `n` is computed from the full reference length (`simulate.rs:352`). These don't match at the edges.
- Het DEL: 0.82× just inside the footprint edge, 0.90× just outside, flank interior 1.04–1.06×.
- SNP at requested AF 0.2: realized 0.211 / 0.217 / 0.222 for 300 / 400 / 550 bp fragments.
- **Fix (tested in scratch):** suppress only fully-contained pairs and use n = cov · v · (L_ref − mean_frag) / mean_frag. That gave edges at 0.96–1.00× and a synthetic/suppressed ratio of 0.984–0.996.

### M2 · Short-insert libraries under-tiled
`simulate.rs:379`: `mean_frag = pool.frag_dist.mean.max(300.0)`.
- With 220 bp inserts (exome ~200–250, cfDNA ~167): requested AF 0.2 → realized 0.160.
- **Fix:** use the real mean; guard only against 0 / NaN.

### M3 · `--flank` below 2000 leaves extra reads
`main.rs:379` hard-codes `hap_flank = 2000`, but suppression only covers event ± `--flank`.
- `--flank 500`, VAF 0.2: depth 1.18× and 1.16× in the unsuppressed zones.
- **Fix:** extract at least event ± (hap_flank + max fragment), or reject `--flank < 2000`.

### M4 · Chromosome-end segments overcount length
`haplotype.rs:81-104` (same pattern at 168-174, 231-236, 283-286): `ref_end` is not clamped, but the fetched sequence is (`reference.rs:171-172`).
- DEL ending 100 bp from the chromosome end: `ref_mapped_len` 4000 vs `total_len` 2100 → 250 fragments instead of 131 (~1.9× depth).
- **Fix:** after fetch, set `ref_end = ref_start + seq.len()`.

### M5 · `merge.sh` loses and duplicates reads; adds a second sample
`merge.sh` removes originals by BED region (`main.rs:1130`), but extraction (`extract.rs:92-100, 176-196`) drops some pairs and pulls in out-of-region mates.
- 10 kb DEL + flank: 1,309 of 9,555 primary records (13.7%) are lost and never replaced: 1,181 duplicates, 64 non-proper pairs, 64 pairs with a low-MAPQ or orphaned mate. 40 out-of-BED mates appear twice.
- `align.sh` tags reads `SM:SIM` (`main.rs:667-671`). The merged BAM had 12 read groups with `SM:NA18488` and 1 with `SM:SIM`. Multi-sample callers will likely split the region into a separate sample (*plausible*).
- **Fix:** pass filtered pairs through unchanged and remove originals by read name (`samtools view -N`). Reuse the original SM in `-R`.

### M6 · Truth VCF unsorted, no `##contig`, `REF=N`
`truth.rs:24-72, 94`; events written in input order (`main.rs:451`).
- `bcftools index` fails ("Unsorted positions"). BCF output fails. `bcftools norm -c e` fails on REF mismatch. truvari only works after sorting and adding contigs.
- **Fix:** sort records, emit `##contig` from the `.fai`, fetch the REF base from the reference.

### M7 · Same `--seed` gives different output
Pairs come out of a `HashMap` in random order (`extract.rs:124`, `285`). `build_read_pool` sorts by `ref_start` only (`extract.rs:376`), so ties stay random, and the suppression draws land on different reads.
- One SNP, two runs: different FASTQ md5s and different read sets.
- **Fix:** sort by `(ref_start, name)`. With that, three runs gave identical md5s.

### M8 · `--region` merged with distant events
`extraction_bounds` (`main.rs:172-180`) takes the min/max of region and event window.
- `--region chr20:30490000-30510000` + fusion partner at chr20:35000001: window 30.49–35.01 Mb, 514,980 pairs written instead of ~4.5k.
- **Fix:** merge only when they overlap.

### M9 · Fusion read pools double-counted
`main.rs:906-911` concatenates A and B pools without removing shared reads. `estimate_coverage_at` (`simulate.rs:471-503`) ignores chromosome.
- Intra-chromosomal breakpoints 3 kb apart: 106 chimeric pairs vs 53 when far apart.
- Different chromosomes with close coordinates: coverage 85 vs 45.
- **Fix:** dedup `pairs_a ∪ pairs_b` by name; filter coverage by chromosome.

### M10 · `validate` coverage check always passes on WGS
`count_depth_in_region` returns reads per bp (`validate.rs:753-756`), ~0.23 at 35x. The guard `flank_depth < 1.0` then reports "no flanking coverage" with `pass: true` (`validate.rs:454-462`).
- Passed on the untouched HG002 BAM and on simulated DEL and DUP BAMs.
- **Fix:** depth = aligned bases / length; "not evaluable" should fail or skip, not pass.

### M11 · `validate` exits 0 when every check errors
`validate.rs:76-80, 136-158` log per-event errors as warnings.
- Truth VCF with `20` vs BAM with `chr20`: "Result: 3/3 PASS", exit 0.
- **Fix:** count errored checks as failures.

### M12 · `validate` split-read check passes with no simulation
`validate.rs:503-517` counts any read with an SA tag in [start − 500, end + 500] and never checks where the SA partner lands.
- The untouched HG002 BAM passed for both DEL and BND.
- **Fix:** require the partner near the expected other breakpoint, and a minimum count scaled by coverage × VAF.

### M13 · All synthetic pairs are F1R2
`synth.rs:497-504, 736-742` always make R1 forward, R2 reverse relative to the haplotype.
- `dup_test/sim.bam`: synthetic R1 forward 659/659; real reads 7305 fwd / 7357 rev.
- Mutect2's read-orientation filter will likely flag simulated SNVs.
- **Fix:** for a random half of fragments, take R1 from the right end on the reverse strand.

### M14 · Missing base qualities → invalid FASTQ
`extract.rs:453-460` does `s.wrapping_add(33)`; missing qualities (0xFF) become byte 32 (space).
- Kept originals are written with spaces as quality characters.
- The quality model learns 32. The new clamp in `sample_quality` turns it into Q0 and the Markov chain sticks there. With 5% of donor pairs missing quality: 4.99% of synthetic reads were mostly Q0 with random bases.
- **Fix:** skip or flag reads with all-0xFF quality at extraction; refuse to write qualities outside 33–126.

### M15 · CRAM extraction ~300× slower than BAM
noodles-cram 0.74 `Query::read_next_container` checks only the reference id, never position, so every container on the chromosome is decoded (`extract.rs:240-243, 314-317`).
- 30 kb region: 0.12 s from BAM, 39.97 s from a CRAM with only 10 Mb of chr20. `samtools view`: 0.015 s.
- **Fix:** filter index entries by position and seek manually, or upgrade noodles after checking the fix.

### M16 · LOH pileup memory ~1 GB per Mb
`loh.rs:580-583` stores one `(u64, u8)` per aligned base per read.
- 3 Mb DEL: 4.58 GB peak RSS with LOH vs 1.37 GB without. Arm-level events will run out of memory.
- **Fix:** call het sites from counts first, then record read alleles only at those sites.

### M17 · `scripts/validate_pipeline.sh` no longer runs
- Step 3 aborts: the filtered truth has overlapping DELs, and the script doesn't pass `--allow-overlap` or filter them.
- `$(grep -vc '^#' f || echo 0)` gives `"0\n0"` on empty input; the `[[ -lt 5 ]]` guard then errors silently and the run continues with 0 events (lines 252, 264, 315, 570).
- Paths at lines 33-35 point to old `data/giab_hg38/` locations. Tool paths are hard-coded to `/home/parlar_ai/dev/sv_caller/.pixi`.

## Low severity

| ID | Problem | Where | Fix |
| --- | --- | --- | --- |
| L1 | Indel error model pads deletions with `N` (Q2): 0.69% of R1 end in N, 1.01% of R2 start with N at rate 0.05 | `synth.rs:697-700, 729-730` | Pass `read_length + 10` bases; generate R2 in sequencing order |
| L2 | CRAM containers with several contigs leak other contigs' reads (synthetic 2-contig CRAM: 5 chrB pairs labelled chrA) | `extract.rs:245-279, 319-362` | Skip records whose ref id or mate ref id differs |
| L3 | bgzipped FASTA read as raw bytes; fails later with misleading "beyond chromosome length" | `reference.rs:33-35` | Use `fasta::io::indexed_reader::Builder` |
| L4 | Final FASTQ flush error ignored (write to `/dev/full` returned `Ok`) | `fastq.rs:47-48` | `r1_gz.finish()?.flush()?` |
| L5 | Panic when mean read length > 1500 (`clamp` with min > max) | `stats.rs:103`; callers `simulate.rs:413`, `synth.rs:534` | Guard min ≤ max |
| L6 | **Fixed.** BND POS is one base past the kept base, vs spec; parser mirrors it so round trips agree | `truth.rs:121`, `vcf_input.rs:103` | Write POS = last kept base |
| L7 | DEL/DUP/INV with no END and no SVLEN silently becomes a 1 bp event | `vcf_input.rs:130-134, 150-154, 169-173` | Error, or derive from REF length |
| L8 | `af=nan` / `--allele-fraction NaN` accepted; suppresses all reads, truth says `SIM_VAF=NaN` | `exon.rs:202`, `main.rs:241` | `if !(v > 0.0 && v <= 1.0)` |
| L9 | REF == ALT checked before uppercasing, so `A:a` passes; VCF path has no check at all | `exon.rs:508`, `vcf_input.rs:226-239` | Uppercase first; add check to VCF path |
| L10 | Sequence INS with multi-base REF duplicates `REF[1..]` (`REF=AT ALT=ATGGG` → 4 bp inserted, not 3) | `vcf_input.rs:197-200` | Strip common prefix |
| L11 | VCF records silently dropped (multi-allelic, unknown SVTYPE, `DUP:TANDEM`, short lines); INFO `AF` (often population frequency) used as VAF | `vcf_input.rs:67-101` | Warn with counts; don't use plain `AF` by default |
| L12 | `align.sh` / `merge.sh` break on relative paths and on `}`, `"`, `$`, backticks in paths | `main.rs:681-683, 1118-1121` | Canonicalize and single-quote paths |
| L13 | Contig names with `:` (525 HLA contigs in this BAM) can't be used | `exon.rs:139`, `main.rs:191`, `vcf_input.rs:314` | `rsplit_once`; match known contigs |
| L14 | Standard BED6: column 5 is score, so every gene is named `"0"` | `exon.rs:79-83` | Use column 4 or detect format |
| L15 | `validate` global checks read only the first 100k records (whole-genome HG002 mean MAPQ 10.0 → FAIL); dup rate always 0 on unmarked BAM | `validate.rs:1031-1140` | Sample across the file / event regions |
| L16 | `validate`: `truncate()` panics on non-ASCII gene names; `escape_json` misses `\t`, `\r` | `validate.rs:1266` | Truncate on char boundary; escape all control chars |
| L17 | R2 quality Markov chain runs backwards (no measurable effect on NovaSeq: mean Q 35.628 vs 35.627) | `synth.rs:413-421, 651-659` | Fixed along with L1 |
| L18 | Reference `N` bases get normal qualities (often Q37); real Illumina N is Q2 | `synth.rs:423-428` | Force Q2 on N |
| L19 | Single-end BAM/CRAM is scanned end to end looking for proper pairs | `bam_stats.rs:81-83, 141-143` | Cap records scanned |

## Uncommitted changes

| File | Change | Assessment |
| --- | --- | --- |
| `extract.rs` | `safe_noodles_position` no longer `expect`s | Correct, no behaviour change. All call sites pass `start + 1`, `end` for 0-based half-open input. |
| `synth.rs` | Clamp sampled quality to Phred+33 range | Correct in itself, but it hides M14 (missing qualities) instead of fixing it. |
| `loh.rs` | Warn when gVCF chromosome names don't match | Wrong in both paths. For `.vcf.gz` it never fires (`bcftools view -r chr17:…` on a `17` file returns nothing, exit 0). For plain VCF it fires whenever the region has no het SNPs and other chromosomes exist. "Falling back to pileup" is only true when the result is empty; if bcftools errors (e.g. unindexed `.gz`), LOH is skipped. |
| `simulate.rs` | `unreachable!` → `bail!`; two new `simulate_event` tests | `bail!` is fine. The tests only assert `> 0` / non-empty and would pass with H5, 100% suppression, or 2 tiled reads. |

## Test gaps

- New and existing `simulate` tests assert "> 0" instead of expected counts. `test_tandem_dup_tiling_count_uses_full_length` asserts `n > 100` where its comment expects 188. `_junction_pairs` is computed and never asserted.
- `test_suppression_*` re-implements the suppression loop instead of calling `simulate_event`. The LOH branch has no test.
- `haplotype.rs` tests never call a real `VariantHaplotype::from_*` constructor (so H4 and M4 went unnoticed), though `SharedReference::from_sequences` exists for this.
- BND parser tests encode the wrong orientation (H3).
- No tests for: `loh.rs`, `truth.rs`, main-level orchestration with several events or `--region`, minus-strand BEDs, script generation, R1/R2 geometry in `generate_read_pair`, `indel_error_rate > 0`, extraction from a real BAM/CRAM, FASTQ round trip.
- A good first regression test for most of H1–H8: simulate on a small synthetic reference and assert the realized VAF / depth / breakpoint position within a tolerance.

## Design notes

- Coordinate conventions are mixed: `del/dup/inv` take VCF POS/END meaning, `snp` is 1-based, `--region` is 1-based inclusive. `del:chr1:0-100` is rejected with "coordinates must be >= 1" although the documented convention is 0-based.
- `del:LDLR:4-8` parses as coordinates on a chromosome named "LDLR". An unknown fusion suffix (`:inverted`, `:rev`) silently gives a forward fusion.
- VCF input ignores FILTER and GT; 0/0 and non-PASS records are simulated.
- Coverage is estimated once, ±1 kb around the first breakpoint, and applied to the whole haplotype. `n.max(2)` emits 2 pairs even at zero coverage. An empty read pool silently gives a constant-Q20 profile.
- The CIGAR walk is duplicated four times in `loh.rs` plus once in `validate.rs`. `generate_read` and `generate_read_from_seq` are ~70 near-duplicate lines that have already diverged (root of L1).
- Quality profile stores every observed quality four times and sorts them unnecessarily; a histogram per bin would do. `SharedReference::load` doubles peak memory via the FIFO cache.
- gVCF sample column is always 9 despite a comment saying it reads the header. Pileup het caller has no base-quality filter and counts overlapping mates twice. A plain-text gVCF is re-read in full per event.

## What is solid

- One haplotype-segment model for all variant types is a sound design. DEL, full and junction DUP, INV (reverse complement and coordinate mapping) and SNP/MNV/indel anchor handling are correct.
- DEL/DUP/INV/INS truth coordinates are right: POS is the preceding base, END is 1-based inclusive, DEL SVLEN is negative.
- Reverse-strand reads are un-reverse-complemented with qualities reversed in both BAM and CRAM paths. Synthetic per-cycle quality matches real reads within ~0.05 Q (R1 35.966 vs 35.962). Errors are applied at P = 10^(−Q/10).
- Expected depth and allele balance for DEL and DUP at VAF 0.5 are right.
- Input intervals are validated before use and flank arithmetic saturates: no reachable panics from normal input were found. REF alleles are checked against the reference.
- The pileup het caller found exactly the 13 truth het SNPs in a 10 kb test region.

## Suggested fix order

1. **H1** (global suppressed-name set) and **M7** (sort tie-break). Small changes, affect every multi-event run.
2. **H2**, **H3**, **H4**: wrong variant placement. Each needs a test using a real constructor.
3. **H5**, **M1**, **M2**, **M3**: VAF math. One shared test that measures realized VAF at 0.05 / 0.2 / 0.5.
4. **H6**, **H8**: haplotype-consistent LOH and carrying sample alleles. Biggest design change; the fixes overlap.
5. **H7**: long insertions.
6. **M6**, **M10**–**M12**: make the truth VCF indexable and make `validate` able to fail.
7. The rest as time allows.
