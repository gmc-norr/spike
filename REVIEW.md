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
| `scripts/validate_pipeline.sh` | **Fixed** (M17): runs end to end, and fails (exit 1) when the spike-in contributed nothing the background does not already carry |

The tests pass, but most would still pass with the high-severity bugs below. See [Test gaps](#test-gaps).

## Summary table

| ID | Severity | Problem | Where |
| --- | --- | --- | --- |
| H1 | High | **Fixed.** Nearby events undo each other's read suppression | `main.rs:431-437`, `main.rs:517-530` |
| H2 | High | **Fixed.** Exon numbers ignore strand | `exon.rs:99-103` |
| H3 | High | **Fixed.** BND orientation misread from VCF input | `vcf_input.rs:325-327` |
| H4 | High | **Fixed.** Inverted fusion simulated ~flank bp from truth position | `haplotype.rs:337-341`, `truth.rs:127-131` |
| H5 | High | **Fixed.** Fusion / junction-DUP get 2× junction reads | `simulate.rs:337-341` |
| H6 | High | **Fixed.** LOH allele chosen at random per SNP | `loh.rs:112-121` |
| H7 | High | **Fixed.** Insertions ≥ ~500 bp yield almost no insert-carrying reads | `synth.rs:745-746` |
| H8 | High | **Fixed.** Synthetic reads lack the sample's own SNPs | `simulate.rs:77-80` |
| M1 | Medium | **Fixed.** Allele fraction drifts at haplotype edges | `simulate.rs:123-130, 352, 437` |
| M2 | Medium | **Fixed.** Short-insert libraries under-tiled | `simulate.rs:379` |
| M3 | Medium | **Fixed** (rejected). `--flank` < 2000 leaves extra reads | `main.rs:379` |
| M4 | Medium | **Fixed.** Chromosome-end segments overcount length | `haplotype.rs:81-104` |
| M5 | Medium | **Fixed.** `merge.sh` loses/duplicates reads; adds sample `SIM` | `main.rs:1130`, `main.rs:667-671` |
| M6 | Medium | **Fixed.** Truth VCF unsorted, no `##contig`, `REF=N` | `truth.rs:24-72, 94`; `main.rs:451` |
| M7 | Medium | **Fixed.** Same `--seed` gives different output | `extract.rs:376` |
| M8 | Medium | **Fixed.** `--region` merged with distant events | `main.rs:172-180` |
| M9 | Medium | **Fixed.** Fusion read pools double-counted (dedup half landed with M8 — see note below) | `main.rs:906-911` (superseded, see note); `simulate.rs:196-203, 566-598` |
| M10 | Medium | **Fixed.** `validate` coverage check always passes on WGS | `validate.rs:753-756` |
| M11 | Medium | **Fixed.** `validate` exits 0 when every check errors | `validate.rs:76-80, 136-158` |
| M12 | Medium | **Fixed.** `validate` split-read check passes with no simulation | `validate.rs:503-517` |
| M13 | Medium | **Fixed.** All synthetic pairs are F1R2 | `synth.rs:497-504, 736-742` |
| M14 | Medium | **Fixed.** Missing base qualities → invalid FASTQ | `extract.rs:453-460` |
| M15 | Medium | **Fixed.** CRAM extraction ~300× slower than BAM | `extract.rs:240-243, 314-317` |
| M16 | Medium | **Fixed.** LOH pileup memory ~1 GB per Mb | `loh.rs:580-583` |
| M17 | Medium | **Fixed.** `validate_pipeline.sh` no longer runs | `scripts/validate_pipeline.sh` |
| L1–L19 | Low | Parsing edge cases, robustness, minor I/O | see [Low](#low-severity) |
| N1 | Medium | **Not fixed** (found during the fix run). `spike validate` scores a cross-sample spike-in against a confounded background, and `split_reads` looks for a signal spike does not emit | `validate.rs:440-520`; `scripts/validate_pipeline.sh` |

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

**Fixed:** every event reads the sample's het and hom-alt SNPs over its footprint, removes reads by copy, and gives synthetic reads their copy's alleles. Same runs (seed 1): 38411227 → 23/17, 38411969 → 12/15, 38427232 → 18/27. Hom-alt SNPs in INV, INS, DEL (VAF 0.8) and fusion footprints get 0 REF reads (before: e.g. 12/32, 18/18, 43/7).

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
- DEL ending 100 bp from the chromosome end: `ref_mapped_len` 4000 vs `total_len` 2100 → 250 fragments instead of 131 (~1.9× depth). **Stale:** the length mismatch is real and reproduces, but the fragment count predates `c07f9d2`; `compute_tiling_count` has scaled by `total_len` ever since, so the DEL tiled-pair count is unchanged. See the measured result below.
- **Fix:** after fetch, set `ref_end = ref_start + seq.len()`.

**Fixed:** every constructor now takes a segment's `ref_end` from the sequence the fetch returned (`haplotype.rs`: `from_deletion`, `from_duplication`, `from_tandem_duplication`, `from_inversion`, `from_insertion`, `from_small_variant`, and the `from_fusion` piece). HG002 on `chr14_KI270723v1_random` (38,115 bp, ~118x to its last base), seed 1: a 10 kb DEL ending 100 bp from the contig end (`del:...:28015-38015`) had `ref_mapped_len` 4000 against `total_len` 2100 → 2100/2100. Its **read count does not move** (229 tiled pairs, byte-identical FASTQ, before and after): `compute_tiling_count` has scaled by `total_len`, not `ref_mapped_len`, since `c07f9d2`, so the 250-vs-131 above belongs to the older `n = cov·VAF·ref_mapped_len / frag` — only its ratio survives (4000/2100 = 1.905 = 250/131). What the unclamped end still broke is *reversed* segments, whose `hap_to_ref` counts down from `ref_end`: a right-right fusion cut 100 bp from the contig end (`fusion:GENEA:exon1:GENEB:exon2`, bp_a 38,015) placed its junction at 39,915, 1,800 bp past the contig, where coverage reads 0.0 and tiling falls to its 2-pair floor — **2 → 88 chimeric pairs** (cov 0.0 → 87.6). A contig-end INV and `dup:chr20:38423496-38427196` are byte-identical before and after.

### M5 · `merge.sh` loses and duplicates reads; adds a second sample
`merge.sh` removes originals by BED region (`main.rs:1130`), but extraction (`extract.rs:92-100, 176-196`) drops some pairs and pulls in out-of-region mates.
- 10 kb DEL + flank: 1,309 of 9,555 primary records (13.7%) are lost and never replaced: 1,181 duplicates, 64 non-proper pairs, 64 pairs with a low-MAPQ or orphaned mate. 40 out-of-BED mates appear twice.
- `align.sh` tags reads `SM:SIM` (`main.rs:667-671`). The merged BAM had 12 read groups with `SM:NA18488` and 1 with `SM:SIM`. Multi-sample callers will likely split the region into a separate sample (*plausible*).
- **Fix:** pass filtered pairs through unchanged and remove originals by read name (`samtools view -N`). Reuse the original SM in `-R`.

**Fixed:** spike writes `replaced_reads.txt` (the names it extracted) and `merge.sh` removes exactly those with `samtools view -N`; `align.sh` reuses the BAM's first `@RG SM`. HG002 chr20:37.5–41.5 Mb slice, `del:chr20:38412500-38422500`, seed 1, flank 10 kb: in-BED primary records lost and never replaced 1,436/10,501 (13.7%) → 0/10,501 (0.0%); originals present twice 53 → 0; merged-BAM `@RG` samples `SM:HG002` + `SM:SIM` → `SM:HG002` only. The 989 suppressed pairs (the deleted copy) stay gone in both.

**Fix pass 2:** README overstated how much of the recovered fraction reaches `spike validate`. `validate.rs` has no proper-pair or mate-unmapped filter (confirmed at `validate.rs:723-733, 755-765, 837-847, 869-879, 964-975`: only `is_unmapped/is_secondary/is_supplementary/is_duplicate/is_qc_fail` and `mq < min_mapq` are skipped), so the 87 non-proper-pair and 39 orphaned-mate records reach measurement, not just the 87 — corrected README to 126/1,436 (1.2% of the 10,501 in-BED records), down from a wrong 87 (0.8%). Also: reworded the shortfall-guard's comment and abort message in the generated `merge.sh` (`main.rs:1176-1188`) to match the hard-floor code (`ACTUAL < NAMES`, not "far fewer") and to name ACTUAL/NAMES for what they count; `outside.bam` is now removed on the abort path too; the `fake_samtools.sh` test stub's `view -c` branch now requires `removed.bam` to actually exist (`[ -f "$3" ] || exit 1`) so a future regression that stops materialising it can't hide behind the stub's path-only dispatch.

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

**Fixed:** `extraction_bounds` returns a list of windows and merges the region with event ± flank only when the two overlap or touch; `extract_pool_for_event` extracts every window and `extract::dedup_pairs_by_name` keeps each fragment once. HG002, `--region chr20:30490000-30510000` with a fusion partner at chr20:35000001 (exon BED, seed 1): one 30,489,999–35,010,000 window, 710,997 donor pairs, **711,130 pairs written → 6,237 donor pairs from two windows (30,489,999–30,510,000 and 34,990,000–35,010,000), 6,304 pairs written** (113×; the review's 514,980/~4.5k used its own breakpoints). A region that genuinely overlaps its event is byte-identical before and after (`del:chr20:38412500-38422500` with `--region chr20:38400000-38440000`, seed 1: same 6,117-pair pool, same R1 md5). The dedup also covers the first half of **M9**: an intra-chromosomal fusion 5 kb across (no `--region`) had its overlapping windows double-counted into a 6,051-pair pool, inflating local coverage and tiling 133 chimeric pairs — now 3,835 pairs and 67 chimeric. M9's chromosome-aware `estimate_coverage_at` is untouched.

### M9 · Fusion read pools double-counted
`main.rs:906-911` concatenates A and B pools without removing shared reads. `estimate_coverage_at` (`simulate.rs:471-503`) ignores chromosome.
- Intra-chromosomal breakpoints 3 kb apart: 106 chimeric pairs vs 53 when far apart.
- Different chromosomes with close coordinates: coverage 85 vs 45.
- **Fix:** dedup `pairs_a ∪ pairs_b` by name; filter coverage by chromosome.

**Stale note (added while fixing M8, `92ad2bd`; this entry is still not marked Fixed):**
this entry's first bullet and first `Fix:` clause — dedup `pairs_a ∪ pairs_b` by
name — landed in `92ad2bd`, but only as a required consequence of the M8 fix
(merge-on-overlap makes a same-chromosome fusion query `--region` once per
side, so the dedup was not optional there), not as separate M9 work. The
`main.rs:906-911` concatenation this entry cites no longer exists: M8 replaced
it with accumulate-then-`extract::dedup_pairs_by_name`, called at `main.rs:968`.
The "Intra-chromosomal breakpoints 3 kb apart: 106 vs 53" figures above no
longer reproduce post-`92ad2bd` — measured on M9's own scenario (an
intra-chromosomal fusion with close breakpoints, no `--region`), the dedup
alone took chimeric pairs 133 → 67, estimated coverage 133.0 → 66.5, and the
donor pool 6,051 → 3,835 (task-5 report, "Overlap with M9"). **Still open:**
the second bullet and second `Fix:` clause only. `estimate_coverage_at` is
chromosome-blind, at `simulate.rs:564-596` as of `92ad2bd` (not `471-503` —
the file has grown since this entry was written; verify the current range
before citing it again). Two breakpoints at similar coordinates on *different*
chromosomes still pool their coverage; untouched, and next in the queue.
(That second half has since landed — see **Fixed** below.)

**Fixed:** `estimate_coverage_at` takes the breakpoint's chromosome and counts
only pool pairs on it; the call site keeps the chromosome `hap_to_ref` already
returns for the first breakpoint. HG002, cross-chromosome fusion
chr20:40001200 >> chr21:40001000 (exon BED, `af=0.5`, seed 1) against the same
chr20 breakpoint with the partner moved to chr21:10001000: estimated coverage
**120.8 vs 61.5 → 61.5 vs 61.5**, chimeric pairs **121 vs 61 → 61 vs 61** (the
entry's "85 vs 45" is the same 2× shape at a different locus; measured
independently from the BAM, chr20 fragment depth over that window is 70.3 and
chr21's is 67.9, so the near side was carrying ~1.97× its own depth). The
distant-partner run is byte-identical before and after, as are a
same-chromosome fusion (`chr20:38421200 >> chr20:38424000`, coverage 60.1) and
`del:chr20:38412500-38422500` (coverage 68.0) — H5's `coverage * v / (1 - v)`
zone arithmetic is untouched, and a single-region event's pool is
single-chromosome so the filter is a no-op there. Two runs at `--seed 1` still
give byte-identical FASTQ (M7).

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
- **Fixed:** `generate_read_pair` and `generate_haplotype_read_pair` now draw one `rng.gen::<bool>()` per fragment, before either read, and it decides which mate is R1. The fragment's left end is still always read forward and its right end reverse, so a flipped pair covers exactly the same interval; only `read_num` (the quality model) and the final R1/R2 assignment swap, so the mate that is R1 keeps the R1 model. HG002 (`dup:chr20:38423496-38427196`, seed 1, bwa-mem2, primary alignments): synthetic R1 forward 808/808 → 402 fwd / 406 rev of 808 (a fair coin gives sd 14.2, so this is 0.1 sd off centre); the real reads in the same BAM are 1544 fwd / 1491 rev. Synthetic proper-pair rate 97.5% → 97.4% and mean |TLEN| 490.2 → 489.4 bp, i.e. the flipped pairs are still FR over the same fragments. Two runs at `--seed 1` still give byte-identical `R1.fq.gz`/`R2.fq.gz`/`truth.vcf`. (The review's 659/659 came from a BAM not available here; the counts differ, the 100%-forward symptom reproduced exactly.)

### M14 · Missing base qualities → invalid FASTQ
`extract.rs:453-460` does `s.wrapping_add(33)`; missing qualities (0xFF) become byte 32 (space).
- Kept originals are written with spaces as quality characters.
- The quality model learns 32. The new clamp in `sample_quality` turns it into Q0 and the Markov chain sticks there. With 5% of donor pairs missing quality: 4.99% of synthetic reads were mostly Q0 with random bases.
- **Fix:** skip or flag reads with all-0xFF quality at extraction; refuse to write qualities outside 33–126.
- **Fixed:** `parse_partial_from_bam_record`/`parse_partial_from_record_buf` now skip (not encode) a record whose raw quality is unusable, counted and logged (`log::warn!`) per extraction call; `write_paired_fastq` refuses (returns `Err`) a quality whose length does not match SEQ or that carries a byte outside 33-126, validating before it creates either file. The `synth.rs` clamp is kept as defense in depth (not proven unreachable — see "Judgement calls" in the task report). HG002 chr20 slice with 5% of donor pairs' quality stripped to `*` (seed 1, `del:chr20:38412500-38422500`): kept-original quality lines containing byte 32 (space) 188/3862 → 0/3667; synthetic (chimeric) reads "mostly Q0" 4.11% (R1) / 4.79% (R2) of 292 → 0.00% of 273; R1 mean quality 34.2 → 36.0 (no longer dragged down by poisoned Q0 runs). Absolute counts differ from the review's numbers because the donor BAM here is a synthetic 5%-stripped HG002 slice, not the original review's BAM.
- **Fix pass 1** (three defects in the fix itself, same slice re-encoded as CRAM with `samtools view -C`): (a) the all-0xFF test never fired on CRAM, because noodles-cram clears the buffer — CRAM donors still yielded 4559 pairs with 188 R1 + 188 R2 records carrying a 151 bp SEQ line and a zero-length QUAL line, and `samtools import` aborted with `truncated file`; now 4331 pairs, 0 zero-length QUAL lines, `samtools import` exits 0. (b) Skipped pairs were in neither `kept_originals` nor `suppressed_names`, so `merge.sh` merged them back untouched: 147 of 2868 eligible records inside `chr20:38412500-38422500` (5.13%) were `*`-quality and un-suppressible → 0.00%, by listing them in `replaced_reads.txt` (4331 → 4562 names; 231 pairs removed without replacement). (c) `missing_qual_count` double-counted any record whose mate pass 1 kept: on a slice with only read1 stripped, 464 reported against a true 230 distinct records → 235 (the tally is now keyed by (name, segment); the 5 extra come from pass 2's wider mate-recovery window, which the message no longer claims to cover). Extraction now also skips any record carrying a raw quality byte above Q93, so "partial 0xFF cannot occur" is enforced rather than assumed.

### M15 · CRAM extraction ~300× slower than BAM
noodles-cram 0.74 `Query::read_next_container` checks only the reference id, never position, so every container on the chromosome is decoded (`extract.rs:240-243, 314-317`).
- 30 kb region: 0.12 s from BAM, 39.97 s from a CRAM with only 10 Mb of chr20. `samtools view`: 0.015 s.
- **Fix:** filter index entries by position and seek manually, or upgrade noodles after checking the fix.
- **Fixed**, and measured: the CRAM reader is now built with a `.crai` pruned to the slices whose `alignment_start`/`alignment_span` can overlap the query interval, so noodles seeks only to containers the region needs; per-record filtering is untouched. Same 30 kb window (`del:chr20:38412500-38422500`, `--seed 1`) on a 10 Mb chr20 CRAM of 338 slices: extraction 42.79 s → 0.47 s (91x), against 0.04 s from the equivalent BAM and 0.016 s for `samtools view -c`; whole run 91.0 s → 49.4 s. Output is unchanged — R1/R2 FASTQ and truth VCF are byte-identical to the slow path on four windows, including two straddling a slice boundary and the first and last slices of the file. Two notes on the entry above: the "~300x" ratio came from a 0.12 s BAM floor, and measured here the CRAM/BAM ratio was 42.79/0.040 ≈ 1070x before and ≈12x after; and the 49 s that remain are two LOH pileup queries (`loh.rs:507, 910`), which open CRAM the same way and are outside M15's scope.

### M16 · LOH pileup memory ~1 GB per Mb
`loh.rs:580-583` stores one `(u64, u8)` per aligned base per read.
- 3 Mb DEL: 4.58 GB peak RSS with LOH vs 1.37 GB without. Arm-level events will run out of memory.
- **Fix:** call het sites from counts first, then record read alleles only at those sites.
- **Fixed** that way: 3 Mb DEL at 1.27 GB peak RSS, with SNPs read over the whole footprint.

### M17 · `scripts/validate_pipeline.sh` no longer runs
- Step 3 aborts: the filtered truth has overlapping DELs, and the script doesn't pass `--allow-overlap` or filter them.
- `$(grep -vc '^#' f || echo 0)` gives `"0\n0"` on empty input; the `[[ -lt 5 ]]` guard then errors silently and the run continues with 0 events (lines 252, 264, 315, 570).
- Paths at lines 33-35 point to old `data/giab_hg38/` locations. Tool paths are hard-coded to `/home/parlar_ai/dev/sv_caller/.pixi`.
- **Fixed**, and measured: the harness now runs to completion. chr20:61.9-64.2 Mb, 8 het DELs, 3 VAFs, NA18488 background — 5 min wall, exit 0, recall 0.500 / 0.375 / 0.500 at VAF 0.5 / 0.25 / 0.1. Before: exit 1 at step 0 (reference not found); with the paths corrected, exit 1 at step 3 ("overlapping events detected", 2 overlapping pairs of 26).
- The `"0\n0"` guards are fixed and now stop the run: with 5 events and `--min-events 6` the script exits 1 (it previously printed `[[: 0\n0: syntax error in expression` and carried on).
- Step 4 no longer re-implements the merge: it runs spike's own `align.sh` and `merge.sh`, so the second half of **M5** (`-L events.bed -U`, hard-coded `SM:SPIKE`) is gone from the harness too, and Delly reports `Sample:NA18488`.
- Correction to the read counts first published here ("600 608 reads from a 601 429-read background — no read lost or duplicated"): those two totals count supplementary records differently, so they were not comparable. Re-measured at VAF 0.5 — background 601 429 records (599 048 primary, 2 381 supplementary) → merged 600 608 (598 142 primary, 2 466 supplementary): **906 primary records fewer**, not zero. Nothing is lost in the merge: `merge.sh` removes 15 068 original pairs by name and spike writes 14 615 back (13 800 originals unchanged, 1 268 suppressed, 815 new ALT-junction pairs `ev*`), and 1 268 − 815 = 453 pairs = 906 records, exactly the difference. The suppressed pairs are the deleted haplotype's own coverage — that is what injecting a het DEL does — not reads the harness dropped: spike's log reads "15068 read pairs (replaced_reads.txt); 0 of them were dropped for unusable quality".
- **Fix pass: the verdict still could not fail.** The only unconditional gates were "truvari recovered 0 at the highest VAF" and "`spike validate` passed 0 checks", and the unspiked background satisfies both on its own — Delly recovers 3 of the 8 truth DELs from the background alone, and `spike validate` run on the background BAM with the same truth VCF passes **7/19**, one *more* than the spiked BAM's 6/19. A run in which spike contributed nothing therefore printed `VALIDATION PASSED` and exited 0, while `README.md` claimed the opposite. Step 7b now runs the same `delly call` and `truvari bench` on the background BAM, writes it as the `background` row of the summary, and step 8 fails the run unless the highest VAF recovers at least `--min-gain` (default 1) more truth events than that control. Measured: legitimate run TP-base 4 vs control 3 → exit 0, and the one event it adds is `sim_del_6` (chr20:63636171), exactly the one truth DEL NA18488 does not carry; the same harness with the merged BAM replaced by the background itself → TP-base 3 vs control 3, `VALIDATION FAILED`, exit 1. Known limit of the gate: it counts events, so on one replicate it cannot separate a very weak spike-in from re-alignment noise — a `--vafs 0.001` run still scores 4 vs 3, but the extra event is `sim_del_8`, a DEL the background carries that re-alignment flipped from FN to TP. Catching that needs replicates or a titration.
- Also in the fix pass: step 3 now skips only when the directory holds everything step 4 needs (`align.sh`, `merge.sh`, `replaced_reads.txt` as well as the FASTQs and the truth VCF) — `data/validation/spike_vaf_0.5` ships only the latter, so `--outdir data/validation` used to skip step 3 and then abort in step 4 with "re-run step 3", a loop; the default `--outdir` moved off the tracked `data/validation` to `validation_run/` and the script now refuses an `--outdir` whose contents git tracks; step 2 no longer needs gawk (POSIX `match()` + `RSTART`/`RLENGTH`, so mawk no longer reports "Only 0 events after filtering"); the summary reports base-side `TP-base` next to the base-side `FN` and `recall`; and `align.sh`'s `grep -c "SA:Z:" ... || echo "0"` — the same `"0\n0"` idiom, which fires on any sim.bam with no supplementary alignments — is fixed in `main.rs:742`.
- **Fix pass 2: the verdict gate was still bypassable.** The gate lived inside `if [[ -f "${truvari_out}/summary.json" ]]`, so a highest VAF with *no* truvari output never reached it: `--skip-to 8` over an outdir whose `spike_vaf_0.5/truvari` is absent printed `0.5  8  28  0  0  0  N/A  N/A  N/A  6/19` and then `VALIDATION PASSED`, exit 0 — step 7's complaint about a missing summary only covers the case where step 7 ran. Step 8 now calls `note_failure` for the highest VAF whenever it has no truvari summary, whichever steps ran: same command, same outdir, now exit 1 ("no truvari summary ... so this run measured nothing that could be attributed to the spike-in"). `--skip-to` is also range-checked now (`--skip-to 9` used to skip every step and reach the same verdict; it exits 1 with a message). A complete outdir still passes: the same `--skip-to 8` on the untouched run gives TP 4 vs control 3, exit 0, and a fresh end-to-end chr20:61.9-64.2 Mb run (1 min 23 s, `--vafs 0.5`) gives background 3/8, VAF 0.5 4/8, `VALIDATION PASSED`, exit 0.
- Also in fix pass 2: `check_outdir_not_tracked` resolved a relative `--outdir` against `$PROJECT_DIR` while everything else resolved it against `$PWD`, and it was wrong in both directions — from `scripts/`, `--outdir ../data/validation` was **accepted** (git answered "outside repository", the error was swallowed by `2>/dev/null`, and the run wrote into the tracked fixture directory), while from `/tmp`, `--outdir data/validation` was **refused** although the target was `/tmp/data/validation`. `OUTDIR` is now canonicalised with `realpath -m` once, before the guard and before `mkdir -p`; the guard asks git only about paths inside the worktree and treats a git failure as unsafe rather than as "not tracked". Measured after: from `wt/scripts`, `--outdir ../data/validation` → exit 1, "holds files git tracks", nothing written; from `/tmp`, `--outdir data/validation` → runs and exits 0. Also: `--min-gain 0` never disabled the gate (`tp >= control_tp + 0` still requires a tie, verified `tp=0 ctl=3 gain=0` → rc=1), so `README.md` and `--help` were corrected to say it is the weakest setting rather than an off switch; `OUTDIR` is documented as an environment override; the slice-detection regex is anchored on `$CHROM`, so a user BAM named `sample.run_1_2.bam` is no longer taken for one of this script's own slices (before: `discover_background_bam` returned 2, "Several BAMs ... pass --background-bam"; after: rc 0, it picks the whole BAM); the tracked-paths list in the refusal message is quoted; and step 7b's cached control is keyed on the background BAM, reference, `--region` and truth VCF, so a re-run in the same outdir with a changed window no longer scores a fresh spiked number against a stale floor.

## Found during the fix run

### N1 · `spike validate` 6/19 is the harness's input, not the checks

*Found while fixing M17, not part of the original review. Measured on
`/home/parlar_ai/spike-review-run/scratch/work/task-07/final`: chr20:61.9–64.2 Mb,
NA18488 background, 8 het HG002 DELs, VAF 0.5. Recorded, not fixed — `validate.rs`
is not the defect, so nothing in it was changed.*

The obvious reading of "6/19 checks passed" is that the checks are too strict or
the spike-in too weak. Neither is true. Running the same `spike validate`, same
truth VCF, on the **unspiked background BAM** passes **7/19** — one *more* than
the spiked BAM:

| event | check | background | spiked | expected |
| --- | --- | --- | --- | --- |
| chr20:61943513-61945040 | coverage_ratio | 0.04 | 0.04 | 0.50 |
| chr20:61946250-61947090 | coverage_ratio | 0.03 | 0.09 | 0.50 |
| chr20:62057603-62058413 | coverage_ratio | 0.08 | 0.09 | 0.50 |
| chr20:63093345-63094243 | coverage_ratio | 0.00 | 0.00 | 0.50 |
| chr20:63134604-63135240 | coverage_ratio | **0.45 pass** | **0.28 pass** | 0.50 |
| chr20:63636171-63636676 | coverage_ratio | 1.00 | **0.60 pass** | 0.50 |
| chr20:63964828-63965924 | coverage_ratio | 0.12 | 0.07 | 0.50 |
| chr20:64127245-64127815 | coverage_ratio | **0.20 pass** | 0.17 | 0.50 |
| chr20:63093345-63094243 | split_reads | **24 pass** | **25 pass** | >=2 |
| chr20:63134604-63135240 | split_reads | **2 pass** | 0 | >=2 |
| the other six | split_reads | 0 | 0 | >=2 |
| [global] | insert_size / dup_rate / mean_mapq | **455±111 / 10.1% / 58.7, all pass** | **identical, all pass** | — |

- **`coverage_ratio`.** Seven of the eight truth DELs are *already depleted in
  the background* — NA18488 carries them too. Spike then correctly multiplies
  what is there by (1 − VAF), so the post-spike ratio tracks
  `background_ratio × 0.5`, while the M10-strengthened check expects
  1 − VAF = 0.50 ± 0.30 against a clean background. The one event the background
  does **not** carry (chr20:63636171, background ratio 1.00) is exactly the one
  whose `coverage_ratio` passes cleanly, at 0.60. The check is behaving
  correctly; the harness's truth-event selection is confounded.
- **`split_reads`.** A different cause. At chr20:63636171 spike emitted 188 ALT
  fragment pairs (`ev0006*`, 376 records): 0 unmapped, 0 supplementary, **0 with
  an `SA:Z` tag**, 11 soft-clipped records with a maximum clip of 26 bp, and 16
  pairs whose TLEN carries the 505 bp deletion. `check_split_reads`
  (`validate.rs:493-520`) requires ≥ 2 SA alignments landing at the partner
  breakpoint, so it reads 0. Of the ~8 013 SA-tagged records in the merged BAM
  only **16 / 4 / 2** are synthetic (`ev*`) at VAF 0.5 / 0.25 / 0.1 — flat in
  VAF, i.e. the SA signal is background-derived. The single passing
  `split_reads` (25 reads at chr20:63093345) is at the event NA18488 carries
  homozygously; the background BAM alone scores 24 there. For sub-2 kb DELs the
  M12 check looks for supplementary alignments that spike plus bwa-mem2 defaults
  do not produce.
- **The three global checks** are byte-identical between the two BAMs: they are
  whole-BAM statistics of a 600 608-record BAM of which 29 381 records are
  synthetic, so the spike-in cannot move them.

**What would settle the rest:** (1) make the per-check background baseline above
part of the harness and read every check as a delta rather than as an absolute;
(2) re-run the harness restricted to truth DELs whose background depth ratio is
≈ 1.0, which isolates the check from the confound.

Two related decisions are **open, and deliberately not taken here**:

- The M17 verdict gate counts recovered truth events, so a `--vafs 0.001` run
  can still pass because re-alignment flipped a background-carried DEL from FN
  to TP (`sim_del_8`, measured above). No gate over TP counts or TP *sets* can
  close that while the truth events are DELs the background carries: it needs
  (2) — selecting truth events the background does not carry — not a stricter
  threshold.
- `spike validate`'s contribution to the verdict ("the report parsed and at
  least one check passed") is satisfiable by the background alone: 7/19 for the
  background against 6/19 for the spiked BAM. Turning it into a delta gate
  would fail legitimate runs, since the spiked BAM legitimately scores *lower*
  on the confounded checks. The honest alternative is to drop the `spike
  validate` check count from the verdict altogether rather than leave an inert
  gate in it — a judgement call for a human, recorded here rather than made.

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
