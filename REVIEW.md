# spike — code review

Date: 2026-09-24 · Reviewed: working tree on `master` at `f5428ce` plus uncommitted edits to `extract.rs`, `loh.rs`, `simulate.rs`, `synth.rs`.

**Verdict:** the code is tidy and hard to crash, but in several common cases it produces quietly wrong output: the realized allele fraction is not what was asked for, or the variant lands in the wrong place. For a tool whose job is to produce truth sets, these matter most, because nothing downstream will flag them. Eight high-severity issues below; fix those before trusting any benchmark made with spike.

## How this was checked

- The code was read in four parts: CLI and input parsing; haplotype and simulation core; read synthesis and I/O; LOH, truth VCF and `validate`.
- Every finding was reproduced in a scratch copy: a unit test, a harness calling the real functions, or a real run on chr20 data (HG002 NovaSeq 35x, NA18488, GRCh38). The numbers below are measured, not estimated.
- Findings marked *plausible* are believed but not fully proven.
- Apart from adding this file, nothing in the repo was changed.

## Health check

This table is a record: the "At review" column is the original measurement
and is never edited. The two right columns were added later (bookkeeping fix,
2026-09-25) once the original numbers had gone stale from the fix run's own
commits — measure both at the row's own command, not by scaling the "At
review" figure.

| Check | At review (`master@f5428ce` + uncommitted) | At branch base (`8d1beba`) | Current (`review-fixes-2` @ `6e0aa57`+, the whole-branch fix pass) |
| --- | --- | --- | --- |
| `cargo build --release` | OK, 1 warning (unused `primary_chrom`, `is_within_single_segment` in `haplotype.rs`) | not re-measured | OK, **1** warning (unused `is_within_single_segment` in `haplotype.rs`) |
| `cargo test` | 128 passed, 0 failed | **173** passed, 0 failed | **394** passed, 0 failed |
| `cargo clippy --all-targets` | Style only: 6× `is_multiple_of`, 4× too many arguments, 2× use `?`, 1× no-effect op, 1× range loop, 1× manual `contains` | **13** (bin) / **14** (test target, 12 duplicates) | **13** (bin) / **14** (test target, 12 duplicates) — unchanged from base; every fix in this run and in the whole-branch pass held the line here |
| `scripts/validate_pipeline.sh` | Broken (see M17) | Broken: exit 1 at step 0, reference not found (M17 fix `29ec590` had not landed yet — `8d1beba` is its ancestor) | **Fixed** (`29ec590` M17/M5; hardened by `3e85a0d`, then `61af374`): runs end to end; fails (exit 1) when the spike-in contributed nothing the background does not already carry; and fails (exit 1) rather than printing `VALIDATION PASSED` when the highest VAF has no truvari summary to grade at all |

The tests pass, but most would still pass with the high-severity bugs below. See [Test gaps](#test-gaps).

## Summary table

| ID | Severity | Problem | Where |
| --- | --- | --- | --- |
| H1 | High | **Fixed.** Nearby events undo each other's read suppression | `main.rs:431-437`, `main.rs:517-530` |
| H2 | High | **Fixed.** Exon numbers ignore strand | `exon.rs:93-103` at `master@f5428ce`; the numbering is now `number_exons`, `exon.rs:144-192` |
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
| M14 | Medium | **Fixed.** Missing base qualities → invalid FASTQ. Fix pass 2: the dropped pairs are now counted in the generated run README, per event and in total, instead of only on stderr | `extract.rs:453-460`; `main.rs` `write_readme` |
| M15 | Medium | **Fixed.** CRAM extraction ~300× slower than BAM | `extract.rs:240-243, 314-317` |
| M16 | Medium | **Fixed.** LOH pileup memory ~1 GB per Mb | `loh.rs:580-583` |
| M17 | Medium | **Fixed.** `validate_pipeline.sh` no longer runs | `scripts/validate_pipeline.sh` |
| L1–L19 | Low | Parsing edge cases, robustness, minor I/O | see [Low](#low-severity) for which are fixed |
| N1 | Medium | **Not fixed** (found during the fix run). `spike validate` scores a cross-sample spike-in against a confounded background, and `split_reads` looks for a signal spike does not emit | `validate.rs:440-520`; `scripts/validate_pipeline.sh` |
| N3 | Medium | **Fixed** (found during the fix run). Five more CRAM query sites walked the whole chromosome's index | `loh.rs:507, 910`; `validate.rs:708, 822, 950` |
| N4 | Medium | **Fixed** (found during the fix run). The same five CRAM query sites also read another contig's records out of a shared container | `count_alleles_cram`, `collect_snp_alleles_cram` (`loh.rs`); `count_depth_in_region`, `split_reads_to_partner`, `pileup_region` (`validate.rs`) -- function names, because the line numbers this row first carried have drifted twice |
| N5 | High | **Fixed** (found during the fix run). An empty or near-empty donor pool was simulated from anyway: exit 0 with a truth VCF and 2 invented read pairs beside it. The pool-size guard alone left the same symptom reachable through a second door (aggregate pool vs. coverage at the breakpoint); now closed where the coverage is measured, at every breakpoint side rather than the first only (N12) | `main.rs:492`, `extract.rs:497`, `simulate.rs:409, 429` (at `66b45a5`); `simulate.rs:203-210, 379-495, 519-521, 573-587` (now) |
| N6 | Medium | **Not fixed** (found during the fix run). Four of `BamStats`'s five fields are read nowhere but its own log line, and one of them, `mean_coverage`, is wrong by ~7000x -- every real BAM prints `est_coverage=0.0x` | `bam_stats.rs:6-17, 258-275`; `main.rs:375` |
| N7 | Medium | **Not fixed** (found during the fix run). A quality profile with 0/1208 usable base-conditioned bins is used without a warning | `synth.rs:92, 199-222` |
| N8 | Medium | **Fixed** (found during the fix run). No `validate` check covered INS, and an uncovered event is a *failed* result, so any truth VCF holding an INS could never report all-PASS -- spike's own round trip, broken for insertions. `ins_reads` now counts reads whose alignment leaves the reference at POS | `validate.rs:133-180` (at `39d9773`); `validate.rs:137-190, 631-686, 1068-1101, 1417-1500` (now) |
| N9 | High | **Fixed** (found by the whole-branch review). `validate`'s per-event `allele_freq` answered `pass: true` on three questions it had not asked -- any indel or MNV, a pileup depth below 5, a non-ACGT alt -- and `load_truth_events` routed unrecognised SVTYPEs into the same arm, so `<CNV>` passed as an indel. A truth record with `END <= POS` PASSed `coverage_ratio` over a region no query read | `validate.rs:601-609, 630-639, 645-655, 397, 1083-1085` (at `39d9773`) |
| N10 | Medium | **Fixed** (found by the whole-branch review). No `validate` check measured a small indel's or an MNV's allele fraction, so once N9 stopped calling them a pass a truth VCF holding one could not report all-PASS -- the same shape as N8, for `snp:` events with multi-base REF or ALT. `allele_freq` now picks a counting rule from the REF/ALT shape: a del/ins/MNV run on the chr20 slice goes from **3/6 PASS, exit 1** to **6/6 PASS, exit 0** | `validate.rs:753-762` (at `6e0c49a`) |
| N12 | Medium | **Fixed** (found while closing N10). N5's donor-coverage refusal measured the **first** breakpoint only, so the same fusion was refused or accepted depending on which partner was named first. Now every breakpoint side is measured, scoped to the loci the pool was extracted from: both sides for a fusion, at least one for a single-locus event | `simulate.rs:196-220` (at `ad9881e`); `simulate.rs:203-210, 379-495` (now) |
| N11 | Low | **Fixed** (found by the whole-branch review). Two `--help` strings contradicted the code (`--allele-fraction (0.0-1.0)` where 0 is refused; `--flank` silent about its 2000 minimum), and spike's refusals were scattered across nine README locations with four not documented at all | `main.rs:93, 120-123` (at `39d9773`) |
| N13 | Critical | **Fixed** (found by the verification review of the fix wave). `cigar_indel_vote`'s deletion **dead zone**: a `D` operation shifted 1..=`indel_len` bases from the junction swallows one of the two reference bases the vote was anchored on, so the read entered **neither** count. `INDEL_POS_PAD = 10` promised a tolerance the code did not deliver, and the same physical 2 bp deletion spelled one repeat unit off left-alignment read **0.04** where the left-aligned spelling read **0.38** -- at `SIM_VAF=0.10` the wrong spelling PASSes and the right one FAILs | `validate.rs:971-1023` (at `99f1a8e`) |

## High severity

### H1 · Nearby events undo each other's read suppression

`main.rs:431-437` collects every event's `kept_originals`, then `dedup_by_name` (`main.rs:517-530`) keeps the last copy of each read name. Each event returns *all* unsuppressed reads in its extraction window (event ± `--flank`, default 10 kb), not just its own footprint. A read suppressed by event 1 comes back through event 2's list.

- DEL chr1:[50000,52000) + SNP at 58000 (no overlap, so allowed): depth in the deleted region **1.000×** (expected 0.5×). The SNP comes out at AF 0.37 instead of ~0.5.
- Two SNPs 5 kb apart on chr20: the first drops from AF **0.556 to 0.355**. 78 of the 83 reads it suppressed reappear.
- With `--region`, every event on that chromosome interacts. Two SNPs 40 kb apart: AF 0.355 and 0.429.

**Fix:** gather suppressed read names from all events into one set and remove them from the final output, instead of keep-last dedup. Or run events in sequence over one shared pool.

### H2 · Exon numbers ignore strand

`exon.rs:99-103` numbers exons 1..n by genomic start. (That range is `master@f5428ce`'s. **Do not follow it into the current tree:** `exon.rs:99-103` there is L14's gene-symbol fix, an unrelated change this document cites separately under L14. The numbering now lives in `number_exons`, `exon.rs:144-192`.) The BED has no strand column, and the exon number in the name (`TP53_exon1`) is ignored. Duplicate transcript lines also shift the numbering. Fusion breakpoints (`exon.rs:392-397`) always take gene A's genomic-left side.

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
- **Fix pass 2 (whole-branch review): a dropped pair was invisible in every
  output file.** Fix pass 1 (b) put the dropped names in `replaced_reads.txt`
  so `merge.sh` removes them, but nothing there or anywhere else distinguishes
  them from replaced pairs, and the only surface was one stderr INFO line
  (`main.rs:579-584`) -- `write_readme` was never given the count. A run whose
  log has scrolled past leaves a localized depth dip with no record, and a
  localized dip is what a depth-based CNV caller reads as signal. Measured on
  the chr20 slice with `QUAL` set to `*` on every 20th record -- **5.0%** of
  records (709 of 14,184) -- over `del:chr20:38412500-38422500 --flank 2000
  --seed 1`: the pool holds 1855 pairs and **201** are dropped, **9.8%** of the
  2056 pairs the window would otherwise have held, because a pair is lost when
  *either* of its records is stripped. The generated run README now carries a
  `Dropped (unusable quality)` column per event and a total line saying how
  many of them are removed with nothing put back. `R1.fq.gz`, `R2.fq.gz`,
  `truth.vcf` and `replaced_reads.txt` are byte-identical before and after the
  column existed (`a052cfcb`, `1c924957`, `5a73bfd8`, `9a23f32f`).
- **N5's message misdiagnosed the failure it now catches, and contradicted
  itself.** A CRAM storing qualities via read features has every record dropped
  by `quality_is_missing`, the pool is then 0, and N5's guard fires -- the right
  outcome -- but the message named three causes (coverage, `--region`,
  `--min-mapq`) and never quality. It also read "has **no** usable donor reads:
  12 read pair(s) extracted", which is false whenever the pool is non-empty.
  `extract_pool_for_event` now measures the per-event delta of
  `unusable_qual_names` and hands it to `finish_donor_pool`, and the message
  opens "has **too few** usable donor reads". Measured on the same slice with
  every `QUAL` set to `*`: before, `has no usable donor reads: 0 read pair(s)
  extracted from chr20:38410500-38424500, fewer than the 30 spike needs. Check
  that the event lies in a covered region ...`; after, `has too few usable donor
  reads: 0 read pair(s) extracted from chr20:38410500-38424500, fewer than the
  30 spike needs (**2097** record(s) in those windows were dropped for unusable
  base qualities and are not in that count)`, and the remedy list ends with the
  CRAM read-feature case.

### M15 · CRAM extraction ~300× slower than BAM
noodles-cram 0.74 `Query::read_next_container` checks only the reference id, never position, so every container on the chromosome is decoded (`extract.rs:240-243, 314-317`).
- 30 kb region: 0.12 s from BAM, 39.97 s from a CRAM with only 10 Mb of chr20. `samtools view`: 0.015 s.
- **Fix:** filter index entries by position and seek manually, or upgrade noodles after checking the fix.
- **Fixed**, and measured: the CRAM reader is now built with a `.crai` pruned to the slices whose `alignment_start`/`alignment_span` can overlap the query interval, so noodles seeks only to containers the region needs; per-record filtering is untouched. Same 30 kb window (`del:chr20:38412500-38422500`, `--seed 1`) on a 10 Mb chr20 CRAM of 338 slices: extraction 42.79 s → 0.47 s (91x), against 0.04 s from the equivalent BAM and 0.016 s for `samtools view -c`; whole run 91.0 s → 49.4 s. Output is unchanged — R1/R2 FASTQ and truth VCF are byte-identical to the slow path on four windows, including two straddling a slice boundary and the first and last slices of the file. Two notes on the entry above: the "~300x" ratio came from a 0.12 s BAM floor, and measured here the CRAM/BAM ratio was 42.79/0.040 ≈ 1070x before and ≈12x after; and the 49 s that remain are two LOH pileup queries (`loh.rs:507, 910`), which open CRAM the same way and are outside M15's scope (tracked and fixed as N3).
- Caveat on that equality proof: it is **single-contig only** — all 338 slices of the test CRAM are chr20. On a multi-contig CRAM the output *can* change, because pruning decodes fewer foreign containers and so leaks fewer of their reads (L2). That direction is an improvement, but it is unproven here.
- **That last sentence is now measured, and it was wrong** (see L2). Pruning changes the leak by nothing at all: htslib writes one `.crai` line per contig for a multi-reference container, all at the same offset, so the queried contig's own entry survives pruning and drags the whole container in with it. Real chr20+chr21 CRAM (`multi_seq_per_slice=1`), `del:chr20:38412500-38422500 --seed 1`: 7875 pairs extracted at `8d1beba`, **7875** at `3a782cb` with pruning in place, 4132 after L2.

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
  **Update after L15 (task 23), measured on a rebuild of the same two BAMs:**
  they are no longer whole-BAM statistics — they are read from the eight truth
  events' own windows — and no longer identical. Background → spiked:
  insert_size 447+/-130 → 449+/-129, dup_rate 9.2% → 9.7%, mean_mapq
  53.1 → 52.7 (the same two BAMs under the pre-L15 binary: 455+/-111 / 10.1% /
  58.7, identical to the digit). So the checks now respond to the spike-in, but
  the inertness this entry is about is **reduced, not removed**: all three
  still pass on the unspiked background, and the run still scores 7/19
  background against 6/19 spiked. A check that moves by 0.4 MAPQ is a
  measurement, not a gate; making the globals part of a verdict still needs the
  per-check background baseline of (1) below.

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

### N3 · Five more CRAM query sites walk the whole chromosome's index

*Found while fixing M15, which names `extract.rs:240-243, 314-317`; those lines
are fixed. These five are the same five-line pattern —
`indexed_reader::Builder::default()` → `read_header()` → a bounded `Region` →
`query` — with no index pruning, and nothing tracked them.*

`loh.rs:507` (`count_alleles_cram`) and `loh.rs:910` (`collect_snp_alleles_cram`)
run on **every event**, with or without `--gvcf`. Measured on the same 10 Mb
chr20 CRAM (338 slices), `del:chr20:38412500-38422500`, `--seed 1`, with M15
already fixed: the two pileup passes take **23.67 s and 23.65 s** of a 50.4 s
run (the same two passes cost 0.027 s and 0.018 s from the equivalent BAM), so
94% of what is left is this defect and a CRAM run is still ~30× the 1.70 s BAM
run. M15's own cited lines are genuinely fixed, but its stated symptom — CRAM
far slower than BAM — survives at most of its original cost.

`validate.rs:708` (`count_depth_in_region`), `validate.rs:822`
(`split_reads_to_partner`) and `validate.rs:950` (`pileup_region`) are the same
defect in the `spike validate` subcommand.

- **Fix:** the same index pruning, through `extract::open_cram_reader_for_region`
  (made `pub(crate)`; each site builds its `Region` before opening instead of after).
- **Fixed**, and measured on the same 10 Mb chr20 CRAM (338 slices), `--seed 1`:
  `del:chr20:38412500-38422500` whole run **50.4 s / 50.8 s → 3.1 s / 3.6 s**
  (16×), against a BAM floor of 1.7 s / 2.1 s that the fix leaves untouched — a
  CRAM run is now 1.7× a BAM run, not 28×. The two LOH pileup passes go
  **23.67 s → 0.18 s** and **23.65 s → 0.17 s**. `spike validate` on the same
  CRAM: **97.3 s / 99.0 s → 3.6 s / 3.6 s** with a DEL truth VCF
  (`count_depth_in_region` + `split_reads_to_partner`), **40.6 s → 3.2 s** with a
  SNP truth VCF (`pileup_region`). Output is byte-identical everywhere: R1, R2,
  `truth.vcf` and `replaced_reads.txt` over seven windows, six of which put a
  query edge exactly on a slice boundary (region end = a slice's first base;
  end = first base − 1; region start = a slice's last base; and the same on the
  file's first and last slices), plus `spike validate`'s stdout on both truth
  VCFs, with the two `allele_freq` pileups aimed at single bases that are
  themselves slice starts. The 1 bp windows discriminate: shifting the event by
  one base changes the extraction (4561 vs 4560 pairs), and both match their own
  unpruned run.

### N4 · The same five CRAM query sites read another contig's records

*Found while fixing L2, which names `extract.rs` only; those two loops are
fixed. `count_alleles_cram` and `collect_snp_alleles_cram` in `loh.rs`, and
`count_depth_in_region`, `split_reads_to_partner` and `pileup_region` in
`validate.rs`, iterate a CRAM `Query` the same way and had the same hole.
They are the same five sites N3 covers. (The line numbers this row first
carried -- `loh.rs:656, 1056`; `validate.rs:712, 823, 948` -- were right when
it was written and are not any more: at `39d9773` `validate.rs:712` is
`&mut self,`. **Follow the function names, not the numbers.** At this commit
the five are `loh.rs:673, 1087` and `validate.rs:1318, 1439, 1682`; L15 added
a sixth guarded site, `sample_region` at `validate.rs:1045`, and N8 a seventh,
`reads_with_inserted_sequence` at `validate.rs:1532`.)*

This is worse than L2 rather than a milder copy of it. L2 diluted the donor
pool; these decide the **truth set**. `loh.rs`'s two passes call the
heterozygous SNPs a region is phased on and then assign each fragment to a
haplotype, which is what decides whose originals are suppressed — and
`validate` is the check that would otherwise notice, so a leak there hides
the leak here.

- **Fix:** `extract::record_is_on_queried_reference` made `pub(crate)` and
  applied at all five, ahead of the flag and MAPQ filters.
- **Fixed**, and measured. HG002 chr20:38.40–38.43 Mb + chr21:38.40–38.43 Mb
  written as one multi-reference CRAM (`multi_seq_per_slice=1`, both contigs in
  one container), `del:chr20:38412500-38422500 --region
  chr20:38405000-38425000 --seed 1`. LOH called **10 510 het SNPs → 15**, in
  **1 phase block → 4**, classifying **4152 fragments → 543**, and wrote
  **4424 read pairs → 3435**. chr21's bases, laid over chr20's reference, had
  made almost every position look heterozygous. Every one of those "after"
  numbers is exactly what the same reads give from a **chr20-only** CRAM, and
  R1, R2, `truth.vcf` and `replaced_reads.txt` are byte-identical to that
  control — the foreign contribution is zero, not merely smaller. Only R1/R2
  actually discriminate here, though: verified directly against the task's
  fixtures (`loh_before_both` vs `loh_before_chr20`) that `truth.vcf` and
  `replaced_reads.txt` were **already** byte-identical to the chr20-only
  control *before* this fix existed — they don't encode per-record
  reference-id information, so this leak could never show up in them.
  Decompressed R1 differs before the fix (multi-ref md5 `4dbe610a…` vs
  control `2fef3416…`) and matches after (both `2fef3416…`).
- **`spike validate` gave a false PASS**, measured on a CRAM whose chr21 reads
  sit only inside the deletion: `coverage_ratio` **1.77 → 0.91** (control
  0.91), so a region with no coverage drop looked like one with a 77% *rise*;
  and `allele_freq` at chr20:38415000, where chr20's reads read G and chr21's
  read C, **0.54 PASS → 0.00 FAIL** (control 0.00 FAIL) — the check passed on
  an unspiked CRAM entirely on the other chromosome's bases.
- **The mate's reference is deliberately not checked at these five**, unlike
  the two extraction loops. They judge records one at a time, and a read on the
  queried contig whose mate lies elsewhere is a genuine record of this contig —
  **27** of 6865 in the test window under spike's own filter set (`-q 20 -F
  0xF04`; re-verified directly against `real_chr20.cram`, `chr20:38405000-
  38425000` — 28 is what `-F 0x904` gives instead, i.e. leaving
  duplicate-flagged records in, which those five sites do drop) — that every
  BAM query returns. Measured: with
  the mate clause added, the same **single-contig** CRAM run classifies **542**
  fragments where both the BAM run and the record-only guard classify **543**.
  So `record_is_on_queried_reference` compares the record's reference id only,
  and `pair_is_on_queried_reference` adds the mate clause for extraction.
- Each guard is independently observable: removing any one of the five reddens
  its own test and only that one. The mate clause is the exception, and is
  reported under L2.

### N5 · An empty donor pool is simulated from anyway

*Found while reviewing the L19 fix, not part of the original review. The
[Design notes](#design-notes) already carry the same observation untracked --
"`n.max(2)` emits 2 pairs even at zero coverage. An empty read pool silently
gives a constant-Q20 profile" -- with no numbered row, no task and no
measurement. This row is that note, numbered, measured and fixed.*

`extract.rs:497` logs `Built read pool: 0 pairs` and continues. `main.rs:492`
then trains `QualityProfile::from_read_pairs` on the empty slice, `stats.rs`
substitutes its default 400/80 fragment distribution, and `simulate.rs:429`
scales the tiling count by a coverage of 0 -- which `n.max(2)` turns into 2
fragments. Nothing between extraction and the FASTQ writer asked whether there
were any donor reads at all: there was no `pairs.is_empty()` check anywhere in
`main.rs` or `extract.rs`. This is reachable on ordinary input -- an event in a
zero-coverage region, an off-target panel BAM, a mistyped `--region`.

Measured on branch HEAD `66b45a5`, `del:chr20:30000000-30010000 --seed 1`
against the chr20 37.5-41.5 Mb HG002 slice (the event lies outside the slice,
so the window has no coverage):

```
INFO spike::extract] Extracted 0 complete read pairs from chr20:29990000-30020000
WARN spike::stats] No valid insert sizes found, using default distribution (mean=400, sd=80)
INFO spike::extract] Built read pool: 0 pairs
INFO spike::synth] Quality profile: 0 pairs, 151 cycles. R1 mean Q: start=0.0 mid=0.0 end=0.0,
     R2: start=0.0 end=0.0. Base-conditioned bins: 0/1208 usable.
     Markov bins: base 0/4832, cycle 0/1208 usable
INFO spike::simulate] Tiling 2 synthetic reads across 4000bp haplotype (cov=0.0, vaf=0.50, bp_only=false)
```

**exit 0**, with `truth.vcf`, `R1.fq.gz`, `R2.fq.gz`, `events.bed`, `align.sh`,
`merge.sh` and `README.md` all written. The two pairs are real reference
sequence at a locus the BAM has no read over, and all 604 of their quality
bytes are the single character `5`.

**On the Design note's "constant-Q20": that phrase names two different
measurable numbers, and they disagree.** The *profile's* mean is Q0 --
`R1 mean Q: start=0.0 mid=0.0 end=0.0` above, the mean of zero observations --
while the *emitted* quality is Q20: with every bin empty, `sample_quality_inner`
falls past all four levels of the hierarchy to `synth.rs:338`,
`b'!' + 20 // last resort: Q20, a fixed byte with no draw at all`. So the note
is right about the FASTQ and wrong about the profile, and reading the log line
alone ("Q0") is right about the profile and wrong about the FASTQ. Neither
value is caught downstream: `5` (byte 53) and `!` (byte 33) both sit inside the
printable Phred+33 range `!`-`~` that M14's `write_paired_fastq` guard checks,
so that guard cannot see either.

Task 26's `test_paired_bam_without_proper_pairs_is_not_mistaken_for_single_end`
deliberately asserts this case is `Ok`. That is right for its own question --
the library *is* paired, and capping the scan there would truncate a paired
BAM's fragment-length estimate (L15) -- but it means the downstream run is the
only thing left that could notice, and it did not.

- **Fix:** a donor pool of fewer than `MIN_DONOR_PAIRS` = **30** read pairs is a
  hard error naming the event and every window it searched. The count is taken
  after `dedup_pairs_by_name`, so a fragment two overlapping windows both hand
  in cannot lift a pool over the floor, and before `FragmentDist::from_read_pairs`,
  so the run fails instead of first warning about a substituted distribution.
  **Why 30:** it is `synth.rs`'s own `MIN_BASE_OBS` (`synth.rs:26`), the
  observation count the quality model requires before it will sample from a
  bin, and a pool of *n* pairs puts exactly *n* observations in each cycle-only
  bin -- level 4, the model's final fallback and the only level an ordinary
  pool always reaches. 30 pairs is therefore the smallest pool at which any
  level of the model is trained to the threshold the model itself sets. "Empty"
  is the wrong floor because one pair trains no bin either and its single
  insert size becomes the whole fragment distribution; 30 is a floor on
  "measured from this library at all", not a coverage requirement -- the same
  30 kb window yields 4559 pairs on the 35x HG002 BAM, so it would have to fall
  to about 0.2x before the floor bound.
- **Fixed**, and measured: the same command and seed now exits **1** with
  `event DEL  chr20:30000001-30010000 (10000bp) has too few usable donor reads:
  0 read pair(s) extracted from chr20:29990000-30020000, fewer than the 30
  spike needs (0 record(s) in those windows were dropped for unusable base
  qualities and are not in that count). ...`  (the message was reworded in the
  whole-branch pass -- see M14's fix pass 2 -- and read "has **no** usable
  donor reads" when this was first measured), and leaves `--output` empty -- **7 files written and exit 0
  before, 0 files and exit 1 after**. A working event is untouched:
  `del:chr20:38412500-38422500 --seed 1` on the same slice still extracts 4559
  pairs and gives byte-identical output before and after (`R1.fq.gz`
  `8cf7964f`, `R2.fq.gz` `1bd652f8`, `truth.vcf` `5a73bfd8`,
  `replaced_reads.txt` `6bc6392a`).
- **The residual `n.max(2)` half, now also fixed.** The pool-size guard is an
  *aggregate*: `extract_pool_for_event` sums `all_pairs` over every extraction
  window and checks the total, while the tiling count is scaled by coverage at
  the **first breakpoint** (`estimate_coverage_at`, a 2 kb window). So a pool
  can clear 30 pairs and still measure coverage 0 where it matters, and
  `n.max(2)` emits 2 invented pairs beside a truth VCF at exit 0 -- N5's
  headline symptom, reached through a second door. Two routes measured on the
  chr20 37.5-41.5 Mb HG002 slice at branch HEAD `39d9773`:

  - `del:chr20:30000000-30010000 --region chr20:38400000-38440000 --seed 1`:
    0 pairs from the event window, **6117** from the region window, guard
    passes on 6117, `Tiling 2 synthetic reads ... (cov=0.0, vaf=0.50)`,
    7 files and **exit 0**.
  - `fusion:GENEA:exon1:GENEB:exon2;af=0.2` with side A at chr20:30.00 Mb (no
    coverage) and side B at chr20:38.42 Mb: 0 pairs from A, **3023** from B,
    `Tiling 2 synthetic reads ... (cov=0.0, vaf=0.20, bp_only=true)`,
    **exit 0**. M8's window split makes the `--region` door cheaper to hit.

  **Fixed** by moving the question to where the coverage is measured:
  `compute_tiling_count` returns 0 for a non-positive (or NaN) coverage instead
  of falling through to the floor, and `simulate_event_with_copies` refuses the
  event right after `estimate_coverage_at`. Measured: both commands above now
  exit **1** with `event chr20:30000000-30010000 has no donor coverage at any
  of its breakpoints (chr20:29999999, chr20:30010000): the pool holds 6117 read
  pair(s) but none of them cover that. ...` and leave `--output` empty --
  **7 files and exit 0 before, 0 files and exit 1 after**. (The refusal first
  keyed on the first breakpoint only; [N12](#n12--n5s-donor-coverage-refusal-measured-the-first-breakpoint-only)
  widened it to every breakpoint side, which is what makes the fusion route
  above refused whichever partner is named first.)

- **The floor's other half: it understated the planted fraction.** At
  low-but-nonzero coverage the requested count rounds below 2 and the floor
  raises it, so the reads planted carry a *higher* fraction than `SIM_VAF`
  records -- the direction that flatters a caller. Measured on a 1%-subsampled
  copy of the same slice (`samtools view -s 1.01`) at branch HEAD `39d9773`,
  `del:chr20:38412500-38422500 --allele-fraction 0.05 --seed 1`: pool **56**
  pairs (clears the 30 guard), `cov=0.7`, mean fragment 442.7, so the formula
  asks for `round(0.7 x 0.05 x (4000 - 442.7) / 442.7) = round(0.28) = 0`
  fragments; the floor emitted **2**, suppression removed **1** original, and
  `truth.vcf` recorded `SIM_VAF=0.050`. Exit 0, no warning.

  **The rule now:** the floor of 2 stands wherever the coverage is *real*,
  because a haplotype shorter than one fragment asks for 0 however good the
  coverage is and planting nothing would leave a truth VCF with no reads behind
  it -- that is the case `test_tiling_count_minimum` pins, and it is kept. What
  changes is that the floor no longer applies silently, and no longer applies
  at all when the coverage is zero:

  - coverage 0 (or NaN): the event is refused, exit non-zero, nothing written.
  - coverage > 0 and the request rounds below 2: emit 2, and warn with the
    request, the emission and the `SIM_VAF` the truth VCF will carry.

  Measured after: the same 1%-subsample command still exits 0 and still plants
  the event, now preceded by `WARN spike::simulate] coverage 0.7x at VAF 0.050
  asks for 0 tiled fragment(s); spike emits the 2 it needs to plant the event at
  all, so the realized allele fraction will be above the 0.050 recorded as
  SIM_VAF in the truth VCF`.

  H5 and M1/M2 are untouched by this: the additive formula
  `n = round(cov x v / (1 - v) x breakpoints.len())` and the interior formula's
  scaling by `total_len` with `pool.frag_dist.mean` used directly are the same
  expressions as before -- only the `.max(2)` tail moved into a named helper.
  Re-verified by `test_tiling_count_breakpoint_only_gives_requested_fraction`,
  `test_tiling_count_dup_breakpoint_only` and
  `test_tiling_uses_actual_fragment_length_for_short_inserts`, all still green.

  `test_simulate_event_keeps_pairs_straddling_footprint_edge` was corrected, not
  weakened: its pool was 500 pairs at `[3800,4200)` and nothing else, so the
  breakpoint at `chr1:999` had coverage 0 and the new guard refused it. It now
  also holds 500 `in_` pairs over the breakpoint, and asserts that no `edge_`
  pair is among the suppressed names -- the same invariant, on a pool that can
  actually be simulated from.

### N6 · `BamStats` carries four dead fields, one of which spike prints wrong

*Found while reviewing the L19 fix. Recorded, not fixed: removing public
struct fields, the `info!` line the user reads them from, and the reason
recorded in `scan_is_complete`'s doc comment is a change of its own, and it
interacts with L19's decision (below).*

`main.rs:375` reads `bam_stats.read_length` and nothing else. `insert_mean`,
`insert_stddev`, `mean_coverage` and `records_sampled` are read nowhere in the
crate outside `bam_stats.rs`'s own `info!` line and its tests' assertions.

`mean_coverage` (`bam_stats.rs:258-261`) is `total_records * read_length /
genome_size`, and the two ends of that fraction do not belong together. The
numerator is the *sample*: the scan stops as soon as it holds `sample_size`
insert sizes. The denominator is the whole header genome -- all 195 GRCh38
contigs, 3,099,922,541 bp -- even when the file is a 4 Mb chr20 slice.
Measured, `--seed 1`, branch HEAD:

| input | records sampled | read_len | header genome | printed | actual depth |
| --- | --- | --- | --- | --- | --- |
| HG002 NovaSeq 35x, whole BAM | 101,798 | 151 | 3,099,922,541 (195 contigs) | `est_coverage=0.0x` (0.00496) | ~35x |
| chr20 37.5-41.5 Mb slice | 100,923 | 151 | 3,099,922,541 (same header) | `est_coverage=0.0x` (0.00492) | ~35x |
| that slice subsampled to 1.2% | 185 | 151 | 3,099,922,541 | `est_coverage=0.0x` | ~0.4x |

Because the scan bounds the numerator, the figure is near-constant for any
151 bp paired WGS BAM carrying a GRCh38 header -- about 0.005x whatever the
file's real depth, 7000x low on the two 35x inputs and 80x low on the 0.4x
one. It is not a poor estimate; it carries no information about the file.

**This is also L19's residual.** `scan_is_complete`'s doc gives, as the reason
for not capping the *paired* record scan, that "a cap on the record count
would truncate that sample on a BAM whose proper pairs are sparse" -- that is,
it protects `insert_mean` and `insert_stddev`. The distribution the simulation
samples from is not those: it is built in `stats.rs` by
`FragmentDist::from_read_pairs` over the extracted donor pool
(`main.rs:1142`), and `FragmentDist::from_stats` is called only by
`default_dist()` (a hard-coded 400/80) and by tests, so `BamStats`'s insert
numbers never reach a read. Removing the four dead fields would leave nothing
for the uncapped paired scan to protect, and L19's exception with it.

### N7 · A degenerate quality profile is used without a warning

*Found while fixing N5. Recorded, not fixed: unlike N5's floor this is a
matter of degree, and choosing the census that deserves a warning needs a
measurement of how far the realised quality distribution drifts before it
matters.*

`QualityProfile::from_read_pairs` (`synth.rs:92`) logs its bin census at
`info!` (`synth.rs:199-222`) and never warns, whatever the census says.
Measured on branch HEAD with N5's 30-pair floor already in place -- a real
run, exit 0, nothing on stderr -- against the chr20 38.40-38.44 Mb slice
subsampled to 1.2%, `del:chr20:38412500-38422500 --seed 1 --flank 2000`:

```
INFO spike::extract] Built read pool: 32 pairs
INFO spike::synth] Quality profile: 32 pairs, 151 cycles. R1 mean Q: start=37.0 mid=35.4 end=35.0,
     R2: start=35.8 end=34.2. Base-conditioned bins: 0/1208 usable.
     Markov bins: base 0/4832, cycle 146/1208 usable
```

**0 of 1208** base-conditioned bins and **0 of 4832** Markov base bins reach
`MIN_BASE_OBS`, and 1062 of 1208 Markov cycle bins do not either, so almost
every quality byte is drawn from the level-4 cycle-only marginal: the model
has collapsed to "the per-cycle distribution of 32 pairs" and says nothing
about having done so. N5 closes only the far end of this band (0 pairs, where
nothing is drawn at all); everything between the floor and a well-trained
profile is silent. Even a healthy run is partly so -- the 4559-pair pool of
`del:chr20:38412500-38422500` on the full HG002 BAM still leaves 1877/4832
Markov base bins and 308/1208 Markov cycle bins unusable.

### N8 · Any truth VCF holding an INS loses one check to "not evaluable"

*Found while reviewing the L15 fix pass (`af9d9fa`). Recorded unfixed at the
time -- an INS-specific check is new work rather than a correction -- and
fixed in the whole-branch pass, because the alternative is shipping a `spike
validate` that cannot pass on spike's own output.*

`check_event` (`validate.rs:133-180`) dispatches on `sv_type` alone:
`coverage_ratio` for DEL/DUP, `split_reads` for DEL/DUP/INV/BND, `allele_freq`
for SNP. INS matches none of them, so it falls to the "no check applies"
branch added by `af9d9fa`, which records an `event_checked` result with
`pass: false`. Nothing about the BAM can change that -- the dispatch never
looks at one -- so **one INS in a truth VCF is one permanent FAIL**, and
`spike validate` on that VCF can never report all-PASS or exit 0, however
good the spike-in is.

Measured on branch HEAD, truth VCF from
`--event del:chr20:38412500-38422500 --event ins:chr20:38430000:500 --seed 1`,
validated against the chr20 slice:

```
DEL chr20:38412500-38422500 (un...  coverage_ratio     0.50   0.91           FAIL
DEL chr20:38412500-38422500 (un...  split_reads        >=2    0              FAIL
INS chr20:38430000 (unknown)        event_checked      a check applies  none for INS  FAIL
[global] insert_size / dup_rate / mean_mapq                                   3x PASS

Result: 3/6 PASS          exit 1
WARN spike::validate] no check applies to INS chr20:38430000 (unknown), so it was not evaluated
```

(The DEL rows fail because this run validates against the *unspiked* slice;
the INS row is the finding, and it is independent of the BAM. Re-measured in
the whole-branch pass against a properly spiked BAM -- `del:chr20:38412500-
38422500` + `ins:chr20:39000000:300 --seed 1`, aligned with the generated
`align.sh` and merged with `merge.sh` -- the two DEL rows PASS and the score is
**5/6 PASS, exit 1**, with the INS row the only failure on a perfect run.)
`af9d9fa`'s branch is right for what it was for -- M11, where an INS-only truth
VCF scored 3/3 PASS on the three global checks alone -- but the cost landed on
every mixed truth set: the honest "we did not check this" was
indistinguishable, in the verdict and in the exit status, from "this spike-in
is wrong". README's step 4 ("verify the spike-in looks correct before running
the caller") was therefore broken for insertions, and this is behaviour the
branch introduced, not a pre-existing gap.

- **Note on the suggested fix.** "Reads carrying the inserted sequence" is not
  available to `validate`: the truth VCF writes INS as a symbolic `<INS>` with
  `SVLEN` only (`truth.rs:264`), so the inserted bases spike generated are not
  in the file `validate` reads. What *is* available is the shape the aligner
  gives them.
- **Fixed:** a new `ins_reads` check counts reads whose alignment **leaves the
  reference at POS** -- an `I` CIGAR operation when the insertion fits inside a
  read that anchors on both sides, a soft clip at the insertion point when it
  does not. A read counts when the operation is at least `min(SVLEN, 50)` bases
  and its reference boundary is within 100 bp of POS; the check passes at 2
  such reads, the same threshold and the same reasoning as `check_split_reads`
  (one clipped read is background anywhere, two at one point are not). An INS
  record with no usable `SVLEN` has no length to look for and is a failed
  check, so the "uncheckable is a failure" principle is not weakened.
- **Measured, same run as above:** the INS row now reads
  `ins_reads  >=2 reads with >=50bp inserted at chr20:39000001  21  PASS`, and
  the run scores **6/6 PASS, exit 0** -- README step 4 works for insertions
  again. Specificity, on the same merged BAM: the same check run at five
  control positions where nothing was planted (chr20:38600000, 38900000,
  39500000, 39750000, 40100000) finds **0, 1, 0, 0 and 0** reads and FAILs all
  five, so the threshold of 2 separates the planted event from the background
  by a factor of 21.
- `spike validate --help` now prints the whole checks-by-event-type table and
  states that a check that cannot run is a failed check, so which types are
  covered is visible without reading the source. README carries the same table.
- `test_event_no_check_applies_to_is_not_evaluable` was kept and re-pointed: it
  used an INS to stand for "an event no check covers", which INS no longer is,
  so it now uses `SVTYPE=CNV`. M11's invariant is unchanged and still tested.
- **Follow-up: the threshold was weak below `SVLEN` ~10, and is now narrowed.**
  `min(SVLEN, 50)` bases of inserted *or clipped* sequence is background once
  SVLEN is small. Measured on the merged HG002 chr20 slice with a truth record
  of `SVTYPE=INS;SVLEN=3` at the same five control positions: **1, 3, 0, 0 and
  1** reads, and 3 is over the two-read pass mark -- chr20:38900000 reported
  `ins_reads 3 PASS` for an insertion that was never planted. A read that
  anchors both sides of a short insertion writes an `I` operation, so a soft
  clip now counts only once the insertion reaches the 50 bp evidence cap, which
  is exactly where a read can no longer anchor both sides. **Measured after:**
  the same five positions give **0, 0, 0, 0 and 0** and all FAIL, while
  genuinely planted 3 bp and 12 bp insertions (`ins:chr20:39300000:3`,
  `ins:chr20:39400000:12`, aligned and merged) still find **12** reads each and
  PASS, and the 300 bp INS above is untouched at **21 PASS, 6/6, exit 0** --
  its evidence is clipped reads, which still count at the cap. The five control
  numbers at the cap are **0, 1, 0, 0, 0**, reproducing the row above exactly.


### N9 · `validate`'s coverage gate is asymmetric in the wrong direction

*Found by the whole-branch review of `review-fixes-2`. The three `pass: true`
returns are older than this branch, but L15 (`af9d9fa`) is the task that
installed `not_evaluable` in the three **global** checks and left the 300 lines
above them untouched, so the branch ends up applying its own stated principle
two opposite ways inside one file: an INS failed loudly (N8) while an indel
passed silently.*

`check_allele_freq` is a single-position pileup of A/C/G/T. Three of its
returns reported `pass: true` for a question it had not asked, and because a
result row *was* pushed, `results.len() != n_before` and L15's "a check
applies" gate counted every one of them as **covered**:

- any REF or ALT longer than one base -- every small indel **and every MNV**
  (`validate.rs:601-609` at `39d9773`);
- a pileup depth below 5 (`validate.rs:630-639`, comment literally
  `pass: true, // not enough data`);
- a non-ACGT alt base (`validate.rs:645-655`).

And `load_truth_events`'s `_ =>` fallback (`validate.rs:397`) routed an
**unrecognised SVTYPE** (`CNV`, `DEL:ME`) into the same SNP arm, where the
symbolic ALT `<CNV>` is longer than one base and took the first exit above.

Measured at branch HEAD `39d9773` against a real spiked+merged chr20 BAM
(`del:chr20:38412500-38422500` + `ins:chr20:39000000:300`, aligned and merged
with the generated scripts), on a hand-written truth VCF of the four cases:

```
SNP chr20:39000499-39000501  allele_freq  0.50  N/A (indel)      PASS   (AC>A)
SNP chr20:39000599-39000601  allele_freq  0.50  N/A (indel)      PASS   (TG>AC)
SNP chr20:39000699-39000700  allele_freq  0.50  unknown alt...   PASS   (T>N)
SNP chr20:39099999-39100000  allele_freq  0.50  N/A (indel)      PASS   (<CNV>)
Result: 7/8 PASS
```

and on a 1%-subsampled BAM, `A>C` at chr20:38410800 where the pileup depth is
1: `allele_freq 0.50 / low depth (1) / **PASS**`, `Result: 4/4 PASS`.

**Separately, a truth record whose `END` is at or before its `POS` PASSed a
coverage ratio taken over a region nothing read.** `count_depth_in_region`
returns a hard-coded `Ok(0.0)` for `start >= end` (`validate.rs:1083-1085`) and
`load_truth_events` never validated END. Measured at `39d9773`, same BAM,
`POS=39200000; END=39199000; SIM_VAF=0.900`:
`DEL chr20:39200000-39199000  coverage_ratio  expected 0.10  observed 0.00
**PASS**`.

- **Fixed:** the three returns now go through `event_not_evaluable`, the
  per-event sibling of L15's `not_evaluable`: the row is still pushed (so the
  event stays covered by the "a check applies" gate) but `pass` is `false` and
  a `WARN` names the reason. `load_truth_events` keeps an unrecognised
  SVTYPE's own type instead of calling it a small variant, so it reaches
  `check_event`'s `event_checked` failure and a `WARN`. A record with
  `END <= POS` is refused as the truth VCF is read, naming the record, and the
  run exits non-zero without grading anything.
- **Measured after,** same BAM and same files: the four-case truth VCF goes
  from **7/8 PASS to 3/8 PASS** with all four rows FAIL and the `<CNV>` row now
  reading `CNV chr20:39100000-39110000  event_checked  none for CNV  FAIL`; the
  depth-1 case goes from `PASS` to `FAIL` (`4/4` to `3/4`); and the `END<POS`
  record now stops the run with `truth record badend at chr20:39200000 has
  END=39199000, at or before its own POS: ...`, **exit 1**, before any check.
- `count_depth_in_region`'s `Ok(0.0)` for a zero-length region is left as it
  is: it is also the answer for a legitimately empty *flank* (an event at the
  very start of a contig), which `check_coverage_ratio` already excludes from
  the flank average by length. With the load-time refusal above, the *event*
  region can no longer be zero-length.

### N10 · No `validate` check measured a small indel's or an MNV's allele fraction

*Found by the whole-branch review, recorded there, and fixed in this pass: it
is the same shape as N8 -- spike's own round trip broken for a variant class
spike itself writes -- and N8's fix left the machinery in place.*

N9 stopped `check_allele_freq` calling an indel or an MNV a pass. What it did
not do is measure one. So a truth VCF holding a record from
`--event "snp:chr20:30000000:ACG:A"` (spike writes it with an explicit
multi-base REF) carried a permanent `allele_freq FAIL`, exactly the position N8
described for INS. Measured before, on a real run of
`snp:chr20:39000000:TGG:T` + `snp:chr20:39100000:T:TCCGG` +
`snp:chr20:39200000:AT:GC --seed 1` on the chr20 37.5-41.5 Mb HG002 slice,
aligned with the generated `align.sh` and merged with `merge.sh`: all three
rows read `N/A (indel or MNV) FAIL`, **3/6 PASS, exit 1** -- and the verdict
was reached before the BAM was opened, so nothing about the data could change
it.

- **Fixed:** `check_allele_freq` now reads the REF/ALT pair's *shape* and picks
  a counting rule for it. A **deletion** (`ACG` > `A`) and an **insertion**
  (`A` > `ACCGG`) are counted off the CIGAR -- an indel is not a column in a
  pileup -- as reads whose alignment carries a `D` or `I` operation of the
  allele's own length at the junction just past the anchor base, against the
  reads that span the same junction without one. The operation may sit within
  10 bp of POS, because an aligner left-aligns an indel to the start of the
  repeat it sits in. An **MNV** (`AT` > `GC`) is counted from the pileup, but
  **jointly**: a read votes only if its bases are the whole alt run or the whole
  ref run, since a fraction per base answers a different question at each offset
  and a read carrying one of the two substitutions is not this variant. All
  three then go through the same depth floor (5) and the same 0.15 tolerance as
  a substitution.
- **Measured after,** same run and same merged BAM: **6/6 PASS, exit 0**, the
  three rows reading 0.40, 0.41 and 0.46 against `SIM_VAF=0.500`. The reads
  behind them are **17, 19 and 13** carrying the variant, against 25, 27 and 15
  spanning without it -- an independent re-implementation of the counting rule
  over `samtools view` reproduces all six numbers exactly. `samtools mpileup`
  sees the same two indels at 16 `-2` and 17 `+4` carriers, a base or two
  fewer because of its own base-quality and BAQ filters, and **0** of either
  at every control position.
- **Specificity, both directions.** The same three truth records against the
  **unspiked** slice BAM read 0.00, 0.00 and 0.00 and all FAIL (**3/6 PASS,
  exit 1**). Against the spiked BAM at five positions where nothing was planted
  (chr20:38600000, 38900000, 39500000, 39750000, 40100000), with the same three
  allele shapes built from each position's own reference bases, all fifteen
  checks read 0.00 and FAIL (**3/18 PASS, exit 1**); the supporting counts there
  are **0, 0, 0, 0 and 0** for each of the three shapes, against 17, 19 and 13
  at the planted sites. Each rule is also silent at the other two events' sites
  (the deletion rule finds 0 at the insertion and MNV positions, and so on).
- **Nothing is left `not_evaluable` that used to be measured, and one case is
  left deliberately.** A **complex** allele -- one that changes length *and*
  rewrites the anchor base, such as `AC` > `GTT`, `A` > `CG` or `AC` > `AGT` --
  has neither a single CIGAR operation nor a single allele run to count, so it
  reports `N/A (complex allele)` and FAILs. spike cannot produce one: its own
  small-variant haplotype is built as `left | ALT | right`, and the truth record
  it writes is whatever REF/ALT the user gave. The depth floor and the non-ACGT
  allele exits are unchanged, and both still FAIL.
- `test_allele_freq_it_cannot_measure_is_not_a_pass` was **corrected, not
  weakened**: three of its four cases (`AC`>`A`, `A`>`ACGT`, `TG`>`AC`) were
  the very shapes this entry gives a check, so they now reach the BAM and error
  on a nonexistent path instead of returning a row. The invariant it pins --
  an unmeasured allele fraction is never a pass -- is unchanged and is now
  asserted over *five* cases (`AC`>`GTT`, `A`>`CG`, `AC`>`AGT`, `T`>`N`,
  `TG`>`AN`). The three that moved are covered by two new tests:
  `test_indel_and_mnv_allele_fractions_are_read_from_the_bam` (the verdict must
  not be reached without opening a BAM) and
  `test_small_indel_and_mnv_allele_fractions_are_measured` (a CRAM fixture with
  real `50M2D50M`, `50M4I46M` and substituted-base CIGARs at half the reads
  each, which must read 0.50 and PASS, and 0.00 and FAIL at a site where
  nothing was planted). N9's `event_not_evaluable` path is untouched.


### N11 · Two `--help` strings contradicted the code, and the refusals had no single home

*Found by the whole-branch review.*

- `main.rs:93` said `Target allele fraction (0.0-1.0)`, but `main.rs:199-206`
  refuses `0` and `0.0` (and NaN, and the infinities). A user who reads the
  help and passes `--allele-fraction 0` is refused by a tool that told them 0
  was in range.
- `main.rs:120-123`'s `--flank` help never mentioned the **2000** minimum
  `main.rs:210-219` enforces -- the likeliest accidental refusal of the lot,
  since 10000 is the default and a user narrowing the window has no reason to
  expect a floor.

**Fixed:** both help strings now state the bound they are checked against, and
`test_help_states_the_bounds_the_code_enforces` renders the long help and
asserts it names `(0.0, 1.0]` and `HAP_FLANK`, next to the validator calls that
enforce them -- so the two cannot drift apart again without a red test. No flag
was removed, renamed or re-defaulted; only the descriptions changed. README's
CLI-reference block was regenerated from the new `--help`.

**Also fixed (the review's Important 10):** the refusals were documented in
nine separate README places and four were not documented at all --
`--allele-fraction` outside `(0,1]` including NaN, `--flank` below 2000, a
gzip-compressed FASTA under a non-`.gz` name (`reference.rs:98-105`), and a
bcftools failure on `--gvcf` (`loh.rs:418`). README now has one **"What spike
refuses"** section listing every command-line-reachable refusal with the exact
message text, each one measured by running it rather than read off the format
string.

- **Correction to the review's list:** a bcftools failure on `--gvcf` is **not**
  a refusal. `loh.rs:418` does bail, but `simulate.rs:94`'s `unwrap_or_else`
  catches it, logs a `WARN` and falls back to `SampleCopies::default()`, so the
  run continues and exits **0**. Measured with an unindexed `.vcf.gz`:
  `WARN spike::simulate] could not read the sample's SNPs in
  chr20:38410500-38424500: bcftools exited with status exit status: 255 on gVCF
  'bad.vcf.gz': Failed to open bad.vcf.gz: not compressed with bgzip. LOH is
  skipped for this region: original reads are suppressed at random.` -- exit 0,
  FASTQ and truth VCF written. The README section records it under **"Not a
  refusal"** with that text, since a user looking for it will expect it there.


### N12 · N5's donor-coverage refusal measured the first breakpoint only

*Found while closing N10, recorded there, **fixed here**.*

`48be0c8` closed N5's second door by refusing an event whose breakpoint has no
donor coverage. It measured **one** position: `simulate.rs:196-220` took
`haplotype.breakpoints().first()`, mapped it back with `hap_to_ref(bp - 1)` --
the last reference base *before* the junction -- and refused only if
`estimate_coverage_at` was 0 or NaN there. A junction has two sides and a
haplotype can have several, and none of the others was ever looked at.

The answer to "is it refused, or accepted, on the wrong evidence?" is **both,
depending on the order the event names its parts**. Measured on the chr20
37.5-41.5 Mb HG002 slice, one BED with `GENEA` at chr20:38.42 Mb (covered) and
`GENEB` at chr20:30.00 Mb (outside the slice), `--seed 1`, `af=0.2`:

```
--event "fusion:GENEB:exon1:GENEA:exon2"   0 pairs from B, 3032 from A
  Error: ... has no donor coverage at its first breakpoint chr20:30000199 ...
  exit 1, --output left empty
--event "fusion:GENEA:exon1:GENEB:exon2"   3022 pairs from A, 0 pairs from B
  INFO Tiling 15 synthetic reads across 4000bp haplotype (cov=58.0, vaf=0.20)
  exit 0, 7 files, truth VCF written
```

Same two loci, same BAM, same seed: naming the covered partner first is enough
to get 15 synthetic reads planted across a junction whose far half has no donor
read behind it -- N5's headline symptom, through a third door. The same holds
without a fusion: `del:chr20:41490000-41600000` on a slice whose reads stop at
41,500,000 has a covered left breakpoint (`cov=57.7`) and **no reads at all** at
its right one, and spike tiles 249 reads across it at exit 0.

**Fixed.** `donor_coverage_for_tiling` now maps back *both* sides of *every*
breakpoint -- `hap_to_ref(bp - 1)` and `hap_to_ref(bp)` for each -- and asks
`estimate_coverage_at` about all of them. Refusing on the first uncovered one
was the obvious next step and is wrong: it turns a DEL whose distal junction
side falls outside the covered window into an error, and that shape is ordinary
input on a sliced or panel BAM.

**The rule this lands on:** *the positions that must carry donor coverage are
the breakpoint sides of each locus the pool was extracted from.*
`extract_pool_for_event` searches **two** windows for a fusion, one per
partner, and **one** for everything else, so:

- **Fusion:** every breakpoint side must be covered. Every fragment spike
  plants for a fusion spans the junction (`breakpoint_only`), so a partner with
  no donor reads makes half of every planted read invention. Testing all sides
  is also what makes the verdict symmetric in the naming order.
- **Every other event:** refused only when **no** side of any of its
  breakpoints is covered -- "no donor coverage anywhere near this event", which
  is what N5 measured. One uncovered side is a thin spot or the far edge of a
  slice, not a refusal. The tiling count is then scaled by the first *covered*
  side; for every event that already ran, that is the first side exactly as
  before, so no working run changes its arithmetic or its bytes (M7).

Measured on the chr20 37.5-41.5 Mb HG002 slice, `--seed 1`, branch HEAD
`ad9881e` -> this commit:

| command | before | after |
| --- | --- | --- |
| `fusion:GENEA:exon1:GENEB:exon2;af=0.2` (covered partner first) | `Tiling 15 synthetic reads ... (cov=58.0)`, 8 files, **exit 0** | **exit 1**, 0 files, `... no donor coverage on one side of its junction ... (chr20:30005000)` |
| `fusion:GENEB:exon1:GENEA:exon2;af=0.2` (same two loci, named the other way) | **exit 1**, 0 files | **exit 1**, 0 files -- the two orders now agree |
| `del:chr20:41490000-41600000` (slice ends 41,500,000; distal side uncovered) | `Tiling 249 synthetic reads ... (cov=57.7)`, **exit 0** | unchanged: `Tiling 249 ... (cov=57.7)`, **exit 0** |
| `del:chr20:37400000-37510000` (slice starts 37,499,851; *near* side uncovered) | **exit 1**, 0 files -- a false refusal of the same shape from the other side | `Tiling 258 synthetic reads ... (cov=60.7)`, **exit 0** |
| `del:chr20:30000000-30010000 --region chr20:38400000-38440000` (N5's route) | **exit 1** | **exit 1**, now naming both sides: `(chr20:29999999, chr20:30010000)` |
| `fusion` with side A at chr20:30.00 Mb, side B at 38.42 Mb (N5's route) | **exit 1** | **exit 1** |

`del:chr20:38412500-38422500 --seed 1` is byte-identical to `ad9881e`'s output
(`R1.fq.gz` `8cf7964f`, `R2.fq.gz` `1bd652f8`, `truth.vcf` `5a73bfd8`,
`replaced_reads.txt` `6bc6392a`) and byte-identical between two runs of the new
binary (M7). `compute_tiling_count`, `floor_tiling_count` and
`estimate_coverage_at` are byte-identical to `ad9881e` -- H5's
`n = round(cov*v/(1-v)*breakpoints.len())` and M1/M2's scaling by `total_len`
with `pool.frag_dist.mean` used directly are untouched; the change is only
*which position* `cov` is read at when the first one has none.

**A correction to what this entry claimed.** It said
`test_simulate_event_keeps_pairs_straddling_footprint_edge` "would newly fail"
under the both-sides check. Measured: it would not. Its pool holds pairs at
[3800,4200), and `estimate_coverage_at` samples a **2 kb** window, so reference
position 3000 measures `cov = 50`, not 0. The entry was reasoning about the
exact base rather than the window the function actually uses. The test passes
unchanged here, and it is the single-locus rule above that keeps the shape it
encodes legitimate.

**Two tests were corrected, not weakened**, both for the same reason: their
synthetic pools held no donor reads at all for one fusion partner, which is the
shape the new rule refuses and which a real extraction never produces.
`test_simulate_event_fusion_is_additive` now holds 100 pairs around each
partner instead of only gene A's; its three assertions (nothing suppressed, all
originals kept, chimeric pairs produced) are unchanged.
`test_cross_chromosome_fusion_ignores_the_other_chromosomes_depth` (M9) now
contrasts **1** chr2 pair against **100** instead of none against 100; its
assertion -- that the chr1 breakpoint tiles the same number of chimeric pairs
either way -- is unchanged, and the contrast still fails loudly if chr2's depth
leaks into chr1's window.

### N13 · `cigar_indel_vote` drops every read whose deletion is spelled off the junction

*Found by the verification review of the fix wave, **fixed here**.*

`ad9881e` gave `allele_freq` a counting rule for small indels: a read carries
the allele if its CIGAR holds an operation of the allele's own kind and length
within `INDEL_POS_PAD = 10` bases of the junction just past the anchor base.
The pad is there because an aligner left-aligns an indel to the start of the
repeat it sits in, so a truth record spelled a few bases the other way still
has to match.

It did not. The vote also required two `M`-covered reference bases -- `pos` and
`pos + ref_len` -- and a `D` operation consumes reference. Shifted left by *s*
it swallows `pos`; shifted right by *s* it swallows `pos + ref_len`. So for
every *s* in `[1, indel_len]` **neither** anchor is covered by an `M` block,
`cigar_indel_vote` returned `None`, and the read entered neither the numerator
nor the denominator. Insertions consume no reference and were never affected.

Measured at unit level (`ACG` > `A` at 0-based 1000, read 950..1050), the vote
was **non-monotonic in the shift**, which is what makes it a bug rather than an
exact-placement rule:

```
shift  -6 -4 -3 | -2 -1 | 0 | +1 +2 | +3 +4 +6
before  C  C  C |  N  N | C |  N  N |  C  C  C      (C = Carries, N = dropped)
after   C  C  C |  C  C | C |  C  C |  C  C  C
```

End to end on the chr20 37.5-41.5 Mb HG002 slice, the **same physical 2 bp
deletion** in an `(AC)n` repeat at chr20:38549586, three legal VCF spellings
differing only by repeat unit, `SIM_VAF=0.10`:

| spelling | before (`99f1a8e`) | after |
| --- | --- | --- |
| `38549585 TAC>T` (left-aligned, where bwa put the `D`) | carries 14, spans 23, **0.38 FAIL** | unchanged: **0.38 FAIL** |
| `38549587 CAC>C` (one `AC` unit right, shift -2) | carries **1**, spans 26, **0.04 PASS** | carries 14, spans 26, **0.35 FAIL** |
| `38549589 CAC>C` (two `AC` units right, shift -4) | carries 14, spans 27, **0.34 FAIL** | unchanged: **0.34 FAIL** |

13 of the 14 carrying reads vanished from *both* numerator and denominator at
shift -2 and came back at shift -4, and the verdict was inverted by a spelling
change that does not change the variant. The N10 entry's fail-safe rationale --
that such reads "would be counted as reference support ... a FAIL, not a false
PASS" -- was **wrong in mechanism**: the denominator shrank with the numerator.

The trigger is a **non-normalized REF/ALT** in the truth record, which spike
accepts without complaint and writes verbatim. A `bcftools norm`-ed record is
safe, because bwa left-aligns too.

**Fixed.** An operation of the allele's own kind and length inside the pad *is*
the junction, so it is now enough on its own: `cigar_indel_vote` returns
`Carries` before it looks at the anchors. The anchors still gate the other
verdict -- they are what separates "spans this junction without the indel" from
"never reached it", and a read carrying a *different* indel over one of them
still votes neither way. Beyond the pad the operation is somebody else's indel
and the read votes `Spans`, as before.

`cigar_indel_vote` had **no direct unit test at all**; its only coverage was a
CRAM fixture whose `D`/`I` ops sit exactly at `POS+1`, the one shift the bug
does not touch -- the weak-test shape [Test gaps](#test-gaps) names. Four
direct tests now pin the whole shift range: every shift in +/-10 votes
`Carries` for a deletion and for an insertion, +/-11 votes `Spans`, and a
clipped read or a different-length `D` over the anchor still votes neither way.

## Low severity

| ID | Problem | Where | Fix |
| --- | --- | --- | --- |
| L1 | **Fixed.** Indel error model pads deletions with `N` (Q2): 0.69% of R1 end in N, 1.01% of R2 start with N at rate 0.05 | `synth.rs:697-700, 729-730` | Pass `read_length + 10` bases; generate R2 in sequencing order |
| L2 | **Fixed.** CRAM containers with several contigs leak other contigs' reads (synthetic 2-contig CRAM: 5 chrB pairs labelled chrA) | `extract.rs:245-279, 319-362` | Skip records whose ref id or mate ref id differs (the same hole at `loh.rs`/`validate.rs` is tracked as N4) |
| L3 | **Fixed.** bgzipped FASTA read as raw bytes; fails later with misleading "beyond chromosome length" (reproduced: real 791 MB bgzipped GRCh38 loaded as a 0-byte chromosome, then `DEL event start on chr20 is at or beyond chromosome length (38412500 >= 0)`) | `reference.rs:33-35` | Used `fasta::io::indexed_reader::Builder`, which picks a bgzf- or plain-file reader by extension; a missing `.gzi` now fails at open with a message naming the `.gzi` index instead of surfacing downstream as chromosome length. **Fixed** that way: same command against the same file now loads chr20 at 61 MB; the truth VCF matches the run against the uncompressed FASTA apart from the `##reference=` path line (`truth.rs:75`), and the FASTQ pair's decompressed content (`zcat \| md5sum`, the correct comparison for a `.gz` pair) matches exactly. See Fix pass 1 below: the detection was extension-only and so still missed a bgzip file under a non-`.gz`/`.bgz` name. |
| L4 | **Fixed.** Final FASTQ flush error ignored (write to `/dev/full` returned `Ok`). Reproduces when the whole gzip output is small enough to sit unflushed inside the `BufWriter`'s own buffer after `finish()` — under 8 KiB by default (measured: a single-pair `write_paired_fastq` call with `R1.fq.gz` symlinked to `/dev/full` returned `Ok(...)`). A large event-scale run (3862-pair `del:chr20:38412500-38422500`, HG002 chr20 slice, `--seed 1`) already failed correctly before this fix, because the loop's own writes overflow that buffer first and hit `/dev/full` mid-stream — but plenty of real spike runs are small enough to stay under that buffer: a single small SNP or short DEL with low local coverage, a tight `--region`, or a demo-scale run. This fix's window tracks the buffer size, not run size in general, so it still covers those. See Fix pass 1 below | `fastq.rs:97-124` | `.finish()?.flush()` on both streams, computed unconditionally so R1 failing first can't leave `r2_gz` to be dropped unfinished and discard its own error the same way — now returns `Err("No space left on device (os error 28)")` for the single-pair case above, and names both streams when both fail |
| L5 | **Fixed.** Panic when mean read length > 1500 (`clamp` with min > max) | `stats.rs:103` (now ~110 after the task 10 rewrite); callers `tile_haplotype_reads` (`simulate.rs:557` at this commit, `413` when the row was written) and `generate_depth_pair` (`synth.rs:616` at this commit, `534` then) | Guard min ≤ max |
| L6 | **Fixed.** BND POS is one base past the kept base, vs spec; parser mirrors it so round trips agree | `truth.rs:121`, `vcf_input.rs:103` | Write POS = last kept base |
| L7 | **Fixed.** DEL/DUP/INV with no END and no SVLEN silently becomes a 1 bp event | `vcf_input.rs:130-132, 148-150, 165-167` (shared via `resolve_sv_end_or_warn`) | Derive the span from the alleles only where they give one unambiguously: a multi-base `REF` whose `ALT` is symbolic or the anchor base alone (`REF=ACGT ALT=A`), or — for DUP — a single-base `REF` anchor whose `ALT` is anchor + duplicated copy (the sequence-resolved INS form). Anything else is rejected with a `log::warn!` (stderr) naming type and `chrom:pos`, instead of guessing 1 bp. Measured on real data (HG002 chr20 slice, `--seed 1`), with the ba50bf8 binary vs. the fix: (a) `POS=38412500 REF=GTTAAAGTTTATCAGAAAATT ALT=GTTAAAG SVTYPE=DEL` (a 14 bp deletion at 38412507-38412520) wrote truth `END=38412520;SVLEN=-20` at `POS=38412500` — 20 bp starting 6 bases too early, exit 0 — and is now skipped with a warning (exit 1, no truth record); (b) `REF`=that 21 bp span, `ALT`=its reverse complement, `SVTYPE=INV` wrote `END=38412520;SVLEN=20` (20 bp at 38412501-38412520, one base short and one right of the true 21 bp span) and is now skipped the same way; (c) the well-formed `POS=38412499 REF=T ALT=T+21 bp SVTYPE=DUP` went from skipped ("single-base REF (no length information)" — false: the length is in ALT) to `POS=38412499;END=38412520;SVLEN=21`, the 21 duplicated bases 38412500-38412520 matching ALT exactly and round-tripping unchanged when the truth VCF is fed back in. The REF-only path for `REF=21 bp ALT=G` is unchanged (`END=38412520;SVLEN=-20`). Stripping a shared REF/ALT prefix, which would let (a) be decoded rather than skipped, was left to L10 and is now done there: (a) and the equal-length INV decode, and the only L7 shape L10 narrowed is a DUP whose single-base REF is not ALT's first base (`REF=T ALT=GGGT`), which no longer takes `ALT[1..]` as the copy but is rejected with the same warning |
| L8 | **Fixed.** `af=nan` / `--allele-fraction NaN` was accepted — `v <= 0.0 \|\| v > 1.0` is false on both sides for NaN — reaching simulate as `vaf=NaN` and landing in the truth VCF as `SIM_VAF=NaN`. Both entry points now use the negated form, which also correctly keeps rejecting `inf`/`-inf` (already did) and `0.0`, and keeps accepting `1.0`. The VCF-input path (`vcf_input.rs::extract_af`) already used the negated form (`v > 0.0 && v <= 1.0`) and needed no change — feeding it `SIM_VAF=nan` falls through to the next INFO key, then to the CLI default, never to NaN; verified below. **Fix pass 2, correcting fix pass 1's note:** "only reads inside the event window were suppressed... not the whole read pool" was right that the pool is bounded but wrong about the mechanism inside it. `copy_rate()` (`synth.rs:794-799`) is `Some(true) -> (2·vaf).min(1.0)`, `Some(false) -> (2·vaf-1.0).max(0.0)`, `None -> vaf`; Rust's NaN-propagating `min`/`max` resolve those first two to a fixed `1.0` and `0.0` for any NaN vaf, and `x < NaN` is always `false`. So within the replaceable pool (pairs entirely inside the haplotype footprint; pairs outside it are untouched at *any* af, valid or not), NaN suppresses event-copy pairs at a hardcoded 100%, and *never* suppresses other-copy or unphased pairs — instead of the ~vaf rate a valid af gives them. Measured, pre-fix binary (built from `1a10346`), `del:chr20:38412500-38422500`, seed 1, chr20 slice BAM: of the 1,971 replaceable pairs (276 event-copy + 265 other-copy + 1,430 unphased) inside the 4,559-pair pool, af=NaN suppressed 276/276 event-copy (100%), 0/265 other-copy (0%), 0/1,430 unphased (0%) = 276 total — matching the earlier note's count, but for the wrong reason (a deterministic per-copy split, not a proportional window rate). A real af=0.7 on the identical pool/seed suppressed 276/276 event-copy (100% — this branch isn't actually wrong at af ≥ 0.5), 113/265 other-copy (42.6%, vs. NaN's 0%), 1,003/1,430 unphased (70.1%, vs. NaN's 0%) = 1,392 total. (af=0.5 alone can't show the other-copy gap: `max(0, 2·0.5−1) = 0` too, a boundary coincidence that also reproduces 276 total for an unrelated reason — af=0.7 was added to break it.) So the unphased majority — 1,430 of 1,971 replaceable pairs, 72.6% — was the largest group silently getting zero suppression instead of a valid af's rate, not a smaller but still-proportional slice of the window. | `exon.rs:281` (`parse_af_value`), `main.rs:297` (extracted into `validate_allele_fraction`, matching the `validate_flank`/`validate_read_length` pattern already in `main.rs`) | `if !(v > 0.0 && v <= 1.0)` |
| L9 | **Fixed.** REF == ALT was checked before uppercasing in the `snp:` CLI spec, so `A:a` (same base, different case) passed as if it were a real variant; the VCF ingest path had no identical-alleles check at all. Either way the "variant" was a silent no-op written straight to the truth VCF. `exon.rs::parse_snp_spec` now uppercases both alleles before comparing, so the existing "REF and ALT alleles are identical" error also catches case-only differences. `vcf_input.rs`'s `SmallVar` arm gained the missing check (`record.ref_allele.eq_ignore_ascii_case(&record.alt)`), matching the exon.rs comparison. Measured on real data (HG002 chr20 slice, `--seed 1`): base binary (`8d1beba`) fed a VCF with `chr20:40000100 REF=A ALT=A` and `chr20:40000200 REF=A ALT=a` wrote both verbatim to `truth.vcf` (`A\tA`, `A\ta`, `SIM_VAF=0.500`, exit 0); the fixed binary rejects the first record immediately with `Error: record same_case has identical REF and ALT alleles: 'A'` (exit 1), and a record differing only by case (`A`/`a`) alone is rejected the same way. A legitimate differing-base record with a lowercase REF (`a`/`c`) still succeeds, confirming the check does not reject genuine soft-masked input. **Follow-up (carried over from the L9 review, now closed):** the accepted record's alleles were written to the truth VCF in the case they arrived in, so `a`/`c` via `--vcf` wrote `a\tc` while the same variant via `--event "snp:chr20:40000100:a:c"` wrote `A\tC` — `types.rs:72-73` documents the alleles as uppercase and `exon.rs::parse_snp_spec` upholds it, the VCF arm did not. The VCF arm now uppercases too. Measured (HG002 chr20 slice, `--seed 1`, `chr20:40000100 REF=a ALT=c`): base binary wrote `a\tc`, the fixed binary writes `A\tC`, the same text the `--event` path writes; the effect is confined to the truth VCF, since the base binary's R1/R2 FASTQs for `a`/`c` and for `A`/`C` are byte-identical (`94ed47c5…`, `48bf4046…`) — `haplotype.rs` uppercases before synthesis and `main.rs::validate_ref_allele` compares case-insensitively | `exon.rs:572-622` (`parse_snp_spec`), `vcf_input.rs:220-232` (`SvTypeTag::SmallVar` arm of `records_to_events`) | Uppercase first; add check to VCF path |
| L10 | **Fixed.** Sequence INS with multi-base REF duplicated `REF[1..]` (`REF=AT ALT=ATGGG` → 4 bp inserted, not 3), and `resolve_sv_end` lacked the same stripping, so a DEL/DUP/INV whose REF and ALT both carried sequence was rejected rather than decoded. Both arms now strip the flanks the alleles share — the common prefix, then any common suffix the leftovers still share (`trim_shared_flanks`) — and start the event after the last base both alleles keep (`event_start`); with no shared prefix at all, an equal-length substitution carries no padding base, so POS is the first affected base. What is left of the alleles must then spell the event: REF keeping bases ALT drops (the span is those bases, for all three types); for DUP, a bare single-base REF anchor whose ALT is anchor + copy (unchanged from L7); for INV, an equal-length pair whose ALT really is REF reverse-complemented. Anything else keeps L7's loud rejection, including a DUP whose ALT is *longer* and carries sequence in both alleles — its copy could equally be REF's own span, which the stripped form would silently relocate. An INS whose REF keeps bases its ALT drops is a complex record, not an insertion, and is skipped with a warning naming `chrom:pos` (`not_an_insertion_warning`), the same shape as L7's. Measured on real data (HG002 chr20 slice, `--seed 1`): (a) `POS=38412500 REF=GTTAAAGTTTATCAGAAAATT ALT=GTTAAAG SVTYPE=DEL` — base binary (`8d1beba`) wrote `POS=38412500;END=38412501;SVLEN=-1` (a 1 bp deletion, the L7 symptom), `b115fed` skipped it (exit 1), and the fix writes `POS=38412506;END=38412520;SVLEN=-14`, the real 14 bp deletion 38412507-38412520; 13 of the synthetic reads carry the exact 40 bp junction `chr20:38412487-38412506 + 38412521-38412540`, which the base binary's output has 0 of; (b) the same 21 bp span with `ALT` its reverse complement and `SVTYPE=INV` — base wrote `END=38412501;SVLEN=1`, `b115fed` skipped it, the fix writes `POS=38412499;END=38412520;SVLEN=21`, the full 21 bp span 38412500-38412520; (c) `POS=40000101 REF=AT ALT=ATGGG SVTYPE=INS` — base and `b115fed` both wrote `POS=40000101;SVLEN=4` (the 4 bp `TGGG`), the fix writes `POS=40000102;SVLEN=3` (the 3 bp `GGG` after the T). All three truth records round-trip through `--vcf` unchanged. **Follow-up (carried over from the review of this commit, now closed):** three shapes still decoded or rejected wrongly. (i) A single-base ALT skipped the stripping entirely, so `REF=GTTA ALT=A SVTYPE=DEL` (the alleles share their *suffix*) wrote `POS=38412500;END=38412503;SVLEN=-3` — deleting `TTA` and leaving `G`, where the record says the result is `A`, exit 0, no warning; only a symbolic `<DEL>` now takes that path, and the same record writes `POS=38412499;END=38412502;SVLEN=-3`, deleting `GTT` and leaving `A`. (ii) Stripping can move the start left of POS, and at `POS=1` that is base 0: `chr20:1 REF=AT ALT=GGGAT SVTYPE=INS` wrote `chr20 0 N <INS>` (exit 0) and `chr20:1 REF=ACGT ALT=T SVTYPE=DEL` wrote `POS=1;END=4`; both are now skipped with `before_first_base_warning` naming `chrom:pos` (exit 1, no truth record). (iii) An INV is no longer trimmed at all (`inv_span`): an inversion's span is stated by the record, and trimming narrowed it whenever the inverted region's ends were their own complements — `chr20:38412503 REF=AAAGTTTAT ALT=ATAAACTTT` (its reverse complement) wrote a 7 bp `POS=38412503;END=38412510` for the 9 bp it spells out, and the padded spelling of the same event (`chr20:38412502 REF=TAAAGTTTAT ALT=TATAAACTTT`) wrote the same 7 bp; both now write `POS=38412502;END=38412511;SVLEN=9`, the 9 inverted bases 38412503-38412511, and round-trip through `--vcf` unchanged. A DUP whose ALT is a strict prefix of REF (`REF=TGTT ALT=TG`) does decode, as the bases REF keeps — the row's "a DUP with sequence in both alleles is rejected" meant the ALT-longer shape, and README and the sentence above now say so | `vcf_input.rs:190-231` (INS arm), `vcf_input.rs:432-576` (`trim_shared_flanks`, `event_start`, `is_reverse_complement`, `resolve_sv_span`, `span_from_alleles`, `not_an_insertion_warning`) | Strip common prefix. The same stripping is what `resolve_sv_end` (L7) lacks: DEL/DUP/INV records whose REF and ALT both carry sequence are currently rejected rather than decoded, so this fix should cover both arms |
| L11 | **Fixed.** VCF records silently dropped (multi-allelic, unknown SVTYPE, `DUP:TANDEM`, short lines); INFO `AF` (often population frequency) used as VAF. Two defects, both fixed. **(a) Silent drops.** Every reason the ingest passes a record over is now counted per reason and reported on stderr (`VcfIngestStats`), including the four L7-L10 warn paths that were loud but uncounted (`no_length_warning`, `not_an_insertion_warning`, `before_first_base_warning` from both the span and the INS arm). The parser's own silent `continue`s each gained a reason: short line, `POS` not a positive integer, multi-allelic `ALT` (its own reason, since taking the first of `A,T` would silently simulate half the record), an `SVTYPE` spike has no model for, and a no-`SVTYPE` record whose alleles are not plain DNA. `DUP:TANDEM` was *not* something to merely report: VCF v4.3 spells subtypes with a colon, so `sv_type_tag` now takes the base type before the first one and `DUP:TANDEM` is simulated as a duplication, `DEL:ME:ALU` as a deletion (byte-identical output to the plain type, verified). **(b) INFO `AF`.** Plain `AF` is no longer read as a VAF by default -- in a population VCF it is the population allele frequency, not this sample's read fraction. Only `SIM_VAF` and `VAF` are read; a record with only `AF` falls back to `--allele-fraction` and is counted. The new `--vcf-info-af` flag restores the old behaviour for a VCF that really does state a VAF there (no existing flag removed, renamed or re-defaulted; documented in README). Folded in as the same quietly-wrong-output shape: an unusable VAF value (`SIM_VAF=nan`, an unparseable number) used to fall through silently to `.unwrap_or(config.allele_fraction)`; it now warns naming the record and `chrom:pos`, and is counted. Also one line, since the counting made it one: spike acts on neither `FILTER` nor `GT`, and now says how many records it read were non-PASS or hom-ref instead of leaving it to be assumed. **Measured** (HG002 chr20 slice, `--seed 1`, a 6-record VCF), `0939ea94` vs. the fix: the base binary loaded 2 events and dropped 4 records with **no output at all** (exit 0); the fix loads 3 and reports `skipped 3 VCF record(s): 1 short line (fewer than 8 columns); 1 multi-allelic ALT; 1 SVTYPE spike does not simulate`. The extra event is the `DUP:TANDEM` record, which now writes truth `POS=38423496;END=38427196;SVLEN=3700` and takes the R1 read count from 2687 to 6595; its truth record round-trips through `--vcf` unchanged, and it is byte-identical (truth VCF and R1 md5) to the same record written `SVTYPE=DUP`. The `AF=0.001` record went from truth `SIM_VAF=0.001` (the population frequency simulated as a VAF) to `SIM_VAF=0.500` plus a warning naming `--vcf-info-af`; with `--vcf-info-af` it is `SIM_VAF=0.001` again. The `SIM_VAF=nan` record writes `SIM_VAF=0.500` both ways, but silently before and with `record badvaf at chr20:40000399 has INFO SIM_VAF=nan, which is not a fraction in (0, 1]` after. On the real phased HG002 chr20 truth VCF (208,757 records) 83,058 (39.8%) are non-PASS and 0 are multi-allelic or hom-ref; a 3-record non-PASS slice of it against the full HG002 BAM produces an identical truth VCF before and after, with the fix additionally reporting `3 record(s) it read are not PASS`. A clean VCF's FASTQ output is byte-identical to `0939ea94`'s. **Fix pass 2** (review of `d8af0ff` itself), four defects: (a) splitting `SVTYPE` on the first colon accepted *any* subtype, including `DUP:DISPERSED`/`DUP:INT` — an interspersed duplication, not the local/tandem copy spike's model produces — so it was simulated as `SVTYPE=DUP` at the wrong span with no warning, a silent drop turned into silently *wrong* truth, exactly what this row exists to close; `sv_type_tag` now allowlists DEL/INS with any suffix (the subtype names what was deleted/inserted, not the event's shape) and DUP only for the bare type or `:TANDEM`, dropping and counting anything else as an unsimulated SVTYPE. (b) `af_unusable` counted unusable *keys*, not records: `SIM_VAF=nan;VAF=0.3` incremented it even though the record found a usable VAF from `VAF` and never touched `--allele-fraction`, and `SIM_VAF=nan;VAF=nan` counted 2 for one record; it's now incremented once per record, only when every key was exhausted. (c) `log_summary` ran only after `ingest_vcf` returned `Ok`, so a `bail!` (L9's identical-allele check, a bad BND) discarded every count from records read before it; `ingest_vcf` now always returns its `VcfIngestStats` alongside the `Result`, so the counts survive and are logged before the error propagates. (d) the skipped-records summary only fired when `dropped > 0`, so "0 skipped" and "the counting silently broke" both logged nothing; a new unconditional `totals_line` ("read R VCF record(s), simulated E, skipped N") logs after every ingest that reaches it. Also: a blank line in the VCF body was counted as a short line (fewer than 8 columns) — it is not a record at all and is no longer counted. **Measured**, HG002 chr20 slice, `--seed 1`: a `SVTYPE=DUP:DISPERSED;END=38427196` record against `d8af0ff` loaded 1 event and wrote truth `SVTYPE=DUP;END=38427196;SVLEN=3700` with no warning at all (exit 0); against the fix it is dropped, `skipped 1 VCF record(s): 1 SVTYPE spike does not simulate`, exit 2 ("no events specified"). A 3-record file (a CNV, `SIM_VAF=nan;VAF=0.3`, a blank line, then a fatal `REF=A ALT=a`) against `d8af0ff` logged only the `SIM_VAF=nan` warning before aborting — the CNV drop vanished with the error; against the fix it additionally logs `skipped 1 VCF record(s): 1 SVTYPE spike does not simulate` before the same abort, and never logs the false "state a VAF that is not a fraction" line for the `af1` record, which found its VAF from `VAF=0.3` | `vcf_input.rs` (`VcfIngestStats`, `ingest_vcf`, `sv_type_tag`, `is_hom_ref_gt`, `extract_af`, `unusable_af_warning`), `main.rs` (`--vcf-info-af`) | Warn with counts; don't use plain `AF` by default |
| L12 | **Fixed.** `align.sh` / `merge.sh` break on relative paths and on `}`, `"`, `$`, backticks in paths. The baked-in defaults were interpolated raw into `REF="${1:-PATH}"` / `ORIGINAL="${1:-PATH}"` / `SAMTOOLS="PATH"`, where a `}` ends the expansion early, a `"` ends the string, and `$` and backticks expand. They are now made absolute (`script_path`: canonicalize, falling back to a lexically absolute path when the file is not there to canonicalize) and single-quoted (`sh_quote`, splicing an embedded `'` back in as `'\''`), emitted as `REF=${1:-'PATH'}` — safe unquoted because an assignment right-hand side is never word-split or globbed. `--samtools` is resolved only when it holds a `/`, so the default bare `samtools` is still looked up on `$PATH`. A path containing a single quote is **handled, not rejected**. Carried over from the M5 review: the `SM` sanitiser rewrote every character outside `[A-Za-z0-9._+@:-]`, so `SM:Patient 123` became `Patient_123`, which no longer matches the original read groups — re-introducing the two-sample merged BAM that M5 exists to prevent. The aligner read-group arguments are now quoted too (bwa-mem2's and minimap2's `-R` already was; bowtie2's `--rg SM:…` was not), so `pick_sample_name` rewrites control characters and a backslash — a tab ends the `SM` field, a newline ends the `@RG` line, and a backslash is unescaped by bwa-mem2/minimap2 inside the `-R` string itself (`SM:LAB\tech01` → a truncated `SM` plus a bogus extra field), none of which any quoting can carry through; a correctly-quoted custom `--aligner` also reached the `echo` banner unquoted and broke the script's own syntax there (fixed: `sh_quote`d for the echo only, the invocation position stays verbatim by design). **Measured** on real data (HG002 chr20 slice hard-linked to a file whose name is `we ird}"$HOME` + a backticked `touch pwned` + `'x.bam`, `del:chr20:38412500-38422500`, `--seed 1`): against `349bd90` both generated scripts fail to parse — `bash -n` and a plain `bash merge.sh` on its own defaults both exit 2 with `syntax error near unexpected token '('`, and no `merged.bam` is produced; with the fix `merge.sh` exits 0 and writes `merged.bam` holding 1,328,169 records out of the 1,337,290-record original (9,121 removed for the 4,559 replaced read names), with `$HOME` unexpanded and the backtick never run. With a relative `-b slice.bam -r ref.fasta`, `merge.sh` run from another directory fails against `349bd90` (`Failed to open file "slice.bam" : No such file or directory`, exit 1, no `merged.bam`) and succeeds with the fix. For a well-behaved sample and path the generated scripts differ from `349bd90`'s only in those three default lines. `scripts/validate_pipeline.sh` is unaffected: it always passes `ORIGINAL`/`REF`/`THREADS` positionally, so the defaults are never used — and positional arguments holding a space and a `$` were measured to still work, now pinned by an executing test for each script. Real minimap2 (2.30) with the exact `-R` string `write_align_script` emits: `SM:LAB\tech01` produced `@RG ID:sim SM:LAB <TAB> ech01 PL:ILLUMINA` before this fix pass (SM truncated, plus a bogus field) and `@RG ID:sim SM:LAB_tech01 PL:ILLUMINA` after | `main.rs` (`sh_quote`, `script_path`, `script_command`, `write_align_script`, `write_merge_script`), `bam_stats.rs` (`pick_sample_name`) | Canonicalize and single-quote paths |
| L13 | **Fixed.** Contig names containing `:` could not be used. `--event` split the spec on every `:` and `--region` split at the first one, so `HLA-A*01:01:01:01` was torn into a contig `HLA-A*01` plus coordinate fields -- `del:HLA-A*01:01:01:01:1000-2000` failed as `gene 'HLA-A*01' not found`, `--region HLA-A*01:01:01:01:1000-2000` as `invalid start in region`. Both now resolve the leading name against the reference's `.fai` contigs first (`reference::split_contig`), longest match winning so `HLA-DRB1*01:01:01:02` beats the `HLA-DRB1*01:01:01` it contains, and fall back to the old first-`:` split when nothing matches -- a gene name, or a contig the reference does not list, so a typo still fails as before. `main.rs` reads the contig list once, before the specs are parsed, and reuses it for the truth VCF's `##contig` lines instead of re-reading the `.fai`. **REVIEW.md is wrong twice in this row.** (a) The third call site, `vcf_input.rs`'s BND ALT decode, already used `rsplit_once` at `8d1beba` and needs no contig list at all: a BND ALT is always `chrom:pos` with a numeric last field, so the last `:` ends the name however many it holds. It is left untouched (H3 owns it) and is now pinned by `test_decode_bnd_partner_contig_name_containing_colons`, which passed on first run. (b) "525 HLA contigs in this BAM" does not hold for this run's data: the HG002 BAM and `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` each carry **195** contigs and **none** contains `:` -- the no_alt analysis set has no HLA contigs -- so the symptom was reproduced on a real HG002 chr20 slice reheadered to the contig name `HLA-A*01:01:01:01` against a matching renamed FASTA. **Measured** (real HG002 reads, `--seed 1`, `--event "del:HLA-A*01:01:01:01:38412500-38422500" --region "HLA-A*01:01:01:01:38400000-38430000"`): against `5237abb` the run dies with `Error: gene 'HLA-A*01' not found. Available:` (exit 1, no output directory contents); with the fix it exits 0, logs `Extraction region: HLA-A*01:01:01:01:38399999-38430000`, writes 4242 read pairs and the truth record `HLA-A*01:01:01:01 38412500 sim_del_1 G <DEL> ... SVTYPE=DEL;END=38422500;SVLEN=-10000` under `##contig=<ID=HLA-A*01:01:01:01,length=64444167>`. On an ordinary contig the change is inert: the same event/region on plain `chr20` against the real reference is byte-identical before and after (`truth.vcf` 5a73bfd8, `R1.fq.gz` cc876711, `R2.fq.gz` 41d5fa9a, `replaced_reads.txt` bdbd3b3b) | `reference.rs` (`split_contig`), `exon.rs` (`split_event_parts`, `parse_event_spec`), `main.rs` (`parse_region`, contig list read before spec parsing); `vcf_input.rs` unchanged | `rsplit_once`; match known contigs |
| L14 | **Fixed.** Standard BED6: column 5 is score, so every gene is named `"0"`. Column 5 is now read as the gene symbol only when it is not a BED score (`.` or an integer 0-1000) -- neither shape is ever a gene symbol, so nothing is guessed and no format flag was added; column 4 keeps supplying the gene when column 5 does not (`LDLR_exon1` -> `LDLR`), which is what standard BED files fall back to. **Measured** on real data (full HG002 BAM, `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`, `--seed 1`, `del:LDLR:exon4-exon6`, `data/ldlr_deletions/ldlr_exons_hg38.bed` rewritten as BED6 with score `0` and strand `+`): the branch tip (`bcf6b1f`) dies with `Error: gene 'LDLR' not found. Available: 0` (exit 1, no output written), and the fix exits 0 and writes `chr19 11105219 sim_del_1 G <DEL> ... SVTYPE=DEL;END=11107514;SVLEN=-2295;SIM_GENE=LDLR;SIM_EXONS=LDLR_exon4,LDLR_exon5,LDLR_exon6` -- byte-identical (`truth.vcf` ab39f3c6, `R1.fq.gz` 0035ebf3, `R2.fq.gz` 1a2c8ba5) to the same run against the original 5-column BED. On that 5-column BED the change is inert: all five output files are byte-identical to `bcf6b1f`'s. The symptom was loud rather than silent -- a single-gene BED6 fails as `gene not found`, and a multi-gene one usually fails earlier still, as `gene 0: exon number 1 appears more than once`, because every gene collapses into one named `0` **Fix pass (task 27):** the range was pinned only by column-5 values `0`, `0`, `.` and `LDLR`, every one of which the mutant `fn is_bed_score(v) -> bool { v == "." \|\| v == "0" }` also satisfies; `\t1000\t` -> gene from the name and `\t1001\t` -> gene `1001` now pin both edges, and each kills a different mutant (`== "0"` kills the first, `score <= 10000` kills the second -- both run). README's "a BED score is never anything else" was false in practice -- bedtools and MACS routinely emit scores above 1000 and UCSC tolerates them, so spike reads `5000` as a gene symbol -- and now says so, states the remedy (put the symbol before the first `_` of column 4; measured: `100_exon1` + score `0` -> gene `100`, but `LDLR_exon1` + `5000` -> gene `5000`, because a non-score column 5 wins), and closes two doc-vs-code gaps in the table: an empty column 5 also falls back to the name (`!g.is_empty()`), and fewer than 4 columns is a hard `bail!`, not a BED3 fallback | `exon.rs:97-103` (gene column), `exon.rs:137-142` (`is_bed_score`) -- the original row cited `93-103` for both, which covers the gene column but not `is_bed_score` | Use column 4 or detect format |
| L15 | **Fixed.** All three global checks read the head of the file — the first 100k records for `mean_mapq`, 200k for `dup_rate`, the first 50k properly-paired for `insert_size` — and the head of a whole-genome BAM is chr1's telomere. They now read one sample taken from the truth events' own windows (event +/- `--flank`) through the same indexed query the per-event checks use, capped at 200,000 records with a per-event share so one long event cannot spend the budget. Event regions rather than a spread across the file: they are what `validate` is judging, and they are already indexed queries, so the sample costs no more than the coverage check on the same event. A file whose duplicates were never marked also stops getting a fabricated `0.0%` — zero duplicate flags is no evidence about a duplicate rate — so `dup_rate` reports `no dup flags` and FAILs (M10's "not evaluable should fail or skip, not pass"); the same for `insert_size` when no properly-paired record is in the windows (it used to inherit `bam_stats`' 350+/-50 fallback, which passes) and for a sample with no records at all. **Measured** on the whole-genome HG002 BAM (54 GB, one truth DEL at `chr20:38412500-38422500`, `--flank 5000`), base binary `8d1beba` vs. the fix: `mean_mapq` **10.0 FAIL → 59.7 PASS** (`samtools` over the same window: 59.74 from 6861 primary records), `dup_rate` 7.1% → 11.7% (`samtools`: 11.7%), `insert_size` 409+/-179 → 416+/-178, whole run **0.52 s → 0.23 s**. On the 10 Mb chr20 CRAM of M15/N3 the same run goes **3.43 s → 0.87 s**, so M15's speedup is kept, not handed back. On that BAM slice with `samtools view --remove-flags 1024`: `dup_rate` **0.0% PASS → `no dup flags` FAIL**. REVIEW.md's "dup rate always 0" is right only for an unmarked file; the HG002 BAM is marked, and there the first-100k defect showed up as 7.1% (chr1's head) against 11.7% in the event's own window **Fix pass 1.** The per-event share had a 1,000-record floor, so a truth VCF with more than 200 events spent the whole 200,000 on the first ~200 in VCF order and events 201+ were never read, with nothing in the log or the output to say so. Measured on `chr20_slice.bam` with 250 DEL events: the pre-fix binary gives the same three numbers for the 250-event VCF as for a VCF holding only its first 200 (`418+/-178` / `12.4%` / `59.6`), and reversing the event order moves them (`418+/-178` / `12.3%` / `59.5`) — the sample depended on where the list started. The floor is gone; every event gets `200,000 / n` records (at least 1), and forward and reversed now agree exactly (`419+/-178` / `12.3%` / `59.6`, log: `[global] sample: 200000 records over 250 of 250 event regions (up to 800 each)`). The sample's size and each region's contribution are logged at INFO, because a check over 117 records prints exactly like one over 200,000. `dup_rate` no longer reads `dups == 0` as proof of an unmarked file: with a duplicate-marking `@PG` (matched on command *tokens*, so a path named `...dedup.bam` is not a marker) zero flags is a measured `0.0%` that passes — on `unmarked.bam`, whose header keeps HG002's `samtools markdup` record, `no dup flags` FAIL → **`0.0%` PASS**; below 1,000 sampled records it is `too few reads` (a 117-record window: `no dup flags` FAIL → `too few reads` FAIL), since at a 1% duplicate rate a 41-record window holds none about two times in three; a file with neither signal still reports `no dup flags` (the same BAM reheadered without the markdup record: unchanged). A truth event **no check applies to** is now a failed `event_checked` result, so an INS-only truth VCF goes from **3/3 PASS, exit 0** to **3/4 PASS, exit 1** — M11's door, one step over. A region that cannot be queried is a WARN and a skip instead of three failed globals. No regression: `scripts/validate_pipeline.sh`'s two BAMs still score 7/19 background and 6/19 spiked with the same three global numbers, and the 10 Mb chr20 CRAM run is 0.885 s → 0.867 s, so M15's speedup is kept. | `validate.rs` (`GlobalSample`, `sample_event_regions`, `sample_region`, `not_evaluable`, the three `check_*`) | Sample across the file / event regions |
| L16 | **Fixed.** `validate`: `truncate()` panics on non-ASCII gene names (it byte-slices at a fixed offset regardless of UTF-8 boundaries); `escape_json` only escaped `\`, `"` and `\n`, so a `\t`, `\r`, or any other C0 control character below `0x20` from a gene name reached the `--json` output unescaped and broke JSON parsing. `escape_json` is the only place spike writes JSON -- `main.rs`'s `summary.json` fixtures are test-only stand-ins for an external tool's (`truvari`) output, not something spike itself serialises. `truncate` now walks back from the byte cut point to the nearest `char_boundary` before slicing; `escape_json` now matches every `char`, keeping the short escapes (`\"`, `\\`, `\n`, `\r`, `\t`, `\b`, `\f`) and emitting `\u00XX` for any other codepoint under `0x20`. **Measured** on real data (full HG002 BAM, full GRCh38 reference, a hand-built truth VCF for `chr20:38412500-38422500` with `SIM_GENE=Lγ2DLR`, a name chosen so the 34-column event-label cut lands inside the 2-byte `γ`): base binary (`8d1beba`) panics -- `byte index 31 is not a char boundary; it is inside 'γ' (bytes 30..32) of \`DEL chr20:38412500-38422500 (Lγ2DLR)\`` (thread `main`, exit unreached); the fix runs to completion (`Result: 3/5 PASS`, exit 1) with the row correctly truncated to `DEL chr20:38412500-38422500 (L...`. Separately, with `SIM_GENE=LDLR<CR>TAIL` (a bare `\r`, e.g. from a hand-edited or Windows-mangled truth VCF -- BED's tab delimiter rules out an embedded `\t` reaching this field, and `BufRead::lines()` already strips a proper `\r\n`) and `--json`: base binary's output fails `json.loads` (`Invalid control character at: line 4 column 50`); the fix's output parses, with the field decoding back to `'DEL chr20:38412500-38422500 (LDLR\rTAIL)'`. No output-file format or CLI change, so README is unchanged **Fix pass (task 27):** the new tests did not reach the `\b` / `\f` short-escape arms or a 4-byte codepoint cut mid-character; both now have one (`escape_json("a\u{8}b\u{c}c")` and `truncate("AB\u{1f600}CDEFGHIJ", 8)` -> `"AB..."`, a 3-byte walk-back against the existing test's 1). Removing the two arms, or the boundary loop, fails them -- both mutants run | `validate.rs:1637-1667` (the original row's `1637-1665` stopped 2 lines short of `escape_json`'s end) | Truncate on char boundary; escape all control chars |
| L17 | **Fixed.** R2 quality Markov chain runs backwards (no measurable effect on NovaSeq: mean Q 35.628 vs 35.627) | `synth.rs:413-421, 651-659` | Fixed along with L1 |
| L18 | **Fixed.** The row is right, and the fix is stated as one rule rather than one branch: **every `N` a synthetic read emits reports Q2**, whatever put it there — a reference `N`, padding past a contig end, or padding after a deletion sequencing error exhausted the template. The three had drifted apart: the template-exhaustion pad hard-coded Q2, while a reference `N` *and* (since L1's fix pass routed contig-end padding through the same "N in template" path) contig-end padding both took a profile-sampled score. `generate_from_template` now pushes a single `N_QUAL` = Q2 constant on both paths; the profile draw is still made, so the *quality* draw stays one per template base, but an `N` skips the `P(error)` draw a called base makes and so consumes strictly fewer random numbers — the stream is not unchanged. Q2 is also what the Markov chain carries forward (the base after an `N` is sampled from the after-a-Q2 bin — previously the chain was conditioned on a quality the read never reported). **Premise confirmed on real data:** all **682** `N` bases in **173,338** HG002 NovaSeq reads over `chr20:60001-500000` are Q2 (100%). **Measured**, `del:chr20:60101-62100` (the left breakpoint sits 100 bp inside chr20's p-telomere, whose first 60,000 bases are `N`), `--flank 2000 --seed 1`, full HG002 BAM + GRCh38: branch HEAD `8763c55` emits **31,440** `N` bases across 230 reads at **Q37 26,140 (83.1%) / Q25 2,379 / Q11 2,919 / Q2 2 (0.006%)**; after the fix, **29,975** `N` bases across 221 reads, **all Q2 (100.0%)**. (The `N` total shifts slightly because the corrected Markov state changes a few downstream error draws. **Isolating control:** a build that pushes `N_QUAL` but leaves `prev_qual` on the discarded draw reproduces HEAD's **31,440** `N` across **230** reads exactly, and is byte-identical to HEAD outside the `N` quality bytes — so the whole −4.7% is RNG-stream divergence from the `prev_qual` line, not a change in how many `N`s are produced. The 2 pre-fix Q2 `N` bases are in **kept donor reads**, written verbatim: the synthetic reads in that run carry 31,438 `N`, none at Q2.) **Carrying Q2 forward was measured, not argued:** 40 seeds at a 1 bp reference gap (`del:chr2:31593000-31594000 --flank 4000`), where an `N`-containing read is 98.3% real sequence, give **mean Q 35.621** over the non-`N` bases of `N`-containing reads when the chain carries Q2 versus **35.615** when it does not (**+0.006 Q, SE 0.022, z = +0.26**; simulated `p_err` 0.2723% vs 0.2763%), and no offset from +1 to +6 after an `N` moves by more than 1.9 SE in either direction. The feared elevated substitution error just past a reference gap does not occur, because real Illumina Q2 is rare enough (**6** of **389,429** donor bases in that window, **4** after-Q2 transitions in the whole profile against `MIN_MARKOV_OBS` = 30) that the after-Q2 bin is never usable and sampling falls through to the non-Markov levels. `test_quality_chain_carries_the_q2_an_n_reported` now pins the line: reverting `prev_qual = Some(N_QUAL)` alone leaves the other 365 tests green and turns only that one red. Q2 = byte 35, inside the FASTQ writer's `!`-`~` check (M14). README's "Quality profile" section now states the rule | `synth.rs:423-428` | Force Q2 on N |
| L19 | **Fixed.** The row reproduces: the scan's only stop condition was `insert_sizes.len() >= sample_size`, which a single-end library never reaches, so it ran to EOF. **Measured** on a 1.24 GB / 9,126,937-record single-end BAM (HG002 chr20, read 1 only with the pairing flags stripped), warm cache: branch HEAD `c19abdc` scanned all **9,126,937 records in 13.55 s** and then went on to write 2 synthetic pairs off an empty quality profile (0 donor pairs extracted, all-zero base qualities); after the fix it stops at **50,000 records in 0.09 s** — **150x** — and exits with `holds no paired reads: SAM flag 0x1 is unset on all 50000 primary records examined`. The scan is capped only for the single-end case, and single-end is *detected* from SAM flag `0x1` (set by a paired library on every read, proper pair or not) rather than inferred from an empty insert-size sample — inferring it would cap a paired BAM's fragment-length estimate at an unrepresentative sample, which is exactly L15's complaint. A paired BAM whose reads never aligned as proper pairs is therefore still scanned in full and still gets the 350/50 fallback. Spike cannot simulate from single-end input at all (extraction keeps only `is_properly_segmented` records), so the cap ends in a clear failure rather than a bad estimate. Paired input is byte-identical before and after (`R1.fq.gz`, `truth.vcf` md5 match on `del:chr20:38412500-38422500 --seed 1`). **Fix pass (task 27):** (a) the message said `(scan capped at 50000)` even when the loop reached end of file -- a 100-record single-end BAM reported `all 100 primary records examined (scan capped at 50000)`. It now says `all 100 primary records in the file` there, and names the cap only when the cap bound. (b) A capped verdict now says what it cannot rule out: 50,000 consecutive primary records with no `0x1`, in a file that does hold pairs, would read as "holds no paired reads" too. Bounded, not impossible -- extraction needs a coordinate-sorted indexed file, where a mixed library's SE and PE reads interleave at every locus -- so the message says so rather than claiming certainty. (c) `compute_stats_cram` had **no test at all**, before or after this fix, on a format M15 and L2 both showed behaving differently from BAM; it now has both paths, through a `flags` parameter on the shared `write_two_contig_cram` fixture: accept (paired, 10 primary records, insert_mean 300, read_length 100) and reject (single-end, 10 records, under the cap). (d) README had cause and effect inverted -- "stops ... *once it has established* the library is single-end"; reaching 50,000 with no `0x1` **is** the establishment | `bam_stats.rs:179-224` (message), `extract.rs:923-935` (fixture) | Cap records scanned |

### L1 + L17 · Reverse mate generated in sequencing order

The Low table has no per-entry prose, so the measured result for these two goes
here; they were one fix, as the L17 row says.

Both defects lived in the `reverse_cycles` branch of `generate_read` /
`generate_read_from_seq`. The reverse mate was generated along the reference and
reverse-complemented afterwards, so (L17) its quality Markov chain conditioned
each cycle on the cycle *after* it, and (L1) a deletion sequencing error ran the
read off the far end of its template — the haplotype path was handed exactly
`read_length` bases, with no slack at all — and the shortfall was padded with
`N` at Q2, landing at the read's 5' end after the reverse-complement. Both mates
now build a template in sequencing order (complemented and walked right to left
for a reverse read, with `INDEL_SLACK` = 10 bases past the 3' end) and share one
generation core, `generate_from_template`.

**The two rows are worded for the wrong mate.** They predate M13, which split
synthetic pairs 50/50 between F1R2 and F2R1; `read_num` and the reverse-strand
flag are independent, so since M13 both defects follow the *reverse mate*, which
is R1 half the time.

- **Measured**, 40 seeds × (`dup:chr20:38423496-38427196` + `del:chr20:38412500-38422500`),
  `--indel-error-rate 0.05`, 44 000 synthetic pairs per arm. At `8d1beba` (pre-M13)
  the original wording holds: **0.777% of synthetic R1 end in `N`, 1.152% of
  synthetic R2 start with `N`**, and nothing at the other end. At branch HEAD
  `65a76b8` (post-M13) it is split across both mates and both ends: R1
  **0.500%** start / **0.393%** end, R2 **0.564%** start / **0.530%** end. After
  the fix: **0.000%** everywhere — no synthetic read contains an `N` at all
  (0/88 000). Real donor reads keep their own `N`s, unchanged (R1 0.140%).
- **L17, measured on the same reads:** mean Q of synthetic reads 35.9094 → 35.9144
  (R1) and 35.5279 → 35.5274 (R2) — no measurable effect on NovaSeq, as the row
  says. The visible part is at the read ends, where the backwards chain had been
  distorting the per-cycle marginal: R2 cycle 0 mean Q 35.576 → 35.744, R1 cycle 150
  34.606 → 34.788, and the RMS deviation of the synthetic per-cycle profile from
  the donor reads it was learned from falls 0.0555 → 0.0501 Q (R1) and
  0.0700 → 0.0651 Q (R2).

**Fix pass 1** (review of `e69af5d` itself). Two real issues in the fix, plus
one claimed issue that did not reproduce:

- **The slack side was pinned by no test.** `generate_read` and
  `generate_haplotype_read_pair` put `INDEL_SLACK` on the correct side (past
  the 3' end — left for a reverse read), but every test that could tell used
  `indel_error_rate = 0.0`, so `slack == 0` on both the correct and a
  wrong-side build and the assertions couldn't distinguish them.
  `test_reversed_read_pair_covers_the_same_fragment` and
  `test_reversed_haplotype_pair_covers_the_same_fragment` now use
  `mock_gen_over_with_indels(..., 93, 93, 1.0)` (Q93: still no sequencing
  errors, but slack is on). **Measured**: flipping the slack side in both
  `generate_read` (`fetch_start`/`fetch_end`) and
  `generate_haplotype_read_pair` (`right_seq`'s start) now fails both tests
  ("F1R2 pair at 1000 has the wrong R2", "R1 at hap 100 matches neither end of
  the fragment"); reverted, both pass again.
- **Contig-end reverse read was shifted instead of `N`-padded** —
  `synth.rs:504-511`. `SharedReference::fetch_sequence` clamps `fetch_end` to
  the contig length, and the reverse template was built by reversing that
  clamped-short window, so its 5' end landed on the contig's last real base
  instead of on the claimed `ref_start + rl - 1` — a full-length, `N`-free
  read silently shifted left of its claimed span, while `ref_end` still
  pointed past the contig. This is a regression from the pre-`e69af5d` code,
  which N-padded the same case at the correct coordinates. **Fixed**: the
  is-reverse branch now computes how many bases the clamp dropped and
  prepends that many `N`s to the template before reversing, so the shortfall
  lands at the read's 5' start (the same "N in template" path already used
  for real reference `N`s), not a shifted window of real bases. **Measured**
  (new test `test_reverse_read_n_pads_past_the_contig_end_instead_of_shifting`,
  rl=150, slack 10, a 5 bp contig-end overhang): before, 0 `N`s and the read
  equals `revcomp(ref[contig_len-150..contig_len])`; after, exactly 5 `N`s at
  the read's start and the remaining 145 bases equal
  `revcomp(ref[ref_start..contig_len])`. Mutation check (drop the `N`-prepend):
  test fails with "expected 5 N bases ... got 0"; restored, passes.
- **Claimed: lowercase survives uncomplemented in `generate_read_from_seq`.**
  Did not reproduce. The claim was that `extract::reverse_complement` relies
  on `complement_base` (whose fallback arm is `other => other`, so it would
  leave lowercase alone), but `reverse_complement` has its own independent
  match arms that already handle `a/c/g/t` (`extract.rs:820-833`) — it never
  calls `complement_base`. A hand-verified test (revcomp of `"acgtacgtac"` is
  the literal `"GTACGTACGT"`, not derived by calling `reverse_complement`
  again) passes on the unmodified code; mutating `reverse_complement` to drop
  its lowercase arms — the literal bug as described — turns it red
  (`left: [67, 65, 84, ...]` i.e. `CATG...` vs expected `GTAC...`), then green
  again once restored. No code change made for this one; kept the test as a
  regression guard. `generate_read`'s own is-reverse path was never at risk —
  it uppercases before calling `complement()`.
- **Tightened the contract**: `generate_read` and `generate_read_from_seq` had
  no caller outside `synth.rs` (checked with `grep -rn` across `src/`), so
  both are now private (`fn`, not `pub fn`) instead of relying on a doc
  comment to keep a caller from reversing the reverse mate a second time.
- Determinism (M7) still holds: two `--seed 1` runs of
  `dup:chr20:38423496-38427196 --indel-error-rate 0.05` (the slack path) on
  the real HG002 chr20 slice give byte-identical `R1.fq.gz`, `R2.fq.gz`,
  `truth.vcf` and `replaced_reads.txt`.
- Full suite: 241 passed, 0 failed (was 239; +2 tests — the third new test
  passed against unmodified code, so it added coverage without pinning a
  fix). Clippy unchanged: 13 (bin) / 14 (test target, 12 duplicates).

### L2 · A multi-contig CRAM container leaks its other contigs' reads

noodles-cram 0.74's `Query` filters returned records on **coordinates only** —
it never compares a record's reference id with the queried one
(`noodles-cram-0.74.0/src/io/reader/query.rs`, `Iterator::next`). A container
written with several contigs in one slice (htslib's `multi_seq_per_slice=1`) is
therefore decoded whole, and every record in it whose position happens to fall
in the queried window is handed to spike. `extract_read_pairs_cram` took them,
paired them and labelled the pairs with the queried contig's name, so another
chromosome's reads entered the donor pool as if they were the event's own
background. noodles-**bam** has no such hole: its `Query` compares
`id == reference_sequence_id` before the interval
(`noodles-bam-0.73.0/src/io/reader/query.rs:72`), which is why only the CRAM
path needed a fix. Both CRAM passes now skip a record whose reference id or
mate reference id is not the queried contig's.

- **Measured, synthetic** (the review's own shape: a 2-contig CRAM built by
  `samtools view -C --output-fmt-option multi_seq_per_slice=1`, 40 chrA pairs +
  5 chrB pairs at chrA-overlapping coordinates, one container, two `.crai`
  lines at the same offset): **45 pairs extracted → 40**, i.e. exactly the 5
  chrB pairs the row names. At `--seed 2` three of those chrB pairs were
  written into the chrA spike-in FASTQ as kept originals (**3 → 0**); at
  `--seed 1` all five were suppressed instead, which is worse in a different
  way — a suppressed foreign read is donor material spike then re-tiles.
- **Measured, real reads**: HG002 chr20:38.40–38.43 Mb + chr21:38.40–38.43 Mb
  merged and written as one multi-reference CRAM,
  `del:chr20:38412500-38422500 --region chr20:38405000-38425000 --seed 1`:
  **7875 pairs → 4132**, i.e. 3743 chr21 pairs removed. 4132 is exactly what
  the same window yields from a **chr20-only** CRAM of the same reads, before
  and after the fix — so the filter removes the foreign reads and nothing else.
- **Task 8's index pruning (M15) does not help here**, measured rather than
  assumed: base `8d1beba` 7875, branch HEAD `3a782cb` (pruning in place) 7875.
  The queried contig's own `.crai` entry points at the shared container, so
  pruning keeps it and the container's other contigs come along. The M15 entry
  above is corrected accordingly.
- Scope note, now closed: `open_cram_reader_for_region` has five more callers
  that iterate a CRAM `Query` the same way and had the same hole — `loh.rs:656`
  (`count_alleles_cram`), `loh.rs:1056`, and `validate.rs:712, 823, 948`. They
  are **fixed and measured under N4**, which is where their numbers live; the
  L2 numbers below are extraction's alone.
- The equality with the single-contig control below says the filter removes the
  foreign reads and nothing else, but it says it by *count*, and a count would
  come out equal even if the mate clause were too broad — a pair can only form
  from two records the same query returned, so both its mates are on the queried
  contig whatever the mate clause does. The strong evidence is the one actually
  taken: R1, R2, `truth.vcf` and `replaced_reads.txt` are byte-identical to the
  control, not merely equinumerous. Separately, and reported rather than hidden:
  the mate clause in `pair_is_on_queried_reference` reddens **no** test even
  with the whole suite run, because `is_properly_segmented` already implies both
  mates are on one reference. It is kept for the reason that rule is another
  tool's flag, not spike's invariant — but it is not, and on well-formed input
  cannot be, independently observable. It is confined to extraction: propagating
  it to N4's five sites would drop genuine records (see N4).
- The pass-2 guard is defence in depth and is **not** independently observable:
  pass 2 only completes a pair whose other mate pass 1 already stored, and the
  pass-1 guard keeps every foreign record out of those maps. Measured —
  removing only the pass-2 guard leaves the new test green; removing the pass-1
  guard, or the comparison inside `record_is_on_queried_reference`, turns it red.

### L3 · bgzipped FASTA read as raw bytes

**Fix pass 1** (review of `7abf5ef` itself). noodles'
`indexed_reader::Builder::build_from_path` picks bgzf vs. raw bytes by
**extension alone** (`.gz`/`.bgz`; `noodles-fasta-0.47.0/src/io/indexed_reader/builder.rs:58-65`).
A bgzip-compressed FASTA under any other name — `ref.fa`, `ref.fna`,
`ref.bgzf` — still took the raw-bytes branch and still died with the exact
"beyond chromosome length" message L3 exists to remove; htslib/samtools sniff
the gzip magic instead of the extension, so spike diverged from the tool that
produced the file. Closed by peeking the file's first two bytes ourselves: if
they are the gzip magic (`1f 8b`) and the name isn't `.gz`/`.bgz`,
`ReferenceReader::open` now fails immediately, naming gzip-compressed content
as the cause. **Failed, not handled**: handling it would need a `.gzi` index,
and one built for a `.gz`/`.bgz` name almost certainly isn't sitting next to a
file that was never given that name.

- **Measured** (`reference::tests::open_bgzip_compressed_fasta_under_a_fa_extension_fails_with_clear_message`,
  a genuinely bgzf-compressed one-contig FASTA named `ref.fa`): before, `open`
  returns `Ok` (silently reading compressed bytes as sequence); after, `open`
  returns `Err` naming "gzip-compressed" and never mentions "beyond
  chromosome length". Mutation check (drop the magic-byte guard): the test's
  `Ok(_) => panic!(...)` arm fires again; restored, passes.
- **The `.gzi` hint in the open-failure message was unconditional** — a
  plain FASTA whose `.fai` is readable but whose FASTA is not (permissions,
  moved file) got the irrelevant "if this is a bgzipped FASTA..." parenthetical
  on the one path every existing run uses. Now gated on the same
  `.gz`/`.bgz` extension check the magic-byte guard uses, so it only appears
  when the file is actually named as bgzipped.
- **Test 1's fixture never exercised the `.gzi` seek itself**: it wrote
  `gzi::Index::default()` (empty), which is only correct because that
  fixture's whole record fits in one bgzf block — for any file with more
  than one block, an empty index gives the *wrong* answer past the first
  block (`gzi::Index::query` degenerates to compressed offset 0 for every
  position). Added
  `reference::tests::fetch_sequence_seeks_correctly_past_a_real_block_boundary`:
  a genuine two-block bgzf fixture with a real, non-empty `.gzi`, fetching
  non-repeating bases (`ACGTACGTAC` / `TGCATGCATG`, not the original
  all-A/all-C style, which would hide an off-by-N seek behind a repeated
  character) from 2 bases into the second block. Verified this actually
  discriminates a correct seek: swapping in an empty `.gzi` (the bug an
  empty-index fixture would hide) makes the test read `TGCATGCA` (the
  sequential-fallback answer — landing at block 2's start instead of 2 bases
  in) instead of the correct `CATGCATG`; restored, passes. The misleading
  "falls back to offset 0 for any position ... exactly right here" comment on
  the original single-block fixture is corrected to say that's true only
  because that fixture is single-block, not a general property.
- **Real-data, re-measured** (previous run's outputs were deleted; redone and
  kept under `scratch/work/task-12-L3-fixpass1/`): same command as the
  original L3 entry (`del:chr20:38412500-38422500 --seed 1` on the chr20
  slice BAM) against the real 791,434,568-byte
  `reference.fna.bgz` (`.bgz`, same code branch as `.gz`) — **base `8d1beba`
  binary**: `Loaded 1 chromosome(s) into shared reference (0 MB)` then `Error:
  DEL event start on chr20 is at or beyond chromosome length (38412500 >=
  0)`; **this branch's binary**: `Loaded 1 chromosome(s) into shared
  reference (61 MB)`, completes, `3570 kept + 292 chimeric, 989 suppressed`.
  Correctness against the same command on the uncompressed FASTA: `truth.vcf`
  identical apart from the `##reference=` line; `R1.fq.gz`/`R2.fq.gz`
  decompressed content (`zcat | md5sum`) matches exactly
  (`dc4c533dfa701a00e81a4b1d713a41d0`, `f666dfa764babe3cf07fd2fcbe0b359a`) —
  in this run the raw `.gz` bytes happened to match too
  (`8cf7964f2f02df1e289cc9c44c195ccb`, `1bd652f8ab40e4d26aeb84f64da3b0fd`),
  which is incidental (gzip's own metadata isn't guaranteed equal run to
  run), not the claim being tested.
- **README now notes the `--align` caveat** — a bgzipped `--reference` loads
  fine but `--align`/`align.sh` never runs `bwa-mem2 index`, only `bwa-mem2
  mem` against `--reference` as the index prefix, so it fails unless an
  index already exists under that exact bgzipped name. The finding as
  handed to this task said "bwa-mem2 cannot use a `.bgz`/`.gz` prefix" —
  **measured and corrected**: `bwa-mem2 index ref.fa.gz` and `bwa-mem2 mem
  ref.fa.gz ...` both succeed against a real bgzipped FASTA when the index
  was built under that name; the actual failure (confirmed: `bwa-mem2 mem`
  against a `.gz` path with no index built at that exact name exits 1,
  `ERROR! Unable to open the file: ref.fa.gz.bwt.2bit.64`) is the missing
  index at that path, not the compression itself.
- **The uncompressed path is provably unchanged**: built `7abf5ef` (the
  commit immediately before this fix pass) in an isolated `git worktree` (no
  stash, nothing uncommitted touched) and ran it against the same BAM/event
  and the real uncompressed GRCh38 FASTA. Against this branch's binary on the
  same inputs: `truth.vcf`, `replaced_reads.txt` and the raw `.gz` bytes of
  both `R1.fq.gz`/`R2.fq.gz` are byte-identical, not just the decompressed
  content.

### L4 · Final FASTQ flush error ignored

**Fix pass 1** (review of `b2d4820` itself). `r1_gz.finish()?.flush()?;
r2_gz.finish()?.flush()?;` still had the same "write error discarded" defect
this task exists to remove, just moved one level up: the `?` on R1's line
short-circuits before R2's ever runs, so a genuine R1 failure drops `r2_gz`
unfinished. `GzEncoder::drop` (flate2 1.1.9, `src/gz/write.rs:161-167`) calls
`try_finish()` and discards its error, and whatever that pushes into R2's
`BufWriter` then hits `BufWriter::drop`'s own swallowed flush — R2's error
(if it has one) never surfaces, even though the function still (correctly,
but only by luck of ordering) returns `Err` because R1's own error already
propagates. Closed by computing both `finish().and_then(flush)` results
unconditionally before propagating either, each tagged with which output
path it came from, and combining both messages when both fail
(`fastq.rs:97-124`).

- **Measured**, two new tests, both against the real (unmocked) `/dev/full`
  device:
  - `fastq::tests::test_write_paired_fastq_r2_flush_error_is_attributed_to_r2`
    (only `R2.fq.gz` symlinked): before, `Err("No space left on device (os
    error 28)")` with no indication of which stream failed; after, the same
    call's `Err` names `R2.fq.gz` and not `R1.fq.gz`.
  - `fastq::tests::test_write_paired_fastq_r1_failure_does_not_swallow_r2_failure`
    (both `R1.fq.gz` and `R2.fq.gz` symlinked): before, `Err("No space left
    on device (os error 28)")` naming neither stream — R1's error masks
    R2's entirely, the defect above, reproduced; after, the same call's
    `Err` names both `R1.fq.gz` and `R2.fq.gz`.
  Both assert on the message content rather than `is_err()` alone: a bare
  `is_err()` check already passes against the masking code in both cases
  (R1's own error still propagates when R1 fails, and R2's own error still
  propagates when only R2 fails), so it would prove nothing about the
  masking defect. Both fail for the right reason (an assertion on that
  content, not a compile error) against `b2d4820`, then pass. Mutation
  check: reverted the whole block to the pre-fix-pass sequential
  `r1_gz.finish()?.flush()?; r2_gz.finish()?.flush()?;`, reran both new
  tests — both failed on the same assertions, for the same reason, against
  the real device; restored from a `scratch/work` backup (never
  `git checkout --`, since the work was uncommitted), reverified green
  (6/6 in `fastq::tests`).
- Regression: the existing `test_write_paired_fastq_reports_error_writing_to_dev_full`
  (R1-only) and the ordinary-path real-data run below are both unaffected.
- **Real-data**: reproducing the masking scenario itself through the CLI
  isn't possible, for the same reason the original L4 real-data check
  found — a real BAM-derived run's own writes overflow the `BufWriter`
  well before the final flush, so `/dev/full` fails loudly mid-stream
  regardless of how many streams are symlinked to it. The narrow
  post-`finish()` window this fix (and this fix pass) target only exists
  at unit/small-run scale, which is exactly why the unit tests above use
  the real device rather than a mock. Confirmed the ordinary path is
  unaffected at event scale instead: same command as the original L4 entry
  (`del:chr20:38412500-38422500 --seed 1`, HG002 chr20 slice, this
  fix pass's release binary) — exit 0, 3862 pairs, `R1.fq.gz`/`R2.fq.gz`
  both 15448 lines decompressed, matching the original L4 measurement
  exactly.
- **Test-gap closed**: the original L4 test only symlinked `R1.fq.gz`;
  nothing exercised R2's own flush or the masking case, so a regression
  that dropped R2's `.flush()?` or reintroduced the short-circuit would
  have passed silently. Closed by the two tests above.
- **Style**: added `#[cfg(unix)]` to all three `/dev/full` tests (the
  original L4 test included), matching the existing convention
  (`main.rs:753, 1271, 1710`) for code using `std::os::unix::fs`. Moot in
  practice — the crate shells out to `samtools`, which isn't available off
  Unix either — but now consistent.

### L5 · Panic when mean read length > 1500

`FragmentDist::sample_in_range(rng, min, max)` (`stats.rs`) rejection-samples a
fragment length in `[min, max]`; after 1000 failures it falls back to
`self.sample(rng).clamp(min, max)`. Both call sites pass the input BAM's mean
read length as `min` and a hard-coded `1500` as `max` (`simulate.rs`'s
per-event tiling, `synth.rs`'s per-duplicate-copy fragment draw). If the BAM's
mean read length exceeds 1500bp — a long-read (PacBio/ONT) library — `min >
max`, every one of the 1000 attempts fails by construction, and the fallback's
`.clamp()` panics with Rust's own "min > max" assertion. The line numbers in
the REVIEW.md row above are stale: `synth.rs` was rewritten by task 10
(commits `e69af5d` + `3a782cb`) after this row was written; the defect and its
fix are unchanged, only the surrounding code moved. The correction was itself
stale by the end of the run: `synth.rs:594` is `ref_end:` inside
`generate_pair_from_sequence`, and the `sample_in_range` call is
`synth.rs:616`, still inside `generate_depth_pair`. **Follow the function
name.**

**Decision**: fail loudly, don't silently clamp. A long-read library can't
produce the fixed-length paired-end reads this tool generates at all — there
is no "close enough" fragment length to substitute — so clamping `min` down
to `max` (or any other silent adjustment) would produce reads shorter than
the mean of the real library and mislabel that as a spike-in run. Instead:
- `main.rs`'s new `validate_read_length` rejects a read length above
  `stats::MAX_FRAGMENT_LEN` (1500, the same constant both callers now use)
  right after it's computed from the BAM, before any per-event extraction,
  quality learning, or haplotype work starts — so the user gets one clear
  `bail!` naming the measured read length and the max, instead of a panic
  raised deep inside whichever event happens to hit the tiling loop first.
- `sample_in_range` itself also gained an `assert!(min <= max, ...)` with a
  message naming both values, replacing std's opaque `clamp` panic, as a
  second line of defense for any future caller that doesn't go through
  `validate_read_length` first.
- Both call sites (`simulate.rs:500`, `synth.rs:594`) were checked
  individually: they hit the identical failure (same `min`, same hard-coded
  `max = 1500`), so one upstream guard in `main.rs` covers both; no
  caller-specific handling was needed.

- **Measured**, three new tests:
  - `stats::tests::test_sample_in_range_min_gt_max_panics_with_clear_message`
    (`#[should_panic(expected = "sample_in_range: min (2000) > max (1500)")]`):
    before the fix, panics with std's own message,
    `"min > max. min = 2000, max = 1500"`, which doesn't contain the expected
    string, so the test fails (not a compile error); after, panics with the
    new message and passes.
  - `tests::test_validate_read_length_rejects_long_read_length` (main.rs):
    `validate_read_length(2000)` — before the fix existed (stub returning
    `Ok(())`), `.expect_err(...)` panics because the call succeeded; after,
    returns an `Err` whose message contains both `"2000"` and `"1500"`.
  - `tests::test_validate_read_length_accepts_normal_illumina_length`:
    `validate_read_length(151)` stays `Ok(())`.
  Mutation check on both fixes: reverted `sample_in_range`'s `assert!` and
  separately `validate_read_length`'s body (each backed up to
  `scratch/work/task-14-L5/` first, restored with `cat`, never
  `git checkout --`) — both new tests failed the same way as the original
  "before" runs above; restored, reverified green.
- **Real-data**: not possible with the resources available. The only BAM on
  hand (`HG002.novaseq.pcr-free.35x...bam`) is 151bp Illumina — checked
  directly (`samtools view | awk '{print length($10)}' | sort -u` over the
  chr20 slice: every read is 151bp) — so no real BAM here can drive the mean
  read length over 1500 to reproduce the panic, and fabricating a synthetic
  long-read BAM was out of scope for this fix. Confirmed instead that the
  ordinary short-read path is unaffected: `del:chr20:38412500-38422500
  --seed 1` against the HG002 chr20 slice (`read_len=151` in the log) runs to
  completion (exit 0, 3862 pairs written), identical in shape to the L4
  real-data run.

**Fix pass 1** (review of `158356e` itself). The two new `validate_read_length`
tests exercised only far-below (151) and far-above (2000) values; nothing
pinned the exact threshold, so a future off-by-one (`>` becoming `>=`) would
pass silently. Added `test_validate_read_length_accepts_exact_max`
(`validate_read_length(1500)` must stay `Ok`) and
`test_validate_read_length_rejects_one_above_max` (`validate_read_length(1501)`
must be rejected, message citing both `1501` and `1500`). Mutation-checked:
flipping `read_length as i64 > max` to `>= max` in `main.rs` (backed up to
`scratch/work/task-14-L5-and-bookkeeping/` first, restored with `cat`, never
`git checkout --`) fails only the 1500-accepts test
(`assertion failed: validate_read_length(1500).is_ok()`); restored, both new
tests pass. Full suite 260/260 (258 + 2 new); clippy unchanged, 13 (bin) / 14
(test, 12 duplicates).

Two Minors raised alongside this finding were checked and did not need a
change: the `bail!` message already names the implication — "spike simulates
fixed-length paired-end reads and does not support long-read (PacBio/ONT)
libraries" is already the full text of the message added in `158356e`, so
there is nothing to add. `MAX_FRAGMENT_LEN`'s name describes its
fragment-length role but the constant also gates read length; the doc
comments on the constant (`stats.rs`) and on `validate_read_length`
(`main.rs`) already say so explicitly, and no name was found that captures
"cap on both fragment length and the read length that must fit inside a
fragment" more clearly than the existing name plus those comments, so it was
left as is.

## Uncommitted changes

| File | Change | Assessment |
| --- | --- | --- |
| `extract.rs` | `safe_noodles_position` no longer `expect`s | Correct, no behaviour change. All call sites pass `start + 1`, `end` for 0-based half-open input. |
| `synth.rs` | Clamp sampled quality to Phred+33 range | Correct in itself, but it hides M14 (missing qualities) instead of fixing it. |
| `loh.rs` | Warn when gVCF chromosome names don't match | **Fixed.** Was wrong in both paths. For `.vcf.gz` it never fired (`bcftools view -r chr17:…` on a `17` file returns nothing, exit 0). For plain VCF it fired whenever the region had no het SNPs and other chromosomes existed. "Falling back to pileup" was only true when the result was empty; if bcftools errors (e.g. unindexed `.gz`), LOH is skipped. |
| `simulate.rs` | `unreachable!` → `bail!`; two new `simulate_event` tests | `bail!` is fine. The tests only assert `> 0` / non-empty and would pass with H5, 100% suppression, or 2 tiled reads. |

- `loh.rs` mismatch warning **Fixed**, and measured: the verdict now comes from the gVCF's own `##contig` names, read straight from the header (plain or bgzipped, no index and no bcftools call), instead of from an empty result, and every outcome names what spike does next. Measured on `1d61e99` + fix, `del:chr20:38412500-38422500 --seed 1`, chr20 37.5–41.5 Mb HG002 slice: a `.vcf.gz` naming the chromosome `20` while `chr20` is asked for gave **0 warnings before, 1 after** ("names chromosome '20', not 'chr20' … Falling back to pileup-based het SNP detection", and it does fall back — 15 het SNPs from pileup); a matching plain VCF whose queried window simply holds no SNPs while `chr21` records exist gave **1 spurious warning before, 0 after**; an unindexed `.vcf.gz` said `bcftools exited with status exit status: 255` before and now says `… 'noindex.vcf.gz': Failed to open …: could not load index. LOH is skipped for this region: original reads are suppressed at random.` On the working path (matching, indexed HG002 chr20 gVCF) R1/R2 FASTQ and `truth.vcf` are byte-identical to `1d61e99`.

## Test gaps

- New and existing `simulate` tests assert "> 0" instead of expected counts. `test_tandem_dup_tiling_count_uses_full_length` asserts `n > 100` where its comment expects 188. `_junction_pairs` is computed and never asserted.
- `test_suppression_*` re-implements the suppression loop instead of calling `simulate_event`. The LOH branch has no test.
- `haplotype.rs` tests never call a real `VariantHaplotype::from_*` constructor (so H4 and M4 went unnoticed), though `SharedReference::from_sequences` exists for this.
- BND parser tests encode the wrong orientation (H3).
- No tests for: `truth.rs`, main-level orchestration with several events or `--region`, minus-strand BEDs, script generation, R1/R2 geometry in `generate_read_pair`, `indel_error_rate > 0`, extraction from a real BAM/CRAM, FASTQ round trip.
- A good first regression test for most of H1–H8: simulate on a small synthetic reference and assert the realized VAF / depth / breakpoint position within a tolerance.

## Design notes

- Coordinate conventions are mixed: `del/dup/inv` take VCF POS/END meaning, `snp` is 1-based, `--region` is 1-based inclusive. `del:chr1:0-100` is rejected with "coordinates must be >= 1" although the documented convention is 0-based.
- `del:LDLR:4-8` parses as coordinates on a chromosome named "LDLR". An unknown fusion suffix (`:inverted`, `:rev`) silently gives a forward fusion.
- VCF input ignores FILTER and GT; 0/0 and non-PASS records are simulated.
- Coverage is estimated once, ±1 kb around the first covered breakpoint side, and applied to the whole haplotype. `n.max(2)` emits 2 pairs even at zero coverage. An empty read pool silently gives a constant-Q20 profile. (Now tracked and fixed as [N5](#n5--an-empty-donor-pool-is-simulated-from-anyway), which also corrects "constant-Q20": the *profile* measures Q0, the *emitted* reads are Q20. `n.max(2)` no longer applies at zero coverage at all -- the event is refused --
and where the coverage is real but the request rounds below 2 it still applies,
now with a warning naming the realized-vs-recorded fraction.)
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
