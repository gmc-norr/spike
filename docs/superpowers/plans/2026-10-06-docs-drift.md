# The README says what master's code does

**Asked 2026-10-06.** The user asked whether the documentation is up to date. Three read-only audits compared README.md, in three parts, with the code at master `97aa3b4`, and I re-read the code behind every finding. About 310 statements match. The 30 below do not. The user picked "fix all of them, plan → fix → re-check".

## Gate A

1. **Principle.** The README describes what the code on master does. A number it gives as measured comes from a run of the current code.
2. **What would kill it.** After the edit, a stale phrase below is still in the README. Or a new sentence disagrees with the code line it cites. Or a re-measured number is not in that run's output. Or the code behaves differently.
3. **Refuted before?** No. `2eb9309` and `7304113` were earlier "README up to date" passes. This one finds what changed since: RF14 (`54c238f`) made `split_reads` advisory for a DEL, and the review fixes, the read-length fix and the align-flags fix all added behaviour.
4. **Simplest thing.** Edit the text, and only the text. The `src/` changes are two comments, one doc comment and the `--aligner` help text. The one `scripts/` change is an error message.
5. **Inputs, and how each was checked.** Each item names the code line I read on `97aa3b4`. The README's CLI block (lines 1443-1558) was diffed against `spike --help` from master's binary and is identical.

Untouched: the "about 60 bp" wording (deferred by the user), the dated plan and review docs, and line 762-764's format example, which is self-consistent and names no run.

## The items

| # | README line | Stale text (K1 greps for it) | Fix | Code |
|---|---|---|---|---|
| D1 | 32 | `By default a non-additive event removes only reads from its donor pool;` | also the duplicates of the pairs it removes | `main.rs:902-919` |
| D2 | 55 | `**Rust** (edition 2021 or later)` | the oldest toolchain K5 builds and tests with | `exon.rs:204` (`Option::is_none_or`, 1.82); `indexmap 2.13.0` declares `rust-version = 1.82`, the highest in Cargo.lock |
| D3 | 76-77 | `for the tests that run the generated \`align.sh\` / \`merge.sh\` and` | also `fastq.sh`, which needs `gzip` or `pigz`, `awk`, `mkfifo`, `realpath` and `mktemp` | `main.rs:3356`; the script text in `main.rs:2236-2400` |
| D4 | 80 | `667 passed; 2 failed; 1 ignored` | M1's line | measured |
| D5 | 258 | `So before any other work, spike counts` | before it extracts any reads | `main.rs:651`, after the stats, header and reference reads |
| D6 | 486 | `` `SM` of the original BAM's first `@RG` line`` | the first `@RG` line that has an `SM` | `bam_stats.rs:282` |
| D7 | 583-585 | `real checks alone: the advisory rows below are described here and in the` | the help lists the advisory rows too | `validate.rs:831-852` |
| D8 | 593-594 | `so no event is left out however many the truth VCF holds` | unless it holds more than 200,000 events; then the cap runs out and spike warns | `validate.rs:2463, 2561, 2567, 2603-2610` |
| D9 | 598-600 | `` `[global] sample: 117 records over 1 of 1 event regions` `` | M2b's two log lines, which end `(up to N each)` | `validate.rs:2626`; measured |
| D10 | 606-607 | `the sample fails only if no region could be read at all` | it stops with an error when no record was sampled and at least one region could not be queried | `validate.rs:2612-2617` |
| D11 | 776-779 | `` `coverage_any_mapq` and `split_reads_each_end` and the line reads `Advisory: 2 `` | M2's rows and line | `validate.rs:622` (RF14); measured |
| D12 | 785-788 | `` `scripts/validate_pipeline.sh` is such a`` | the script recounts the non-advisory rows from each check's `advisory` flag; the `validate.rs` comment says the same | `validate_pipeline.sh:329-337, 888-893`; `validate.rs:3726-3729` |
| D13 | 1069, 1147 | `The three rows that measure nothing` | add the `no coverage` and `too shallow (…)` outcomes | `validate.rs:1885-1923, 2822-2838` |
| D14 | 1246-1253 | `` `13/13 PASS` and exits 0`` | M3's numbers | measured |
| D15 | 1264-1266 | `now scores **12/12 PASS, exit 0** at **0.50, 0.53 and 0.42**` | M4's numbers | measured |
| D16 | 1355-1367 | (missing row) | `sim.bam.bai` | `main.rs:1505` |
| D17 | 1421 | `The truth VCF contains one record per simulated event with:` | two for a fusion, `sim_fus_N` and `sim_fus_N_mate` | `truth.rs:293-306` |
| D18 | 1421-1432 | (missing fields) | `SIM_EXONS` (DEL records), `MATEID` (BND), FORMAT `GT` (`1/1` when the requested fraction is ≥ 0.9) | `truth.rs:192, 196, 201, 232, 267, 306, 468-473` |
| D19 | 1579 | `` `validate_interval` / `validate_point` `` | `validate_range` / `validate_point` | `main.rs:1646, 1683` |
| D20 | 1580 | `check your interval\` \| \`main.rs\`` | `exon.rs` (`parse_region_spec`) | `exon.rs:340, 361` |
| D21 | 1583 | `` Remove these events from the input` \| `carried.rs` `` | `carried.rs` (the rule) and `main.rs` (`refuse_carried_alleles`) | `main.rs:1021, 1056` |
| D22 | 1634 | `read-name shape and single-end check` | also the `@RG` sample name and the input's bwa `@PG` options | `bam_stats.rs` `sample_name`, `bwa_options` |
| D23 | 1663, 1686 | ``are suppressed at rate `P = VAF` `` and `are suppressed at the VAF rate` | the per-copy rate, as Read suppression details states it | `simulate.rs:259`; `synth.rs:997-1003` |
| D24 | 1834, 1846 | `` `[read_length, 1500]` `` with no trimmed case | `[1, 1500]` for an adapter-trimmed library | `types.rs:211-216`; `main.rs:1916`; `simulate.rs:987` |
| D25 | 1894 | `Each level requires at least 30 observations before it is used` | levels 1-3 need 30; level 4 takes any cycle with one | `synth.rs:28, 35, 388, 397, 407, 414` |
| D26 | 1946 | `so it is the smallest pool at which any level of the model is trained to its own threshold` | 30 matches levels 1-3's per-bin count; it does not guarantee any bin reaches it; same fix in `MIN_DONOR_PAIRS`'s doc comment | `synth.rs:414`; `main.rs:1850-1857` |
| D27 | 2188 | `can be overridden with an` | each listed variable, then `PATH`, then `EXTRA_TOOL_DIRS`; `python3` from `PATH` only | `validate_pipeline.sh:41-77, 297, 330, 383` |
| D28 | 2203-2205 | `step N (1-8); a value outside that range is refused` | 0, the default, runs every step; the script's own error message says the same | `validate_pipeline.sh:112, 182-184` |
| D29 | 2248-2249 | `the report parsed and at least one check passed` | at least one non-advisory check | `validate_pipeline.sh:888-893` |
| D30 | 1443-1558 | (help text) | `--aligner`'s help says the bwa-mem2 preset copies the input's `@PG` options; the README block is regenerated from `--help` | `main.rs:146-148` |

## Measurements (locked recipes)

All runs use master `97aa3b4`'s release binary, copied to the scratchpad. The slice is `chr20:38500000-40200000` of the 35x HG002 BAM, indexed, the same slice `7304113` used. `--seed` is left at its default.

- **M1.** `cargo test --release` with every `.pixi` directory removed from `PATH`. That takes samtools, bwa-mem2, bcftools and minimap2 off it; bgzip, tabix, delly, truvari and bowtie2 are not on PATH at all.
- **M3.** `--event del:chr20:38900000-38910000 --event ins:chr20:39000000:300` on the slice, then `align.sh`, `merge.sh` and `spike validate` on `merged.bam`. Record the `Result:` line, the exit status, the `Advisory:` line and the `ins_planted` row.
- **M2.** M3's `truth.vcf` with only the DEL record kept, and both census INFO fields and both census header lines removed, through `spike validate` on M3's `merged.bam`. Record the advisory rows and the `Advisory:` line.
- **M2b.** M2 again with `--flank 100`. Record the two `[global]` log lines that D9 quotes.
- **M4.** `snp:chr20:39000000:TGG:T`, `snp:chr20:39100000:T:TCCGG` and `snp:chr20:39200000:AT:GC` on the slice, aligned, merged and validated the same way. Record the `Result:` line, the exit status, the three `allele_freq` observed values and the `Advisory:` line.

The controls that D14 and D15's paragraphs quote from earlier runs are not re-run. Their text already says they come from an earlier run.

## Checks (locked)

- **K1, the stale text is gone.** Each stale string in the appendix has 0 hits in its file, and each new string at least 1. **Control:** on master every stale string has 1 hit and every new string 0 (run before this commit, to prove the strings).
- **K2, the CLI block is `--help`.** The README's CLI block equals the new binary's `spike --help`, byte for byte. **Control:** master's block differs from the new `--help`, because D30 changes the help.
- **K3, measured numbers are fresh.** Every number D4, D9, D11, D14 and D15 put in the README appears in the saved M1-M4 output.
- **K4, the code is unchanged.** Every changed line in `src/` is a comment, a doc comment or the `--aligner` help. The `scripts/` diff is the one error-message line. `cargo test --release` passes 730. A run of `del:chr20:14530000-14531000 --seed 1` on the 35x BAM gives `R1.fq.gz`, `R2.fq.gz`, `truth.vcf`, `replaced_reads.txt`, `fastq_removed_reads.txt` and `align.sh` byte-identical to master's binary.
- **K5, the Rust floor.** `cargo +1.82 build --release --locked` succeeds and `cargo +1.82 test --release --locked` passes. If 1.82 fails, D2 names the first of 1.85 and current stable that passes.
- **K6, every citation was read.** The result's table names the code line behind each item, read on this branch.

## Appendix: K1's search strings (locked)

Tab-separated: item, file, the exact string, searched with `grep -cF`. Before the fix, every stale string was found once, and every new string zero times, on `97aa3b4`. That run is the control, done before this commit so the strings themselves are known to be right.

Stale (0 hits each after the fix):

```
D1	README.md	By default a non-additive event removes only reads from its donor pool;
D2	README.md	**Rust** (edition 2021 or later)
D3	README.md	for the tests that run the generated `align.sh` / `merge.sh` and
D4	README.md	667 passed; 2 failed; 1 ignored
D5	README.md	So before any other work, spike counts
D6	README.md	`SM` of the original BAM's first `@RG` line
D7	README.md	real checks alone: the advisory rows below are described here and in the
D8	README.md	so no event is left out however many the truth VCF holds
D9	README.md	`[global] sample: 117 records over 1 of 1 event regions`
D10	README.md	the sample fails only if no region could be read at all
D11	README.md	`coverage_any_mapq` and `split_reads_each_end` and the line reads `Advisory: 2
D12	README.md	`scripts/validate_pipeline.sh` is such a
D12	src/validate.rs	`scripts/validate_pipeline.sh`'s step-5 guard is that consumer
D13	README.md	The three rows that measure nothing
D14	README.md	`13/13 PASS` and exits 0
D15	README.md	now scores **12/12 PASS, exit 0** at **0.50, 0.53 and 0.42**
D17	README.md	The truth VCF contains one record per simulated event with:
D19	README.md	`validate_interval` / `validate_point`
D20	README.md	check your interval` | `main.rs`
D21	README.md	Remove these events from the input` | `carried.rs`
D22	README.md	read-name shape and single-end check
D23	README.md	are suppressed at rate `P = VAF`
D23	README.md	are suppressed at the VAF rate
D24	README.md	`[read_length, 1500]` -- so the count
D24	README.md	the observed insert sizes in `[read_length, 1500]` and nothing else
D25	README.md	Each level requires at least 30 observations before it is used
D26	README.md	so it is the smallest pool at which any level of the model is trained to its own threshold
D26	src/main.rs	the only level an ordinary pool always reaches
D27	README.md	can be overridden with an
D28	README.md	step N (1-8); a value outside that range is refused
D28	scripts/validate_pipeline.sh	must be a step number from 1 to 8
D29	README.md	the report parsed and at least one check passed
```

New (at least 1 hit each after the fix):

```
D16	README.md	| `sim.bam.bai` |
D18	README.md	`SIM_EXONS`
D18	README.md	`MATEID`
D18	README.md	FORMAT `GT`
D30	README.md	of the input BAM's bwa @PG line
D30	src/main.rs	of the input BAM's bwa @PG line
```

## Result (2026-10-06): supported

The fix is `0c4e406`. Master `97aa3b4`'s release binary is `4265ba92`; the new one is `e21d282b`. Scratch: `scratchpad/ddp/`.

**K1: PASS.** All 32 stale strings have 0 hits. All 6 new strings have at least 1 hit. Before the plan commit, on `97aa3b4`, every stale string had exactly 1 hit and every new string 0, so the check can fail.

**K2: PASS.** The README's CLI block (116 lines) equals the new `spike --help` byte for byte. **Control:** master's block differs from the new help in one line, `--aligner`'s.

**K3: PASS.** Every number written is in the run's saved output:

| Item | Written | From |
|---|---|---|
| D4 | `728 passed; 2 failed; 1 ignored`, the two being the bcftools tests | M1: the suite with `PATH` = `~/.cargo/bin:/usr/sbin:/usr/bin:/sbin:/bin`. On that PATH, samtools, bwa-mem2, bcftools, minimap2, bgzip, tabix, delly, truvari and bowtie2 are all missing. `/usr/local/bin` holds a samtools too, so dropping `.pixi` alone would not have done it |
| D9 | `[global] sampled 1976 records from chr20:38899901-38910100`, `[global] sample: 1976 records over 1 of 1 event regions (up to 200000 each)` | M2b |
| D11 | advisory `coverage_any_mapq`, `split_reads`, `split_reads_each_end`; `Advisory: 3 checks, 3 PASS, 0 FAIL` | M2 (`8/8 PASS`, exit 0) |
| D14 | `ins_planted` 20, `ins_reads` 16, `15/15 PASS`, exit 0, `Advisory: 9 checks, 9 PASS, 0 FAIL` | M3 |
| D15 | `12/12 PASS`, exit 0, `allele_freq` 0.49, 0.36 and 0.58, `Advisory: 6 checks, 6 PASS, 0 FAIL` | M4 |

**K4: PASS.**
- Every changed line in `src/` starts with `//`.
- The one changed line in `scripts/` is `--skip-to`'s error message.
- `cargo test --release`: `730 passed; 0 failed; 1 ignored`.
- `del:chr20:14530000-14531000 --seed 1` on the 35x BAM: `R1.fq.gz`, `R2.fq.gz`, `truth.vcf`, `replaced_reads.txt`, `fastq_removed_reads.txt` and `align.sh` are byte-identical between the two binaries.

**K5: PASS.** With rustc `1.82.0 (f6e511eec 2024-10-15)`, `cargo +1.82 build --release --locked` exits 0. `cargo +1.82 test --release --locked` gives `730 passed; 0 failed; 1 ignored`, both on `97aa3b4`'s code and again on `0c4e406`. The first build was started before the plan commit; its logs were read after it.

**K6: PASS.** Every code line in the items table was read on this branch, either while the plan was written or while the fix was made.

**Beyond the plan.**
- The output-file table also gains `align.log` and `merged.bam.bai`, because a real run directory (M3's) holds both.
- Line 32's "suppressed at the target VAF rate" stays. It is an overview, and the per-copy rates average to the VAF (`synth.rs:995-996`).

**Side effects on this machine.** rustup installed toolchain 1.82 and updated itself from 1.29.0 to 1.29.1.
