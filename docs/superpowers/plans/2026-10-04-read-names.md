# Read names: spike's reads named like the input's, and marked SPIKE

**Asked 2026-10-04.** The user wants two things:
- spike's reads can be told apart from real ones;
- every tool in the raredisease pipeline still works on them.

Real names differ from one dataset to the next, so spike must copy the shape of the input's names rather than use one fixed format. The user picked this option: put `SPIKE` in the machine-name part of the name.

## Measured before this plan

- **What raredisease runs.** The from_rv HG001 30x run (`pipeline_info/nf_core_raredisease_software_mqc_versions.yml`) uses Picard 3.3.0, FastQC 0.12.1, bwa-mem2 2.2.1 and GATK 4.5.0.0. Picard MarkDuplicates runs with `--READ_NAME_REGEX <optimized ...>` (the default) and `--OPTICAL_DUPLICATE_PIXEL_DISTANCE 100`.
- **Picard leans on the names.** Its metrics count 10,361,338 optical duplicate pairs out of 17,385,340 duplicate pairs.
- **The read group does not come from the names.** It is `ID:<fastq file name>`, set by bwa-mem2 `-R` from the sample sheet.
- **How Picard 3.3.0 reads a name.** Read with `javap` from `picard/sam/util/ReadNameParser.class`:
  - it splits the name on `:` (byte 58) and accepts 5 or 7 fields;
  - the last three fields are the tile, x and y;
  - any other count logs, once per run: `Default READ_NAME_REGEX '%s' did not match read name '%s'. You may need to specify a READ_NAME_REGEX in order to correctly identify optical duplicates. ...`
  - such a read gets no position, so it can never count as an optical duplicate.
- **The input names differ by machine.** The first read of each local BAM, and of the hospital BAM above:
  - `A00744:46:HV3C3DSXX:2:1221:8775:9361` (GIAB HG002 NovaSeq 6000; 4 BAMs);
  - `D00360:96:H2YLYBCXX:1:2105:5916:51581` (GIAB HG002 HiSeq 2500; 2 BAMs);
  - `LH00352:…` (the hospital's NovaSeq X).

  All have 7 fields. 200,000 hospital chr20 reads span 247 lane:tile pairs.
- **spike's names today:**
  - `ev{:04}_hap_{:06}` (`simulate.rs` `tile_haplotype_reads`);
  - `ev{:04}_dup_depth_{:06}` (`synth.rs` `generate_dup_depth_copies`, `--dup-model junction` only).

  They are written to FASTQ as `@<name>/1` and `/2`.
- **Code that recognises spike's names:**
  - `validate.rs` `planted_read_prefix` (the `ins_planted` and `del_planted` rows), pinned to the namer by `simulate.rs:3004` `test_an_insertions_tiled_reads_carry_the_prefix_validate_looks_for`;
  - `scripts/review_sv_model.py:169`, `startswith("ev")`;
  - `scripts/rf13_planted.py`, `rf13b_planted.py` and `rf14_planted.py`. These read runs that are already done, so they stay as they are.
- **spike never edits a real read in place.** Grepping `src/` for writes to `seq1`/`seq2` finds hits only in test modules. Every read spike changes is a new read, so a mark on spike's reads covers all of them, and nothing else.

## Design (locked)

**Learning the shape.** `bam_stats` already reads the input's first 50,000 primary, mapped, non-duplicate, non-QC-fail records. It now also keeps their names.

- Each name falls into one class:
  - **7-part:** it splits on `:` into exactly 7 fields, and the last three are each 1-9 ASCII digits;
  - **5-part:** the same, with exactly 5 fields;
  - **other:** anything else.
- **The shape** is the class most sampled names fall into. A tie, or no names at all, gives **other**.
- For a 7-part or 5-part shape, among the names of that class:
  - **middle:** the most common run of middle fields. That is fields 2-4 (`RUN:FLOWCELL:LANE`) for 7-part, and field 2 (`LANE`) for 5-part. A tie goes to the smallest string, so the choice is deterministic.
  - **tiles:** the distinct tiles of the names with that middle.
  - **x and y ranges:** the min and max of those names' x, and of their y.

**The name.** spike's internal names stay as they are (`ev0001_hap_000123`, `ev0001_dup_depth_000123`). Each is written as:
- `SPIKE_<internal>:<middle>:<tile>:<x>:<y>` for a 7-part or 5-part shape;
- `SPIKE_<internal>` for **other**.

  In the **other** case Picard cannot read the real names either, so a position on spike's reads would give them nothing.

The tile, x and y come from a hash of the internal name: FNV-1a 64, then splitmix64 steps, one per value. They never come from the run's random stream. So every base, quality, position and choice in a run stays exactly as before, and only the names change.

**Prefixes.** `planted_read_prefix(n)` becomes `SPIKE_ev{:04}_hap_`. `review_sv_model.py` tests `startswith("SPIKE_")`.

**What the user sees:**
- a log line naming the learned shape and one example of spike's names;
- a README section on read names.

**Tests, written first and seen red:**
- the class rule for 7-part, 5-part, 6-part, non-digit tails, SRA-style names and empty input;
- the majority and tie rules for the shape and for the middle;
- tiles taken only from names with the chosen middle;
- x and y within range;
- the name for each shape, and that it parses back under Picard's rule (5 or 7 fields, last three digits);
- the same internal name always gives the same name;
- naming draws nothing from the random stream (same seed, with and without a shape: identical bases);
- both namers (`_hap_` and `_dup_depth_`) use it;
- `main` hands the learned shape to the generator;
- the `planted_read_prefix` pin.

**Mutation checks (each must turn a test red):**
1. a 6-part name counts as 7-part;
2. the digits check is dropped;
3. the middle comes from the first name, not the most common;
4. a tie goes to 7-part instead of other;
5. the tile is fixed at the first tile seen;
6. x is not clamped to the range;
7. the `SPIKE_` prefix is dropped;
8. the hap namer skips the shape;
9. the dup-depth namer skips the shape;
10. `main` does not hand over the shape;
11. the tile is drawn from the run's random stream;
12. `planted_read_prefix` goes back to `ev`.

## Checks (locked before running)

Two binaries:
- **old** = master `3684f74`;
- **new** = this branch.

Both are built with a separate `CARGO_TARGET_DIR`, and both are md5'd. Runs use 16 threads. Each run's stderr goes to a log file that is read in full, and its exit status is checked.

**K0: only the names change.** Every comparison runs old against new.
- **K0a:** run 0 of the read-length K2 rerun. That is the first command in `readlen/full/spike.log`: HG002 hospital 30x, 700 events (654 small variants, 46 DUPs), `--seed 1 --threads 16 --align`. It is run with each binary from its own directory, using a relative `-o out`.
- **K0b:** the same BAM, reference and seed, with `--dup-model junction --align` and these events:
  - `del:chr20:10500000-10500300;af=0.5`
  - `ins:chr20:10600000:100;af=0.5`
  - `dup:chr20:10700000-10700200;af=0.5`

  If spike refuses one (RF8), all three are moved +100,000 bp. The move is the same for both binaries, and is tried at most 3 times. Then `spike validate --bam out/sim.bam --truth out/truth.vcf --reference REF`.
- **Pass, in both runs:**
  - **R1.fq.gz and R2.fq.gz:** the same record count, in the same order.
    - Every record whose old name starts with `ev` has the new name `SPIKE_<old name>` followed by either nothing or `:<middle>:<tile>:<x>:<y>`. The middle is the learned one, the tile is in the learned set, and x and y are in the learned ranges.
    - Every other name is unchanged.
    - Every sequence and quality is identical.
  - **truth.vcf, events.bed, replaced_reads.txt and merge.sh:** byte-identical.
  - **sim.bam:** record by record, with names mapped as above. Flag, contig, position, MAPQ, CIGAR, mate contig, mate position, TLEN, sequence and quality are all identical.
  - **K0b:** the `spike validate` stdout is identical, including the `ins_planted` and `del_planted` rows.

**K1: the tools read the names.** Run on K0a's and K0b's outputs, old and new alike.
- **Picard 3.3.0 MarkDuplicates.** The jar is md5 `63ed3f5d6da8934d4199e06b1ac3c176`. It runs on `sim.bam` with the defaults raredisease uses (optimized `READ_NAME_REGEX`, pixel distance 100).
  - **Pass:** new exits 0, and its stderr does not contain `did not match read name`.
  - **Control:** old's stderr *does* contain it. If old's doesn't, the Picard check is no control, and K1-Picard is reported inconclusive.
  - **Reported, not judged:** `READ_PAIR_OPTICAL_DUPLICATES`, old against new.
- **FastQC 0.12.1** on `R1.fq.gz`.
  - **Pass:** new's `fastqc_data.txt` has a `>>Per tile sequence quality` module.
  - **Control:** old's lacks it. If old has it too, FastQC never needed the names, and that is reported. It is not a fail.
- **The mark finds exactly spike's reads.**
  - **Pass:** in new `sim.bam`, the distinct names starting `SPIKE_` number the same as the distinct names starting `ev` in old `sim.bam`, and the other names are the same set in both.

**K2: the shape is learned right.**
- **Inputs:** new spike, with one event `snp:chr20:10400000:<ref>:<alt>` (no `--align`), on three inputs:
  - the hospital HG002 BAM (taken from K0a's log);
  - `data/validation/hg002_novaseq_chr20.bam` (A00744);
  - `data/giab_hg38/HG002/HG002.GRCh38.chr20.bwamem2.bam` (D00360).
- **Pass:** for each input, the logged class, middle, tile set and x/y ranges equal what an independent `LC_ALL=C awk` gives over the same 50,000 records. Those records are the first 50,000 primary, mapped, non-duplicate, non-QC-fail records from `samtools view -F 0xF04`.

**Outcomes:**
- **K0, K1 and K2 pass:** supported. Merge on the user's word.
- **K0 fails:** stop. A name reaches something it should not, and the cause is found before anything else.
- **K1 Picard pass fails:** refuted for Picard. Look at why.
- **K2 fails:** the learner has a bug. Fix it with a test first, then rerun K2.

**Not run here:** the full raredisease pipeline, which needs the hospital. Also GATK's mitochondrial steps (RevertSam, SamToFastq, MergeBamAlignment): they pair mates by name, and old and new names are both unique, so they cannot tell old from new.

## Result (2026-10-04)

**Verdict.**
- K0 passes, and so does K2.
- K1:
  - **FastQC: passes**, and its control works;
  - **the mark: passes**;
  - **Picard: inconclusive by the locked rule.** The control did not fire, because the plan locked the wrong warning text. Evidence found afterwards is below. It points one way, but it is not the locked result.

**Setup:**
- binaries: old `918e92fc` (master `3684f74`, the same md5 as the read-length binary) and new `b558405f` (`03debe6`);
- Picard jar: `63ed3f5d`.

Every run exited 0. `spike validate` exited 1 for both binaries on K0b's `sim.bam`, and its output is what K0 compares. K0b needed no shift.

**K0: only the names change. Passes.** Old and new, record by record:

|  | K0a (700 events) | K0b (DEL, INS, DUP junction) |
|---|---|---|
| R1 / R2 records | 1,667,537 each | 7,193 each |
| spike's reads (renamed as the plan says) | 146,837 | 457 (4 of them `_dup_depth_`) |
| bad names, or changed bases or qualities | 0 | 0 |
| truth.vcf, events.bed, replaced_reads.txt, merge.sh | identical | identical |
| sim.bam records (11 SAM fields, names mapped) | 3,338,735, 0 differ | 14,425, 0 differ |
| `spike validate` stdout | | identical; `del_planted` 8 PASS, `ins_planted` 20 PASS |

The learned shape for the hospital BAM is `45:227NC2LT1:2`, with 247 tiles, x 1000-52208 and y 1016-29759. An example name is `SPIKE_ev0001_hap_000000:45:227NC2LT1:2:2218:22853:29387`.

The name checker was tested against bad names: a wrong middle, an unseen tile, x out of range, no tail, and the old name. It rejected each one.

**K2: the shape is learned right. Passes.** spike's log equals `shape.awk` over the same 50,000 records on all three inputs (classes, middle, the full tile list, and the x and y ranges):

| input | 7-part / 5-part / other | middle | tiles | x | y |
|---|---|---|---|---|---|
| hospital NovaSeq X | 50,000 / 0 / 0 | `45:227NC2LT1:2` | 247 | 1000-52208 | 1016-29759 |
| `hg002_novaseq_chr20.bam` (NovaSeq 6000) | 50,000 / 0 / 0 | `46:HV3C3DSXX:2` | 936 | 1027-32922 | 1000-37043 |
| `HG002.GRCh38.chr20.bwamem2.bam` (HiSeq 2500) | 50,000 / 0 / 0 | `93:H2YHMBCXX:2` | 64 | 1043-21291 | 2058-101292 |

The HiSeq BAM mixes several runs. Its first read is `D00360:96:H2YLYBCXX:...`, but the most common middle is `93:H2YHMBCXX:2`.

`shape.awk` gave the unit tests' answers on their six cases. The judge failed when given mismatched sides: the hospital awk against the NovaSeq log, and one tile dropped.

**K1: the tools read the names.**

| | K0a old | K0a new | K0b old | K0b new |
|---|---|---|---|---|
| Picard exit | 0 | 0 | 0 | 0 |
| locked text `did not match read name` | no | no | no | no |
| any `WARNING` line | 1 | **0** | 1 | **0** |
| optical / duplicate / examined pairs | 0 / 29 / 1,667,537 | 0 / 29 / 1,667,537 | 0 / 0 / 7,193 | 0 / 0 / 7,193 |
| FastQC `>>Per tile sequence quality` | **absent** (10 modules) | **present** (11) | **absent** (10) | **present** (11) |

- **Picard: inconclusive by the locked rule.** New never prints the locked text, but old doesn't either, so the control does not fire.

  The text was read from `ReadNameParser.class`. Picard 3.3.0 MarkDuplicates prints a different message for these names, from `AbstractOpticalDuplicateFinderCommandLineProgram`: `A field field parsed out of a read name was expected to contain an integer and did not. Read name: ev0291_hap_000133.`

  **Found after looking (not locked):**
  - old prints that message once in each run, and new prints no warning of any kind.
  - A diagnostic in `scratchpad/names/diag` copied one of spike's pairs in K0b's `sim.bam` under a second name, with the same tile and x/y. Picard called that pair an **optical** duplicate with the new names (1), and could not with the old (0). So Picard reads the new names' positions.
- **FastQC: passes, and the control works.** With the old names, FastQC 0.12.1 drops the per-tile quality module for the whole file. With the new names it keeps it.

  The first K1 run reported FastQC as failing in all four cells. That was a bug in `k1.py`, which looked in `R1.fq_fastqc/` instead of FastQC's `R1_fastqc/`, so it never found the report. The rule was not changed; with the path fixed, K1 was rerun in full. The first log is kept as `k1-pathbug.log`.
- **The mark: passes.** The number of distinct `SPIKE_` names in new `sim.bam` equals the number of `ev` names in old: 146,837 in K0a and 457 in K0b. All other names are the same set: 1,520,700 and 6,736.

**Left as it is.** `spike validate`'s expected column still reads `>=1 ev0001 read carrying`. That is only the row's label: the prefix it matches is now `SPIKE_ev0001_hap_`, and the rows found the reads.

**Not run:** the full raredisease pipeline. Also how often spike's reads land in a duplicate set in a merged whole-genome BAM, which is when Picard reads their names.

## K1b: Picard, on fresh data (locked before running, 2026-10-04)

This settles K1-Picard, which was inconclusive. The user asked for it ("1": firm up Picard first). Nothing below has been run on the K1b inputs.

**The real warning.** Copied from a run, as the case file now requires: Picard 3.3.0 MarkDuplicates on K0b's old `sim.bam` (`scratchpad/names/k1/k0b/old/picard.log`) printed:

```
WARNING ... AbstractOpticalDuplicateFinderCommandLineProgram  A field field parsed out of a read name was expected to contain an integer and did not. Read name: ev0001_hap_000178. ...
```

No new log so far has a line containing `read name` (case-insensitive).

**Inputs.**
- Binaries: old `918e92fc` and new `b558405f`.
- Picard jar: `63ed3f5d`, with defaults.
- Seed **2**; it was 1 before.
- Events, all new: `del:chr20:11500000-11500300`, `ins:chr20:11600000:100` and `dup:chr20:11700000-11700200`, each `af=0.5`, under the default `--dup-model`, with `--align`. If spike refuses one (RF8), all three move +100,000 bp, the same for both binaries, at most 3 times.
- Run on two BAMs:
  - **H:** the hospital HG002 30x (NovaSeq X);
  - **N:** `data/validation/hg002_novaseq_chr20.bam` (NovaSeq 6000).

**Checks, on each of H and N:**
1. **The warning.** Picard on `out/sim.bam`.
   - **Pass:** new exits 0, and no log line contains `read name` (case-insensitive).
   - **Control:** old's log has a line containing `expected to contain an integer`.
2. **The position is read.** Take the first of spike's pairs in `sim.bam`, in file order, with both mates present as primary, proper-pair records. Add a copy of it under a second name that differs only in spike's part (an `x` after the internal name, so the tile and x/y stay the same). Then sort, and run Picard on that file and on `sim.bam` alone.
   - **Pass:** for new, `READ_PAIR_OPTICAL_DUPLICATES` rises by exactly 1.
   - **Control:** for old, it rises by 0.
   - **Sanity check for both:** `READ_PAIR_DUPLICATES` rises by exactly 1.

**Outcomes:**
- **Supported:** both checks pass on H and N, with every control firing and the sanity check holding.
- **Refuted:** new warns, or new's optical count does not rise by 1, while the controls fire.
- **Inconclusive:** any control or sanity check fails.

### K1b result (2026-10-04): supported

Both inputs ran with no shift. Every spike run exited 0, and so did every Picard run. Log: `scratchpad/names/k1b.log`.

| | H old | H new | N old | N new |
|---|---|---|---|---|
| `sim.bam` alone: log lines containing `read name` | 1 (`expected to contain an integer`) | **0** | 1 (`expected to contain an integer`) | **0** |
| with the copied pair: `READ_PAIR_DUPLICATES` | 0 -> 1 | 0 -> 1 | 0 -> 1 | 0 -> 1 |
| with the copied pair: `READ_PAIR_OPTICAL_DUPLICATES` | 0 -> **0** | 0 -> **1** | 0 -> **0** | 0 -> **1** |

The copied pairs were:
- H: `ev0001_hap_000175` (old), and `SPIKE_ev0001_hap_000175:45:227NC2LT1:2:1265:9538:12252` copied as `SPIKE_ev0001_hap_000175x:...` (new);
- N: `ev0001_hap_000044` (old), and `SPIKE_ev0001_hap_000044:46:HV3C3DSXX:2:2219:27166:13336` (new).

On both inputs the old names make Picard warn, and Picard cannot place an old-named read on the flowcell. With the new names Picard prints no name warning, reads the tile and x/y, and finds the optical duplicate. Every control fired, and the sanity check held.

**K1-Picard: supported.** With K0, K1 (FastQC and the mark) and K2, the read-names change is **supported**.
