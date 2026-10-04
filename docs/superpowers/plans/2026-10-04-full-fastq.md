# Full FASTQ: spike's reads into the sample's raw FASTQ

**Asked 2026-10-04.** The user will spike variants into demo data at the hospital and run the raredisease pipeline on the result, which starts from raw FASTQ. The user picked option 1: build the full FASTQ from the **raw** FASTQ.
- Drop only the originals spike removed.
- Keep every other raw read exactly as it was.
- Add spike's own reads.

The other route was rejected. It builds the FASTQ from `merged.bam` with `samtools collate | fastq`, which would hand raredisease reads that fastp has already trimmed, corrected and filtered.

## Measured before this plan

- **The raw FASTQ is not on this machine.** It is at the hospital, in the FASTQ directory the cohort's `generate3.py` names.
  - Its FastQC report (`fastqc/*/..._1_fastqc.zip`) shows 380,000,000 reads per mate, with lengths 35-151. So the sequencer's adapter trimming is already in it.
  - The sample is a single FASTQ pair (`LNUMBER1`).
  - raredisease aligns it with `bwa-mem2 mem -M -K 100000000` (the BAM's `@PG`).
- **What spike writes**, from K0b of the read-names plan (`scratchpad/names/k0b/new/out`):
  - `replaced_reads.txt` lists all 7,171 originals spike extracted.
  - `R1.fq.gz` holds 7,193 pairs: 6,736 of those originals, which spike kept unchanged, and 457 of spike's own pairs, named `SPIKE_...`.
  - So 435 originals were removed and not kept. `merge.sh` removes all 7,171 from the BAM and adds back `sim.bam` (the 6,736 plus the 457).
- **What a BAM name is.** bwa-mem2 (without `-C`) names a record by the FASTQ header's first word, with a trailing `/1` or `/2` dropped. So a BAM name equals the raw header's first word, minus `@` and that suffix.
- **Every read spike makes starts with `SPIKE_`, and no other read does.** This was measured in read-names K1.
- **Tools:** `pigz` is present here, and `gzip` is the fallback.

## Design (locked)

**The new list.** spike also writes `fastq_removed_reads.txt`: the names in `replaced_reads.txt` minus the names of the pairs spike writes to `R1.fq.gz`/`R2.fq.gz`. These are the originals it removed and did not keep: the suppressed pairs, the pairs dropped for unusable quality, and under `--edit-model origin` the reads `origin` removed. Sorted, one per line.

**The new script.** spike also writes `fastq.sh`:

```
bash fastq.sh RAW_R1 RAW_R2 OUT_R1 OUT_R2 [THREADS]
```

For each mate, the output is two parts, in this order:

1. **Every raw record whose name is not in `fastq_removed_reads.txt`.** The name is the header's first word, with `@` and a trailing `/1` or `/2` removed. Records are kept byte for byte, in order.
2. **spike's own records** (names starting `SPIKE_`) from `R<mate>.fq.gz`, in order. Each header becomes `@<name><style>`, where `<style>` copies the raw file's first header:
   - its text from the first space on (for example ` 1:N:0:ACGTACGT+TGCATGCA`);
   - else `/<mate>`, if that header ends in `/1` or `/2`;
   - else nothing.

**Compression.** The output is gzip, written with `pigz -p THREADS` when `pigz` is on the PATH and with `gzip` otherwise. Input is read the same way.

**Guards.** The script exits 1, deletes both outputs and says why:
- if either mate drops a different number of records than `fastq_removed_reads.txt` has names;
- if R1 and R2 add different numbers of spike's reads.

The message: `RAW_R1/RAW_R2 must be the sample's full raw FASTQ (every lane, concatenated), from the run whose BAM spike was given`.

**Paths.** All paths reach `awk` through the environment, never through `-v` or the program text, so a path holding a backslash, quote, `$` or backtick works. The THREADS default is spike's `--threads`.

**What the user sees.** README: an output-file row for each new file, and a section on building the full FASTQ. The "Next steps" log line names `fastq.sh`.

**Tests, written first and seen red:**
- the removed list equals replaced minus written;
- `fastq.sh` on small gzip FASTQs:
  - with comment-style headers: listed names are dropped, other records are byte-identical and in order, spike's records are appended with the raw comment (R1 ` 1:...`, R2 ` 2:...`), and the originals spike kept are not appended a second time;
  - with `/1`-style headers: appended names end in `/1` and `/2`;
  - with a listed name missing from the raw FASTQ: the guard fires, exit 1, no outputs;
  - with `HOSTILE_NAME` in the output directory and output paths: it runs, no command is injected;
- `main` writes `fastq_removed_reads.txt` and `fastq.sh`.

**Mutation checks (each must turn a test red):**
1. the removed list is all of `replaced_reads.txt`;
2. the removed list is empty;
3. the `/1` suffix is not stripped before matching;
4. the header's whole line is matched instead of its first word;
5. the originals spike kept are appended too;
6. the raw comment is not copied;
7. the `/1` style is dropped;
8. the per-mate count guard is removed;
9. the R1/R2 count guard is removed;
10. outputs are kept on failure;
11. a path is passed with `awk -v`.

## Checks (locked before running)

**The input.**
- **The stand-in raw FASTQ:** every primary record (`-F 0x900`, duplicates included, as raw FASTQ has them) of the hospital HG002 30x BAM over chr20:9,000,000-13,000,000, through `samtools collate | samtools fastq -n` (singletons dropped). Every header is given `" 1:N:0:ACGTACGT+TGCATGCA"` (R1) or `" 2:N:0:ACGTACGT+TGCATGCA"` (R2), the bcl-convert style, and the files are gzipped.
- **The spike run:** the new binary on the full hospital BAM, `--seed 1 --threads 16 --align`, with these events: `del:chr20:10500000-10500300;af=0.5`, `ins:chr20:10600000:100;af=0.5` and `dup:chr20:10700000-10700200;af=0.5` (default `--dup-model`).
- **Fallbacks:**
  - If the guard fires because a removed pair lies outside the window, the window becomes chr20:8,000,000-14,000,000, once.
  - If spike refuses an event (RF8), all three move +100,000 bp, at most 3 times.

**F1: the FASTQ is exactly what was meant.** `fastq.sh` exits 0, and an independent Python checker finds, for each mate:
- the output equals: the stand-in's records whose names are not in `fastq_removed_reads.txt` (byte-identical, in order), then spike's `SPIKE_` records (sequence and quality identical, header `@<name> <mate>:N:0:ACGTACGT+TGCATGCA`);
- each mate dropped exactly as many records as the list has names;
- R1 and R2 have the same names in the same order;
- no name appears twice;
- the list equals `replaced_reads.txt` minus the names in spike's `R1.fq.gz`, computed with `comm`.

**F1 controls.**
- The checker, fed the stand-in plus spike's reads with nothing dropped, must fail.
- `fastq.sh` on a stand-in from chr20:20,000,000-21,000,000, which lacks the event's reads, must exit non-zero and leave no output.

**F2: raredisease's aligner sees spike's reads as `sim.bam` does.** Align the F1 output with `bwa-mem2 mem -M -K 100000000 -t 16 -R '@RG\tID:sim\tPL:ILLUMINA\tSM:sim'`, then `samtools sort`. Then:
- **Pass:** at least 99% of spike's primary records (names starting `SPIKE_`) have the same contig, position, CIGAR and strand as in `sim.bam`. MAPQ agreement is reported, not judged.
- **Pass:** `spike validate --bam <that BAM> --truth truth.vcf` gives `del_planted` and `ins_planted` PASS. Every other row is reported, not judged.

**F3: speed.** `fastq.sh` is timed on the stand-in. Records per second are reported, with a projection for 2 x 380,000,000 labelled `PREDICTED (not run)`. Not judged.

**Outcomes:**
- **Supported:** F1 and F2 pass and both F1 controls behave. Merge on the user's word.
- **F1 fails:** the script is wrong. Fix it with a test first.
- **F2 fails:** look at why before anything else, for example bwa's per-batch insert-size estimate with spike's reads all at the end.

**Not run here:**
- raredisease itself, and the real raw FASTQ (both at the hospital);
- the gzip fallback (`pigz` is installed here);
- samples with several FASTQ pairs (the user concatenates them first, as the guard message says).

## Result (2026-10-04): supported

**Setup.**
- Binary `3b36bfbe` (code `4033041`).
- 659 tests. 11 of the 11 planned mutants are caught. A 12th, `main` not calling `write_fastq_route`, survives the unit tests and is caught here: without it the run has no `fastq.sh`.
- Logs: `scratchpad/fastq/{f_run,f1,f1_control,f2}.log`.

**The input.**
- **The stand-in raw FASTQ:** 478,815 pairs over chr20:9-13 Mb. samtools dropped 7,828 singletons.
- **The control stand-in:** 126,141 pairs over chr20:20-21 Mb.
- **The spike run:** exit 0, with no shift.
  - `replaced_reads.txt` holds 7,171 pairs, of which spike kept 6,517 in R1/R2.
  - `fastq_removed_reads.txt` holds 654 pairs.
  - spike's own reads: 648 pairs.

**F1: passes.**
- `fastq.sh` exited 0 and printed `removed 654 original pairs, added 648 of spike's pairs`.
- The checker, for each mate:
  - 478,809 records, as wanted (478,815 - 654 + 648);
  - byte-identical to the raw records minus the listed ones, in order, then spike's 648 with the header `@<name> <mate>:N:0:ACGTACGT+TGCATGCA`;
  - 654 of the 654 listed names dropped;
  - every name unique.
- R1 and R2 hold the same names in the same order.
- The list equals `replaced_reads.txt` minus the names in spike's R1, computed with `comm`.

**F1 controls: both behave.**
- The checker, fed the stand-in plus spike's reads with nothing dropped, says FAIL: 479,463 records, first difference at record 1,334.
- `fastq.sh` on the control stand-in exited 1 with `R1: found 0 of the 654 originals listed`, and left no output.

**F2: passes.** The full FASTQ was aligned with raredisease's `bwa-mem2 mem -M -K 100000000` and sorted. Then:
- All **1,296 of 1,296** of spike's primary records have the same contig, position, CIGAR and strand as in `sim.bam`, and the same MAPQ.
- `spike validate` on that BAM: `del_planted` PASS (8 reads) and `ins_planted` PASS (20 reads).
- Reported, not judged: 8 rows pass on both the full-FASTQ BAM and `sim.bam`. The only failing row in both is the global `dup_rate` (`no dup flags`), because neither BAM went through duplicate marking. All 13 advisory rows pass.

**F3: speed** (reported). 478,815 pairs took 2.46 s wall, both mates, `pigz -p 8`, with a 4 MB maximum resident size. That is about 194,600 pairs/s.

`PREDICTED (not run)`: the hospital sample's 380,000,000 pairs would take about 33 minutes at that rate. Real files are read from disk, not page cache, so it may be slower.

**Not run here:**
- raredisease itself, and the real raw FASTQ;
- the `gzip` fallback;
- a real bcl-convert header. The stand-in's ` 1:N:0:ACGTACGT+TGCATGCA` was written by `f_run.sh`, in the bcl-convert form.

## `--raw-fastq`: spike runs fastq.sh itself (locked before running, 2026-10-04)

**Asked.** The user picked "1" on 2026-10-04: one command instead of two, the way `--align` runs `align.sh`.

**Design.**
- **The option:** `--raw-fastq RAW_R1 RAW_R2`, exactly two values.
- **Up front:** spike checks that both paths are files before it reads the BAM, so a typo costs seconds, not a run.
- **At the end:** after `--align` if that is given too, spike runs the run's own `fastq.sh` with `RAW_R1 RAW_R2 <output>/spiked_R1.fastq.gz <output>/spiked_R2.fastq.gz <threads>`.
  - If the script exits non-zero, spike exits non-zero and names its exit code. The script's own message reaches stderr, and it leaves no output.
- **Docs:** the README and the run README say so.

**Tests, written first and seen red:**
- the option parses two values, and refuses one;
- the up-front check refuses a missing file and accepts two present ones;
- `run_fastq` on the fastq.sh test files writes `spiked_R*.fastq.gz` equal to `fastq_expected`;
- `run_fastq` returns an error when fastq.sh refuses.

**Mutation checks (each must turn a test red):**
1. the up-front check passes a missing file;
2. `run_fastq` ignores the script's exit status;
3. `run_fastq` passes R1 as both raw files;
4. the clap option takes one value.

**Check G (on the F1 stand-in and spike run).** The new binary is run with the F1 run's exact command, plus `--raw-fastq <stand-in raw R1> <stand-in raw R2>`.
- **Pass:**
  - spike exits 0;
  - its `spiked_R1/R2.fastq.gz`, decompressed, equal F1's verified `standin/main/spiked_R1/R2.fastq.gz`, decompressed, byte for byte.
- **Controls:**
  - with the F1 control stand-in (chr20:20-21 Mb), spike exits non-zero, the log holds fastq.sh's `found 0 of the`, and no `spiked_R*.fastq.gz` is left;
  - with a raw path that does not exist, spike exits non-zero before any `Extracting read pairs` log line.

### `--raw-fastq` result (2026-10-04): supported

**Setup.**
- Binary `8ae6e449` (code `58a1ed8`).
- 663 tests; clippy has master's 14 warnings.
- The 4 planned mutants are caught. Two more mutants, where `main` skips the up-front check or never runs `fastq.sh`, survive the unit tests and are caught below.
- Log: `scratchpad/rawfq/g.log`.

**Check G.**
- **Main run:** the F1 run's command plus `--raw-fastq` on the F1 stand-in. spike exits 0, and `spiked_R1.fastq.gz` and `spiked_R2.fastq.gz`, decompressed, are byte-identical to F1's verified output.
- **Control stand-in (chr20:20-21 Mb):** spike exits 1, the log holds fastq.sh's `found 0 of the ...`, and no `spiked_R*.fastq.gz` is left.
- **A raw R2 that does not exist:** spike exits 1 with `--raw-fastq ... is not a file`, and there is no `Extracting read pairs` line in the log.

## `--raw-fastq` becomes `--into-fastq`; which FASTQ to use is said plainly (plan, 2026-10-04)

**Asked.** The user asked whether it is obvious that `--raw-fastq` makes spiked FASTQ files. It is not:
- The name says only what goes in.
- What comes out is in the help text's second sentence and in the run's last log line.
- An output directory holds two FASTQ pairs:
  - `R1.fq.gz`/`R2.fq.gz`, spike's reads around the events (100,513 pairs in the duplicates after-check);
  - `spiked_R1.fastq.gz`/`spiked_R2.fastq.gz`, the whole sample (478,809 pairs on the F1 stand-in).

  The README's file table calls the first pair "Forward reads" and "Reverse reads", and the run's summary prints `FASTQ: out/R1.fq.gz and out/R2.fq.gz`. Nothing says which pair a pipeline needs.

The user picked option 1: rename, and say plainly which FASTQ is which.

**Design (locked).**
- **The option.** `--into-fastq RAW_R1 RAW_R2`, read as "spike the events into these FASTQ files".
  - The help's first sentence says what comes out: the whole sample, spiked, in `<output>/spiked_R1.fastq.gz` and `spiked_R2.fastq.gz`, the pair to give a pipeline that starts from FASTQ.
  - `--raw-fastq` is no longer accepted. It was public for a few hours, and an alias would keep the unclear name alive.
  - Nothing else about the option changes: the up-front file check (its message names `--into-fastq`), `fastq.sh`, and the output paths.
- **What the run says.**
  - The summary line for R1/R2 says they are spike's reads around the events, not the whole sample.
  - With `--into-fastq`, the last line names both `spiked_` files and says they are the whole sample, spiked.
  - Without it, "Next steps" names `--into-fastq` beside `fastq.sh`.
- **The run's README.**
  - The R1/R2 row says the same as the summary line.
  - A new row for `spiked_R1.fastq.gz`/`spiked_R2.fastq.gz` says they are the whole sample, spiked: written at the end of this run when `--into-fastq` was given; otherwise how to make them.
  - The workflow block names `--into-fastq`.
- **The README.** The same wording in the R1/R2 rows and in the Full FASTQ section. The `--help` copy is regenerated from the binary. `scripts/duplicates/run.sh` uses the new name. Done plans keep the old name, since that is what they ran.

**Tests, written first and seen red:**
- `--into-fastq a b` parses into the pair;
- `--raw-fastq a b` is refused;
- the up-front check names `--into-fastq`;
- the run README's R1/R2 row says "not the whole sample", and its `spiked_` row differs with and without `--into-fastq`;
- the two log messages, built by functions, say what is planned above.

**Check H (locked).** The duplicates after-check's spike step (same 43 events, seed and stand-in, `scratchpad/dups/run-fix/spike`) is run with the new binary and `--into-fastq` in place of `--raw-fastq`.
- **H1:** spike exits 0. `spiked_R1.fastq.gz` and `spiked_R2.fastq.gz` decompress byte-identical to the after-check's. `R1.fq.gz`, `R2.fq.gz` (decompressed), `replaced_reads.txt`, `fastq_removed_reads.txt` and `truth.vcf` are identical too.
  - **Control:** the same comparison of `spiked_R1.fastq.gz` against the raw stand-in R1 must say "differs".
- **H2:** the same command with `--raw-fastq` exits non-zero before any work, and no output directory is made. clap's message is recorded.
- **H3:**
  - the log's last line names both `spiked_` files and the words "whole sample";
  - the run README holds the new rows;
  - the README's `--help` copy equals `spike --help`.

**Outcomes.**
- **Supported:** H1, H2 and H3 pass. Merge on the user's word.
- **Otherwise:** fix it, with a test first.

### Addendum, before any code: `--fastq-prefix` (locked)

**Asked.** The user, mid-plan: the output FASTQ files must take a prefix, because they will run many samples and need to tell the output files apart.

**Design.**
- **The option.** `--fastq-prefix NAME`, default `spiked`. The whole-sample pair is `<output>/NAME_R1.fastq.gz` and `<output>/NAME_R2.fastq.gz`, so the default keeps today's names.
- **The folder** is `-o`'s. NAME is a plain file name: empty, `.`, `..` and anything holding `/` are refused before any work, with a message naming `--fastq-prefix`.
- **It needs `--into-fastq`.** Given alone, it is refused before any work: there would be nothing to name.
- **Everywhere the plan above names `spiked_R1.fastq.gz`, the run uses the prefixed names:**
  - the run README's row and workflow line;
  - the last log line;
  - `--into-fastq`'s help (`<output>/<--fastq-prefix>_R1.fastq.gz`).
- `R1.fq.gz`/`R2.fq.gz` and every other output keep their names. `-o` already separates runs.

**Tests, written first and seen red:**
- the default is `spiked`;
- `--fastq-prefix S1` names `S1_R1.fastq.gz`/`S1_R2.fastq.gz` in the output folder;
- `a/b`, `..`, `.` and an empty name are refused;
- `--fastq-prefix` without `--into-fastq` is refused;
- the run README and the last log line carry the prefixed names.

**Check H4 (locked), added to H.** H1's command plus `--fastq-prefix S1`, into a fresh output folder.
- spike exits 0.
- `S1_R1.fastq.gz` and `S1_R2.fastq.gz` decompress byte-identical to the after-check's `spiked_` pair.
- No `spiked_` file is written.
- The last log line and the run README name the `S1_` files.
- `--fastq-prefix a/b` and `--fastq-prefix S1` without `--into-fastq` each exit non-zero before any work.
