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
