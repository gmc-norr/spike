# Plan: the quality sample stays on the reference's contigs, and LDLR exon 1 is where MANE puts it

Date: 2026-10-07. Base: master `4a5672d`. Two small defects found by the improvement-roadmap workflow (ideas M32 and M29) and reproduced by hand before this plan.

## The defects, as measured

**D1: v2's startup sample aborts on an input aligned to a larger reference.**
- Input: `data/giab_hg38/HG001/NA12878.final.cram` (3,366 `@SQ`, decoys included) with the no-alt FASTA (195 contigs). Event: `del:chr19:11100000-11102000 --threads 8`.
  - Master `4a5672d` exits 1: `Error: invalid reference sequence name: chrUn_JTFH01000277v1_decoy`.
  - The pre-v2 binary (`e21d282b`) exits 0 on the same command.
- **Cause.** `quality::sample_input` places its 20 blocks over every indexed window of every `@SQ` contig (`indexed_windows`). A block on a contig the FASTA lacks fails to read, and any block failure is fatal.
- Regression from v2 (`7f28bb0`). v2's K4 ran only on inputs whose `@SQ` set equals the FASTA's.

**D2: the bundled LDLR exon BED puts exon 1 in intron 1.**
- `data/ldlr_deletions/ldlr_exons_hg38.bed` has been unchanged since `a63115e`. Its exon 1 is `chr19:11090578-11090919`.
- MANE v1.0 `NM_000527.5` exon 1 is `chr19:11089462-11089615` (0-based, half-open). The two do not overlap.
- The CDS start (ATG) is at chr19:11,089,549 (1-based), and the reference reads `ATG` there. It lies inside MANE exon 1.
- Exons 2-17 match MANE exactly. Exon 18 ends at 11133816 against MANE's 11133820.
- Any exon-level event naming LDLR exon 1 therefore lands about 1 kb inside intron 1.

## Gate A

1. **Principle.** spike learns its sample from the reads it can read against the reference it was given. An exon spec means the exon the clinic means: the MANE transcript.
2. **What would kill this.**
   - D1's fix fails if the HG001 CRAM still aborts.
   - It also fails if any input whose contigs are all in the FASTA gets different output, because placement must not move.
   - D2's fix fails if exon 1 does not hold the ATG, or if exons 2-17 change.
3. **Refuted here before?** No. Neither appears in the case file or the plans. These are bugs, not mechanisms.
4. **Simplest thing.**
   - D1: drop the indexed windows whose contig is not in the FASTA's `.fai`, before the blocks are placed. Log how many were dropped. Nothing else changes; block failures stay fatal, since a failure on a contig the FASTA has is a real error.
   - D2: replace exon 1's and exon 18's coordinates with MANE's.
5. **Input assumptions, and how each was checked.**
   - The 35x BAM's 195 `@SQ` are all in the FASTA (`comm` of the header against the `.fai`: 0 missing).
   - The 31-value chr20 BAM lists 2,385 `@SQ` the FASTA lacks, but holds reads on chr20 only. If its index has a chunk on any such contig, its output can move. C2 measures this rather than assumes it.

## Locked checks

- **C1 (D1), deciding.** The HG001 CRAM command above exits 0 with the new binary, and its log names the dropped windows.
- **C2, byte identity.** The master binary (`4a5672d`) against the new one. For each case, the md5 of the decompressed `R1.fq.gz`, `R2.fq.gz`, `truth.vcf`, `replaced_reads.txt` and `fastq_removed_reads.txt`:
  - (a) 35x BAM, `--event dup:chr20:14550000-15550000 --event "snp:chr20:38550000:T:C;af=0.1" --seed 1 --threads 8`;
  - (b) 31-value BAM, `--event "snp:chr20:38550000:T:C;af=0.1" --event del:chr20:38900000-38910000 --threads 8`.

  PASS only if all 10 match and every R1 holds SPIKE_ reads. The new exon BED cannot affect these runs: they name no exons.
- **C3 (D2), deciding.**
  - The bundled BED's exon 1 is `11089462-11089615` and holds 0-based position 11089548.
  - Exon 18 is `11131280-11133820`.
  - Exons 2-17 are byte-identical to before.
  - A unit test reads the bundled file and checks all three.
- **C4, tests.** Every test passes: 750 plus the new ones, nothing removed. Clippy is no worse than on `4a5672d`.
- **C5, mutants.** Each must make a test FAIL:
  - remove the contig filter;
  - put the old exon 1 back.
- **Reported, not gated.** The HG001 CRAM run's startup sample (blocks, pairs) and its time.

## Out of scope

The bundled `ldlr_known_deletions_hg38.vcf`, whose round-number breakpoints are named like published alleles (M29 part 2). MANE GFF input. Decoy contigs present in a FASTA, which are readable and so are kept.

## Result (2026-10-07): every locked check passes

Code: `1103835` (the fix and its two tests) and `9964b38` (an end-to-end test of the wiring). Binaries: master `4a5672d` is `df547601`, and the fix is `e533ec11` (built from `1103835`; `9964b38` adds a test only).

**C1: PASS.** The HG001 CRAM (`NA12878.final.cram`, no-alt FASTA, `del:chr19:11100000-11102000 --threads 8`) exits 0 in 5 s.
- Master exits 1 with `invalid reference sequence name: chrUn_JTFH01000277v1_decoy`.
- The log names the dropped contigs: `Quality sample: 2620 contig(s) with reads are not in the reference FASTA and are left out (e.g. HLA-A*01:01:01:01, HLA-A*01:01:01:02N, HLA-A*01:01:38L)`.
- The sample: 20 blocks from chr1:73995837 to chrX:144483718, 110,789 pairs. Read in 2.5 s, learned in 0.9 s.
- This CRAM has eight quality values (2, 3, 4, 5, 6, 10, 20, 30).

**C2: PASS.** Master's binary against the fix, 10 of 10 files md5-identical:

| case | R1 | R2 | truth.vcf | replaced_reads | fastq_removed_reads | SPIKE_ reads in R1 |
|---|---|---|---|---|---|---|
| (a) 35x dup + snp, seed 1 | `0bcd4052fb` | `27e7b36a4b` | `3c9571b20c` | `36db938278` | `71ff026371` | 120,099 |
| (b) 31-value snp + del | `fd005167f5` | `2a563a9898` | `43ecf37fa7` | `06b4805f71` | `33397ca685` | 325 |

On both, the new "left out" line is absent: no window was dropped, despite the 31-value BAM's 2,385 extra `@SQ`.

**C3: PASS.** `test_the_bundled_ldlr_exons_are_manes` checks all 18 exons against MANE v1.0 `NM_000527.5` and that exon 1 holds the ATG at 0-based 11089548.
- It FAILED on the old file: `left: [(1, 11090578, 11090919), …] right: [(1, 11089462, 11089615), …]`.
- It passes on the new one.
- The file diff is two lines: exon 1, and exon 18's end. Exons 2-17 are byte-identical.

**C4: PASS.**
- 753 passed, 0 failed, 3 ignored. That is 750 plus 3 new tests, with nothing removed.
- Clippy is identical to `4a5672d`: 16 lines, the same warnings.

**C5: PASS.** Each mutant turned a test FAILED:
- **The old exon 1 back:** `test_the_bundled_ldlr_exons_are_manes` FAILED. This is the pre-fix run above, on the same file content.
- **The contig filter removed** (`windows_on_reference` keeps everything): `test_sample_windows_on_contigs_the_reference_lacks_are_left_out` FAILED.
- **The filter left out of `sample_input`'s wiring:** this mutant survived the unit test above, so an end-to-end test was added (`9964b38`). With the mutant it FAILED at `blocks.expect("the startup sample aborted on a contig the FASTA lacks")`, the same abort as the bug. Without the mutant it passes.

**Not done here.**
- The bundled `ldlr_known_deletions_hg38.vcf` still has round-number breakpoints named like published alleles (M29 part 2).
- A FASTA that holds decoys is readable, so decoy blocks are still sampled from it.

## Follow-up after review (2026-10-07)

Two independent reviewers, one on the code and one on the LDLR data, found no critical or important problem.
- The data reviewer rebuilt the BED from the MANE GFF: byte-identical, 18 lines. They checked the ATG and the stop codon (TGA at chr19:11131314-11131316, inside exon 18).

Their minors, fixed in `70e25c9`:
1. **A contig-name mismatch blamed the wrong cause** ("20" against "chr20", reachable with an SV event). Every window was dropped and the sample came back empty, so `from_input` bailed with "check that the file is indexed and that --min-mapq…". Now `sample_input` refuses such an input itself, naming the cause. Measured with a reheadered no-chr BAM against a chr20 FASTA: `Error: none of the contigs … lists (1 of them, e.g. 20) are in the reference FASTA …: the input and the FASTA name their contigs differently (e.g. "20" against "chr20") or come from different builds`. Master's error was `chromosome '20' not found in FASTA index`.
2. **The end-to-end test passed vacuously.** A filter dropping every window passed all 753 tests. The test is now two:
   - a chrA BAM against a chrA FASTA must sample all 20 pairs;
   - a decoy-only BAM must be refused, with the mismatch named.
3. **The log line said "contig(s) with reads"**, which was wrong on the no-index path. It now says "contig(s) of <input>". `sample_input`'s doc now says the windows are filtered.
4. **The exon test also checks** that there is one LDLR, on chr19.

**Mutants, each FAILED a test:**
- drop-all wiring → `test_the_startup_sample_reads_the_contigs_the_fasta_holds`;
- filter not wired → `test_the_startup_sample_names_a_contig_the_fasta_lacks`;
- no mismatch error → the same test;
- LDLR rows on chr1 → `test_the_bundled_ldlr_exons_are_manes`.

**Tests and clippy.** 754 passed, 0 failed, 3 ignored. Clippy is identical to `4a5672d`. A `.err().expect()` first added a warning and was replaced by `expect_err` before the commit was final.

**Re-run with the new binary (`23560464`):**
- C1: the HG001 CRAM exits 0.
- C2: 10 of 10 files are md5-identical to master on the 35x and 31-value runs (120,099 and 325 SPIKE_ reads).

**Left open, both older than this change:**
- A CRAM multi-reference slice can mix FASTA and non-FASTA contigs. On HG001, 3 of 76,890 slices do, and the sample cannot reach them. An event on chrEBV aborts there, as it did before.
- The README's gVCF example `del:chr19:11090000-11133000` misses MANE exon 1. It does not claim to be the whole gene.
