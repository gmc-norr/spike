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
