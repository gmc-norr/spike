# `--indel-error-rate` must be a probability

**Asked 2026-10-04.** This is review finding 8 (P2, `docs/review/2026-10-04-independent-review.md`), the user's third pick. We reproduced it on master `2764dbc`: `NaN`, `2.0` and `-0.5` all exited 0 and wrote truth and output.

## Read in the code

- `synth.rs:565`: after a sequencing error is drawn, it is an indel when `rng.gen::<f64>() < indel_error_rate`. So the value is a probability: the share of sequencing errors that are indels.
  - NaN and negative values turn indels off without a word.
  - A value above 1 behaves like 1.
- The README already says "the fraction of sequencing errors that are indels". The `--help` text says "Indel error rate per base", then "(fraction of total error ...)". It contradicts itself.

## Design (locked)

- `validate_indel_error_rate` refuses a value that is not finite or lies outside `[0.0, 1.0]`. The message is: `--indel-error-rate must be in [0.0, 1.0]: it is the fraction of sequencing errors that are indels; got <value>`.
- It runs with the other input checks, before the output folder is made. A refused run therefore leaves nothing behind.
- The `--help` text becomes: "Fraction of the sequencing errors in synthetic reads that are indels rather than substitutions, in [0, 1]. Default 0.0: substitutions only. Typical Illumina: 0.0 to 0.05."
- The README gets the regenerated `--help` copy and a row in "What spike refuses".

**Tests, written first and seen red:**
- 0, 0.05 and 1 are accepted;
- NaN, inf, -inf, -0.5, 2.0 and 1.0000001 are refused, with `--indel-error-rate` in the message.

**Mutation checks** (each must turn the test red; the unmutated tests are green first):
1. written as `rate < 0.0 || rate > 1.0`, which lets NaN through;
2. the upper bound dropped;
3. the lower bound dropped;
4. `> 1.0` made `>= 1.0`, which refuses 1.

## Checks (locked)

**I1.** The reviewer's script (the scratch copy from the carried-allele result) against the new debug binary: `indel_rate_NaN`, `indel_rate_2.0` and `indel_rate_-0.5` each exit non-zero without writing `truth.vcf`. Every other result is as before.

**I2.** On the hospital BAM, the new and master binaries give identical `truth.vcf`, `R1.fq.gz`, `R2.fq.gz`, `replaced_reads.txt` and `fastq_removed_reads.txt` for F1's three structural events with `--indel-error-rate 0.05 --seed 1 --threads 16`.

**I3.** `spike --indel-error-rate nan` and `2` (with valid other arguments) exit 1 within 5 s, print the message, and do not create the `-o` folder.

**I4.** The README's `--help` copy equals `spike --help`.
