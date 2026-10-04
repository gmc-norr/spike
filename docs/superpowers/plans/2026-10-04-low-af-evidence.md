# truth.vcf says how many planted pairs show each event

**Asked 2026-10-04.** This is the review's low-count note (`docs/review/2026-10-04-independent-review.md:127`). The user picked "fix the low-AF note" after finding 2 was merged. The review reads: "Record informative fragment counts and distinguish expected from observed AF; do not condition on positive evidence silently, since that would bias sensitivity estimates."

**Reproduced.** The reviewer's case is SNP `af=0.001`. Spike emitted 2 pairs, recorded `SIM_VAF=0.014`, and neither pair carried the SNP.

**Measured with master `2ab9fc7`'s release binary** (`920d2b7c`), using `scratchpad/ocp/lowaf/probe.py`. These are clean-mode runs, and finding 2's merge (`6c39e54`) left clean output byte-identical. It used the hospital BAM, clean mode and `--seed 1 --threads 16`, at the first 20 sites of the origin-conservation R draw that ran. It counts the spike reads holding the ALT 21-mer (10 bases each side):

| AF | tiled pairs (mean) | spike reads with the ALT (mean) | sites with none |
|---|---|---|---|
| 0.001 | 2.0 | 0.20 | 16 of 20 |
| 0.01 | 4.3 | 0.25 | 15 of 20 |
| 0.02 | 8.8 | 0.40 | 13 of 20 |
| 0.05 | 22.1 | 2.20 | 2 of 20 |
| 0.1 | 44.0 | 4.00 | 0 of 20 |

So this is not only the two-pair floor. Below about 5% at this depth, most events have no planted read that shows them. Every one is in `truth.vcf`, and a caller scored against it is charged with a miss it could not make.

## Read in the code

- **`SIM_VAF`** is `compute_tiling_count`'s fraction (`src/simulate.rs:723-805`): the request, or the cap's or floor's realized fraction. It comes from the count of pairs, not from where they land. Its header says "the fraction of the depth the fragments spike planted make up" (`src/truth.rs:134-140`).
- **The floor's warning** says spike emits "the 2 it needs to plant the event at all" (`src/simulate.rs:842-853`). Two pairs drawn uniformly over a ~4 kb haplotype usually miss a 1 bp event.
- **The pair's ends.** `tile_haplotype_reads` (`:944-1067`) draws each pair's haplotype start and fragment length. `generate_haplotype_read_pair` (`src/synth.rs:874-951`) cuts the forward mate from the start and the reverse mate from the end, and trims each mate's 3' end separately (`trim_adapter_start`, `:498-504`). So the two mates can differ in length, and which one is R1 is a coin flip.
- **Small variants** are three segments: left flank, the ALT bases, right flank (`src/haplotype.rs:404-464`).
- **Structural haplotypes** join segments. Not every boundary is a junction: a tandem DUP's left flank runs straight into its first copy (`:172-232`).
- **`truth.vcf`'s readers** (`src/validate.rs`, `src/vcf_input.rs`, `scripts/validate_pipeline.sh` and the scoring scripts) all look `SIM_*` fields up by key. A new key changes none of them.

## Gate A

1. **Principle.** `truth.vcf` says what the planted reads show, beside what spike aimed for. An event that no planted pair shows is labelled as such. Spike still never redraws to make it show, because that would bias sensitivity.
2. **What would kill it.** The new count disagrees with an independent count, made from the FASTQ, of the same runs' spike pairs that show the event.
3. **Refuted before?** No. `git log` has `SIM_VAF`'s cap and floor work, but nothing that counts what the reads show.
4. **Simplest thing.** Count at tiling time, from positions spike already knows. An observed allele fraction (ALT over all reads at the site) is not done. It needs the original reads kept there, and those are only settled when `merge.sh` runs. The count is what the review asks for first.
5. **Inputs, and how each was checked.**
   - The mates' haplotype spans are taken from their final lengths: forward `[start, start + len)`, reverse `[end - len, end)`. This was read in `src/synth.rs:874-951`.
   - A read with an indel error covers one base more or less than its length. **Not measured**; the independent count will show where it matters.
   - Both copies share one layout (`src/simulate.rs:1049`), so a pair from the other copy is placed the same way.

## Design (locked)

**A pair shows its event:**
- **For a small variant (`SmallVariant`)**, when one of its reads covers the changed bases and one base on each side. The changed bases are REF and ALT with their common prefix and then suffix trimmed, as in `carried::changed_span`. That is the span rule the carried-allele check uses. In haplotype coordinates, with `L` the left flank's length, the changed ALT bases are `[a, b) = [L + p, L + len(ALT) - q)`. A read `[s, e)` shows them when `s <= a - 1` and `e >= b + 1`. A pure deletion has `a == b`, so the read must cover the bases on both sides of the join.
- **For a structural event**, when its fragment `[s, e)` covers both bases around one of the haplotype's junctions (`s <= j - 1` and `e >= j + 1`). That makes it a split read or a discordant pair. A junction is a segment boundary where the reference does not simply continue: one side is novel sequence, the chromosome or strand changes, or the reference positions do not meet. A junction DUP's interior depth copies are not counted; they carry dosage, not a junction.

**Output.**
- `truth.vcf` gets `SIM_ALT_FRAGS=<n>` on every record, after `SIM_REQ_VAF`. Its header line explains it.
- `SIM_VAF`'s header says it is the fraction expected from the number of pairs, and that `SIM_ALT_FRAGS` says how many show the event.
- When `n` is 0, spike warns that no planted pair shows the event. It names the event, the request and the number of pairs tiled.
- The floor's warning no longer says two pairs plant the event.

Reads, their order and the random stream are unchanged: the count reads positions and draws nothing.

**Tests, written first and seen red:**
- the small-variant rule on an SNV, a pure deletion and an insertion: a read ending on the changed base does not count; a read ending one base past it does;
- the junction rule: a fragment straddling the junction with neither read crossing it counts; a fragment ending at the junction does not;
- junctions of a tandem DUP: only the copy-to-copy boundary; of a deletion: one; of an inversion: two;
- the counted mates: unequal mate lengths are matched to the right end.

**Mutation checks** (each must turn a test red; the unmutated tests are green first):
1. the flank base dropped (`s <= a`, `e >= b`);
2. the structural rule on reads, not the fragment;
3. every segment boundary counted as a junction;
4. every tiled pair counted;
5. the common prefix not trimmed.

## Checks (locked)

**The independent count.** A Python locator, `scripts/low_af/locate.py`, rebuilds the event's haplotype from the reference: left flank, ALT, right flank, with 2,000 bases a side, or the two flanks of a DEL. It places each spike read by 20-mer seeds every 20 bases, on either strand, taking the start most seeds agree on. It then applies the same rule. It does not use spike's code.

**K1: SNVs.** The 100 runs of the table above (20 sites x 5 AFs), with the new release binary:
- `SIM_ALT_FRAGS` equals the locator's count in at least 98 of 100 runs, and differs by at most 1 where it does not;
- every run with `SIM_ALT_FRAGS=0` has 0 from the locator and 0 ALT 21-mers.

**K2: small deletions.** 20 sites from the same draw, each `snp:chr20:POS:R:A` with REF the 5 reference bases from POS and ALT their first base (a 4 bp deletion), at AF 0.02 and 0.05. Equal in at least 39 of 40 runs, and off by at most 1 elsewhere.

**K3: structural deletions.** 20 sites from the same draw, each `del:chr20:POS-(POS+1000)`, at AF 0.02 and 0.05. Equal in at least 39 of 40 runs, and off by at most 1 elsewhere. A run refused by the same rule under master and new is skipped, and listed.

**P: the reviewer's probe.** `low_fraction` gives `SIM_ALT_FRAGS=0`, the warning is logged, and every other probe result is as before.

**B: the reads do not change.** In every K1 run, the new and master binaries give identical `R1.fq.gz`, `R2.fq.gz`, `replaced_reads.txt` and `fastq_removed_reads.txt`. Their `truth.vcf` files are identical once the `SIM_ALT_FRAGS` header line, the `SIM_VAF` header line and the `;SIM_ALT_FRAGS=n` field are removed.

**R (reported).** From K1-K3: the share of runs with `SIM_ALT_FRAGS=0` per AF and event type.

No speed check. The count is one comparison per pair.
