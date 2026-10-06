# align.sh aligns spike's reads with the sample's own bwa-mem2 options

**Asked 2026-10-06.** The user picked "fix 1 and 2" after a read-level look at 20 transplanted GIAB deletions. That look is in the session scratchpad, `inspect/`. This plan is fix 2. Fix 1 (low-quality 3' tails) gets its own plan.

**Found.** spike's default `align.sh` runs `bwa-mem2 mem -t N -R <spike's RG> REF R1 R2`. The hospital BAM was aligned with `bwa-mem2 mem -M -K 100000000 ...`, as its `@PG` line says. Without `-M`, bwa-mem2 marks the shorter part of a split read SUPPLEMENTARY. With `-M`, it marks it SECONDARY. So in a BAM built with `merge.sh`, spike's split reads carry a flag the sample's own reads never carry.

Measured on the hospital BAM, LDLR `del:chr19:11106493-11114129;af=het`, seed 1, master `bdd568a`. The two runs have the same reads and the same truth record, and differ only in the aligner command:

| aligner command | records | SUPPLEMENTARY | SECONDARY | SECONDARY with SA |
|---|---|---|---|---|
| master's default (`slice-safety` run, `data/B2/spike/sim.bam`) | 6,144 | 22 | 0 | 0 |
| `--aligner "bwa-mem2 mem -M -K 100000000 -t 4 -R ..."` (`inspect/ld2/run/sim.bam`) | 6,144 | 0 | 22 | 22 |

This touches only the BAM route (`align.sh`, then `merge.sh`), which is how the preliminary caller tests run. The FASTQ route is not touched, because the pipeline realigns every read itself.

## Gate A

1. **Principle.** spike's reads go through the same aligner command as the sample's own reads, so their flags and placements are made the same way.
2. **What would kill it.** On the hospital BAM, the new default `align.sh` still writes a SUPPLEMENTARY record. Or its `sim.bam` differs from the hospital command's `sim.bam` in anything other than the read group. Or a BAM whose `@PG` holds none of these options gets a different `align.sh` than master's.
3. **Refuted before?** No. `git log` has quoting fixes to `align.sh` (L12) and the read group's SM, but nothing that copies aligner options.
4. **Simplest thing.**
   - Read the input's `@PG` line for bwa-mem2 (or bwa), and copy the options that shape alignment onto the `bwa-mem2` preset's command line.
   - Always adding `-M` would be simpler, but wrong for a sample aligned without it. The 35x GIAB HG002 BAM is one: its `@PG` has no `-M`.
5. **Inputs, and how each was checked.**
   - The hospital BAM's `@PG ID:bwa-mem2 ... CL:bwa-mem2 mem -M -K 100000000 -R @RG\t...\tSM:D24-14230_Seq25-7600_30x -t 12 ./bwamem2Index/genome.fa ...`: read with `samtools view -H`. The `-R` value has no space in it.
   - The 35x HG002 BAM's `@PG ID:bwa-mem2 ... CL:<path>/bwa-mem2 mem -t 36 -R @RG\t...\tSM:HG002 <ref> <fastqs>`: read the same way. Note the program is given as a path.
   - bwa-mem2 2.2.1's options and which take a value: read from `bwa-mem2 mem` with no arguments.

## Design (locked)

**Which `@PG`.** The first `@PG` whose `PN` or `ID` is `bwa-mem2` or `bwa`, or starts with one of them followed by `.` (samtools numbers repeats that way). Its `CL` is split on whitespace. The token whose file name (after the last `/`) is `bwa-mem2` or `bwa` must be followed by `mem`. The options come after that.

**Which options.**
- Copied, with their value: `-k -w -d -r -y -c -D -W -m -A -B -O -E -L -U -T -h -I -K -x`.
- Copied, without a value: `-S -P -j -5 -q -a -V -Y -M`.
- Not copied:
  - `-t` (spike passes its own threads);
  - `-R` (spike's own read group);
  - `-o` and `-v` (output and log only);
  - `-H` (adds header lines);
  - `-p` (smart pairing would ignore spike's R2);
  - `-C` (it would turn spike's FASTQ comments into SAM tags).
- Any other token that starts with `-` is not copied, and spike warns, naming it.
- Tokens that do not start with `-` (the index, the FASTQs) are not options and are skipped.

**Where it goes.**
- Only the `bwa-mem2` preset's command line gets the copied options: `bwa-mem2 mem <copied> -t "$THREADS" -R ...`.
- `minimap2`, `bowtie2` and a custom `--aligner` are unchanged.
- Each copied token is shell-quoted, because header text is not trusted. A plain word such as `-M` or `100000000` stays as it is.
- When options are copied, spike logs them at info level, and `align.sh` gets a comment line naming them and the `@PG` they came from.
- With no bwa `@PG`, or none of these options, `align.sh` is byte-identical to master's.

**Tests, written first and seen red:**
- The hospital `CL` gives `-M -K 100000000`.
- The 35x BAM's `CL`, with a path to the program, gives nothing.
- `-t`, `-R` (with its value), `-o`, `-v`, `-H`, `-p` and `-C` are not copied, and a value after a value-taking option is not read as an option.
- An unknown option is reported and not copied.
- A `bwa mem` (not mem2) `CL` is read too.
- A header with no bwa `@PG` gives nothing.
- `align.sh` for the `bwa-mem2` preset holds the copied options, quoted. The other presets do not.
- A value with shell metacharacters reaches bwa-mem2 as one word.

**Mutation checks.** Each must turn a test red, and the unmutated suite is green first:
1. `-t` is copied.
2. A value-taking option loses its value.
3. The `@PG` is ignored (nothing is ever copied).
4. `-p` is copied.
5. Copied tokens are not quoted.

## Checks (locked)

**K1, the hospital BAM.** `del:chr19:11106493-11114129;af=het`, `--seed 1`, default aligner, new binary.
- `align.sh`'s bwa-mem2 line holds `-M -K 100000000`.
- `sim.bam` has 0 SUPPLEMENTARY records, and at least 1 SECONDARY record, every one carrying SA.

**K2, the same command as the hospital's.** K1's `sim.bam` equals `inspect/ld2/run/sim.bam` record for record, once the `RG` tag is removed. That file came from the custom `--aligner "bwa-mem2 mem -M -K 100000000 -t 4 -R ..."` on the same event and seed. `-K` makes bwa-mem2's output independent of its thread count.

**K3, a BAM without these options.** The 35x HG002 BAM, `del:chr20:14530000-14531000`, seed 1. The new binary's `align.sh` is byte-identical to master's.

**B, the reads do not change.** In K1 and K3, `R1.fq.gz`, `R2.fq.gz`, `truth.vcf` and `replaced_reads.txt` are byte-identical to master's. `truth.vcf` is compared without its `##reference` temp-path line, if that differs.
