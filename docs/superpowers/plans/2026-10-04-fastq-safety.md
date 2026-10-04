# fastq.sh: never touch a raw input, and check that the mates pair up

**Asked 2026-10-04.** The independent review (`docs/review/2026-10-04-independent-review.md`) found two faults in the full-FASTQ route. We reproduced both on master `2764dbc`. The user picked "fix the hospital-path ones first"; these two come first because both are in `fastq.sh`, which the hospital run uses.

- **Finding 1 (P1).** If an output path is the same file as a raw input, `fastq.sh` truncates the input (`> OUT`) before it reads it. Its failure trap then deletes both outputs, which is the raw input. Measured: `alias_R1.fq.gz` given as RAW_R1 and OUT_R1 gave exit 1, and the raw file was gone.
- **Finding 7 (P2).** `fastq.sh` reads R1 and R2 one after the other, never together. Measured: RAW_R1 names `[a,b]` and RAW_R2 names `[b,a]` gave exit 0, and the outputs kept the wrong order. A pipeline would then align mispaired reads.

## Measured before this plan

- **bash `-ef`** (bash 5.2.21) says "same file" for: the same name, `./a` and `a`, a hard link, and a symlink. It says "not" when either file is missing.
- **`realpath -m`** (`/usr/bin/realpath`) gives one spelling for a path that does not exist yet (`nope/../x` gives `<cwd>/x`).
- **awk:** the default `awk` is gawk (`/etc/alternatives/awk -> /usr/bin/gawk`). mawk is also here, and the hospital's awk is not known.
- **The old script's time** on the full-FASTQ F1 stand-in (478,809 pairs, 8 threads) was 2.46 s wall (`scratchpad/fastq/standin/main/fastq.log`).
- **The files spike writes into `-o`:** `R1.fq.gz`, `R2.fq.gz`, `truth.vcf`, `events.bed`, `replaced_reads.txt`, `fastq_removed_reads.txt`, `fastq.sh`, `align.sh`, `merge.sh`, `README.md`, and with `--align` `sim.bam` and `sim.bam.bai`. Under `--into-fastq` it also writes `<prefix>_R1.fastq.gz` and `<prefix>_R2.fastq.gz`. Found by grepping `src/main.rs` and `src/fastq.rs` for paths joined to the output folder.

## Design (locked)

### In `fastq.sh`

**1. Refuse a clash before any work.** The script exits 1 and writes nothing if:
- OUT_R1 or OUT_R2 is the same file (`-ef`) as RAW_R1, RAW_R2, or `R1.fq.gz`, `R2.fq.gz` or `fastq_removed_reads.txt` in the script's folder;
- OUT_R1 and OUT_R2 are the same path (`-ef`, or the same `realpath -m`);
- RAW_R1 and RAW_R2 are the same file (`-ef`).

The message names both paths and says they are the same file.

**2. Write beside, then publish.** Each output is first written to a temporary file in the output's own folder (`mktemp` there, a name starting `.`). Only after every check below has passed, for both mates, is each moved onto its output name with `mv -f`.
- On any failure, the trap deletes only those temporary files and the script's own temporary folder.
- A file that already sits at OUT_R1 or OUT_R2 is never deleted or changed by a failed run.

**3. Read both mates together.** One awk pass reads RAW_R1 on its input and RAW_R2 through a named pipe, one line of each at a time. It refuses (exit 1, with the reason) if:
- at some record, the two names differ (name = the header's first word, without `@` and a trailing `/1` or `/2`). The message gives the record number and both names, and says the mates are not in the same order;
- one mate has more records than the other;
- a header does not start with `@`, a third line does not start with `+`, or a sequence and its quality differ in length;
- the file ends inside a record;
- a name in `fastq_removed_reads.txt` is not found (as now), or is found more than once.

The records it keeps go out unchanged and in order, as now, and the header style is read from each mate's first header, as now. awk runs as `LC_ALL=C awk` (this repo's rule).

**4. spike's own reads** are added at the end of each mate's temporary file, as a second gzip member, with the same awk and the same R1/R2 count guard as now. (A file made of several gzip members is ordinary gzip. `cat` of gzip lanes makes the same kind, and the script's own comment already relies on it.)

### In spike (`--into-fastq`)

`validate_into_fastq` also refuses, before the output folder is made and before any work:
- RAW_R1 and RAW_R2 being the same file;
- either being the same file as one that the run writes into `-o` (the list measured above, with the `--fastq-prefix` pair).

"Same file" means the same device and inode, after following symlinks. A path that does not exist yet cannot be the same file as an input that does exist. `--fastq-prefix` is checked first, because the pair's names depend on it.

### Tests, written first and seen red

- `fastq.sh` refuses, exit 1, with the raw file's bytes unchanged:
  - OUT_R1 = RAW_R1 (same name), a hard link to it, and a symlink to it;
  - OUT_R2 = OUT_R1;
  - OUT_R1 = the folder's `R1.fq.gz` (spike's R1.fq.gz unchanged);
  - RAW_R2 = RAW_R1.
- A failed run leaves files already at OUT_R1 and OUT_R2 as they were, and leaves no temporary file in their folder.
- `fastq.sh` refuses, exit 1, no outputs:
  - R2 in a different order (the review's `[a,b]`/`[b,a]`);
  - R2 one record short, and R2 one record longer;
  - the last record cut short;
  - a sequence and quality of different lengths;
  - a header without `@`;
  - a listed name present twice in the raw pair.
- The present `fastq.sh` tests still pass unchanged (same output, both header styles, the hostile path).
- `validate_into_fastq` refuses: RAW_R1 = RAW_R2 (as a name, and as a symlink); RAW_R1 at `<out>/S1_R1.fastq.gz` with prefix S1; RAW_R1 at `<out>/R1.fq.gz`. It accepts two different files.

### Mutation checks (each must turn a named test red; the runner first checks the unmutated suite is green)

1. the output-input clash check is removed;
2. the OUT_R1 = OUT_R2 check is removed;
3. the RAW_R1 = RAW_R2 check is removed;
4. outputs are written in place, not beside;
5. the trap deletes the outputs again (`rm -f "${OUT[@]}"`);
6. the name comparison is removed;
7. the "one mate has more records" check is removed;
8. the sequence/quality length check is removed;
9. the "ends inside a record" check is removed;
10. the "listed name found twice" check is removed;
11. `validate_into_fastq`'s same-file check is removed;
12. `validate_into_fastq` leaves out `R1.fq.gz` from its list.

## Checks (locked before running)

**The input.** The full-FASTQ F1 stand-in (`scratchpad/fastq/standin/main/raw_R{1,2}.fastq.gz`, 478,809 pairs from the hospital HG002 30x BAM over chr20:9-13 Mb). The spike run is F1's (`scratchpad/fastq/f_run.sh`: same BAM, reference, events, `--seed 1 --threads 16`), with the new binary and `--into-fastq` on the stand-in. F1's run used shift 0 (`scratchpad/fastq/run/shift`).

**G1: same output as before.** F1's own output cannot be the reference: the duplicates fix (`2eb9309`) came after F1 and changed the removed list. The reference is master's script instead. `old.sh` is the `fastq.sh` that master `2764dbc`'s binary writes for the same arguments; it is run by hand on the new run's output folder and the stand-in.
- The new run's `spiked_R1.fastq.gz` and `spiked_R2.fastq.gz`, decompressed, equal `old.sh`'s outputs, decompressed, byte for byte.
- The new `fastq.sh`, run by hand on the same input with `awk` pointing to mawk, gives the same bytes again.
- **Control:** the same comparison against the stand-in's raw R1 must say "differs".

**G2: the reviewer's probes.** `scripts/review_20261004.py` run against the new debug binary:
- `input_output_alias`: exit not 0, and `raw_input_still_exists` is true;
- `mismatched_mates`: exit not 0.

The other findings' numbers are reported, not judged; nothing here fixes them.

**G3: refusals on real-size data.** On copies of the stand-in, each must exit not 0, leave no file at the output names and no temporary file in their folder, and leave the raw files' md5 unchanged:
- (a) R2 with records 200,000 and 200,001 swapped;
- (b) R2 without its last record;
- (c) OUT_R1 a hard link to RAW_R1;
- (d) `spike --into-fastq <out>/S1_R1.fastq.gz <out>/S1_R2.fastq.gz --fastq-prefix S1 -o <out>` (raw pair copied there first): refused within 5 s, with "same file" in the message.

**G4: time.** The new script's wall time on the stand-in, 8 threads, best of 3, is at most 1.5 times the old script's, best of 3, measured in the same session. Above 1.5 times, the change is not ready: at the hospital's 380,000,000 reads per mate, the old script is about 800 times the stand-in.

**Verdict.** Supported if G1, G2, G3 and G4 all pass. Any failure is reported as measured, and nothing is merged.
