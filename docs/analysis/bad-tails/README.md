# What real bad tails look like

October 7, 2026. An offline analysis (Python, outside spike) of the reads whose 3′ end "crashes" (10 or more of the last 20 qualities below Q15), done before the second quality model was built. It follows [the soft-clip analysis](../soft-clips/README.md), which found bad ends to be the largest cause of soft clips.

## Data and scripts

- **Reads:** the same chr20 slice (38.5–40.2 Mb) of the GIAB HG002 NovaSeq PCR-free 35x BAM as the soft-clip analysis. Everything is in sequencing orientation: reverse reads are flipped and complemented.
- `collect.py SLICE_BAM REF OUT.pkl [all]` collects per-base data: quality, called base, reference base, clip state, and whether the base is counted. It needs `cause_per_read.pkl` from `../soft-clips/clips/summary.py`, in `clips/`. Without `all` it keeps crashed and bad-end-clipped reads plus 1 in 10 good reads; with `all` it keeps 1 in 8 of every read (an unbiased sample).
- `analyze.py` to `analyze5.py` are steps 2 to 6; each says in its docstring what it measures. Steps 2–4 use the grouped sample; steps 5–6 (`analyze4.py`, `analyze5.py`) the unbiased one.
- `slice.out*` are their printed outputs. Each file starts with the section it prints.
- The pickles (`slice.pkl`, 142 MB; `slice_all.pkl`, 157 MB) are not kept; `collect.py` rebuilds them.

## Findings

- **Q2 only marks N.** The slice's qualities are 2, 11, 25 and 37. "Low" means Q11, mixed with Q25 and Q37 right to the read's end. A run of lows at the end is a median of 1 base, and no read stays low once low. So count lows over a window, not in a run.
- **One-letter runs in the template set off crashes.** After a run, the strand reading into it has more low qualities and more errors for the rest of the read; 60 cycles later it is still 7x for runs of 12 or more. The other strand is flat.
  - Runs of 9–11: low quality 2–2.5x, errors 2.7–3.5x.
  - Runs of 12+: low quality 7–19x, errors 14–43x.
  - Poly-C (in sequencing direction) is worst: runs of 7–8 already give 6x. Poly-G hardly matters.
  - Reads with a run of 12+ crash 15.0% of the time, against 0.46% with runs of 4 or less. About 38% of crashes are set off by runs of 7 or more.
  - At shared clip spots, same-strand reads go from 1.1% to 8.2% errors after the spot; the other strand stays at 0.8%.
- **How bad a tail is belongs partly to the read.** Q11 bases err 12% of the time in mild crashed tails and 51% in bad ones; even Q37 bases err 7% in bad tails. The per-read error rate has SD 0.16, against 0.108 if every read erred alike. The qualities predict part of it (tails with more Q11 err at 0.22 to 0.43). R1 and R2 tail rates correlate 0.25.
- **Not it.** Tails are not shifted by a base (0.9% of clips). Errors barely clump inside a tail (0.435 against 0.408). Wrong calls lean slightly to A and G (two-colour chemistry).
- **A correction to the soft-clip causes.** About 950 crashed reads were put under "the sample's own variant" by the rule "3 or more reads clip at the same place". But 27% of crashed clipped reads share a boundary on the same strand, against about 3% when shuffled and 2.5% on the other strand. Most of them are bad ends set off by the sequence, so bad ends are more than the 44.7% in the soft-clip table.

## What came of it

spike knows each read's template, so it can use "the longest one-letter run read so far" as a context for quality and errors, which a compressor-style model such as fqzcomp cannot. The second quality model (on master since 4a5672d) does this with its run bins; see `docs/analysis/quality-model-v2/`.
