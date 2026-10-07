# Quality model v2: the big held-out clip test

K2b of `docs/superpowers/plans/2026-10-07-quality-model-v2.md`.

**What it checks.** Whether spike's reads soft-clip at a bad end as often as the sample's own reads, on many more reads than the 22 hospital events of K2. Those events hold about 8,000 reads a side, too few to see a shortfall of 8%.

**How.**
1. The ignored test `synth::tests::measure_clip_share` samples the BAM as a run does (`quality::sample_input`).
2. It splits the pairs by a hash of their name and learns from the even half.
3. It writes the odd half, |TLEN| at least one read and both mates full length, twice:
   - as sequenced;
   - as spike makes them from the reference at each mate's own 5' position and strand, each block's pairs with that block's class mix.
4. `k2b.sh` aligns each set with bwa-mem2 and scores it with `clipscore.py`. A clip is bad-end when it is 1-4 bp, or matches the reference at at least 50% of its placed bases, with no Q100 variant within 10 bp of its boundary.

**Run.** From the repository root:

```bash
BAM=... REF=... Q100_VCF=... OUT=... docs/analysis/quality-model-v2/k2b.sh
```

| Variable | What it is |
|---|---|
| `BAM` | an indexed BAM; the plan used HG002 35x |
| `REF` | GRCh38 no-alt FASTA, bwa-mem2 indexed |
| `Q100_VCF` | the GIAB HG002 T2T-Q100 v1.1 benchmark VCF |
| `OUT` | an empty directory |

It writes the real set once and spike's sets for seeds 21, 22 and 23.

## The two correction scripts (K6 and K3)

Both back the `## Correction (2026-10-07): K3 and K6` section of
`docs/superpowers/plans/2026-10-07-quality-model-v2.md`. Neither is a gate; they exist so the two
withdrawn claims can be re-checked.

### `k6_runbins.py` — K6's crash shares, by template run bin

**What it shows.** That K6's "24.3% against 14.4%" was a composition effect. It splits both sets by
the longest one-letter run in the reference under each read (binned as `quality::run_bin` bins it)
and prints the crashed-read share per bin, per mate, plus what spike's reads *would* crash at if
they took the donors' rate in each bin. Within each bin the two sets agree; spike's reads simply
concentrate over a T13 about 200 bp left of the deletion edge — 23.8% of them sit in bin 3 against
6.7% of the window's reads.

**Run.** From the repository root:

```bash
python3 docs/analysis/quality-model-v2/k6_runbins.py SIM_BAM ORIGINAL_BAM REFERENCE_FASTA
```

| Argument | What it is |
|---|---|
| `SIM_BAM` | spike's aligned reads for the K6 events (`align.sh`'s `sim.bam`); only `SPIKE_` reads are counted |
| `ORIGINAL_BAM` | the 31-value `HG002.GRCh38.chr20.bam` spike ran on; every read in the two windows is counted |
| `REFERENCE_FASTA` | GRCh38 no-alt FASTA, faidx'd |

Needs `pysam` and `numpy`. The two windows are the K6 events' own, hard-coded as `WINDOWS`.

### `mates.py` — the mate crash link from paired FASTQ

**What it shows.** That K3's "4.74 against 11.6" rested on five both-crashed pairs. For each
prefix it prints P(R1 crashes), P(R2 crashes), P(both) and the link P(both) / (P(R1) · P(R2)) — 1
means mates crash independently. On 62,682 held-out 35x pairs spike's link is 10.55-11.91 against
the sample's 7.77, so spike does not under-link there.

**Run.** From the repository root, with one or more prefixes; each needs `<prefix>_R1.fq` and
`<prefix>_R2.fq`, **uncompressed and in the same order**:

```bash
python3 docs/analysis/quality-model-v2/mates.py OUT/real OUT/spike21 OUT/spike22 OUT/spike23
```

`k2b.sh` above writes exactly those files, so its `OUT` is the directory to point this at. Pure
Python; no dependencies. A read counts as crashed on the same rule K1, K2 and K6 use: 10 or more of
its last 20 qualities below Q15.
