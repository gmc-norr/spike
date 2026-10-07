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
