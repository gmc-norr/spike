# How spike works: slide deck

The source of a 27-slide talk on spike's inner workings: what it does, the statistical models behind each random step, and measurements from real runs. October 2026.

## Layout

| Path | What it holds |
|---|---|
| `deck.json` | Slide order, the four sections (Introduction, Methods, Results, Discussion) and the fonts |
| `slides/*.html` | One slide each, a single `<section>` on a 1920 × 1080 canvas; the speaker notes are in its `<aside>` |
| `figures/*.png` | The charts. Slides refer to them as `../figures/<name>.png` |
| `data/` | The numbers behind the charts and the statistics quoted on the slides |
| `scripts/` | How the data and the charts were made |

## Data

- **Sample:** a chr20 slice (38.5–40.2 Mb) of the GIAB HG002 NovaSeq PCR-free 35x BAM, aligned with bwa-mem2 2.2.1.
- **spike version:** master ba9f4e7 (October 7, 2026), after the quality model was rebuilt. The deck first used 97aa3b4; slides 10-11 and 18-20 changed most.
- **Runs:**
  - a heterozygous and a homozygous 10 kb deletion at chr20:38,900,000–38,910,000;
  - 25 SNVs, 60 kb apart, at allele fractions 0.1, 0.25, 0.5, 0.75 and 1.0.
- **Other tools (slides 23-24):** not from runs. Each cell of Table 8 was read from the tool's paper, documentation or code in October 2026: BAMSurgeon 1.4.1 and its unreleased 2026 rewrite, SomatoSim 1.0.0, VarBen 1.0, SVEngine 1.0.0, NEAT 4.7.1, VarSim 0.8.6 and VISOR 1.1.2.1. Re-check them before reusing the slides; these tools change (BAMSurgeon's structural-variant method did in 2026).

## Rebuilding the figures

```bash
REF=GRCh38_no_alt.fasta SLICE=hg002_chr20_slice.bam SPIKE=spike bash scripts/runs.sh
REF=GRCh38_no_alt.fasta SLICE=hg002_chr20_slice.bam SPIKE=spike bash scripts/extract.sh
```

- `runs.sh` makes the three spike runs and validates each one.
- `extract.sh` turns those runs into `data/` and `figures/`.
- The BAMs and the 16 MB read extract are not included.

Every chart is drawn from these runs; none is illustrative. Slides that show a schematic say "not data".
