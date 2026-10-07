"""Figure: crashed 3' ends by the template's longest one-letter run, real vs spike's reads (the reads of Figure 11).

A read is crashed when 10+ of its last 20 qualities (sequencing order) are below Q15, as in Figure 13.
Its run bin is the longest one-letter run in the reference under the read, in sequencing order, binned
as spike's quality::run_bin: A/G/T 7-8 -> 1, 9-11 -> 2, 12+ -> 3; C 5-6 -> 2, 7+ -> 3; else 0.
Needs REF (the FASTA) in the environment and snv_reads.sam from extract.sh. Writes crash_runbin.csv
and chart_runbin.png, and prints the table with Wilson 95% CIs.
"""
import csv, math, os, re
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

SURFACE = "#FDFDFC"; INK = "#17202E"; MUTED = "#4A5363"; GRID = "#E2E2DE"
BLUE = "#2a78d6"; ORANGE = "#eb6834"
plt.rcParams.update({
    "font.family": "Noto Sans", "font.size": 11, "axes.edgecolor": GRID,
    "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.unicode_minus": False,
})
PX = 1 / 200
COMP = str.maketrans("ACGTN", "TGCAN")
fa = pysam.FastaFile(os.environ["REF"])


def run_bin(base, n):
    if base == "C":
        return 3 if n >= 7 else 2 if n >= 5 else 0
    if base in "AGT":
        return 3 if n >= 12 else 2 if n >= 9 else 1 if n >= 7 else 0
    return 0


def template_bin(chrom, pos1, cigar, flag, length):
    start = pos1 - 1
    span = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar))
    if flag & 0x10:
        end = start + span
        t = fa.fetch(chrom, max(0, end - length), end).upper().translate(COMP)[::-1]
    else:
        t = fa.fetch(chrom, start, start + length).upper()
    best, base, run = 0, "N", 0
    for b in t:
        if b == base and b != "N":
            run += 1
        else:
            base, run = b, int(b != "N")
        best = max(best, run_bin(base, run))
    return best


def wilson(k, n):
    if n == 0:
        return (0.0, 0.0)
    z = 1.96; p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


seen = set()
t = {k: np.zeros((4, 2), int) for k in ("real", "spike")}
for line in open("snv_reads.sam"):
    f = line.split("\t"); name, flag, chrom, pos, cigar, seq, qual = f[0], int(f[1]), f[2], int(f[3]), f[5], f[9], f[10]
    key = (name, flag & 0xC0)
    if key in seen or cigar == "*":
        continue
    seen.add(key)
    q = np.frombuffer(qual.encode(), dtype=np.uint8).astype(int) - 33
    if flag & 0x10:
        q = q[::-1]
    crashed = int((q[-20:] < 15).sum() >= 10)
    h = template_bin(chrom, pos, cigar, flag, len(seq))
    k = "spike" if name.startswith("SPIKE_") else "real"
    t[k][h, 0] += 1; t[k][h, 1] += crashed

rows = []
with open("crash_runbin.csv", "w", newline="") as out:
    w = csv.writer(out); w.writerow(["set", "run_bin", "reads", "crashed", "share", "lo95", "hi95"])
    for k in ("real", "spike"):
        for h in range(4):
            n, c = t[k][h]; lo, hi = wilson(c, n)
            w.writerow([k, h, n, c, f"{c / n if n else 0:.5f}", f"{lo:.5f}", f"{hi:.5f}"])
            rows.append((k, h, n, c, lo, hi))
for k, h, n, c, lo, hi in rows:
    print(f"{k:5s} bin {h}: {c:4d} / {n:6d} = {100 * c / max(n, 1):5.2f}% ({100 * lo:.2f}-{100 * hi:.2f}%)")
for k in ("real", "spike"):
    n, c = t[k].sum(axis=0)
    print(f"{k:5s} all: {c} / {n} = {100 * c / n:.2f}%; run bin 3 holds {100 * t[k][3, 0] / n:.1f}% of reads")

fig, ax = plt.subplots(figsize=(888 * PX, 497 * PX))
labels = ["none", "A/G/T 7-8", "A/G/T 9-11\nor C 5-6", "A/G/T 12+\nor C 7+"]
x = np.arange(4); wbar = 0.36
for i, (k, col, lab) in enumerate([("real", BLUE, "Real reads"), ("spike", ORANGE, "spike")]):
    share = np.array([c / n if n else 0 for n, c in t[k]]) * 100
    lo = np.array([wilson(c, n)[0] for n, c in t[k]]) * 100
    hi = np.array([wilson(c, n)[1] for n, c in t[k]]) * 100
    xs = x + (i - 0.5) * wbar
    ax.bar(xs, share, wbar * 0.92, color=col, label=lab, zorder=2)
    ax.errorbar(xs, share, yerr=[share - lo, hi - share], fmt="none", ecolor=INK, elinewidth=1, capsize=3, zorder=3)
ax.set_yscale("log")
ax.set_ylim(0.1, 40)
ax.set_yticks([0.1, 1, 10]); ax.minorticks_off()
ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_xticks(x, labels)
ax.set_xlabel("Longest one-letter run in the read's template")
ax.set_ylabel("Crashed reads (%, log)")
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig("chart_runbin.png", dpi=400)
