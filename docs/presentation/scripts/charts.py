"""Charts for the spike deck, from the runs in this folder (real data only).

Each PNG is drawn at the size it takes on the 1920x1080 slide, at 2x pixels.
Colours: real reads blue, spike's reads orange (dataviz slots 1-2, validated
on the slide surface); a third series aqua (slot 3).
"""
import csv
import statistics
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SURFACE = "#FDFDFC"
INK = "#17202E"
MUTED = "#4A5363"
GRID = "#E2E2DE"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
BAND = "#EEEEEA"

plt.rcParams.update({
    "font.family": "Noto Sans",
    "font.size": 11,
    "axes.edgecolor": GRID,
    "axes.labelcolor": MUTED,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "axes.unicode_minus": False,
})

# slide px -> inches at 200 px/in; saved at dpi 400 = 2x pixels
PX = 1 / 200


def fig(w, h):
    return plt.figure(figsize=(w * PX, h * PX))


def save(f, name):
    f.savefig(name, dpi=400)
    plt.close(f)


# 1. Depth across the planted 10 kb deletion
rows = list(csv.DictReader(open("depth_bins.csv")))
x = [(int(r["start"]) - 38900000) / 1000 for r in rows]
f = fig(1664, 600)
ax = f.add_axes([0.07, 0.19, 0.92, 0.78])
ax.axvspan(0, 10, color=BAND, zorder=0, lw=0)
for key, col, label in [("hom", AQUA, "Hom deletion (af 1.0)"),
                        ("het", ORANGE, "Het deletion (af 0.5)"),
                        ("original", BLUE, "Before spike")]:
    y = [float(r[key]) for r in rows]
    ax.step(x, y, where="post", color=col, lw=2)
ax.set_xlim(x[0], x[-1] + 0.25)
ax.set_ylim(0, 76)
ax.set_xlabel("Position from the deletion's start (kb), chr20, 250 bp bins")
ax.set_ylabel("Read depth")
ax.grid(axis="x", visible=False)
# direct labels inside the deleted span, beside each line
ax.text(5, 62.5, "Before spike: 47x", color="#1c5cab", ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.text(5, 11, "Het (af 0.5): 24.5x", color="#c4501f", ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.text(5, 4.5, "Hom (af 1.0): 0.3x", color="#13805a", ha="center", va="bottom", fontsize=11, fontweight="bold")
ax.text(5, 75, "planted 10 kb deletion (shaded)", color=MUTED, ha="center", va="top", fontsize=10)
save(f, "chart_depth.png")

# 2. Requested vs measured allele fraction, 25 SNVs, with Wilson 95% CIs from each site's own count
rows = list(csv.DictReader(open("snv_af_counts.csv")))
f = fig(1000, 680)
ax = f.add_axes([0.13, 0.14, 0.83, 0.82])
ax.plot([0, 1], [0, 1], color=MUTED, lw=1, ls=(0, (4, 3)))
seen = {}
for r in rows:
    a = float(r["asked"]); o = float(r["frac"]); lo = float(r["lo"]); hi = float(r["hi"])
    k = seen.get(a, 0); seen[a] = k + 1
    xj = a + (k - 2) * 0.012
    ax.plot([xj, xj], [lo, hi], color=ORANGE, lw=1.4, alpha=0.75, zorder=2)
    if r["graded"] == "yes":
        ax.scatter(xj, o, s=40, color=ORANGE, edgecolor=SURFACE, linewidth=1.2, zorder=3)
    else:
        ax.scatter(xj, o, s=40, facecolor=SURFACE, edgecolor=ORANGE, linewidth=1.6, zorder=3)
ax.set_xlim(0, 1.05)
ax.set_ylim(0, 1.05)
ax.set_xticks([0, 0.1, 0.25, 0.5, 0.75, 1.0])
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xlabel("Allele fraction requested")
ax.set_ylabel("Allele fraction measured (alt / all fragments)")
ax.text(0.47, 0.06, "Bars: 95% Wilson CI per site\nHollow: too few fragments to grade", color=MUTED, fontsize=10, ha="left")
save(f, "chart_snv.png")

# 3. Fragment length, real vs spike
rows = list(csv.DictReader(open("insert.csv")))
f = fig(1000, 600)
ax = f.add_axes([0.10, 0.17, 0.86, 0.79])
bins = list(range(0, 1001, 20))
medians = {}
for kind, col in [("real", BLUE), ("spike", ORANGE)]:
    v = [int(r["tlen"]) for r in rows if r["kind"] == kind]
    medians[kind] = statistics.median(v)
    counts = [0] * (len(bins) - 1)
    for t in v:
        if t < 1000:
            counts[t // 20] += 1
    tot = len(v)
    ys = [c / tot * 100 for c in counts]
    ax.step(bins[:-1], ys, where="post", color=col, lw=2)
ax.set_xlim(0, 1000)
ax.set_ylim(0, None)
ax.text(520, ax.get_ylim()[1] * 0.80, f"Real pairs: median {medians['real']:g} bp", color="#1c5cab", fontsize=11, fontweight="bold")
ax.text(520, ax.get_ylim()[1] * 0.68, f"spike's pairs: median {medians['spike']:g} bp", color="#c4501f", fontsize=11, fontweight="bold")
ax.set_xlabel("Fragment length (bp)")
ax.set_ylabel("Share of pairs (%, 20 bp bins)")
ax.grid(axis="x", visible=False)
save(f, "chart_insert.png")

# 4. Mean base quality by cycle, R1 and R2
rows = list(csv.DictReader(open("qual_cycle.csv")))
f = fig(1664, 520)
cyc = [int(r["cycle"]) for r in rows]
for i, (mate, title) in enumerate([("R1", "Read 1"), ("R2", "Read 2")]):
    ax = f.add_axes([0.065 + i * 0.49, 0.22, 0.425, 0.68])
    ax.plot(cyc, [float(r["real_" + mate]) for r in rows], color=BLUE, lw=2)
    ax.plot(cyc, [float(r["spike_" + mate]) for r in rows], color=ORANGE, lw=2)
    ax.set_xlim(1, 151)
    ax.set_ylim(28, 38)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_xlabel("Cycle (base position in the read)")
    if i == 0:
        ax.set_ylabel("Mean base quality")
        ax.text(8, 30.6, "Real reads", color="#1c5cab", fontsize=11, fontweight="bold")
        ax.text(8, 29.2, "spike's reads", color="#c4501f", fontsize=11, fontweight="bold")
    ax.grid(axis="x", visible=False)
save(f, "chart_qual.png")
print("ok")
