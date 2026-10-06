"""Figure: spread of base quality, real vs spike's reads (same read set as Figure 11).
A, B: per-cycle standard deviation of base quality, read 1 and read 2.
C: distribution of each read's mean base quality.
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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

seen = set(); Q = {k: {1: [], 2: []} for k in ("real", "spike")}
for line in open("snv_reads.sam"):
    f = line.split("\t"); name, flag, qual = f[0], int(f[1]), f[10]
    key = (name, flag & 0xC0)
    if key in seen: continue
    seen.add(key)
    kind = "spike" if name.startswith("SPIKE_") else "real"
    q = np.frombuffer(qual.encode(), dtype=np.uint8).astype(float) - 33
    if flag & 0x10: q = q[::-1]
    Q[kind][1 if flag & 0x40 else 2].append(q)
for k in Q:
    for m in (1, 2): Q[k][m] = np.array(Q[k][m])

fig = plt.figure(figsize=(1664 * PX, 440 * PX))
cyc = np.arange(1, 152)
for i, (m, title) in enumerate([(1, "A  Read 1"), (2, "B  Read 2")]):
    ax = fig.add_axes([0.07 + i * 0.265, 0.24, 0.205, 0.64])
    for k, col in (("real", BLUE), ("spike", ORANGE)):
        ax.plot(cyc, Q[k][m].std(axis=0, ddof=1), color=col, lw=2)
    ax.set_xlim(1, 151); ax.set_ylim(0, 9)
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_xlabel("Cycle")
    if i == 0:
        ax.set_ylabel("SD of base quality")
        ax.text(10, 1.6, "Real reads", color="#1c5cab", fontsize=11, fontweight="bold")
        ax.text(10, 0.5, "spike's reads", color="#c4501f", fontsize=11, fontweight="bold")
    ax.grid(axis="x", visible=False)

ax = fig.add_axes([0.63, 0.24, 0.36, 0.64])
edges = np.arange(20, 37.5 + 1e-9, 0.5)
for k, col in (("real", BLUE), ("spike", ORANGE)):
    rm = np.vstack([Q[k][1], Q[k][2]]).mean(axis=1)
    h, _ = np.histogram(np.clip(rm, 20, 37.49), bins=edges)
    y = np.where(h > 0, h / len(rm) * 100, np.nan)   # empty bins break the line
    ax.stairs(y, edges, color=col, lw=2, baseline=None)
ax.set_yscale("log"); ax.set_ylim(0.003, 100); ax.set_xlim(20, 37.5)
ax.set_yticks([0.01, 0.1, 1, 10, 100]); ax.minorticks_off()
ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
ax.set_title("C  Mean quality of each read", loc="left", color=INK, fontsize=12, fontweight="bold")
ax.set_xlabel("Mean base quality of the read")
ax.set_ylabel("Share of reads (%, log)")
ax.grid(axis="x", visible=False)
ax.text(21, 12, "Real: SD 2.08", color="#1c5cab", fontsize=11, fontweight="bold")
ax.text(21, 4, "spike: SD 0.52", color="#c4501f", fontsize=11, fontweight="bold")
fig.savefig("chart_qspread.png", dpi=400)
print("ok")
