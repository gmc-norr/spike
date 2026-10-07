"""Method-illustration charts, from real data (qual_stats.json, snv_accept.csv)."""
import csv, json, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
exec(open("charts.py").read().split("# 1. Depth")[0])  # shared style + helpers

st = json.load(open("qual_stats.json"))
from matplotlib.ticker import FuncFormatter
PLAIN = FuncFormatter(lambda v, _: f"{v:g}")
RAMP = LinearSegmentedColormap.from_list("blue", ["#eef4fc", "#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"])

# A. First-order transitions P(Q_i | bin of Q_{i-1}), real vs spike
f = fig(1000, 560)
labels_prev = ["Q10–19", "Q20–29", "Q30+"]
labels_next = ["Q11", "Q25", "Q37"]
for i, kind in enumerate(["real", "spike"]):
    ax = f.add_axes([0.14 + i * 0.44, 0.17, 0.36, 0.66])
    t = st["trans"][kind]
    M = []
    for a in (1, 2, 3):
        row = [t.get(f"{a}>{b}", 0) for b in (11, 25, 37)]
        s = sum(row); M.append([v / s for v in row])
    ax.imshow(M, cmap=RAMP, vmin=0, vmax=1, aspect="auto")
    for r in range(3):
        for c in range(3):
            v = M[r][c]
            ax.text(c, r, f"{v:.3f}", ha="center", va="center", fontsize=11,
                    color="#ffffff" if v > 0.55 else INK, fontweight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels(labels_next)
    ax.set_yticks(range(3)); ax.set_yticklabels(labels_prev if i == 0 else [])
    ax.grid(False)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title("Real reads" if kind == "real" else "spike's reads", loc="left", fontsize=12, fontweight="bold", color=INK)
    ax.set_xlabel("Quality of base i")
    if i == 0: ax.set_ylabel("Quality bin of base i-1")
save(f, "chart_markov.png")

# B. Low-quality bases among the last 20 cycles, real vs spike (share of reads, log scale)
f = fig(1000, 560)
ax = f.add_axes([0.13, 0.16, 0.84, 0.78])
ax.axvspan(9.5, 20.5, color=BAND, lw=0, zorder=0)
n_reads = {}
for kind, col in [("real", BLUE), ("spike", ORANGE)]:
    c = {int(k): v for k, v in st["lowtail"][kind].items()}
    tot = sum(c.values())
    n_reads[kind] = tot
    xs = list(range(21))
    ax.plot(xs, [c[x] / tot * 100 if c.get(x, 0) else float("nan") for x in xs], color=col, lw=2, marker="o", ms=4)
ax.set_yscale("log")
ax.set_ylim(0.002, 200)
ax.yaxis.set_major_formatter(PLAIN)
ax.set_xlim(-0.5, 20.5)
ax.set_xticks([0, 5, 10, 15, 20])
ax.set_xlabel("Bases below Q15 among a read's last 20 cycles")
ax.set_ylabel("Share of reads (%, log scale)")
ax.grid(axis="x", visible=False)
ax.text(10.2, 70, "crashed end (10 or more)", color=MUTED, fontsize=10)
ax.text(10.4, 14, f"Real reads (n = {n_reads['real']:,})", color="#1c5cab", fontsize=11, fontweight="bold")
ax.text(10.4, 5, f"spike's reads (n = {n_reads['spike']:,})", color="#c4501f", fontsize=11, fontweight="bold")
save(f, "chart_lowtail.png")

# C. Mismatch rate by base quality (aligned bases, variant positions excluded)
def wilson(k, n, z=1.96):
    c = (k + z * z / 2) / (n + z * z); h = z * math.sqrt(k * (n - k) / n + z * z / 4) / (n + z * z)
    return c - h, c + h
f = fig(1000, 560)
ax = f.add_axes([0.13, 0.16, 0.84, 0.78])
qs = [x / 2 for x in range(16, 82)]
ax.plot(qs, [10 ** (-q / 10) * 100 for q in qs], color=MUTED, lw=1.2, ls=(0, (4, 3)))
ax.text(28, 10 ** (-2.4) * 100 * 1.6, "10^(-Q/10)", color=MUTED, fontsize=10)
for kind, col, dx in [("real", BLUE, -0.6), ("spike", ORANGE, 0.6)]:
    for q in (11, 25, 37):
        k, n = st["mm"][kind][str(q)]
        lo, hi = wilson(k, n)
        ax.plot([q + dx, q + dx], [lo * 100, hi * 100], color=col, lw=2)
        ax.scatter(q + dx, k / n * 100, s=44, color=col, edgecolor=SURFACE, linewidth=1.2, zorder=3)
ax.set_yscale("log")
ax.yaxis.set_major_formatter(PLAIN)
ax.set_xlim(8, 40)
ax.set_xticks([11, 25, 37])
ax.set_xlabel("Base quality (Phred)")
ax.set_ylabel("Mismatch rate (%, log scale)")
ax.grid(axis="x", visible=False)
ax.text(27, 3, "Real reads", color="#1c5cab", fontsize=11, fontweight="bold")
ax.text(27, 1.6, "spike's reads", color="#c4501f", fontsize=11, fontweight="bold")
save(f, "chart_errors.png")

# D. validate's allele-fraction rule on the 25 real sites: 99% acceptance range at each site's n
rows = list(csv.DictReader(open("snv_accept.csv")))
f = fig(1664, 500)
ax = f.add_axes([0.07, 0.17, 0.92, 0.78])
groups = {}
for r in rows:
    groups.setdefault(float(r["asked"]), []).append(r)
x = 0; ticks = []
for g, (fa, rs) in enumerate(sorted(groups.items())):
    xs = []
    for r in rs:
        fr = float(r["frac"])
        if r["verdict"] == "too shallow":
            ax.scatter(x, fr, s=46, facecolor=SURFACE, edgecolor=ORANGE, linewidth=1.6, zorder=3)
        else:
            ax.plot([x, x], [float(r["lo"]), float(r["hi"])], color="#B9C2CC", lw=9, solid_capstyle="butt", zorder=1)
            ax.scatter(x, fr, s=46, color=ORANGE, edgecolor=SURFACE, linewidth=1.2, zorder=3)
        xs.append(x); x += 1
    ax.plot([xs[0] - 0.4, xs[-1] + 0.4], [fa, fa], color=INK, lw=1.2, ls=(0, (3, 2)))
    ticks.append(((xs[0] + xs[-1]) / 2, f"requested {fa:g}"))
    x += 1.5
ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks])
ax.set_ylim(0, 1.05); ax.set_xlim(-1, x - 1)
ax.set_ylabel("Allele fraction")
ax.grid(axis="x", visible=False)
save(f, "chart_accept.png")
print("ok")
