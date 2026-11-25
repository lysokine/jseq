#!/usr/bin/env python3
"""Generate the two documentation figures.

    python docs/make_figures.py          # writes docs/pipeline.{png,pdf} and docs/modes.{png,pdf}

Schematics only: every coordinate here comes from config/demo.json, the synthetic project, so the
figures carry no real experiment's sequence or dimensions. Palette is Okabe-Ito, which stays
readable in greyscale and under the common colour-vision deficiencies.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
CFG = json.load(open(os.path.join(HERE, "..", "config", "demo.json")))
A5 = CFG["amplicons"]["Demo_5p"]
A3 = CFG["amplicons"]["Demo_3p"]

# Okabe-Ito
BLUE, ORANGE, GREEN, SKY = "#0072B2", "#E69F00", "#009E73", "#56B4E9"
VERM, PURPLE, GREY, INK = "#D55E00", "#CC79A7", "#999999", "#222222"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "figure.dpi": 150})


def box(ax, x, y, w, h, label, fc="white", ec=INK, fs=9, lw=1.2, weight="normal", tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                fc=fc, ec=ec, lw=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fs,
            zorder=3, weight=weight, color=tc, linespacing=1.45)


def arrow(ax, x1, y1, x2, y2, color=INK, lw=1.3, style="-|>", ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=11,
                                 color=color, lw=lw, linestyle=ls, zorder=1,
                                 shrinkA=2, shrinkB=2))


# --------------------------------------------------------------------------- figure 1: pipeline

def pipeline():
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.965, "What the program does", ha="center", fontsize=13, weight="bold")
    ax.text(0.5, 0.925, "Two text files define a project. Everything below follows from them.",
            ha="center", fontsize=8.5, color="#555555")

    box(ax, 0.035, 0.80, 0.42, 0.085,
        "config.json\namplicons · anchors · landmarks · thresholds", fc="#EAF3FA", ec=BLUE, fs=8.0)
    box(ax, 0.545, 0.80, 0.42, 0.085,
        "samples.tsv\none row per library: files · sample · culture", fc="#EAF3FA", ec=BLUE, fs=8.0)

    arrow(ax, 0.245, 0.80, 0.34, 0.735, BLUE)
    arrow(ax, 0.755, 0.80, 0.66, 0.735, BLUE)

    box(ax, 0.30, 0.655, 0.40, 0.078, "reads  (FASTQ)", fc="white", ec=GREY, fs=9)
    ax.text(0.725, 0.694, "paired reads merged on\ntheir overlap first (fastp)",
            ha="left", va="center", fontsize=7.4, color="#555555", style="italic")

    arrow(ax, 0.50, 0.655, 0.50, 0.585)

    box(ax, 0.165, 0.495, 0.67, 0.088,
        "count each read   →   net change in junction length\n"
        "spacing mode  or  landmark mode  (next figure)",
        fc="#FFF4E3", ec=ORANGE, fs=8.6)

    arrow(ax, 0.50, 0.495, 0.50, 0.425)

    box(ax, 0.235, 0.335, 0.53, 0.085,
        "counts kept separately for every imaging tile\ncounts/<library>.json", fc="white",
        ec=GREY, fs=8.8)

    arrow(ax, 0.50, 0.335, 0.50, 0.265)
    ax.text(0.775, 0.300, "re-summarise without\nre-reading the FASTQs", ha="left", va="center",
            fontsize=7.4, color="#555555", style="italic")

    box(ax, 0.055, 0.175, 0.89, 0.085,
        "tile quality control\ndo the tiles of one library agree?  if not, that is the instrument talking, not the sample",
        fc="#E7F6F0", ec=GREEN, fs=8.2)

    arrow(ax, 0.50, 0.175, 0.50, 0.108)

    box(ax, 0.055, 0.018, 0.89, 0.085,
        "results\nper library · per culture (mean ± SEM) · indel spectrum · heatmap · report.md",
        fc="#F3EAF4", ec=PURPLE, fs=8.2)

    fig.savefig(os.path.join(HERE, "pipeline.png"), dpi=300, bbox_inches="tight",
                facecolor="white")
    fig.savefig(os.path.join(HERE, "pipeline.pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --------------------------------------------------------------------------- figure 2: the modes

def read_track(ax, y, x0, x1, color="#E9E9E9", h=0.030, ec="#BBBBBB"):
    ax.add_patch(Rectangle((x0, y), x1 - x0, h, fc=color, ec=ec, lw=0.8, zorder=2))


def seg(ax, y, x0, x1, color, h=0.030, label=None, fs=7.4, tc="white", lw=0.0, hatch=None):
    ax.add_patch(Rectangle((x0, y), x1 - x0, h, fc=color, ec=INK, lw=lw, zorder=3, hatch=hatch))
    if label:
        ax.text((x0 + x1) / 2, y + h / 2, label, ha="center", va="center", fontsize=fs,
                color=tc, zorder=4, weight="bold")


def modes():
    fig, ax = plt.subplots(figsize=(8.4, 7.4))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.980, "Two ways to measure the same thing", ha="center", fontsize=13,
            weight="bold")
    ax.text(0.5, 0.950, "Both ask: is the junction the length it was designed to be?",
            ha="center", fontsize=8.5, color="#555555")

    # ---------------- A: spacing
    ax.text(0.03, 0.905, "A", fontsize=12, weight="bold")
    ax.text(0.072, 0.905, "spacing  —  CRIS.py's method", fontsize=10.5, weight="bold")
    ax.text(0.072, 0.879, "Two anchors must match exactly. The distance between them is the answer.",
            fontsize=8.3, color="#555555")

    y = 0.820
    read_track(ax, y, 0.075, 0.900)
    seg(ax, y, 0.170, 0.265, BLUE, label="anchor")
    seg(ax, y, 0.680, 0.775, BLUE, label="anchor")
    ax.text(0.915, y + 0.015, "read", ha="left", va="center", fontsize=8, color="#555555")
    ax.annotate("", xy=(0.170, y - 0.020), xytext=(0.775, y - 0.020),
                arrowprops=dict(arrowstyle="<|-|>", color=INK, lw=1.1))
    ax.text(0.4725, y - 0.047, "measured distance", ha="center", fontsize=8)

    y2 = 0.705
    read_track(ax, y2, 0.075, 0.900, color="#F7F7F7")
    seg(ax, y2, 0.170, 0.265, BLUE, label="anchor")
    seg(ax, y2, 0.630, 0.725, BLUE, label="anchor")
    ax.text(0.915, y2 + 0.015, "read", ha="left", va="center", fontsize=8, color="#555555")
    ax.annotate("", xy=(0.170, y2 - 0.020), xytext=(0.725, y2 - 0.020),
                arrowprops=dict(arrowstyle="<|-|>", color=VERM, lw=1.1))
    ax.text(0.4475, y2 - 0.047, "shorter by 9 bp  \u2192  a 9 bp deletion", ha="center",
            fontsize=8, color=VERM)

    ax.text(0.072, 0.632,
            "Simple and exact. It fails when the junction itself makes the bases after it "
            "unreadable,\nbecause then the second anchor is missing from most reads.",
            fontsize=8.2, color="#555555", va="top")

    ax.plot([0.03, 0.97], [0.560, 0.560], color="#DDDDDD", lw=1)

    # ---------------- B: landmark
    ax.text(0.03, 0.520, "B", fontsize=12, weight="bold")
    ax.text(0.072, 0.520, "landmark  —  for junctions that corrupt the read", fontsize=10.5,
            weight="bold")
    ax.text(0.072, 0.494,
            "One exact anchor, then two landmarks placed clear of the damage. Matching tolerates "
            "errors.", fontsize=8.3, color="#555555")

    # the offset annotation sits above the track, with room of its own
    ax.text(0.515, 0.455, "slid over a range of offsets; best match wins", ha="center",
            fontsize=7.6, color="#8A6200")
    ax.annotate("", xy=(0.395, 0.437), xytext=(0.640, 0.437),
                arrowprops=dict(arrowstyle="<|-|>", color=ORANGE, lw=1.0))

    y3 = 0.385
    read_track(ax, y3, 0.075, 0.900)
    seg(ax, y3, 0.115, 0.200, BLUE, label="anchor", fs=7.0)
    seg(ax, y3, 0.225, 0.370, VERM, label="misread stretch", fs=7.0, hatch="///")
    seg(ax, y3, 0.405, 0.590, ORANGE, label="landmark 1", tc=INK, fs=7.4)
    seg(ax, y3, 0.630, 0.855, GREEN, label="landmark 2", fs=7.4)
    ax.text(0.915, y3 + 0.015, "read", ha="left", va="center", fontsize=8, color="#555555")

    ax.text(0.2975, y3 - 0.022, "bases here are wrong,\nand it does not matter", ha="center",
            va="top", fontsize=7.4, color=VERM)

    box(ax, 0.072, 0.190, 0.40, 0.088,
        "landmark 1  \u2192  how long?\noffset 0 = the designed length,\nanything else = a small indel",
        fc="#FFF4E3", ec=ORANGE, fs=8.0)
    box(ax, 0.528, 0.190, 0.40, 0.088,
        "landmark 2  \u2192  is that true?\nconfirms offset 0, and recovers\ndeletions that removed landmark 1",
        fc="#E7F6F0", ec=GREEN, fs=8.0)

    ax.text(0.5, 0.152,
            "Without landmark 2, a tandem duplication of landmark 1 would pass as the designed "
            "length.", ha="center", fontsize=8.0, color="#555555")

    box(ax, 0.072, 0.018, 0.856, 0.102,
        "What each one buys, and what it costs\n"
        "spacing reads the junction base by base, so it can check the sequence as well as the length.\n"
        "landmark measures length only: compensating insertions and deletions look designed.",
        fc="#F4F4F4", ec=GREY, fs=8.0)

    fig.savefig(os.path.join(HERE, "modes.png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(os.path.join(HERE, "modes.pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    pipeline()
    modes()
    print("wrote docs/pipeline.{png,pdf} and docs/modes.{png,pdf}")
