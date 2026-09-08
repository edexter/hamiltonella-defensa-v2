#!/usr/bin/env python3
"""
Disjoint-replicate test of the five delivered APSE genomes.

Each lineage's read pile was split into NON-OVERLAPPING subsets of ~150x and
each subset assembled independently. No read is shared between any two
replicates, so agreement is not an artefact of resampling the same reads.

Panel A  every replicate's assembled length against the delivered reference
Panel B  the one locus where replicates disagreed, and what is there

Reads the measured table results/replicate_test.tsv.

Usage:  python3 scripts/plot_replicate_test.py
Output: results/figures/replicate_test.png
"""

import os
import collections

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREEN = "#1baf7a"
PURPLE = "#7a4fbd"
RED = "#b0402a"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"

SHIPPED = {"S07_RHS": 39306, "S07_CdtB": 38993, "S05_RHS": 39289,
           "S05_hapR": 38980, "S05_hapA": 36248}
ORDER = ["S07_RHS", "S07_CdtB", "S05_RHS", "S05_hapR", "S05_hapA"]
COL = {"S07_RHS": BLUE, "S07_CdtB": ORANGE, "S05_RHS": BLUE,
       "S05_hapR": ORANGE, "S05_hapA": GREEN}
NICE = {"S07_RHS": "S07  APSE-RHS", "S07_CdtB": "S07  APSE-CdtB",
        "S05_RHS": "S05  APSE-RHS", "S05_hapR": "S05  hapR", "S05_hapA": "S05  hapA"}

# measured allele counts at the polymorphic locus, S05 RHS pile,
# reads carrying >=100 diagnostic 31-mers of one allele and <10 of the other
ALLELE1, ALLELE2 = 471, 345


def main():
    rows = collections.defaultdict(list)
    with open("results/replicate_test.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            rows[f["lineage"]].append(f)

    fig = plt.figure(figsize=(12.6, 8.6), facecolor=SURFACE)
    fig.suptitle("Replicate test — every genome rebuilt from non-overlapping sets of reads",
                 fontsize=15, color=INK, y=0.975, fontweight="semibold")
    fig.text(0.5, 0.937,
             "Each pile was cut into disjoint ~150× subsets and each assembled on its own. No read appears in two replicates (checked: 0 in every case).\n"
             "27 replicate assemblies in total. Agreement here cannot come from resampling the same reads.",
             ha="center", va="top", fontsize=9.8, color=INK2, linespacing=1.55)

    axA = fig.add_axes([0.235, 0.500, 0.715, 0.335])
    axA.set_title("A.   Assembled length of every replicate, against the delivered reference",
                  fontsize=11.5, color=INK, loc="left", pad=10, fontweight="medium")
    yt, ylab = [], []
    for i, lin in enumerate(ORDER):
        y = len(ORDER) - 1 - i
        yt.append(y); ylab.append(f"{NICE[lin]}\n{SHIPPED[lin]:,} bp")
        ref = SHIPPED[lin]
        exact = [r for r in rows[lin] if r["delta_len"] == "+0"]
        other = [r for r in rows[lin] if r["delta_len"] != "+0"]
        axA.axhline(y, color="#ecebe6", lw=8, zorder=0)
        if exact:
            axA.scatter([0] * len(exact), [y] * len(exact), s=170, color=COL[lin],
                        zorder=3, edgecolor="white", lw=1.5)
            axA.text(0, y + 0.30, f"{len(exact)}×", ha="center", fontsize=9.5,
                     color=COL[lin], fontweight="bold")
        for r in other:
            d = int(r["delta_len"])
            axA.scatter([d], [y], s=130, color=PURPLE if abs(d) < 100 else RED,
                        zorder=3, edgecolor="white", lw=1.4, marker="D")
    axA.axvline(0, color=INK2, lw=1.2, ls=(0, (3, 2)))
    axA.text(0, len(ORDER) - 0.42, "delivered length", ha="center", fontsize=9,
             color=INK2, fontweight="medium")
    axA.set_yticks(yt); axA.set_yticklabels(ylab, fontsize=9.5, color=INK, linespacing=1.5)
    axA.set_xlim(-3000, 400)
    axA.set_ylim(-0.55, len(ORDER) - 0.25)
    axA.set_xlabel("difference from the delivered reference (bp)", fontsize=9.5, color=INK2)
    axA.tick_params(labelsize=8.5, colors=MUTED, length=3)
    axA.grid(True, axis="x", color="#f0eee9", lw=0.6); axA.set_axisbelow(True)
    for s in ("top", "right", "left"):
        axA.spines[s].set_visible(False)
    axA.spines["bottom"].set_color("#dcdbd5")
    axA.set_facecolor(SURFACE)
    fig.text(0.235, 0.437,
             "21 of 27 replicates reproduce the delivered length exactly, at 100 % coverage and zero differences. The six that do not are the two\n"
             "S05 RHS replicates that found a second allele (purple), and both hapA replicates (red).",
             fontsize=8.8, color=INK2, va="top", linespacing=1.55)

    axB = fig.add_axes([0.075, 0.075, 0.875, 0.275])
    axB.axis("off")
    axB.set_title("B.   The one locus where replicates disagreed — S05's RHS phage is polymorphic there",
                  fontsize=11.5, color=INK, loc="left", pad=24, fontweight="medium")
    L = 39289
    for i, (lab, col, note) in enumerate([
            ("allele 1 — delivered as S05_APSE_RHS_v2", BLUE,
             f"{ALLELE1} reads  ·  {100*ALLELE1/(ALLELE1+ALLELE2):.0f} %  ·  4 of 6 replicates"),
            ("allele 2 — delivered as S05_APSE_RHS_alt_v2", PURPLE,
             f"{ALLELE2} reads  ·  {100*ALLELE2/(ALLELE1+ALLELE2):.0f} %  ·  2 of 6 replicates")]):
        y = 0.60 - i * 0.24
        axB.add_patch(Rectangle((0.30, y), 0.50, 0.13, facecolor="#dcdad4",
                                edgecolor="none"))
        axB.add_patch(Rectangle((0.30 + 0.50 * 29702 / L, y),
                                max(0.50 * 436 / L, 0.006), 0.13, facecolor=col, edgecolor="none"))
        axB.text(0.29, y + 0.065, lab, ha="right", va="center", fontsize=9.5,
                 color=col, fontweight="medium")
        axB.text(0.815, y + 0.065, note, ha="left", va="center", fontsize=9,
                 color=INK2)
    axB.annotate("", xy=(0.30 + 0.50 * 29920 / L, 0.80), xytext=(0.30 + 0.50 * 29920 / L, 0.93),
                 arrowprops=dict(arrowstyle="-", color=INK2, lw=1))
    axB.text(0.30 + 0.50 * 29920 / L, 0.965,
             "436 bp cassette at 29,702–30,138, immediately upstream of attP",
             ha="center", va="bottom", fontsize=9, color=INK, fontweight="medium")
    axB.text(0.02, 0.22,
             "The two alleles share no 31-mer — unrelated sequence, not variants of one another. Allele 1 is also carried by S05's hapR\n"
             "and hapA and by S07's CdtB phage; allele 2 is identical to the same locus in S07's RHS phage. Reads carrying either\n"
             "allele carry integrase A and map best to RHS, so this is variation WITHIN the RHS lineage, not contamination from\n"
             "another lineage. Both versions are now delivered.",
             fontsize=9, color=INK2, va="top", linespacing=1.6)

    legend = [Line2D([], [], marker="o", ls="", ms=10, mfc=BLUE, mec="white",
                     label="reproduces the delivered length exactly"),
              Line2D([], [], marker="D", ls="", ms=9, mfc=PURPLE, mec="white",
                     label="found the alternative allele"),
              Line2D([], [], marker="D", ls="", ms=9, mfc=RED, mec="white",
                     label="shorter — hapA, as expected at half depth")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.005))

    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/replicate_test.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


main()
