#!/usr/bin/env python3
"""
Read depth in 1 kb windows across the five APSE-derived contigs of the
pre-curation S07 assembly, against the coverage Flye reported for each.

The point of the figure: the number in assembly_info.txt is not a measurement
of read depth for a contig that Flye designated an alternative haplotype. Two
contigs reported at 20x carry sequence that hundreds to thousands of reads
cover. That number is what the deduplication rule was applied to.

Depth is measured by mapping all 25,153 S07 reads against the five contigs
together, so reads are shared out between them exactly as an assembler would
have to. Shared backbone is therefore split; lineage-specific stretches are not.

Usage:  python3 scripts/plot_S07_contig_coverage.py <five_cov.tsv>
Output: results/figures/S07_contig_coverage.png
"""

import os
import sys
import collections

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREEN = "#1baf7a"
RED = "#b0402a"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"

# contig -> (Flye reported coverage, Flye multiplicity, alt_group, lineage, marker)
META = {
    "contig_4": (2221, 49, "8",  "RHS",  BLUE,   "RHS toxin"),
    "contig_6": (  20,  1, "8",  "CdtB", ORANGE, "cdtB toxin"),
    "contig_3": (2045, 44, "10", "RHS",  BLUE,   "integrase A"),
    "contig_1": (  20,  1, "10", "CdtB", ORANGE, "integrase B"),
    "contig_5": (1613, 34, "—",  "CdtB", ORANGE, "integrase B + cdtB toxin"),
}
ORDER = ["contig_4", "contig_6", "contig_3", "contig_1", "contig_5"]
CHROM_DEPTH = 34


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "five_cov.tsv"
    cov = collections.defaultdict(list)
    for line in open(path):
        c, w, d = line.split()
        cov[c].append((int(w), float(d)))
    for c in cov:
        cov[c].sort()

    fig, axes = plt.subplots(5, 1, figsize=(11.0, 12.2), facecolor=SURFACE)
    fig.suptitle("Read depth across the five APSE contigs of S07, against the coverage Flye reported",
                 fontsize=13.5, color=INK, y=0.984, fontweight="semibold")
    fig.text(0.5, 0.958,
             "Every window is measured from the sample's own reads. The dashed line is the number in assembly_info.txt — the number the 2024\n"
             "deduplication rule was applied to. For the two contigs Flye called alternative haplotypes it is off by roughly sixtyfold.",
             ha="center", va="top", fontsize=9.6, color=INK2, linespacing=1.55)

    for ax, c in zip(axes, ORDER):
        rep, mult, grp, lineage, colour, marker = META[c]
        xs = [w + 0.5 for w, _ in cov[c]]
        ys = [d for _, d in cov[c]]
        ax.fill_between(xs, 0.5, ys, color=colour, alpha=0.28, zorder=2)
        ax.plot(xs, ys, color=colour, lw=1.6, zorder=3)

        ax.axhline(rep, color=RED, lw=1.5, ls=(0, (4, 2)), zorder=4)
        ax.text(len(xs) * 0.995, rep * 1.35, f"Flye reported {rep:,}×", fontsize=8.6,
                color=RED, ha="right", va="bottom", fontweight="medium")
        ax.axhline(CHROM_DEPTH, color=MUTED, lw=1, ls=(0, (1.5, 2)), zorder=4)
        ax.text(0.15, CHROM_DEPTH * 1.25, "chromosome 34×", fontsize=7.8, color=MUTED,
                va="bottom")

        mean = sum(ys) / len(ys)
        grp_txt = f"alt_group {grp}, multiplicity {mult}" if grp != "—" else "not in an alt_group"
        ax.set_title(f"{c}   ·   {lineage} lineage   ·   {marker}   ·   {grp_txt}",
                     fontsize=10.2, color=INK, loc="left", pad=7, fontweight="medium")
        ax.text(len(xs) + 1.0, mean, f"measured mean {mean:,.0f}×",
                ha="left", va="center", fontsize=9.5, color=colour,
                fontweight="semibold")

        ax.set_yscale("log")
        ax.set_ylim(8, 9000)
        ax.set_xlim(0, 39.5)
        ax.set_ylabel("depth (log)", fontsize=8.8, color=INK2)
        ax.tick_params(labelsize=8, colors=MUTED, length=3)
        ax.grid(True, which="major", color="#ecebe6", lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("bottom", "left"):
            ax.spines[s].set_color("#dcdbd5")
        ax.set_facecolor(SURFACE)

    axes[-1].set_xlabel("position along the contig (kb)", fontsize=9.5, color=INK2)

    legend = [Line2D([], [], color=BLUE, lw=3, label="RHS lineage"),
              Line2D([], [], color=ORANGE, lw=3, label="CdtB lineage"),
              Line2D([], [], color=RED, lw=1.5, ls=(0, (4, 2)),
                     label="coverage reported by Flye")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.004))

    fig.subplots_adjust(left=0.075, right=0.975, top=0.912, bottom=0.062, hspace=0.42)
    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/S07_contig_coverage.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
