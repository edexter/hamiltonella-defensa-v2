#!/usr/bin/env python3
"""
Three-panel figure explaining the APSE architecture of sample S05 (strain H76).

Panel A  where the prophages sit in the chromosome, and the fact that one site is
         occupied by different phages in different cells
Panel B  the mosaic composition of hapA, a recombinant phage
Panel C  the read evidence underlying both claims

Every number plotted is measured, not schematic. Sources are noted per panel.

Usage:  python3 scripts/plot_S05_architecture.py
Output: results/figures/S05_APSE_architecture.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.lines import Line2D

# Validated categorical palette (light mode); CVD-checked, see dataviz reference
RHS = "#2a78d6"      # RHS-carrying phage
HAPR = "#eb6834"     # CdtB phage, integrase B
HAPA = "#1baf7a"     # CdtB phage, integrase A -- the recombinant
CHROM = "#c9c8c1"    # bacterial chromosome
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"

# --- measured values -------------------------------------------------------
# All values below are measured against the final reconstructions (2026-08-03)
# using the locally available filtered read set for S05 (9,828 reads). An
# earlier version of this figure used the full raw read set and the round-4
# reconstructions; where the two disagree it is noted.

# Integration sites: split-read clustering against the chromosome plus the three
# phages. Coordinates refer to the delivered 2024 curated contigs.
SITE1 = dict(contig="contig_02", length=210195, pos=187054, occupants=[("hapR", 26)])
SITE2 = dict(contig="contig_11", length=811192, pos=805721,
             occupants=[("hapA", 7), ("RHS", 5)])
CHROM_DEPTH = {"site1": 11.9, "site2": 14.5}

# hapA parentage in 1 kb windows, minimap2 -x asm20 against each parent. A window
# is assigned to a parent only when the identity margin exceeds 0.3 %; below that
# the parents are themselves identical and the window carries no information
# about its origin.
HAPA_WINDOWS = [(0, 3000, "hapR"), (3000, 8000, "same"), (8000, 14000, "hapR"),
                (14000, 18000, "same"), (18000, 25000, "RHS"),
                (25000, 26000, "same"), (26000, 27000, "RHS"),
                (27000, 34000, "same"), (34000, 36248, "RHS")]
HAPA_LEN = 36248

# feature coordinates, from minimap2 placement of the toxin and integrase panels
FEATURES = {
    "RHS":  dict(length=39289, toxin=("RHS", 9545, 13907), integrase=("A", 30817, 31981)),
    "hapR": dict(length=38980, toxin=("CdtB", 17576, 18554), integrase=("B", 87, 1111)),
    "hapA": dict(length=36248, toxin=("CdtB", 1, 926), integrase=("A", 18462, 19626)),
}

# genomic position of hapR<->RHS junctions, from split alignments of the reads
# against the two parent phages. 185 reads split between the parents; 166 of them
# break at exactly hapR position 5,883.
BREAKPOINTS = {1: 4, 4: 1, 5: 168, 9: 9, 17: 1, 34: 2}   # hapR coordinate (kb) -> reads
HAPA_JUNCTION_KB = 5.9

COL = {"RHS": RHS, "hapR": HAPR, "hapA": HAPA, "same": "#d8d7d0"}


def panel_a(ax):
    """Chromosome with its two integration sites, and who occupies each."""
    ax.set_title("A.  Two integration sites in the S05 chromosome",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")

    rows = [(SITE1, 1.0, "site1"), (SITE2, 0.0, "site2")]
    for site, y, key in rows:
        L = site["length"]
        ax.add_patch(Rectangle((0, y), L, 0.16, facecolor=CHROM, edgecolor="none"))
        ax.text(-L * 0.015, y + 0.08, site["contig"], ha="right", va="center",
                fontsize=8.5, color=INK2)
        ax.text(L * 0.5, y - 0.055, f"{L/1000:.0f} kb chromosome fragment · {CHROM_DEPTH[key]:.1f}× depth",
                ha="center", va="top", fontsize=7.5, color=MUTED)

        x = site["pos"]
        total = sum(n for _, n in site["occupants"])
        # stack the occupants above the att site, height proportional to read support
        base = y + 0.19
        for name, n in site["occupants"]:
            h = 0.30 * n / 26
            ax.add_patch(Rectangle((x - L * 0.055, base), L * 0.11, h,
                                   facecolor=COL[name], edgecolor="white", lw=1.2))
            ax.text(x + L * 0.075, base + h / 2, f"{name}   {n} junction reads",
                    ha="left", va="center", fontsize=8.5, color=INK2)
            base += h + 0.035
        ax.plot([x, x], [y + 0.16, y + 0.19], color=INK2, lw=1)
        ax.text(x, y + 0.055, f"att ~{x/1000:.0f} kb", ha="center", va="center",
                fontsize=7.5, color=INK, fontweight="bold")

    ax.text(SITE1["length"] * 0.30, 1.585,
            "one occupant only", fontsize=8.5, color=INK2, style="italic")
    ax.text(SITE2["length"] * 0.22, 0.60,
            "TWO different phages — different cells carry different ones",
            fontsize=8.5, color=INK2, style="italic")

    ax.set_xlim(-130000, 1230000); ax.set_ylim(-0.22, 1.62)
    ax.axis("off")


def panel_b(ax):
    """hapA is a mosaic: some segments match the CdtB phage, one matches RHS."""
    ax.set_title("B.  hapA is a recombinant — hapR at one end, RHS at the other, one junction between them",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")

    y = 1.05
    for a, b, src in HAPA_WINDOWS:
        ax.add_patch(Rectangle((a, y), b - a, 0.22, facecolor=COL[src],
                               edgecolor="white", lw=1.4))
    ax.text(-1200, y + 0.11, "hapA", ha="right", va="center", fontsize=9, color=INK2)
    ax.annotate("", xy=(18000, y + 0.27), xytext=(36248, y + 0.27),
                arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.2))
    ax.text(27500, y + 0.315, "RHS-derived from here to the end",
            ha="center", va="bottom", fontsize=8, color=INK2)
    ax.plot([14000, 14000], [y - 0.09, y], color=INK, lw=2.4)
    ax.plot([18000, 18000], [y - 0.09, y], color=INK, lw=2.4)
    ax.annotate("", xy=(14000, y - 0.05), xytext=(18000, y - 0.05),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.4))
    ax.text(16000, y - 0.085,
            "the junction is somewhere in this 4 kb window\n"
            "(166 split reads place it here; the parents are\nidentical across it, so it cannot be pinned down)",
            ha="center", va="top", fontsize=8, color=INK, fontweight="medium")
    ax.text(0, -0.11,
            "grey = the two parents are identical here, so those segments cannot be assigned to either",
            ha="left", va="top", fontsize=7.8, color=MUTED, style="italic")

    # the two parents for comparison
    for name, yy in (("hapR", 0.42), ("RHS", 0.02)):
        f = FEATURES[name]
        ax.add_patch(Rectangle((0, yy), f["length"], 0.18,
                               facecolor=COL[name], edgecolor="none", alpha=0.55))
        ax.text(-1200, yy + 0.09, name, ha="right", va="center", fontsize=9, color=INK2)

    # mark toxin and integrase on each element
    for name, yy, h in (("hapA", 1.05, 0.22), ("hapR", 0.42, 0.18), ("RHS", 0.02, 0.18)):
        f = FEATURES[name]
        tl, ts, te = f["toxin"]
        ax.add_patch(Rectangle((ts, yy), te - ts, h, facecolor="none",
                               edgecolor=INK, lw=1.6))
        ax.text((ts + te) / 2, yy + h + 0.03, tl, ha="center", va="bottom",
                fontsize=7.5, color=INK)
        il, is_, ie = f["integrase"]
        ax.add_patch(Rectangle((is_, yy), max(ie - is_, 400), h, facecolor="none",
                               edgecolor=INK, lw=1.6, ls=(0, (2, 1.5))))
        ax.text((is_ + ie) / 2, yy + h + 0.03, f"int{il}", ha="center", va="bottom",
                fontsize=7.5, color=INK)

    ax.set_xlim(-4200, 41500); ax.set_ylim(-0.30, 1.72)
    ax.set_xlabel("position in phage genome (bp)", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=7.5, colors=MUTED, length=3)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#dcdbd5")


def panel_c(ax):
    """The junction sits at ONE position -- a discrete lineage, not a swarm."""
    ax.set_title("C.  The recombination junction is at a single fixed position",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")
    xs = sorted(BREAKPOINTS)
    ys = [BREAKPOINTS[x] for x in xs]
    ax.bar(xs, ys, width=0.85, color=[HAPA if x == 5 else "#b8b7b0" for x in xs])
    for x, yv in zip(xs, ys):
        ax.text(x, yv + 5, str(yv), ha="center", va="bottom", fontsize=9.5,
                color=INK if x == 5 else MUTED,
                fontweight="bold" if x == 5 else "normal")
    ax.set_xlabel("position of the junction in the hapR genome (kb)",
                  fontsize=8.5, color=INK2)
    ax.set_ylabel("reads spanning", fontsize=8.5, color=INK2)
    ax.set_xlim(0, 41); ax.set_ylim(0, 215)
    ax.tick_params(labelsize=7.5, colors=MUTED, length=3)
    ax.grid(True, axis="y", color="#ecebe6", lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_color("#dcdbd5")
    ax.text(24, 150,
            "166 of 185 split reads (90 %) break at the same point.\n"
            "Ongoing recombination would smear this across the genome;\n"
            "a single peak means one stable mosaic lineage.",
            fontsize=8.5, color=INK2, va="top")
    ax.text(24.2, 78, "reads that map partly to hapR and partly to the RHS phage",
            fontsize=7.8, color=MUTED, style="italic", va="top")


def main():
    fig = plt.figure(figsize=(11.5, 12.2), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.90, 1.20, 0.78], hspace=0.42)
    for f, ax in zip((panel_a, panel_b, panel_c),
                     [fig.add_subplot(gs[i]) for i in range(3)]):
        ax.set_facecolor(SURFACE)
        f(ax)

    legend = [Line2D([], [], marker="s", ls="", ms=10, mfc=RHS, mec="none",
                     label="RHS phage (integrase A)"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=HAPR, mec="none",
                     label="hapR — CdtB phage (integrase B)"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=HAPA, mec="none",
                     label="hapA — CdtB phage (integrase A), recombinant"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=COL["same"], mec="none",
                     label="parents identical — origin cannot be assigned")]
    fig.legend(handles=legend, loc="lower center", ncol=2, frameon=False,
               fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.004))

    fig.suptitle("S05 carries three APSE lineages across two integration sites",
                 fontsize=14, color=INK, y=0.988, fontweight="semibold")
    fig.text(0.5, 0.963,
             "Every cell carries hapR at ~187 kb. The ~805 kb site holds either the RHS phage or hapA — a stable mosaic of the other two.\n"
             "So three lineages coexist in the population while any single chromosome carries two prophages, not three.",
             ha="center", va="top", fontsize=10, color=INK2, linespacing=1.5)

    fig.subplots_adjust(left=0.10, right=0.97, top=0.878, bottom=0.075)
    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/S05_APSE_architecture.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
