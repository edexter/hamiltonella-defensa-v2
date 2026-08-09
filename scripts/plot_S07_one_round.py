#!/usr/bin/env python3
"""
One page, one sample, one round: how the two APSE phages of S07 were separated.

Four steps, left to right and top to bottom:
  1  build the reference from S07's own contigs, grouped by two marker genes
  2  sort every read into a pile
  3  assemble each pile from scratch
  4  test the result against the reads that built it

Every number shown is measured. Contig lengths and coverages come from the
pre-curation Flye assembly (assembly_info.txt); marker assignments from
minimap2 placement of data/reference/integrases.fasta and toxins_CDS.fasta;
read counts from the k-mer classification; assembly sizes from Flye; the test
from mapping each pile back onto its own assembly.

Usage:  python3 scripts/plot_S07_one_round.py
Output: results/figures/S07_one_round.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

BLUE = "#2a78d6"      # the RHS-carrying phage
ORANGE = "#eb6834"    # the CdtB-carrying phage
GREY = "#c9c8c1"      # shared / unassignable
GREEN = "#1baf7a"
RED = "#b0402a"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"
BLUEBG = "#e2ecf9"
ORANGEBG = "#fce6dc"
GREYBG = "#eeece6"

# --- step 1: the five APSE-derived contigs of the pre-curation S07 assembly --
# name, length, coverage, marker text, which panel, colour
CONTIGS = [
    ("contig_3",  7158, 2045, "integrase A",            "RHS",  BLUE),
    ("contig_4", 17114, 2221, "RHS toxin",              "RHS",  BLUE),
    ("contig_5", 38642, 1613, "integrase B + cdtB toxin", "CdtB", ORANGE),
    ("contig_1", 10921,   20, "integrase B",            "CdtB", ORANGE),
    ("contig_6", 13248,   20, "cdtB toxin — duplicates contig_5", "unused", MUTED),
]

# --- step 2: where the 15,912 phage-derived reads went ----------------------
SORT = [
    ("only RHS-diagnostic 31-mers", 10006, BLUE,   "→ RHS pile"),
    ("only CdtB-diagnostic 31-mers", 5261, ORANGE, "→ CdtB pile"),
    ("no diagnostic 31-mer at all",   261, GREY,   "→ BOTH piles"),
    ("both, in separate blocks\n(recombinant molecules)", 384, RED, "→ held out"),
]
PILE_RHS, PILE_CDTB = 10267, 5522
DOWN_RHS, DOWN_CDTB = 662, 663
DEPTH_RHS, DEPTH_CDTB = 147, 146


def box(ax, x, y, w, h, fc, ec="none", lw=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.006,rounding_size=0.012",
                                facecolor=fc, edgecolor=ec, lw=lw, zorder=1))


def step_header(ax, x, y, n, title, sub):
    ax.text(x, y, f"{n}", fontsize=17, color=INK, fontweight="bold", va="top")
    ax.text(x + 0.026, y, title, fontsize=12.5, color=INK, fontweight="semibold",
            va="top")
    ax.text(x + 0.026, y - 0.035, sub, fontsize=9, color=INK2, va="top",
            linespacing=1.45)


def main():
    fig = plt.figure(figsize=(15.5, 10.6), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.set_facecolor(SURFACE)

    fig.text(0.5, 0.975,
             "How the two APSE phages in sample S07 were separated — one round",
             ha="center", va="top", fontsize=16.5, color=INK, fontweight="semibold")
    fig.text(0.5, 0.940,
             "The two phages are 97–99 % identical, so an assembler merges them into one. The reference built in step 1 is used only to decide which pile each read belongs in. "
             "It contributes no sequence:\nboth genomes in step 3 are assembled from scratch out of S07's own reads.",
             ha="center", va="top", fontsize=10, color=INK2, linespacing=1.5)

    # ================= STEP 1 =================
    step_header(ax, 0.030, 0.880, "1", "Build the reference from S07's own contigs",
                "The assembly made before curation already contains pieces of both phages — it just doesn't know it.\n"
                "Two genes sort them, and neither is ambiguous: the two integrase types share no 31-mer, and the two toxins are unrelated.")

    y = 0.760
    SC = 0.0000047          # bp -> axis width
    for name, length, cov, marker, panel, col in CONTIGS:
        ax.text(0.045, y, name, fontsize=10, color=INK, ha="right", va="center",
                fontweight="medium")
        ax.add_patch(Rectangle((0.055, y - 0.011), length * SC, 0.022,
                               facecolor=col, edgecolor="none",
                               alpha=1.0 if panel != "unused" else 0.35))
        ax.text(0.345, y, f"{length:,} bp   {cov:,}×",
                fontsize=8.5, color=MUTED, ha="right", va="center")
        ax.text(0.362, y, marker, fontsize=9.5, color=INK if panel != "unused" else MUTED,
                ha="left", va="center", linespacing=1.35,
                fontweight="medium" if panel != "unused" else "normal")
        y -= 0.055

    ax.text(0.030, 0.505,
            "Coverage corroborates but never overrules a marker. The chromosome sits at 34×, so contig_3\n"
            "and contig_4 at ~2,100× are free phage particles, while contig_1 at 20× is a single integrated copy.",
            fontsize=8.8, color=MUTED, va="top", style="italic", linespacing=1.55)

    # the two reference panels
    box(ax, 0.575, 0.690, 0.140, 0.090, BLUEBG)
    ax.text(0.645, 0.753, "RHS reference panel", fontsize=10, color=INK,
            ha="center", fontweight="semibold")
    ax.text(0.645, 0.714, "contig_3 + contig_4\n24,272 bp", fontsize=9, color=INK2,
            ha="center", va="center", linespacing=1.4)
    box(ax, 0.575, 0.565, 0.140, 0.090, ORANGEBG)
    ax.text(0.645, 0.628, "CdtB reference panel", fontsize=10, color=INK,
            ha="center", fontweight="semibold")
    ax.text(0.645, 0.589, "contig_5 + contig_1\n49,563 bp", fontsize=9, color=INK2,
            ha="center", va="center", linespacing=1.4)
    for yy, col in ((0.735, BLUE), (0.610, ORANGE)):
        ax.add_patch(FancyArrowPatch((0.570, yy), (0.540, yy), arrowstyle="->",
                                     color=col, lw=1.6, mutation_scale=14))
    ax.text(0.645, 0.545,
            "Neither panel is a finished genome — the RHS\npanel is two fragments covering 62 % of that\nphage. That is fine: its only job is to tell\nreads apart.",
            fontsize=8.6, color=MUTED, ha="center", va="top", style="italic",
            linespacing=1.55)

    # ================= STEP 2 =================
    step_header(ax, 0.750, 0.880, "2", "Sort every read",
                "Compare each read's 31-mers against\nboth reference panels.")

    y = 0.760
    for label, n, col, dest in SORT:
        ax.add_patch(Rectangle((0.768, y - 0.010), 0.013, 0.020, facecolor=col,
                               edgecolor="none"))
        ax.text(0.788, y + 0.006, label, fontsize=9, color=INK, va="center",
                linespacing=1.3)
        ax.text(0.788, y - 0.022, f"{n:,} reads   {dest}", fontsize=9, color=col,
                va="center", fontweight="medium")
        y -= 0.062

    ax.text(0.768, 0.520,
            "The grey case matters most. About two-thirds\n"
            "of an APSE genome is identical between the\n"
            "two phages, so a read from there cannot say\n"
            "where it came from — and need not. It goes to\n"
            "both piles: both assemblies need those bases.",
            fontsize=8.8, color=INK2, va="top", linespacing=1.55)

    # ================= STEP 3 =================
    step_header(ax, 0.030, 0.430, "3", "Assemble each pile from scratch",
                "No base of the reference is copied. Flye 2.9.6 --pacbio-hifi --meta, on reads only.")

    for i, (name, pile, down, depth, length, col, bg) in enumerate([
            ("RHS pile",  PILE_RHS,  DOWN_RHS,  DEPTH_RHS,  39306, BLUE,   BLUEBG),
            ("CdtB pile", PILE_CDTB, DOWN_CDTB, DEPTH_CDTB, 38993, ORANGE, ORANGEBG)]):
        yy = 0.305 - i * 0.135
        box(ax, 0.045, yy - 0.048, 0.115, 0.096, bg)
        ax.text(0.1025, yy + 0.022, name, fontsize=10.5, color=INK, ha="center",
                fontweight="semibold")
        ax.text(0.1025, yy - 0.016, f"{pile:,} reads\nthinned to {down} ({depth}×)",
                fontsize=9, color=INK2, ha="center", va="center", linespacing=1.4)
        ax.add_patch(FancyArrowPatch((0.205, yy), (0.168, yy), arrowstyle="->",
                                     color=col, lw=1.6, mutation_scale=14))
        ax.text(0.212, yy + 0.020, f"{length:,} bp, circular", fontsize=10.5,
                color=col, va="center", fontweight="semibold")
        ax.text(0.212, yy - 0.018,
                "one toxin, one integrase,\ncorrectly paired", fontsize=8.8,
                color=INK2, va="center", linespacing=1.4)

    ax.text(0.045, 0.062,
            "Thinning is practical, not statistical: the untrimmed pile runs to\nthousands of fold and crashed the assembler.",
            fontsize=8.6, color=MUTED, va="top", style="italic", linespacing=1.5)

    # ================= STEP 4 =================
    step_header(ax, 0.430, 0.430, "4", "Test it against the reads that built it",
                "Map each pile back onto its own assembly and\ncount positions where the reads disagree.")

    box(ax, 0.450, 0.205, 0.245, 0.105, "#eaf5ef")
    ax.text(0.5725, 0.288, "0 positions   ·   0.00 per kb", fontsize=13, color=GREEN,
            ha="center", fontweight="semibold")
    ax.text(0.5725, 0.240,
            "for both phages. A clonal genome should give exactly this,\nand sample S12 — which genuinely carries one phage — does,\nat 1,382× coverage. So zero is a real result, not a threshold.",
            fontsize=9, color=INK2, ha="center", va="center", linespacing=1.5)

    ax.text(0.450, 0.178,
            "For comparison, the single APSE reference delivered in 2024\n"
            "gives 6.94 per kb on the same reads — not polymorphism, but\n"
            "the difference between two phages read as variation within one.",
            fontsize=9.2, color=RED, va="top", linespacing=1.55)

    ax.text(0.450, 0.092,
            "Two further checks passed: each genome carries exactly one\n"
            "toxin and one integrase of the expected type, and the att site\n"
            "where each phage integrates matches the chromosome exactly.",
            fontsize=8.8, color=INK2, va="top", linespacing=1.55)

    # ---- the loop back ----
    ax.add_patch(FancyArrowPatch((0.795, 0.310), (0.795, 0.425),
                                 arrowstyle="->", color=INK2, lw=1.6,
                                 mutation_scale=15))
    ax.text(0.812, 0.375,
            "If the test had NOT returned zero,\nthe new assemblies would replace\nthe step-1 panels and the round\nwould repeat.",
            fontsize=9, color=INK2, va="center", linespacing=1.6)
    ax.text(0.812, 0.250,
            "Here it converged. Re-sorting the\nreads against the new assemblies\nmoved almost nothing: 99.8 % kept\nthe same pile.",
            fontsize=9, color=INK2, va="center", linespacing=1.6)

    # dividers
    ax.plot([0.730, 0.730], [0.500, 0.905], color="#e4e2dc", lw=1)
    ax.plot([0.030, 0.975], [0.470, 0.470], color="#e4e2dc", lw=1)
    ax.plot([0.415, 0.415], [0.045, 0.455], color="#e4e2dc", lw=1)
    ax.plot([0.772, 0.772], [0.045, 0.455], color="#e4e2dc", lw=1)

    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/S07_one_round.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
