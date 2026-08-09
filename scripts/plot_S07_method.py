#!/usr/bin/env python3
"""
Two methods figures for the S07 two-phage reconstruction.

  S07_base_algorithm.png     the entry point — what you run on a new sample
  S07_refinement_rounds.png  what the three later rounds each changed, and why

Round 3 is the true beginning: it is the first round whose inputs are only the
sample's own pre-curation contigs and the sample's own reads. Rounds 1 and 2
were failed experiments and nothing from them feeds into it.

Every number is measured. Sources: assembly_info.txt for lengths and
circularity; pile sizes from the recruitment id lists; heterozygosity from
mapping a common 1,614-read sample back onto each assembly; window ancestry
from minimap2 against the finished pair.

Usage:  python3 scripts/plot_S07_method.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREY = "#c9c8c1"
GREEN = "#1baf7a"
RED = "#b0402a"
AMBER = "#b8860b"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"
BLUEBG = "#e2ecf9"
ORANGEBG = "#fce6dc"
PANEL = "#f4f2ec"

# ---- round 3, the base algorithm -------------------------------------------
CONTIGS = [
    ("contig_3",  7158, 2045, "integrase A",              "RHS",  BLUE),
    ("contig_4", 17114, 2221, "RHS toxin",                "RHS",  BLUE),
    ("contig_5", 38642, 1613, "integrase B + cdtB toxin",  "CdtB", ORANGE),
    ("contig_1", 10921,   20, "integrase B",              "CdtB", ORANGE),
    ("contig_6", 13248,   20, "cdtB toxin — duplicates contig_5", "unused", MUTED),
]

# ---- the three refinement rounds -------------------------------------------
REFINE = [
    dict(n="4", title="Feed the assemblies back",
         change="The reference stops being contig fragments and becomes the two\ngenomes round 3 produced.",
         why="Round 3's RHS reference covered only 62 % of that phage — 15 kb of it\n"
             "had no reference at all, so reads from there had nowhere correct to go.\n"
             "Round 3's assembly covers 98 %.",
         eff=[("RHS", "38,524 bp, not circular", "39,306 bp, CIRCULAR", True),
              ("RHS", "0.00 /kb", "0.00 /kb", True),
              ("CdtB", "4.68 /kb", "3.91 /kb", False)],
         lesson="Iteration cures an INCOMPLETE reference. RHS closes here and never\nchanges again. It does almost nothing for CdtB.",
         ok=None),
    dict(n="5", title="Stop duplicating ambiguous reads",
         change="Reads that match both references equally are dropped instead of being\nput into both piles.",
         why="While the reference still carried chimeric sequence, some reads scored\n"
             "ambiguous because the REFERENCE was wrong there, not because the two\n"
             "phages are identical there. Copying those into both piles spread the error.",
         eff=[("piles", "53 reads in both", "0 reads in both", True),
              ("CdtB", "7,007 reads", "6,372 reads", None),
              ("CdtB", "3.91 /kb", "0.94 /kb", True)],
         lesson="Buys purity at the cost of coverage. Affordable at 150×; it is also what\nshortens a recombinant lineage when coverage is thin.",
         ok=None),
    dict(n="6", title="Patch the reference, sort on k-mers",
         change="The 36 bases the reads contradicted are corrected IN THE\nREFERENCE, and reads are sorted on parent-diagnostic 31-mers\nrather than on alignment score.",
         why="The residual had moved to ~95 % in one 2.3 kb block — the signature of a\n"
             "minority allele taken, not of two lineages mixed. Alignment score was too\n"
             "blunt: 8 kb of shared backbone drowns twenty discriminating positions.",
         eff=[("CdtB", "0.94 /kb", "0.03 /kb", True),
              ("CdtB", "2 foreign windows", "0 foreign windows", True),
              ("both", "—", "99.8 % of reads keep their pile", True)],
         lesson="Chimerism needs an explicit repair; no number of extra rounds removes it.\nWhen the classification stops changing, the loop is done.",
         ok=True),
]

TRAJ_CDTB = [("3", 4.68), ("4", 3.91), ("5", 0.94), ("6", 0.03)]
TRAJ_RHS = [("3", 0.00), ("4", 0.00), ("5", 0.00), ("6", 0.00)]


def box(ax, x, y, w, h, fc, ec="none", lw=1.2, z=1):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.006,rounding_size=0.012",
                                facecolor=fc, edgecolor=ec, lw=lw, zorder=z))


def step(ax, x, y, n, title, sub):
    ax.text(x, y, n, fontsize=16, color=INK, fontweight="bold", va="top")
    ax.text(x + 0.024, y, title, fontsize=12, color=INK, fontweight="semibold", va="top")
    ax.text(x + 0.024, y - 0.038, sub, fontsize=8.8, color=INK2, va="top",
            linespacing=1.5)


# ============================ figure 1 ======================================
def base_algorithm():
    fig = plt.figure(figsize=(15.5, 10.4), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    fig.text(0.5, 0.977, "The base algorithm — what to run on a new sample that may carry two phages",
             ha="center", va="top", fontsize=16.5, color=INK, fontweight="semibold")
    fig.text(0.5, 0.941,
             "Shown as it ran on S07. Its only inputs are the sample's own pre-curation contigs and the sample's own reads: nothing is borrowed from another sample and nothing\n"
             "is carried over from an earlier attempt. The reference decides which pile each read joins and contributes no sequence — both genomes are assembled from scratch.",
             ha="center", va="top", fontsize=9.8, color=INK2, linespacing=1.55)

    # ---------------- step 1
    step(ax, 0.030, 0.878, "1", "Group the sample's own contigs by marker gene",
         "The pre-curation assembly already contains pieces of both phages. Two genes separate them and neither is ambiguous:\n"
         "the two integrase types share no 31-mer, and the two toxins are unrelated genes.")

    y = 0.762
    SC = 0.0000047
    for name, length, cov, marker, panel, col in CONTIGS:
        ax.text(0.045, y, name, fontsize=9.5, color=INK, ha="right", va="center",
                fontweight="medium")
        ax.add_patch(Rectangle((0.055, y - 0.010), length * SC, 0.021, facecolor=col,
                               edgecolor="none", alpha=1.0 if panel != "unused" else 0.35))
        ax.text(0.345, y, f"{length:,} bp   {cov:,}×", fontsize=8.3, color=MUTED,
                ha="right", va="center")
        ax.text(0.362, y, marker, fontsize=9.2,
                color=INK if panel != "unused" else MUTED, ha="left", va="center",
                fontweight="medium" if panel != "unused" else "normal")
        y -= 0.052

    box(ax, 0.575, 0.690, 0.150, 0.086, BLUEBG)
    ax.text(0.650, 0.752, "RHS reference", fontsize=10, color=INK, ha="center",
            fontweight="semibold")
    ax.text(0.650, 0.714, "contig_3 + contig_4\n24,272 bp", fontsize=8.8, color=INK2,
            ha="center", va="center", linespacing=1.4)
    box(ax, 0.575, 0.578, 0.150, 0.086, ORANGEBG)
    ax.text(0.650, 0.640, "CdtB reference", fontsize=10, color=INK, ha="center",
            fontweight="semibold")
    ax.text(0.650, 0.602, "contig_5 + contig_1\n49,563 bp", fontsize=8.8, color=INK2,
            ha="center", va="center", linespacing=1.4)
    for yy, col in ((0.733, BLUE), (0.621, ORANGE)):
        ax.add_patch(FancyArrowPatch((0.570, yy), (0.540, yy), arrowstyle="->",
                                     color=col, lw=1.6, mutation_scale=14))
    ax.text(0.650, 0.556,
            "Neither is a finished genome, and they are\nflawed in DIFFERENT ways:",
            fontsize=8.5, color=MUTED, ha="center", va="top", style="italic",
            linespacing=1.5)
    ax.text(0.650, 0.500, "RHS reference:   covers 62 %,  0 foreign windows",
            fontsize=8.6, color=BLUE, ha="center", va="top", fontweight="medium")
    ax.text(0.650, 0.474, "CdtB reference:  covers 99 %,  6 foreign windows",
            fontsize=8.6, color=ORANGE, ha="center", va="top", fontweight="medium")

    # ---------------- step 2
    step(ax, 0.760, 0.878, "2", "Sort every read",
         "Map all reads against BOTH references at once —\nnever one at a time.")
    rows = [("matches one reference uniquely (MAPQ ≥ 20)", "→ that pile", INK),
            ("matches both equally (MAPQ 0)", "→ BOTH piles", MUTED)]
    yy = 0.762
    for a, b, c in rows:
        ax.text(0.762, yy, a, fontsize=8.8, color=c, va="center")
        ax.text(0.762, yy - 0.026, b, fontsize=9, color=INK2, va="center",
                fontweight="medium")
        yy -= 0.070
    ax.text(0.762, 0.640, "RHS pile    5,547 reads", fontsize=9.6, color=BLUE,
            va="center", fontweight="semibold")
    ax.text(0.762, 0.612, "CdtB pile  10,356 reads", fontsize=9.6, color=ORANGE,
            va="center", fontweight="semibold")
    ax.text(0.762, 0.560,
            "About two-thirds of an APSE genome is identical\n"
            "between the two phages. A read from there cannot\n"
            "say where it came from and does not need to —\n"
            "both assemblies need those same bases.",
            fontsize=8.6, color=INK2, va="top", linespacing=1.55)

    # ---------------- step 3
    step(ax, 0.030, 0.400, "3", "Assemble each pile from scratch",
         "Flye 2.9.6 --pacbio-hifi --meta, thinned to ~150×. No base of the reference is copied.")
    for i, (nm, pile, down, ln, circ, col, bg) in enumerate([
            ("RHS pile", 5547, 700, "38,524 bp", "not circular", BLUE, BLUEBG),
            ("CdtB pile", 10356, 700, "39,308 bp", "circular  (+2 fragments)", ORANGE, ORANGEBG)]):
        yy = 0.288 - i * 0.120
        box(ax, 0.045, yy - 0.043, 0.118, 0.086, bg)
        ax.text(0.104, yy + 0.019, nm, fontsize=10, color=INK, ha="center",
                fontweight="semibold")
        ax.text(0.104, yy - 0.016, f"{pile:,} reads\nthinned to {down}", fontsize=8.6,
                color=INK2, ha="center", va="center", linespacing=1.4)
        ax.add_patch(FancyArrowPatch((0.205, yy), (0.170, yy), arrowstyle="->",
                                     color=col, lw=1.6, mutation_scale=14))
        ax.text(0.212, yy + 0.016, f"{ln}   ·   {circ}", fontsize=10, color=col,
                va="center", fontweight="semibold")

    # ---------------- step 4
    step(ax, 0.420, 0.400, "4", "Test each one against the reads that built it",
         "Map the pile back onto its own assembly and count positions where the reads disagree.")
    for i, (nm, val, txt, col, good) in enumerate([
            ("RHS", "0.00 /kb", "0 positions — meets the single-phage benchmark", BLUE, True),
            ("CdtB", "4.68 /kb", "184 positions at intermediate frequency", ORANGE, False)]):
        yy = 0.288 - i * 0.090
        box(ax, 0.435, yy - 0.033, 0.500, 0.066, "#eaf5ef" if good else "#faeeea")
        ax.text(0.452, yy, f"{nm}", fontsize=10.5, color=col, va="center",
                fontweight="semibold")
        ax.text(0.505, yy, val, fontsize=12, color=GREEN if good else RED, va="center",
                fontweight="bold")
        ax.text(0.585, yy, txt, fontsize=9.2, color=INK2, va="center")
        ax.text(0.925, yy, "✓" if good else "✗", fontsize=14,
                color=GREEN if good else RED, ha="right", va="center", fontweight="bold")

    ax.text(0.435, 0.108,
            "One genome is finished at the first pass. The other is not, and the shape of its residual says why:\n"
            "184 positions at intermediate allele frequency means two versions of the same region are still sharing a pile.",
            fontsize=9.2, color=INK2, va="top", linespacing=1.6)
    ax.text(0.435, 0.045, "→  Iterate: the next figure shows what the following three rounds each changed.",
            fontsize=9.8, color=AMBER, va="top", fontweight="medium")

    ax.plot([0.745, 0.745], [0.44, 0.905], color="#e4e2dc", lw=1)
    ax.plot([0.030, 0.960], [0.430, 0.430], color="#e4e2dc", lw=1)
    ax.plot([0.400, 0.400], [0.030, 0.412], color="#e4e2dc", lw=1)

    out = "results/figures/S07_base_algorithm.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"wrote {out}")


# ============================ figure 2 ======================================
def refinement():
    fig = plt.figure(figsize=(15.5, 10.0), facecolor=SURFACE)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    fig.text(0.5, 0.977, "What each round after the first one changed — and why the loop alone was not enough",
             ha="center", va="top", fontsize=16.5, color=INK, fontweight="semibold")
    fig.text(0.5, 0.941,
             "After the base algorithm, the RHS phage was already correct and only needed closing. The CdtB phage carried two further faults, and each needed a different\n"
             "intervention. Repeating the loop cures an incomplete reference; it does nothing at all for a chimeric one.",
             ha="center", va="top", fontsize=9.8, color=INK2, linespacing=1.55)

    W = 0.290
    for i, R in enumerate(REFINE):
        x = 0.035 + i * (W + 0.023)
        ax.text(x, 0.885, f"ROUND {R['n']}", fontsize=14, color=INK, va="top",
                fontweight="bold")
        ax.text(x, 0.848, R["title"], fontsize=11.5, color=INK, va="top",
                fontweight="medium")

        box(ax, x, 0.700, W, 0.115, PANEL)
        ax.text(x + 0.012, 0.798, "WHAT CHANGED", fontsize=8.4, color=MUTED,
                va="top", fontweight="bold")
        ax.text(x + 0.012, 0.772, R["change"], fontsize=9.0, color=INK, va="top",
                linespacing=1.55)

        ax.text(x, 0.672, "WHY", fontsize=8.4, color=MUTED, va="top", fontweight="bold")
        ax.text(x, 0.648, R["why"], fontsize=8.8, color=INK2, va="top", linespacing=1.6)

        ax.text(x, 0.520, "EFFECT", fontsize=8.4, color=MUTED, va="top",
                fontweight="bold")
        yy = 0.492
        for who, before, after, good in R["eff"]:
            c = BLUE if who == "RHS" else ORANGE if who == "CdtB" else INK2
            ax.text(x, yy, who, fontsize=8.8, color=c, va="center", fontweight="medium")
            ax.text(x + 0.048, yy, before, fontsize=8.8, color=MUTED, va="center")
            ax.text(x + 0.150, yy, "→", fontsize=9, color=MUTED, va="center")
            ax.text(x + 0.170, yy, after, fontsize=8.8,
                    color=GREEN if good else (INK2 if good is None else RED),
                    va="center", fontweight="medium" if good else "normal")
            yy -= 0.030

        box(ax, x, 0.300, W, 0.078, "#eaf5ef" if R["ok"] else "#f7f1e2")
        ax.text(x + 0.012, 0.362, R["lesson"], fontsize=8.9,
                color=GREEN if R["ok"] else AMBER, va="top", linespacing=1.6,
                fontweight="medium")

        if i < 2:
            ax.add_patch(FancyArrowPatch((x + W + 0.019, 0.60), (x + W + 0.004, 0.60),
                                         arrowstyle="->", color=MUTED, lw=1.6,
                                         mutation_scale=14))

    # ---- trajectory
    axT = fig.add_axes([0.075, 0.075, 0.400, 0.175])
    xs = list(range(4))
    axT.plot(xs, [v for _, v in TRAJ_CDTB], color=ORANGE, lw=2.4, marker="o", ms=7,
             mec="white", mew=1.4, zorder=3, label="CdtB")
    axT.plot(xs, [v for _, v in TRAJ_RHS], color=BLUE, lw=2.4, marker="o", ms=7,
             mec="white", mew=1.4, zorder=3, label="RHS")
    for i, (_, v) in enumerate(TRAJ_CDTB):
        axT.text(i, v + 0.28, f"{v:.2f}", ha="center", fontsize=9, color=ORANGE,
                 fontweight="semibold")
    axT.text(0, 0.30, "0.00 from the very first round", fontsize=8.6, color=BLUE,
             fontweight="medium")
    axT.set_xticks(xs); axT.set_xticklabels([f"round {n}" for n, _ in TRAJ_CDTB])
    axT.set_ylim(-0.45, 5.6); axT.set_xlim(-0.35, 3.35)
    axT.set_ylabel("heterozygous sites per kb", fontsize=9, color=INK2)
    axT.set_title("Heterozygosity, round by round", fontsize=10.5, color=INK,
                  loc="left", pad=8, fontweight="medium")
    axT.tick_params(labelsize=8.5, colors=MUTED, length=3)
    axT.grid(True, axis="y", color="#ecebe6", lw=0.6); axT.set_axisbelow(True)
    for s in ("top", "right"):
        axT.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        axT.spines[s].set_color("#dcdbd5")
    axT.set_facecolor(SURFACE)
    axT.legend(frameon=False, fontsize=9, labelcolor=INK2, loc="upper right")

    ax.text(0.520, 0.245,
            "How to know which fault you have, and when to stop",
            fontsize=11, color=INK, va="top", fontweight="semibold")
    ax.text(0.520, 0.205,
            "Watch heterozygosity each round. If it keeps falling, keep iterating. If it PLATEAUS, the loop\n"
            "has done all it can — for S07's CdtB that plateau was 4.68 → 3.91, an 18 % improvement, at\n"
            "the same round where RHS went to zero outright.\n\n"
            "Then read the SHAPE of what is left. Sites at intermediate frequency spread over a region\n"
            "mean two lineages still share a pile: add a lineage or tighten the sort. Sites near 95 %\n"
            "packed into a short block mean the assembly took a minority allele: patch those bases in\n"
            "the reference. Stop when the read classification stops changing.",
            fontsize=9.3, color=INK2, va="top", linespacing=1.65)

    out = "results/figures/S07_refinement_rounds.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE); plt.close(fig)
    print(f"wrote {out}")


def main():
    os.makedirs("results/figures", exist_ok=True)
    base_algorithm()
    refinement()


if __name__ == "__main__":
    main()
