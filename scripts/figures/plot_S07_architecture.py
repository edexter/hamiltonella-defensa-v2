#!/usr/bin/env python3
"""
Four-panel figure explaining the APSE architecture of sample S07 (strain H101).

Panel A  the two prophages and where they sit in the chromosome
Panel B  how the two phage genomes differ: one conserved backbone and two
         variable cassettes, one carrying the toxin and one carrying the
         integrase and its own att site
Panel C  what the 2024 reference actually was -- a chimera of the two
Panel D  recombination between the two phages, measured from raw reads

Every number plotted is measured. Sources are given per panel.

Usage:  python3 scripts/plot_S07_architecture.py
Output: results/figures/S07_APSE_architecture.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

# Validated categorical palette (light mode); CVD-checked, see dataviz reference
RHS = "#2a78d6"       # the RHS-carrying phage
CDTB = "#eb6834"      # the CdtB-carrying phage
CHROM = "#c9c8c1"     # bacterial chromosome
SAME = "#d8d7d0"      # the two phages are identical here
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"
COL = {"RHS": RHS, "CdtB": CDTB, "same": SAME}

# --- measured values --------------------------------------------------------
# Integration sites: split-read clustering of all 25,153 S07 reads against the
# chromosome plus both reconstructed phages. Coordinates refer to the delivered
# 2024 chromosome contig_02 (2,066,358 bp), mean depth 35.6x.
CHR_LEN = 2066358
CHR_DEPTH = 35.6
SITES = [
    dict(pos=757236, core=62, phage="RHS", lysogenic=88, excised=4,
         integrase="A", depth=2251),
    dict(pos=1878390, core=58, phage="CdtB", lysogenic=89, excised=1,
         integrase="B", depth=1184),
]

# Phage genomes, de novo from partitioned reads (Flye 2.9.6, 150x).
PHAGE = {
    "RHS":  dict(length=39306, toxin=(35637, 39306), toxin2=(1, 1356),
                 integrase=(17631, 18795), att=(17562, 17623)),
    "CdtB": dict(length=38993, toxin=(12693, 13671), toxin2=None,
                 integrase=(34197, 35221), att=(34040, 34097)),
}

# Segments of each phage with no counterpart in the other, from a whole-genome
# alignment of the two reconstructions (minimap2 -x asm10). Everything else
# aligns at 97.3-99.96 % identity.
UNIQUE = {
    "RHS":  [(34200, 39306), (0, 2186), (16493, 19256), (21835, 22567),
             (25666, 26099), (11533, 11660)],
    "CdtB": [(30092, 36009), (12340, 15611), (38589, 38993), (0, 1008),
             (4112, 4239), (25132, 25259)],
}
# The two cassettes that matter, named for what they carry.
CASSETTE = {
    "toxin":       {"RHS": (34200, 39306), "CdtB": (12340, 15611)},
    "integration": {"RHS": (16493, 19256), "CdtB": (30092, 36009)},
}

# Ancestry of the single APSE reference delivered in 2024 (contig_05,
# pH101_APSE, 38,642 bp), scored in 1 kb windows against both reconstructions.
# A window is assigned only when the identity margin exceeds 0.3 %.
DELIVERED_LEN = 38642
DELIVERED = [
    (0, "CdtB"), (1000, "CdtB"), (2000, "CdtB"), (3000, "CdtB"), (4000, "CdtB"),
    (5000, "same"), (6000, "same"), (7000, "same"), (8000, "same"),
    (9000, "CdtB"), (10000, "CdtB"), (11000, "CdtB"), (12000, "CdtB"),
    (13000, "same"), (14000, "same"), (15000, "same"), (16000, "same"),
    (17000, "same"), (18000, "same"), (19000, "same"), (20000, "same"),
    (21000, "RHS"), (22000, "same"), (23000, "RHS"), (24000, "RHS"),
    (25000, "RHS"), (26000, "RHS"), (27000, "same"),
    (28000, "CdtB"), (29000, "CdtB"), (30000, "CdtB"), (31000, "CdtB"),
    (32000, "CdtB"), (33000, "RHS"), (34000, "same"), (35000, "same"),
    (36000, "same"), (37000, "same"),
]

# Heterozygous sites per kb measured by mapping the sample's own reads back.
# S12, which genuinely carries one phage, sets the zero point.
DISCORDANCE = [("2024 pH101_APSE\n(single reference)", 6.94, "#b0402a"),
               ("v2 APSE-RHS", 0.0, RHS), ("v2 APSE-CdtB", 0.0, CDTB),
               ("S12 single-phage\ncontrol", 0.0, MUTED)]

# Read classification from parent-diagnostic 31-mers, 15,912 phage-derived
# reads. Alignment-free: k-mers are read from the raw sequence.
CLASSES = [("pure RHS", 9560, RHS), ("pure CdtB", 4913, CDTB),
           ("uninformative", 826, SAME), ("mosaic", 352, "#7a4fbd")]

# Crossover positions of the mosaic reads, in RHS coordinates, 2 kb bins.
# Each transition is the point where a read stops carrying one phage's
# diagnostic 31-mers and starts carrying the other's; 196 transitions in 352
# mosaic reads, falling at 57 distinct coordinate pairs.
CROSSOVERS = {2: 6, 6: 6, 8: 9, 10: 33, 12: 7, 16: 13, 18: 14, 20: 37,
              22: 63, 24: 5, 34: 3}


def panel_a(ax):
    ax.set_title("A.  Two prophages, two chromosomal sites, two integrase types",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")
    y = 0.55
    ax.add_patch(Rectangle((0, y), CHR_LEN, 0.14, facecolor=CHROM, edgecolor="none"))
    ax.text(CHR_LEN / 2, y - 0.06,
            f"Hamiltonella defensa chromosome, {CHR_LEN/1e6:.2f} Mb  ·  {CHR_DEPTH:.1f}× depth",
            ha="center", va="top", fontsize=8, color=MUTED)
    for s in SITES:
        x = s["pos"]
        L = PHAGE[s["phage"]]["length"]
        ax.add_patch(Rectangle((x - L * 1.6, y + 0.17), L * 3.2, 0.15,
                               facecolor=COL[s["phage"]], edgecolor="white", lw=1.2))
        ax.plot([x, x], [y + 0.14, y + 0.17], color=INK2, lw=1)
        ax.text(x, y + 0.345, f"APSE-{s['phage']}   {L:,} bp   integrase {s['integrase']}",
                ha="center", va="bottom", fontsize=9, color=INK, fontweight="medium")
        ax.text(x, y + 0.048, f"attB {x:,}", ha="center", va="center",
                fontsize=7.5, color=INK, fontweight="bold")
        ax.text(x, y - 0.20,
                f"{s['core']} bp att core, in a tRNA gene\n"
                f"{s['lysogenic']} reads cross into the phage\n"
                f"{s['excised']} read{'s' if s['excised']!=1 else ''} cross{'' if s['excised']!=1 else 'es'} the empty site\n"
                f"free phage at {s['depth']:,}× — {s['depth']/CHR_DEPTH:.0f}× the chromosome",
                ha="center", va="top", fontsize=8, color=INK2, linespacing=1.45)
    ax.text(CHR_LEN * 0.5, y + 0.60,
            "Each phage carries its own integrase and its own att site, so the two target different tRNA genes\n"
            "and do not compete. Both are integrated in essentially every cell.",
            ha="center", va="center", fontsize=8.5, color=INK2, style="italic",
            linespacing=1.5)
    ax.set_xlim(-160000, CHR_LEN + 160000)
    ax.set_ylim(-0.52, 1.02)
    ax.axis("off")


def panel_b(ax):
    ax.set_title("B.  One conserved backbone, two variable cassettes",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")
    rows = [("RHS", 0.62), ("CdtB", 0.16)]
    for name, y in rows:
        p = PHAGE[name]
        ax.add_patch(Rectangle((0, y), p["length"], 0.20, facecolor=SAME,
                               edgecolor="none"))
        for a, b in UNIQUE[name]:
            ax.add_patch(Rectangle((a, y), b - a, 0.20, facecolor=COL[name],
                                   edgecolor="none"))
        ax.text(-900, y + 0.10, f"APSE-{name}", ha="right", va="center",
                fontsize=9.5, color=INK2)
        ax.text(p["length"] + 700, y + 0.10, f"{p['length']:,} bp",
                ha="left", va="center", fontsize=8, color=MUTED)
        for lab, span, style in (("toxin", p["toxin"], "-"),
                                 ("toxin", p["toxin2"], "-"),
                                 (f"int {'A' if name=='RHS' else 'B'}",
                                  p["integrase"], (0, (2, 1.5))),
                                 ("attP", p["att"], "-")):
            if span is None:
                continue
            a, b = span
            ax.add_patch(Rectangle((a, y), max(b - a, 260), 0.20, facecolor="none",
                                   edgecolor=INK, lw=1.5, ls=style))
            if lab == "attP":
                ax.text((a + b) / 2, y - 0.035, lab, ha="center", va="top",
                        fontsize=7.5, color=INK)
            elif lab != "toxin" or span is p["toxin"]:
                ax.text((a + b) / 2, y + 0.235, lab, ha="center", va="bottom",
                        fontsize=7.5, color=INK)
    for cname, ends in CASSETTE.items():
        ra, rb = ends["RHS"]
        ca, cb = ends["CdtB"]
        ax.plot([(ra + rb) / 2, (ca + cb) / 2], [0.62, 0.36], color=MUTED,
                lw=0.9, ls=(0, (3, 2)))
        ax.text((ra + rb) / 2, 0.98,
                f"{cname} cassette", ha="center", va="bottom", fontsize=8.5,
                color=INK2, style="italic")
    ax.text(28000, 0.02,
            "grey = the two phages are identical here (28.1 kb of 39 kb, 97.3–99.96 % identity)",
            ha="center", va="top", fontsize=8, color=MUTED, style="italic")
    ax.set_xlim(-5200, 44500)
    ax.set_ylim(-0.22, 1.22)
    ax.set_xlabel("position in phage genome (bp)", fontsize=8.5, color=INK2)
    ax.tick_params(labelsize=7.5, colors=MUTED, length=3)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color("#dcdbd5")


def panel_c(ax):
    ax.set_title("C.  The APSE reference delivered in 2024 is a chimera of both phages",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")
    y = 0.42
    for start, call in DELIVERED:
        ax.add_patch(Rectangle((start, y), 1000, 0.30, facecolor=COL[call],
                               edgecolor="white", lw=0.7))
    ax.text(-800, y + 0.15, "pH101_APSE", ha="right", va="center", fontsize=9,
            color=INK2)
    ax.text(DELIVERED_LEN + 700, y + 0.15, f"{DELIVERED_LEN:,} bp", ha="left",
            va="center", fontsize=8, color=MUTED)
    n_r = sum(1 for _, c in DELIVERED if c == "RHS")
    n_c = sum(1 for _, c in DELIVERED if c == "CdtB")
    n_s = sum(1 for _, c in DELIVERED if c == "same")
    ax.text(DELIVERED_LEN / 2, y - 0.10,
            f"{n_c} of the 38 one-kilobase windows can only have come from the CdtB phage and {n_r} only from the RHS phage;\n"
            f"{n_s} are shared and cannot be told apart. One delivered reference, two sources.",
            ha="center", va="top", fontsize=8.5, color=INK2, linespacing=1.5)
    ax.set_xlim(-4200, 43500)
    ax.set_ylim(-0.40, 0.76)
    ax.axis("off")


def panel_c2(ax):
    labels = [d[0] for d in DISCORDANCE]
    vals = [d[1] for d in DISCORDANCE]
    cols = [d[2] for d in DISCORDANCE]
    ax.barh(range(len(vals))[::-1], vals, color=cols, height=0.62)
    for i, v in enumerate(vals):
        ax.text(v + 0.11, len(vals) - 1 - i, f"{v:.2f}", va="center",
                fontsize=8.5, color=INK if v else MUTED)
    ax.set_yticks(range(len(vals))[::-1])
    ax.set_yticklabels(labels, fontsize=8, color=INK2)
    ax.set_xlabel("heterozygous sites per kb when the sample's own reads are mapped back",
                  fontsize=8.5, color=INK2)
    ax.set_xlim(0, 8.2)
    ax.set_ylim(-0.7, 3.7)
    ax.tick_params(labelsize=7.5, colors=MUTED, length=3)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_color("#dcdbd5")
    ax.text(1.05, 1.85,
            "The delivered reference reports polymorphism that exists in neither phage:\n"
            "it is the difference between the two, read as variation within one.",
            fontsize=8.5, color=INK2, va="top", linespacing=1.5)


def panel_d(ax):
    ax.set_title("D.  The two phages recombine, but rarely — no dominant mosaic lineage",
                 fontsize=11, color=INK, loc="left", pad=14, fontweight="medium")
    tot = sum(n for _, n, _ in CLASSES)
    x = 0
    for lab, n, c in CLASSES:
        ax.add_patch(Rectangle((x, 0.30), n, 0.34, facecolor=c, edgecolor="white",
                               lw=1.2))
        if n / tot > 0.20:
            ax.text(x + n / 2, 0.47, f"{lab}\n{100*n/tot:.1f} %", ha="center",
                    va="center", fontsize=8.5, color="white")
        elif n / tot > 0.03:
            ax.text(x + n / 2, 0.66, f"{lab}  {100*n/tot:.1f} %", ha="center",
                    va="bottom", fontsize=8, color=INK2)
        x += n
    ax.text(0, 0.70, f"{tot:,} phage-derived reads, classified by parent-diagnostic 31-mers",
            ha="left", va="bottom", fontsize=8.5, color=MUTED, style="italic")
    ax.annotate("", xy=(tot * 0.995, 0.28), xytext=(tot * 0.93, 0.08),
                arrowprops=dict(arrowstyle="->", color=INK2, lw=1.1))
    ax.text(tot * 0.925, 0.05, "352 mosaic reads (2.2 %)", ha="right", va="top",
            fontsize=8.5, color=INK)
    ax.set_xlim(-tot * 0.01, tot * 1.01)
    ax.set_ylim(-0.15, 0.92)
    ax.axis("off")


def panel_d2(ax):
    xs = sorted(CROSSOVERS)
    ax.bar([k + 1 for k in xs], [CROSSOVERS[k] for k in xs], width=1.8,
           color="#7a4fbd")
    ax.set_xlabel("position of the crossover in the RHS phage genome (kb)",
                  fontsize=8.5, color=INK2)
    ax.set_ylabel("crossovers", fontsize=8.5, color=INK2)
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 92)
    ax.tick_params(labelsize=7.5, colors=MUTED, length=3)
    ax.grid(True, axis="y", color="#ecebe6", lw=0.6)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    for sp in ("bottom", "left"):
        ax.spines[sp].set_color("#dcdbd5")
    ax.text(0.5, 88,
            "196 crossovers at 57 distinct positions, spread over the shared backbone.\n"
            "In S05 one position carried 95 % of them, because there a recombinant had\n"
            "become a fixed lineage occupying an integration site. Nothing like that here.",
            fontsize=8.5, color=INK2, va="top", linespacing=1.5)


def main():
    fig = plt.figure(figsize=(11.5, 14.2), facecolor=SURFACE)
    gs = fig.add_gridspec(6, 1, height_ratios=[0.98, 1.00, 0.60, 0.58, 0.44, 0.80],
                          hspace=0.55)
    for f, ax in zip((panel_a, panel_b, panel_c, panel_c2, panel_d, panel_d2),
                     [fig.add_subplot(gs[i]) for i in range(6)]):
        ax.set_facecolor(SURFACE)
        f(ax)

    legend = [Line2D([], [], marker="s", ls="", ms=10, mfc=RHS, mec="none",
                     label="APSE-RHS (integrase A, RHS toxin)"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=CDTB, mec="none",
                     label="APSE-CdtB (integrase B, cdtB toxin)"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=SAME, mec="none",
                     label="the two phages are identical — origin cannot be assigned")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.004))

    fig.suptitle("S07 carries two complete APSE phages at two separate integration sites",
                 fontsize=14, color=INK, y=0.991, fontweight="semibold")
    fig.text(0.5, 0.9715,
             "Both are reconstructed here as circular genomes and spliced back into the chromosome. The 2024 delivery contained one APSE reference,\n"
             "assembled from the reads of both, which is why mapping against it produces polymorphism that exists in neither phage.",
             ha="center", va="top", fontsize=10, color=INK2, linespacing=1.5)

    fig.subplots_adjust(left=0.115, right=0.965, top=0.915, bottom=0.070)
    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/S07_APSE_architecture.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
