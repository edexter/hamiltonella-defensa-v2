#!/usr/bin/env python3
"""
Three-panel figure explaining how the iterative read-partitioning method works.

Panel A  the cycle: what the reference is for, where reads go, what is tested
Panel B  the three things the test can say, plotted from real measurements
Panel C  what actually happened, round by round, in S07 and S05

Panels B and C plot measured values. Sources are noted per panel.

Usage:  python3 scripts/plot_method_iteration.py
Output: results/figures/method_iteration.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D

# Validated categorical palette (light mode); CVD-checked
BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREEN = "#1baf7a"
PURPLE = "#7a4fbd"
RED = "#b0402a"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"
BOX = "#eceae4"
SAME = "#d8d7d0"

# --- measured values --------------------------------------------------------
# Panel B, left: a correct assembly. S12 carries one phage genuinely; its APSE
# gives no discordant site at 1,382x. Our finished reconstructions match it.
# Panel B, middle: the S07 CdtB assembly of round 4. 36 sites, all above 95 %,
# confined to 28.9-31.2 kb -- the assembly holds the minority allele there,
# inherited from a chimeric template.
CHIMERIC = [(28869, .958), (28956, .975), (28959, .959), (29040, .959), (29046, .959),
            (29073, .959), (29199, .959), (30537, .960), (30550, .960), (30782, .959),
            (30789, .959), (30797, .957), (30800, .955), (30831, .957), (30833, .957),
            (30835, .957), (30836, .957), (30842, .957), (30851, .956), (30855, .957),
            (30858, .957), (30861, .957), (30863, .957), (30866, .957), (30867, .957),
            (30875, .957), (30884, .956), (30885, .968), (30890, .957), (30893, .957),
            (31152, .968), (31154, .967), (31155, .968), (31166, .978), (31206, .968),
            (31214, .967)]
# Panel B, right: the S05 CdtB assembly when only two lineages were assumed.
# 126 sites at intermediate frequency, drifting smoothly -- two haplotypes in
# one pile, not an assembly error.
TWO_LINEAGES = [(14768, .275), (16835, .581), (16847, .572), (16866, .582), (16871, .570),
                (16874, .584), (17009, .625), (17327, .635), (17496, .643), (17499, .644),
                (17505, .669), (17511, .637), (17514, .639), (17517, .634), (17520, .637),
                (17526, .641), (17532, .644), (17535, .644), (17538, .642), (17552, .640),
                (17555, .638), (17561, .635), (17562, .642), (17571, .656), (17580, .613),
                (17581, .619), (17582, .612), (17584, .626), (17586, .630), (17592, .628),
                (17593, .627), (17594, .631), (17595, .638), (17598, .639), (17601, .637),
                (17606, .635), (17610, .636), (17612, .638), (17616, .639), (17619, .635),
                (17622, .635), (17843, .668), (17862, .659), (17894, .655), (17896, .650),
                (17906, .653), (17907, .652), (17912, .648), (17914, .653), (17915, .651),
                (17920, .648), (17924, .650), (17927, .655), (18117, .649), (18118, .659),
                (18119, .653), (18125, .657), (18159, .690), (18188, .711), (18195, .694),
                (18197, .690), (18200, .697), (18203, .692), (18233, .697), (18236, .695),
                (18252, .697), (18305, .694), (18314, .701), (18320, .707), (18328, .709),
                (18329, .701), (18335, .702), (18516, .719), (18610, .722), (18613, .727),
                (18658, .732), (18772, .743), (18782, .739), (18783, .732), (18787, .741),
                (18808, .744), (18844, .751), (18856, .753), (18892, .757), (18894, .759),
                (18898, .756), (18903, .754), (18907, .753), (18913, .752), (19240, .778),
                (19253, .780), (19765, .805), (19787, .804), (20033, .808), (20412, .816),
                (20540, .822), (21137, .845), (21188, .851), (21206, .850), (21364, .853),
                (21473, .857), (21632, .858), (21665, .855), (21923, .865), (22164, .863),
                (27007, .779), (27059, .772), (27064, .774), (27085, .770), (27088, .768),
                (27120, .768), (27128, .767), (27144, .765), (27145, .763), (27157, .763),
                (27160, .761), (27163, .765), (27166, .763), (27174, .765), (27175, .764),
                (27181, .763), (27193, .755), (27200, .756), (27201, .759), (27208, .763),
                (27225, .760)]

# Panel C: heterozygosity (sites per kb) and assembly length per round.
# "-" means the round produced nothing usable for that lineage.
S07_TRACE = [
    ("1", "clean genomes borrowed\nfrom other samples",
     [("RHS", None, "56,981 bp, BOTH integrases")], "missing lineage"),
    ("2", "one template, permissive\nread recruitment",
     [("RHS", None, "78,792 bp dimer")], "sorting not competitive"),
    ("3-4", "this sample's own contigs,\ncompetitive sorting",
     [("RHS", 0.00, "39,306 bp"), ("CdtB", 2.72, "39,299 bp")], "chimeric patch in CdtB"),
    ("5", "round-4 assemblies,\npatched, k-mer sorting",
     [("RHS", 0.00, "39,306 bp"), ("CdtB", 0.00, "38,993 bp")], "converged, 99.8 %"),
]
S05_TRACE = [
    ("earlier", "two lineages assumed",
     [("RHS", 0.13, "39,287 bp"), ("CdtB", 4.16, "38,979 bp")], "missing lineage -> hapA"),
    ("1", "three lineages,\ncontainment sorting",
     [("RHS", 0.13, "39,288 bp"), ("hapR", 0.71, "39,286 bp"), ("hapA", 0.66, "36,345 bp")],
     "chimeric patches"),
    ("2", "round-1 assemblies",
     [("RHS", 0.13, "39,289 bp"), ("hapR", 0.00, "38,980 bp"), ("hapA", 0.00, "36,248 bp")],
     "SHIPPED"),
    ("3", "round-2 assemblies",
     [("RHS", 0.13, "39,279 bp"), ("hapR", 0.00, "38,980 bp"), ("hapA", 0.00, "33,092 bp")],
     "hapA shrinking - stop"),
]
LCOL = {"RHS": BLUE, "CdtB": ORANGE, "hapR": ORANGE, "hapA": GREEN}


def box(ax, x, y, w, h, text, fc=BOX, ec="none", fs=8.5, weight="normal", tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                                facecolor=fc, edgecolor=ec, lw=1.2))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, linespacing=1.45, fontweight=weight)


def arrow(ax, xy, xytext, label=None, lx=0, ly=0, color=INK2, fs=8, style="->"):
    ax.add_patch(FancyArrowPatch(xytext, xy, arrowstyle=style, color=color, lw=1.3,
                                 mutation_scale=13, shrinkA=2, shrinkB=2))
    if label:
        ax.text((xy[0] + xytext[0]) / 2 + lx, (xy[1] + xytext[1]) / 2 + ly, label,
                ha="center", va="center", fontsize=fs, color=color, linespacing=1.4)


def panel_a(ax):
    ax.set_title("A.  The cycle — the reference sorts reads, it never contributes sequence",
                 fontsize=11.5, color=INK, loc="left", pad=12, fontweight="medium")

    box(ax, 0.01, 0.60, 0.155, 0.22,
        "REFERENCE\nthis sample's own\ncontigs, grouped by\nintegrase and toxin",
        fc="#dfe9f7", fs=8.2, weight="medium")
    arrow(ax, (0.215, 0.71), (0.168, 0.71))

    box(ax, 0.218, 0.60, 0.135, 0.22,
        "compare every\nread against\nboth groups", fs=8.2)

    # three destinations
    arrow(ax, (0.415, 0.90), (0.357, 0.76))
    arrow(ax, (0.415, 0.71), (0.357, 0.71))
    arrow(ax, (0.415, 0.52), (0.357, 0.66))
    box(ax, 0.418, 0.83, 0.20, 0.135, "carries lineage-A markers\n→  pile A",
        fc="#dfe9f7", fs=8.2)
    box(ax, 0.418, 0.645, 0.20, 0.135, "carries lineage-B markers\n→  pile B",
        fc="#fbe3d8", fs=8.2)
    box(ax, 0.418, 0.455, 0.20, 0.135,
        "carries NO markers\n(≈ 2/3 of the genome)\n→  BOTH piles", fc=SAME, fs=8.2,
        weight="medium")

    arrow(ax, (0.665, 0.71), (0.622, 0.71))
    box(ax, 0.668, 0.60, 0.145, 0.22,
        "assemble each\npile DE NOVO\n(no reference base\nreaches the output)",
        fc="#dff0e6", fs=8.2, weight="medium")

    arrow(ax, (0.86, 0.71), (0.816, 0.71))
    box(ax, 0.863, 0.60, 0.13, 0.22,
        "TEST\nmap the pile's own\nreads back and\ncount disagreements",
        fc="#f6e7c9", fs=8.2, weight="medium")

    # verdicts
    arrow(ax, (0.966, 0.505), (0.966, 0.595), color=GREEN)
    ax.text(0.966, 0.462, "0.00 per kb\n→  done", ha="center", va="center", fontsize=8.6,
            color=GREEN, fontweight="medium", linespacing=1.4)

    # feedback loop
    # feedback loop, routed as a right-angle path so it cannot cross any label
    ax.plot([0.882, 0.882], [0.575, 0.30], color=RED, lw=1.6, solid_capstyle="round")
    ax.plot([0.882, 0.088], [0.30, 0.30], color=RED, lw=1.6, solid_capstyle="round")
    ax.add_patch(FancyArrowPatch((0.088, 0.575), (0.088, 0.30), arrowstyle="->",
                                 color=RED, lw=1.6, mutation_scale=15, shrinkA=0, shrinkB=0))
    ax.text(0.47, 0.375,
            "otherwise: repair the REFERENCE and go round again\n"
            "the assembly you just built is a better reference than the one you started with",
            ha="center", va="center", fontsize=9, color=RED, linespacing=1.5,
            fontweight="medium")
    ax.text(0.47, 0.185,
            "Round 1's reference covered only 62 % of one S07 phage (two fragments, 24 kb of 39 kb).\n"
            "Round 1's assembly covers 100 %. That is what each turn of the loop buys.",
            ha="center", va="center", fontsize=8.3, color=INK2, linespacing=1.5,
            style="italic")
    ax.set_xlim(-0.01, 1.02); ax.set_ylim(0.13, 1.0)
    ax.axis("off")


def sig(ax, pts, title, colour, note, ymax=1.0):
    if pts:
        ax.scatter([p[0] / 1000 for p in pts], [p[1] for p in pts], s=13, color=colour,
                   zorder=3, edgecolor="none")
    ax.set_xlim(0, 40); ax.set_ylim(0, 1.03)
    ax.set_xlabel("position in phage (kb)", fontsize=8, color=INK2)
    ax.tick_params(labelsize=7.5, colors=MUTED, length=3)
    ax.grid(True, color="#eceae4", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color("#dcdbd5")
    ax.set_facecolor(SURFACE)
    ax.set_title(title, fontsize=9.2, color=INK, loc="left", pad=8, fontweight="medium")
    ax.text(0.5, -0.40, note, transform=ax.transAxes, ha="center", va="top",
            fontsize=8.2, color=INK2, linespacing=1.5)


def panel_c(ax, trace, sample, x0):
    ax.text(x0 + 0.19, 1.05, sample, fontsize=10.5, color=INK, fontweight="semibold",
            ha="center")
    y = 0.94
    for rnd, what, res, verdict in trace:
        label = rnd if rnd == "earlier" else f"round {rnd}"
        ax.text(x0 - 0.045, y - 0.030, label, fontsize=8.6, color=INK,
                ha="right", va="center", fontweight="medium")
        ax.text(x0, y, what, fontsize=8, color=INK2, ha="left", va="top", linespacing=1.35)
        yy = y - (0.062 if "\n" in what else 0.036)
        for name, het, length in res:
            c = LCOL[name]
            if het is None:
                ax.text(x0 + 0.012, yy, f"✗  {name}: {length}", fontsize=8.3, color=RED,
                        ha="left", va="top")
            else:
                mark = "●" if het == 0 else "○"
                ax.text(x0 + 0.012, yy, f"{mark}  {name}  {length}   {het:.2f} /kb",
                        fontsize=8.3, color=c if het == 0 else INK2, ha="left", va="top",
                        fontweight="medium" if het == 0 else "normal")
            yy -= 0.040
        vc = GREEN if verdict == "SHIPPED" else (RED if ("missing" in verdict or "not" in verdict) else INK2)
        ax.text(x0 + 0.012, yy - 0.002, verdict, fontsize=8.1, color=vc, ha="left",
                va="top", style="italic", fontweight="medium" if verdict == "SHIPPED" else "normal")
        y = yy - 0.070
        if y > 0.02:
            ax.annotate("", xy=(x0 - 0.028, y + 0.030), xytext=(x0 - 0.028, y + 0.062),
                        arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.1))


def main():
    fig = plt.figure(figsize=(12.4, 14.2), facecolor=SURFACE)
    gs = fig.add_gridspec(3, 1, height_ratios=[0.66, 0.80, 1.18], hspace=0.40)

    ax_a = fig.add_subplot(gs[0]); ax_a.set_facecolor(SURFACE); panel_a(ax_a)

    ax_b = fig.add_subplot(gs[1]); ax_b.axis("off")
    ax_b.set_title("B.  What the test can say — three signatures, each calling for a different repair",
                   fontsize=11.5, color=INK, loc="left", pad=12, fontweight="medium")
    subs = [ax_b.inset_axes([0.035, 0.28, 0.27, 0.60]),
            ax_b.inset_axes([0.375, 0.28, 0.27, 0.60]),
            ax_b.inset_axes([0.715, 0.28, 0.27, 0.60])]
    sig(subs[0], [], "correct — nothing to do", GREEN,
        "No position where reads disagree.\nThe single-phage control S12 gives this\nat 1,382× coverage, so zero is real.")
    sig(subs[1], CHIMERIC, "one bad patch — repair and re-sort", PURPLE,
        "S07 CdtB, round 4. 36 sites, ALL above 95 %,\nin one 2.3 kb block: the assembly took the\nminority allele. The reads say which bases to fix.")
    sig(subs[2], TWO_LINEAGES, "two lineages in one pile — add a lineage", RED,
        "S05 CdtB, two lineages assumed. 126 sites at\nINTERMEDIATE frequency, drifting smoothly.\nNo number of rounds fixes this.")
    for s in subs:
        s.set_ylabel("fraction of reads\ndisagreeing", fontsize=7.8, color=INK2)
        s.axhspan(0.30, 0.90, color="#f2efe8", zorder=0)

    ax_c = fig.add_subplot(gs[2]); ax_c.set_facecolor(SURFACE)
    ax_c.set_title("C.  What actually happened, round by round",
                   fontsize=11.5, color=INK, loc="left", pad=12, fontweight="medium")
    ax_c.text(0.5, 1.135,
              "●  meets the single-phage benchmark      ○  measurable disagreement remains      ✗  failed outright",
              ha="center", va="top", fontsize=8.4, color=INK2)
    panel_c(ax_c, S07_TRACE, "S07  (strain H101)", 0.10)
    panel_c(ax_c, S05_TRACE, "S05  (strain H76)", 0.58)
    ax_c.plot([0.505, 0.505], [0.01, 1.06], color="#e4e2dc", lw=1)
    ax_c.set_xlim(0, 1); ax_c.set_ylim(-0.01, 1.16)
    ax_c.axis("off")

    legend = [Line2D([], [], marker="s", ls="", ms=10, mfc="#dfe9f7", mec="none",
                     label="RHS-carrying lineage"),
              Line2D([], [], marker="s", ls="", ms=10, mfc="#fbe3d8", mec="none",
                     label="CdtB-carrying lineage"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=SAME, mec="none",
                     label="shared between lineages — goes to both piles")]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=9.5, labelcolor=INK2, bbox_to_anchor=(0.5, 0.004))

    fig.suptitle("Separating two collapsed phage genomes: sort the reads, assemble each pile, test, repeat",
                 fontsize=14, color=INK, y=0.991, fontweight="semibold")
    fig.text(0.5, 0.9755,
             "Two APSE phages in one sample are 97–99 % identical, so an assembler merges them. The reference below is used only to decide which pile each read\n"
             "belongs in; every base of every delivered genome comes from the sample's own reads. The starting guess is imperfect by design — the test is what catches it.",
             ha="center", va="top", fontsize=9.8, color=INK2, linespacing=1.55)

    fig.subplots_adjust(left=0.055, right=0.975, top=0.916, bottom=0.050)
    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/method_iteration.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
