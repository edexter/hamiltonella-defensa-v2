#!/usr/bin/env python3
"""
Nucleotide dot-plots comparing complete APSE bacteriophage genomes.

Produces a two-panel figure:

  Panel A  two phages carrying the SAME toxin family (CdtB)      -- the control
  Panel B  two phages carrying DIFFERENT toxin families          -- the comparison

Both pairs are colinear along a shared backbone. The difference between the
panels is structural rather than a matter of overall similarity: the same-toxin
pair aligns essentially end to end, whereas the different-toxin pair has a
multi-kilobase segment that does not align at all, and that segment is where the
toxin coding sequence sits.

Note that whole-genome k-mer sharing is a poor summary here -- both pairs share
roughly 70 % of their 25-mers, but for opposite reasons. The same-toxin pair is
colinear with scattered substitutions across its whole length; the
different-toxin pair is near-identical where it aligns but carries a large
unaligned block. Alignment coverage and per-block identity are reported instead.

Usage:  python3 scripts/plot_APSE_dotplot.py
Output: results/figures/APSE_dotplot.png
"""

import collections
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Validated categorical palette (light mode) -- see docs/, dataviz reference
FWD = "#2a78d6"          # forward-strand match
REV = "#eb6834"          # reverse-strand match
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8880"
SURFACE = "#fcfcfb"
TOXIN_BAND = "#eda100"   # slot 4, used only to mark toxin-gene extents

K = 25                   # k-mer length: long enough to be locus-specific
STEP = 2                 # query stride, purely to keep the point count sane

GENOMES = "docs/files_from_giacomo"
TOXINS = "data/reference/toxins_CDS.fasta"


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def dotplot_points(a, b):
    """Return forward and reverse (x, y) match coordinates between two sequences."""
    index = collections.defaultdict(list)
    for i in range(len(b) - K + 1):
        index[b[i:i + K]].append(i)
    fwd_x, fwd_y, rev_x, rev_y = [], [], [], []
    for i in range(0, len(a) - K + 1, STEP):
        mer = a[i:i + K]
        for j in index.get(mer, ()):
            fwd_x.append(i); fwd_y.append(j)
        for j in index.get(revcomp(mer), ()):
            rev_x.append(i); rev_y.append(j)
    return (fwd_x, fwd_y), (rev_x, rev_y)


def locate(subject, query):
    """Find a toxin CDS inside a phage genome, either strand."""
    for seq, flip in ((query, False), (revcomp(query), True)):
        pos = subject.find(seq)
        if pos >= 0:
            return pos, pos + len(seq), flip
    return None


def orient_to(reference, target, anchor=4000):
    """
    Rotate and if necessary reverse-complement a circular genome so that it
    starts at the same locus as the reference.

    These phage genomes are circular, so the coordinate at which each assembly
    was linearised is arbitrary. Left uncorrected, that offset displaces the
    diagonal of a dot-plot and can split it into two segments, which reads as a
    rearrangement when nothing has actually moved.

    The offset is found by taking k-mers from the start of the reference,
    locating them in the target on both strands, and using the modal implied
    offset -- which tolerates the substitutions that make an exact match of a
    long anchor sequence fail.
    """
    index = collections.defaultdict(list)
    for i in range(len(target) - K + 1):
        index[target[i:i + K]].append(i)

    votes = collections.Counter()
    rc_target = revcomp(target)
    rc_index = collections.defaultdict(list)
    for i in range(len(rc_target) - K + 1):
        rc_index[rc_target[i:i + K]].append(i)

    for i in range(0, min(anchor, len(reference) - K)):
        mer = reference[i:i + K]
        for j in index.get(mer, ()):
            votes[(False, (j - i) % len(target))] += 1
        for j in rc_index.get(mer, ()):
            votes[(True, (j - i) % len(target))] += 1

    if not votes:
        return target, False, 0
    (flipped, offset), _ = votes.most_common(1)[0]
    seq = rc_target if flipped else target
    return seq[offset:] + seq[:offset], flipped, offset


def main():
    phage = {
        "S05_RHS": list(read_fasta(f"{GENOMES}/PGAP_S05/S05_APSE.fasta").values())[0],
        "S07_CdtB": list(read_fasta(f"{GENOMES}/PGAP_S07/S07_APSE.fasta").values())[0],
        "S12_CdtB": list(read_fasta(f"{GENOMES}/PGAP_S12/S12_APSE.fasta").values())[0],
    }
    toxin = read_fasta(TOXINS)
    toxin_of = {"S05_RHS": "S05_RHS", "S07_CdtB": "S07_cdtB", "S12_CdtB": "S12_CdtB"}

    # Put every genome on the same circular start coordinate and strand as the
    # reference, so that a straight diagonal means colinearity and a break in it
    # means a real difference.
    reference = "S07_CdtB"
    orientation = {reference: (False, 0)}
    for name in ("S05_RHS", "S12_CdtB"):
        phage[name], flipped, offset = orient_to(phage[reference], phage[name])
        orientation[name] = (flipped, offset)
        print(f"  {name}: rotated by {offset:,} bp"
              f"{', reverse-complemented' if flipped else ''}")

    # Alignment summaries, from: minimap2 -x asm20 --cs -c
    # (aligned fraction of the x-axis genome, number of blocks, identity range)
    panels = [
        ("S12_CdtB", "S07_CdtB", "A", "Same toxin family (CdtB vs CdtB)",
         "aligns end to end — 95 % of the genome in 2 blocks"),
        ("S05_RHS", "S07_CdtB", "B", "Different toxin families (RHS vs CdtB)",
         "75 % aligns in 4 blocks — ~10 kb unaligned at the toxin locus"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.8), facecolor=SURFACE)

    for ax, (xname, yname, letter, subtitle, note) in zip(axes, panels):
        a, b = phage[xname], phage[yname]
        (fx, fy), (rx, ry) = dotplot_points(a, b)

        ax.set_facecolor(SURFACE)
        ax.scatter(fx, fy, s=0.6, c=FWD, marker=".", linewidths=0, rasterized=True)
        ax.scatter(rx, ry, s=0.6, c=REV, marker=".", linewidths=0, rasterized=True)

        # Mark the toxin gene on each axis
        for name, axis in ((xname, "x"), (yname, "y")):
            hit = locate(phage[name], toxin[toxin_of[name]])
            if not hit:
                continue
            start, end, _ = hit
            if axis == "x":
                ax.axvspan(start, end, color=TOXIN_BAND, alpha=0.30, lw=0, zorder=0)
            else:
                ax.axhspan(start, end, color=TOXIN_BAND, alpha=0.30, lw=0, zorder=0)

        ax.set_xlim(0, len(a)); ax.set_ylim(0, len(b))
        rc_note = "  (rev. comp.)" if orientation.get(xname, (False, 0))[0] else ""
        ax.set_xlabel(f"{xname}   ({len(a):,} bp){rc_note}", fontsize=9, color=INK_SECONDARY)
        ax.set_ylabel(f"{yname}   ({len(b):,} bp)", fontsize=9, color=INK_SECONDARY)
        ax.set_title(f"{letter}.  {subtitle}", fontsize=10.5, color=INK,
                     loc="left", pad=9, fontweight="medium")
        ax.text(0.5, -0.155, note, transform=ax.transAxes, ha="center", va="top",
                fontsize=8.8, color=INK_SECONDARY)

        ax.tick_params(labelsize=8, colors=INK_MUTED, length=3)
        for s in ax.spines.values():
            s.set_color("#dcdbd5"); s.set_linewidth(0.8)
        ax.grid(True, color="#ecebe6", lw=0.6, zorder=-1)
        ax.set_axisbelow(True)
        ax.set_aspect("equal", adjustable="box")

    legend = [
        Line2D([], [], marker="o", ls="", ms=6, mfc=FWD, mec=FWD, label="forward match"),
        Line2D([], [], marker="o", ls="", ms=6, mfc=REV, mec=REV, label="reverse match"),
        Line2D([], [], marker="s", ls="", ms=7, mfc=TOXIN_BAND, mec="none",
               alpha=0.5, label="toxin coding sequence"),
    ]
    fig.legend(handles=legend, loc="lower center", ncol=3, frameon=False,
               fontsize=9, labelcolor=INK_SECONDARY, bbox_to_anchor=(0.5, -0.005))

    fig.suptitle("APSE phages share a colinear backbone and differ at the toxin locus",
                 fontsize=12.5, color=INK, x=0.5, y=0.985, fontweight="semibold")
    fig.text(0.5, 0.935,
             "Genomes rotated to a common start coordinate.\n"
             "The diagonal runs unbroken through the toxin band when both phages carry the same toxin (A);\n"
             "it breaks at the toxin locus when they differ (B).",
             ha="center", va="top", fontsize=9.2, color=INK_SECONDARY, linespacing=1.5)

    fig.tight_layout(rect=(0, 0.085, 1, 0.855))
    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/APSE_dotplot.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
