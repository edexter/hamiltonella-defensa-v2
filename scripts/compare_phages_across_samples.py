#!/usr/bin/env python3
"""
Qualitative cross-sample comparison of every APSE genome in the panel.

The question this answers: are the phages in different samples the same lineage,
or should each be treated as its own entity? It is a naming question, not a
phylogenetic one -- no tree is built, because these genomes recombine and a
single topology would describe a history that did not happen.

Two measures, deliberately not one:
  CONTAINMENT   fraction of the shorter genome's canonical 31-mers that are also
                in the longer one. Alignment-free. A shared backbone with
                different cassettes gives an intermediate value rather than the
                misleadingly high one a whole-genome identity would report.
  MODULE TYPE   which toxin and which integrase each genome carries. Modules are
                known to assort independently of one another here, so this is
                shown alongside rather than folded into the similarity score.

Clustering is average-linkage on 1 - containment, implemented directly to avoid
a scipy dependency.

Usage:  python3 scripts/compare_phages_across_samples.py <phages.fa> <mod.tsv>
Output: results/figures/cross_sample_phages.png
        results/cross_sample_containment.tsv
"""

import os
import sys
import collections

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

K = 31
COMP = str.maketrans("ACGT", "TGCA")
INK, INK2, MUTED, SURFACE = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"
TOXCOL = {"RHS": "#2a78d6", "CdtB": "#eb6834", "MAC": "#1baf7a",
          "Shiga": "#7a4fbd", "LeucineRich": "#d4a017", "-": "#d8d7d0"}
INTCOL = {"A": "#2a78d6", "B": "#eb6834", "-": "#d8d7d0"}


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def read_fa(path):
    d, name, order = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            d[name] = []
            order.append(name)
        else:
            d[name].append(line.strip())
    return order, {k: "".join(v).upper() for k, v in d.items()}


def family(t):
    t = t.upper()
    for f in ("RHS", "CDTB", "MAC", "SHIGA", "LEUCINERICH"):
        if f in t:
            return {"CDTB": "CdtB", "SHIGA": "Shiga",
                    "LEUCINERICH": "LeucineRich"}.get(f, f)
    return t


def average_linkage(names, dist):
    """returns leaf order from average-linkage clustering"""
    clusters = {i: [i] for i in range(len(names))}
    d = {(i, j): dist[i][j] for i in range(len(names)) for j in range(len(names)) if i < j}
    while len(clusters) > 1:
        (a, b) = min(d, key=d.get)
        merged = clusters[a] + clusters[b]
        del clusters[a], clusters[b]
        new = max(clusters) + 1 if clusters else a
        for k in list(d):
            if a in k or b in k:
                del d[k]
        for c, members in clusters.items():
            v = sum(dist[x][y] for x in merged for y in members) / (len(merged) * len(members))
            d[(min(c, new), max(c, new))] = v
        clusters[new] = merged
    return list(clusters.values())[0]


def main():
    fa = sys.argv[1] if len(sys.argv) > 1 else "phages.fa"
    modf = sys.argv[2] if len(sys.argv) > 2 else "mod.tsv"
    order, seqs = read_fa(fa)

    mods = collections.defaultdict(lambda: collections.defaultdict(set))
    for line in open(modf):
        g, k, v = line.split()
        mods[g][k].add(family(v) if k == "tox" else v)

    ks = {n: {canon(seqs[n][i:i + K]) for i in range(len(seqs[n]) - K + 1)}
          for n in order}
    N = len(order)
    cont = [[0.0] * N for _ in range(N)]
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            small, big = (ks[a], ks[b]) if len(ks[a]) <= len(ks[b]) else (ks[b], ks[a])
            cont[i][j] = len(small & big) / len(small)

    with open("results/cross_sample_containment.tsv", "w") as fh:
        fh.write("genome\t" + "\t".join(order) + "\n")
        for i, a in enumerate(order):
            fh.write(a + "\t" + "\t".join(f"{cont[i][j]:.4f}" for j in range(N)) + "\n")

    dist = [[1 - cont[i][j] for j in range(N)] for i in range(N)]
    leaf = average_linkage(order, dist)
    lab = [order[i] for i in leaf]

    fig = plt.figure(figsize=(13.2, 10.4), facecolor=SURFACE)
    fig.suptitle("Every APSE genome in the panel, compared with every other",
                 fontsize=15.5, color=INK, y=0.975, fontweight="semibold")
    fig.text(0.5, 0.938,
             "Shared fraction of canonical 31-mers, clustered. No tree is drawn: these phages recombine, so a single topology would describe a history that did not happen.\n"
             "Toxin and integrase are shown separately because they assort independently of each other and of overall similarity.",
             ha="center", va="top", fontsize=9.6, color=INK2, linespacing=1.55)

    ax = fig.add_axes([0.135, 0.175, 0.575, 0.690])
    M = [[cont[leaf[i]][leaf[j]] for j in range(N)] for i in range(N)]
    im = ax.imshow(M, cmap="magma_r", vmin=0.0, vmax=1.0, aspect="equal")
    ax.set_xticks(range(N)); ax.set_yticks(range(N))
    ax.set_xticklabels(lab, rotation=90, fontsize=8.2, color=INK)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
    ax.set_yticklabels(lab, fontsize=8.2, color=INK)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for i in range(N):
        for j in range(N):
            v = M[i][j]
            if i != j and v >= 0.55:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6.2,
                        color="white" if v > 0.75 else INK)

    cax = fig.add_axes([0.905, 0.175, 0.014, 0.690])
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.set_label("fraction of the shorter genome's 31-mers found in the other",
                 fontsize=8.8, color=INK2)
    cb.ax.tick_params(labelsize=8, colors=MUTED, length=3)
    cb.outline.set_visible(False)

    axm = fig.add_axes([0.722, 0.175, 0.075, 0.690])
    axm.set_xlim(0, 2); axm.set_ylim(N - 0.5, -0.5); axm.axis("off")
    axm.text(0.5, -1.1, "toxin", ha="center", fontsize=9, color=INK2, rotation=90,
             va="bottom")
    axm.text(1.5, -1.1, "int", ha="center", fontsize=9, color=INK2, rotation=90,
             va="bottom")
    for i, g in enumerate(lab):
        t = sorted(mods[g]["tox"]) or ["-"]
        n = sorted(mods[g]["int"]) or ["-"]
        axm.add_patch(Rectangle((0.08, i - 0.42), 0.84, 0.84,
                                facecolor=TOXCOL.get(t[0], "#d8d7d0"), edgecolor="white", lw=1))
        axm.add_patch(Rectangle((1.08, i - 0.42), 0.84, 0.84,
                                facecolor=INTCOL.get(n[0], "#d8d7d0"), edgecolor="white", lw=1))
        axm.text(2.15, i, f"{len(seqs[g]):,} bp", fontsize=7.6, color=MUTED, va="center")

    handles = [Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="none")
               for c in ["#2a78d6", "#eb6834", "#1baf7a", "#7a4fbd", "#d4a017", "#d8d7d0"]]
    fig.legend(handles, ["RHS", "CdtB", "MAC", "Shiga", "LeucineRich", "none detected"],
               loc="upper left", bbox_to_anchor=(0.135, 0.115), frameon=False,
               fontsize=9, labelcolor=INK2, ncol=6, title="toxin carried",
               title_fontsize=9)
    print("wrote results/cross_sample_containment.tsv")

    os.makedirs("results/figures", exist_ok=True)
    out = "results/figures/cross_sample_phages.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    print(f"wrote {out}")


main()
