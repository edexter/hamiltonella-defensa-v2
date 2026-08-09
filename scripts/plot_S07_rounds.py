#!/usr/bin/env python3
"""
One page per round of the S07 reconstruction, showing the diagnostics as they
were applied at that round.

Each page has the same three panels, so the pages can be flipped through and
compared:

  A  COMPOSITION  — every 1 kb window of the assembly, coloured by which of the
     two finished phages it matches. A correct assembly is one colour plus grey
     (grey = the phages are identical there, so the window carries no
     information). Windows of the *other* colour are the other phage's sequence
     sitting inside this one. Toxin and integrase positions are marked beneath.
  B  HETEROZYGOSITY — heterozygous sites per 1 kb window. Each candidate assembly
     gets its OWN axes in its OWN coordinates, because at this stage they are two
     different sequences and a shared x axis would be meaningless. Beneath each
     line is the mean alternate-allele fraction of the sites in that window: the
     shaded band means two lineages in one pile, a line at the top means the
     assembly took a minority allele.
  C  VERDICT — what the diagnostics said and what was changed for the next round.

All panels are measured. Every round is tested against the SAME read set (a
1,614-read sample of the S07 phage reads, ~358x) so the rounds are comparable;
the per-round figures quoted elsewhere used each round's own pile and differ
slightly.

Usage:  python3 scripts/plot_S07_rounds.py
Output: results/figures/S07_round1.png ... S07_round5.png
"""

import os
import collections

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D

BLUE = "#2a78d6"
ORANGE = "#eb6834"
SHARED = "#d8d7d0"
NEITHER = "#ffffff"
GREEN = "#1baf7a"
RED = "#b0402a"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
SURFACE = "#fcfcfb"
COL = {"RHS": BLUE, "CdtB": ORANGE, "same": SHARED, "neither": NEITHER}

EXPECTED = 39300          # a complete APSE genome in this sample

ROUNDS = {
 "r1": dict(n="1", ref="Clean single-phage genomes borrowed from other samples\n(S01's RHS phage, S12's CdtB phage)",
            changed="First attempt. The tidiest available references were used.",
            verdict="FAILED — a lineage was missing from the reference",
            why="Both borrowed genomes carry a type-A integrase. The CdtB phages of S05 and S07 carry\n"
                "type B, so the integrase-B region had no matching reference at all and its reads were\n"
                "not merely mis-sorted but sorted into the wrong pile.\n\n"
                "The result is visible in panel A: the RHS assembly is 16 windows of CdtB sequence mixed\n"
                "into 19 windows of RHS. It is a fusion, not a genome — and at 56,981 bp it is 45 % longer\n"
                "than an APSE genome, carrying BOTH integrases (type B at 7.2 kb, type A at 49.5 kb).",
            action="Build the reference from S07's OWN contigs, so every lineage present is represented.",
            ok=False),
 "r2": dict(n="2", ref="A single template, with permissive (non-competitive)\nread recruitment",
            changed="Reference rebuilt from S07's contigs, but reads were recruited to one template at a time.",
            verdict="FAILED — sorting must stay competitive",
            why="With no competing template, the other phage's reads were recruited through the shared\n"
                "backbone. The assembler was handed two phages' worth of reads for one genome and joined\n"
                "them end to end: 78,792 bp, twice the expected length.\n\n"
                "This is the worst round by every diagnostic. 40 of its 78 windows match NEITHER finished\n"
                "phage cleanly, and heterozygosity is 8.76 per kb — the reads disagree with it everywhere.",
            action="Nothing from this round is carried forward. Sort every read against ALL references at\n"
                   "once, every round, without exception.",
            ok=False),
 "r2b": dict(n="3", ref="S07's OWN contigs, grouped by integrase and toxin,\nwith COMPETITIVE recruitment",
            changed="Both faults corrected at once. This is the round where the method starts working.",
            verdict="THE FIX — one lineage clean, the other nearly so",
            why="Two changes, together. The reference is now contig_3 + contig_4 for RHS and contig_5 +\n"
                "contig_1 for CdtB — S07's own sequence, verified identical to its raw contigs — so every\n"
                "lineage in the sample is represented. And every read is now mapped against BOTH groups at\n"
                "once, taking whichever it matches uniquely; reads matching both equally go to both piles.\n\n"
                "The effect is stark. Round 2 put 15,867 reads into ONE pile — essentially the whole phage\n"
                "read set. Competition splits them into 8,752 and 7,007, and the RHS assembly drops from\n"
                "78,792 bp to 38,524 with every heterozygous position gone.",
            action="Iterate: use these two assemblies as the reference, which are far better than the\n"
                   "fragments this round started from.",
            ok=None),
 "r3": dict(n="4", ref="The round-3 assemblies (38,524 bp and 39,308 bp),\nwith competitive recruitment",
            changed="The reference is now two near-complete genomes rather than contig fragments.",
            verdict="PARTLY CORRECT — RHS is done, CdtB is not",
            why="The RHS phage is finished: 39,306 bp, circular, one toxin, one integrase, and NOT ONE\n"
                "position where the reads disagree. It never changes again in any later round.\n\n"
                "The CdtB assembly is close but still carries 5 windows of RHS-derived sequence, and its\n"
                "heterozygosity of 3.91 per kb sits at intermediate frequencies — reads split between two\n"
                "versions of the same region.",
            action="Iterate again — the CdtB assembly is still carrying the other phage in places.",
            ok=None),
 "r4": dict(n="5", ref="The round-4 assemblies, with stricter read assignment",
            changed="Same reference logic, tighter thresholds on which reads count as uniquely assigned.",
            verdict="ALMOST — one bad patch left in CdtB",
            why="CdtB improved from 5 foreign windows to 2, and from 3.91 to 0.94 per kb. But the\n"
                "residual has changed character: the remaining sites sit near the TOP of panel B, at 95 %\n"
                "and above, packed into one short block.\n\n"
                "That is a different fault from the previous round. It means the assembly took the minority allele\n"
                "across a 2.3 kb stretch — sequence inherited from the chimeric contig in the reference.\n"
                "The reads themselves say, unambiguously, what the correct bases are.",
            action="Patch those bases in the REFERENCE (never in the output) and switch to sorting on\n"
                   "parent-diagnostic 31-mers, which alignment score had been too blunt to resolve.",
            ok=None),
 "r5": dict(n="6", ref="The round-5 assemblies, patched, with reads sorted by\nparent-diagnostic 31-mers",
            changed="Sorting is now alignment-free and the reference no longer carries the bad patch.",
            verdict="CONVERGED — both phages meet the benchmark",
            why="Both assemblies are circular, one APSE genome long, carry exactly one toxin and one\n"
                "integrase of the expected type, and contain no window of the other phage's sequence.\n\n"
                "Heterozygosity is 0.00 and 0.03 per kb. The single-phage control sample S12 gives 0.00 at\n"
                "1,382× coverage, so this is the floor of the assay rather than a threshold we chose.\n\n"
                "The RHS assembly is byte-identical to the previous round's. Re-sorting every read against these two\n"
                "genomes moved 0.2 % of them — the loop has nothing left to do.",
            action="Ship. Both genomes go forward to att-site mapping and the v2 assembly.",
            ok=True),
}
ORDER = {"r1": ["r1_RHS", "r1_CdtB"], "r2": ["r2_RHS"], "r2b": ["r2b_RHS", "r2b_CdtB"],
         "r3": ["r3_RHS", "r3_CdtB"],
         "r4": ["r4_RHS", "r4_CdtB"], "r5": ["r5_RHS", "r5_CdtB"]}
NICE = {"RHS": "RHS pile", "CdtB": "CdtB pile"}


def load():
    summ = {}
    for line in open("results/S07_round_summary.tsv"):
        f = line.split()
        summ[f[0]] = dict(length=int(f[1]), circ=f[2], frags=int(f[3]))
    anc = collections.defaultdict(list)
    for i, line in enumerate(open("results/S07_round_ancestry.tsv")):
        if i == 0:
            continue
        a, w, c = line.split()
        anc[a].append((int(w), c))
    het = collections.defaultdict(list)
    for line in open("results/S07_round_het_sites.tsv"):
        a, p, f = line.split()
        het[a].append((int(p), float(f)))
    cont = collections.defaultdict(list)
    for line in open("results/S07_round_content.tsv"):
        f = line.rstrip("\n").split("\t")
        cont[f[0]].append((f[1], f[2], int(f[3]), int(f[4])))
    return summ, anc, het, cont


def page(rk, summ, anc, het, cont):
    R = ROUNDS[rk]
    keys = ORDER[rk]
    fig = plt.figure(figsize=(11.5, 9.6), facecolor=SURFACE)

    fig.text(0.045, 0.972, f"ROUND {R['n']}", fontsize=19, color=INK, fontweight="bold",
             va="top")
    fig.text(0.045, 0.928, R["ref"], fontsize=11, color=INK, va="top", linespacing=1.5,
             fontweight="medium")
    fig.text(0.045, 0.862, R["changed"], fontsize=9.3, color=MUTED, va="top",
             style="italic")

    vc = GREEN if R["ok"] else (RED if R["ok"] is False else "#b8860b")
    fig.text(0.955, 0.972, R["verdict"], fontsize=12.5, color=vc, ha="right", va="top",
             fontweight="semibold")

    # ---------- Panel A : composition ----------
    axA = fig.add_axes([0.075, 0.650, 0.860, 0.155])
    axA.set_title("A.   Composition — each 1 kb window coloured by which finished phage it matches",
                  fontsize=11, color=INK, loc="left", pad=10, fontweight="medium")
    y = 0.55
    for k in keys:
        lin = k.split("_")[1]
        L = summ[k]["length"]
        for w, c in anc[k]:
            axA.add_patch(Rectangle((w / 1000.0, y), 1.0, 0.30, facecolor=COL[c],
                                    edgecolor="#f4f2ec", lw=0.5))
        axA.text(-1.2, y + 0.15, NICE[lin], fontsize=10, color=INK2, ha="right",
                 va="center", fontweight="medium")
        n_other = sum(1 for _, c in anc[k] if c != lin and c in ("RHS", "CdtB"))
        n_nei = sum(1 for _, c in anc[k] if c == "neither")
        bits = [f"{L:,} bp", "circular" if summ[k]["circ"] == "Y" else "NOT circular"]
        if summ[k]["frags"] > 1:
            bits.append(f"+{summ[k]['frags']-1} fragment" + ("s" if summ[k]["frags"] > 2 else ""))
        axA.text(L / 1000.0 + 1.5, y + 0.235, "  ·  ".join(bits), fontsize=9,
                 color=INK2, va="center")
        flag = []
        if n_other:
            flag.append(f"{n_other} foreign window" + ("s" if n_other != 1 else ""))
        if n_nei:
            flag.append(f"{n_nei} unmatched")
        if flag:
            axA.text(L / 1000.0 + 1.5, y + 0.065, ",  ".join(flag), fontsize=9,
                     color=RED, va="center", fontweight="medium")
        seen_feat = []
        for kind, typ, a, b in sorted(cont.get(k, []), key=lambda x: x[2]):
            if any(kind == sk and typ == st and a < sb and sa < b
                   for sk, st, sa, sb in seen_feat):
                continue
            seen_feat.append((kind, typ, a, b))
            axA.add_patch(Rectangle((a / 1000.0, y), max((b - a) / 1000.0, 0.35), 0.30,
                                    facecolor="none", edgecolor=INK, lw=1.5,
                                    ls="-" if kind == "toxin" else (0, (2, 1.5))))
            axA.text((a + b) / 2000.0, y - 0.055,
                     f"{typ} toxin" if kind == "toxin" else f"int {typ}",
                     fontsize=8, color=INK, ha="center", va="top")
        y -= 0.52
    axA.axvline(EXPECTED / 1000.0, color=MUTED, lw=1.2, ls=(0, (3, 2)))
    axA.text(EXPECTED / 1000.0, 1.10, "one APSE genome (39.3 kb)", fontsize=8.5,
             color=MUTED, ha="center", va="bottom")
    axA.set_xlim(-9, 90)
    axA.set_ylim(-0.20 if len(keys) > 1 else 0.32, 1.14)
    axA.tick_params(labelsize=8, colors=MUTED, length=3)
    axA.set_yticks([])
    for s in ("top", "right", "left"):
        axA.spines[s].set_visible(False)
    axA.spines["bottom"].set_color("#dcdbd5")
    axA.set_facecolor(SURFACE)

    # ---------- Panel B : heterozygosity, one plot per candidate assembly ------
    def windows(k):
        """sites per 1 kb window, and the mean alternate-allele fraction in it"""
        L = summ[k]["length"]
        nwin = (L + 999) // 1000
        cnt = [0] * nwin
        tot = [0.0] * nwin
        for pos, frac in het.get(k, []):
            w = min((pos - 1) // 1000, nwin - 1)
            cnt[w] += 1
            tot[w] += frac
        mean = [(tot[i] / cnt[i]) if cnt[i] else float("nan") for i in range(nwin)]
        return list(range(nwin)), cnt, mean

    fig.text(0.075, 0.600,
             "B.   Heterozygosity in 1 kb windows — each candidate assembly on its own axes, in its own coordinates",
             fontsize=11, color=INK, va="top", fontweight="medium")

    ncol = len(keys)
    span = 0.860 if ncol == 1 else 0.405
    gap = 0.050
    for i, k in enumerate(keys):
        x0 = 0.075 + i * (span + gap)
        lin = k.split("_")[1]
        L = summ[k]["length"]
        xs, cnt, mean = windows(k)
        n = len(het.get(k, []))
        rate = 1000.0 * n / L
        xmax = L / 1000.0 * 1.02

        ax1 = fig.add_axes([x0, 0.455, span, 0.105])
        ax1.fill_between([x + 0.5 for x in xs], 0, cnt, color=COL[lin], alpha=0.28,
                         step="mid", zorder=2)
        ax1.plot([x + 0.5 for x in xs], cnt, color=COL[lin], lw=1.8,
                 drawstyle="steps-mid", zorder=3)
        tag = f"{NICE[lin]}   ·   {L:,} bp   ·   {n} sites   ·   {rate:.2f} /kb"
        ax1.set_title(tag + ("   ✓" if n == 0 else ""), fontsize=9.8,
                      color=COL[lin] if n == 0 else INK, loc="left", pad=6,
                      fontweight="semibold" if n == 0 else "medium")
        ax1.set_xlim(0, xmax); ax1.set_ylim(0, max(max(cnt) if cnt else 1, 1) * 1.35)
        ax1.set_ylabel("sites per kb", fontsize=8.6, color=INK2)
        ax1.tick_params(labelsize=7.8, colors=MUTED, length=3, labelbottom=False)
        ax1.grid(True, color="#ecebe6", lw=0.6); ax1.set_axisbelow(True)
        if n == 0:
            ax1.text(xmax / 2, 0.5, "no position where the reads disagree",
                     ha="center", va="center", fontsize=9.5, color=COL[lin],
                     transform=ax1.get_xaxis_transform(), fontweight="medium")

        ax2 = fig.add_axes([x0, 0.362, span, 0.080])
        ax2.axhspan(0.30, 0.90, color="#f2efe8", zorder=0)
        ax2.plot([x + 0.5 for x in xs], mean, color=COL[lin], lw=1.6, marker="o",
                 ms=3.0, mec="none", zorder=3)
        ax2.set_xlim(0, xmax); ax2.set_ylim(0, 1.05)
        ax2.set_xlabel(f"position in this assembly (kb)", fontsize=8.8, color=INK2)
        ax2.set_ylabel("alt-allele\nfraction", fontsize=8.4, color=INK2)
        ax2.set_yticks([0, 0.5, 1.0])
        ax2.tick_params(labelsize=7.8, colors=MUTED, length=3)
        ax2.grid(True, axis="x", color="#ecebe6", lw=0.6); ax2.set_axisbelow(True)
        if i == ncol - 1:
            ax2.text(xmax * 0.99, 0.60, "two lineages\nin one pile", fontsize=7.6,
                     color=MUTED, ha="right", va="center", style="italic",
                     linespacing=1.35)
            ax2.text(xmax * 0.99, 0.985, "minority allele taken", fontsize=7.6,
                     color=MUTED, ha="right", va="top", style="italic")
        for ax in (ax1, ax2):
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
            for sp in ("bottom", "left"):
                ax.spines[sp].set_color("#dcdbd5")
            ax.set_facecolor(SURFACE)

    # ---------- Panel C : verdict ----------
    axC = fig.add_axes([0.045, 0.028, 0.915, 0.242])
    axC.set_xlim(0, 1); axC.set_ylim(0, 1); axC.axis("off")
    axC.add_patch(FancyBboxPatch((0.0, 0.0), 1.0, 1.0,
                                 boxstyle="round,pad=0.008,rounding_size=0.02",
                                 facecolor="#f4f2ec", edgecolor="none"))
    axC.text(0.018, 0.93, "C.   What the diagnostics said", fontsize=11, color=INK,
             va="top", fontweight="medium")
    axC.text(0.018, 0.775, R["why"], fontsize=9.1, color=INK2, va="top", linespacing=1.55)
    nlines = R["why"].count("\n") + 1
    ay = max(0.775 - nlines * 0.083, 0.08)
    axC.text(0.018, ay, "→  " + R["action"], fontsize=9.4,
             color=GREEN if R["ok"] else RED, va="top", linespacing=1.6,
             fontweight="medium")

    legend = [Line2D([], [], marker="s", ls="", ms=10, mfc=BLUE, mec="none",
                     label="matches the RHS phage"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=ORANGE, mec="none",
                     label="matches the CdtB phage"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=SHARED, mec="none",
                     label="phages identical here — uninformative"),
              Line2D([], [], marker="s", ls="", ms=10, mfc=NEITHER, mec="#c9c8c1",
                     label="matches neither")]
    fig.legend(handles=legend, loc="lower center", ncol=4, frameon=False, fontsize=9,
               labelcolor=INK2, bbox_to_anchor=(0.5, 0.292))

    out = f"results/figures/S07_round{R['n']}.png"
    fig.savefig(out, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    os.makedirs("results/figures", exist_ok=True)
    summ, anc, het, cont = load()
    for rk in ["r1", "r2", "r2b", "r3", "r4", "r5"]:
        page(rk, summ, anc, het, cont)


if __name__ == "__main__":
    main()
