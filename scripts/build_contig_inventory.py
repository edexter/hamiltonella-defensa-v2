#!/usr/bin/env python3
"""
Write a one-line description of every contig delivered, for every sample.

WHAT THIS IS FOR

A reader opening one of these FASTA files should be able to find out what each
contig is without running anything. The headers carry that information already,
but they are written for NCBI rather than for a person, and a sample with forty
contigs is tedious to read that way. This produces the same information as a
short list, one line per contig, grouped so that the interesting contigs are not
buried among the routine ones.

Contig numbers are kept exactly as delivered, including where a number is now
missing because a contaminant contig was removed. The gap is deliberate: these
numbers correspond to the version 1.0 delivery and to the public NCBI record, and
renumbering them to close a gap would break that correspondence for the sake of
tidiness.

Runs of unplaced contigs are collapsed into a single line once there are more
than a few, because listing forty of them separately tells the reader nothing
that one line does not.

USAGE

    python3 scripts/build_contig_inventory.py results/replicon_table.tsv \\
            assemblies_v2 results/contig_inventory.txt

DEPENDENCIES

    python3    standard library only
"""

import collections
import os
import re
import sys

# Above this many unplaced contigs in one sample, they are listed as a single
# collapsed line rather than one line each.
COLLAPSE_ABOVE = 5


def read_headers(path):
    """Return [(contig, header, length)] in the order the file delivers them."""
    out, name, header, n = [], None, None, 0
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                out.append((name, header, n))
            header = line.rstrip("\n")
            name = header[1:].split()[0]
            n = 0
        else:
            n += len(line.strip())
    if name is not None:
        out.append((name, header, n))
    return out


def contig_sort_key(name):
    """Sort contig_02 before contig_10 rather than the other way round."""
    m = re.search(r"(\d+)", name)
    return (int(m.group(1)) if m else 0, name)


def describe(label, header, sample, n_chromosome_contigs):
    """One phrase saying what a contig is."""
    if label == "chromosome":
        if n_chromosome_contigs == 1:
            text = "H. defensa chromosome"
        else:
            text = "Segment of H. defensa chromosome"
        for note in re.findall(r"\[note=([^\]]*prophage[^\]]*)\]", header):
            text += f" — {note}"
        return text
    if label == "plasmid":
        m = re.search(r"\[plasmid-name=([^\]]+)\]", header)
        return f"H. defensa plasmid {m.group(1)}" if m else "H. defensa plasmid"
    if label == "APSE phage":
        return "APSE"
    if label == "plasmid, not circularised":
        return "H. defensa plasmid (could not be circularized)"
    if label == "plasmid fragment":
        return "Fragment of an H. defensa plasmid that is complete in another sample"
    if label == "mis-joined plasmids":
        return "Two H. defensa plasmids joined end to end by the assembler in error"
    return "Unplaced H. defensa contig"


def collapse(names):
    """Render [contig_02, contig_03, contig_04, contig_07] as '02-04, 07'."""
    nums = sorted(int(re.search(r"(\d+)", x).group(1)) for x in names)
    runs, start, prev = [], nums[0], nums[0]
    for x in nums[1:]:
        if x == prev + 1:
            prev = x
            continue
        runs.append((start, prev))
        start = prev = x
    runs.append((start, prev))
    return ", ".join(f"{a:02d}" if a == b else
                     (f"{a:02d}, {b:02d}" if b == a + 1 else f"{a:02d}-{b:02d}")
                     for a, b in runs)


def main():
    table_path, asm_dir, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

    labels = collections.defaultdict(dict)
    with open(table_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            labels[f[col["sample"]]][f[col["contig"]]] = f[col["label"]]

    lines = []
    for sample in sorted(labels):
        path = f"{asm_dir}/{sample}_v2.fasta"
        if not os.path.exists(path):
            continue
        entries = read_headers(path)
        entries.sort(key=lambda e: contig_sort_key(e[0]))
        total_bp = sum(n for _, _, n in entries)
        n_chrom = sum(1 for c, _, _ in entries if labels[sample].get(c) == "chromosome")

        lines.append(f"File = {sample}_v2.fasta   "
                     f"({len(entries)} contigs, {total_bp:,} bp)")

        unplaced = [(c, n) for c, _, n in entries
                    if labels[sample].get(c, "unplaced") == "unplaced"]
        collapsing = len(unplaced) > COLLAPSE_ABOVE

        for contig, header, n in entries:
            label = labels[sample].get(contig, "unplaced")
            if label == "unplaced" and collapsing:
                continue
            topology = "circular" if "topology=circular" in header else "linear"
            lines.append(f"  {contig} ({n:,} bp {topology}): "
                         f"{describe(label, header, sample, n_chrom)}")

        if collapsing:
            sizes = sorted(n for _, n in unplaced)
            lines.append(f"  contig_{collapse([c for c, _ in unplaced])} "
                         f"({len(unplaced)} contigs, {sizes[0]:,}-{sizes[-1]:,} bp, "
                         f"all linear): Unplaced H. defensa contigs")

        # The phage reference files belonging to this sample.
        for d in ("results/phage_references", "results/phage_reconstructions"):
            if not os.path.isdir(d):
                continue
            for fn in sorted(os.listdir(d)):
                if not fn.startswith(f"{sample}_") or not fn.endswith(".fasta"):
                    continue
                for contig, header, n in read_headers(f"{d}/{fn}"):
                    topology = "circular" if "topology=circular" in header else "linear"
                    m = re.search(r"\[lineage=([^\]]+)\]", header)
                    complete = "complete" if "completeness=complete" in header else "partial"
                    detail = f"APSE bacteriophage ({complete})"
                    if m:
                        detail += f", lineage {m.group(1)}"
                    lines.append("")
                    lines.append(f"File = {fn}")
                    lines.append(f"  {contig} ({n:,} bp {topology}): {detail}")
        lines.append("")

    with open(out_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  inventory for {len(labels)} samples written to {out_path}")


main()
