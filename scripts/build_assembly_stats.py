#!/usr/bin/env python3
"""
Build the version 2.0 assembly statistics table, in the same shape as the
version 1.0 spreadsheet so that the two can be compared row by row.

The assembly columns are computed from the delivered FASTA files. The BUSCO
columns are read from BUSCO's own short_summary files rather than recomputed, so
that the numbers in the table are exactly the numbers BUSCO reported.

A NOTE ON COMPARABILITY

The version 1.0 figures were produced with BUSCO 5.1.2 against
enterobacterales_odb10, in genome mode with prodigal, searching 440 BUSCO groups.
A different lineage dataset would change the denominator and make the two
versions incomparable, so the script checks the dataset named in each summary
file and refuses to mix datasets silently.

Most samples should be unchanged from version 1.0. Removing an aphid contig
cannot remove a bacterial BUSCO, and moving a prophage from its own contig into
the chromosome does not alter gene content. The samples where a change is
expected are those that gained sequence -- S05 and S07 -- and a change anywhere
else is worth investigating rather than reporting.

USAGE

    python3 scripts/build_assembly_stats.py assemblies_v2 <busco_dir> \\
            docs/curated_assembly_stats_v2.xlsx

<busco_dir> holds one directory per sample containing BUSCO's output. If it does
not exist, or a sample has no summary, the assembly columns are still written and
the BUSCO columns are left empty.

DEPENDENCIES

    python3, openpyxl (only for writing the spreadsheet; a TSV is always written)
"""

import glob
import os
import re
import sys

COLUMNS = ["sampleID", "strainID", "bases", "contigs", "meanContigLength",
           "largestContig", "smallestContig", "N50", "L50",
           "Complete BUSCOs (C)", "Single-copy BUSCOs (S)",
           "Duplicated BUSCOs (D)", "Fragmented BUSCOs (F)",
           "Missing BUSCOs (M)", "Total BUSCO groups (n)", "BUSCO %"]


def read_lengths(path):
    """Contig lengths and the strain name from the FASTA headers."""
    lengths, strain, n = [], "", 0
    for line in open(path):
        if line.startswith(">"):
            if n:
                lengths.append(n)
            n = 0
            m = re.search(r"\[strain=([^\]]+)\]", line)
            if m and not strain:
                strain = m.group(1)
        else:
            n += len(line.strip())
    if n:
        lengths.append(n)
    return sorted(lengths, reverse=True), strain


def n50(lengths):
    """N50 and L50. N50 is the length of the contig at which half the assembly
    is contained in contigs of that length or longer; L50 is how many contigs
    that takes."""
    total = sum(lengths)
    running = 0
    for i, x in enumerate(lengths, 1):
        running += x
        if running >= total / 2:
            return x, i
    return 0, 0


def read_busco(busco_dir, sample):
    """Pull the six counts out of BUSCO's short_summary file."""
    pattern = os.path.join(busco_dir, f"{sample}*", "short_summary*.txt")
    hits = glob.glob(pattern) or glob.glob(
        os.path.join(busco_dir, f"short_summary*{sample}*.txt"))
    if not hits:
        return None, None
    text = open(hits[0]).read()
    dataset = ""
    m = re.search(r"The lineage dataset is:\s*(\S+)", text)
    if m:
        dataset = m.group(1)
    counts = {}
    for label, key in (("Complete BUSCOs", "C"),
                       ("Complete and single-copy BUSCOs", "S"),
                       ("Complete and duplicated BUSCOs", "D"),
                       ("Fragmented BUSCOs", "F"),
                       ("Missing BUSCOs", "M"),
                       ("Total BUSCO groups searched", "n")):
        m = re.search(rf"(\d+)\s+{re.escape(label)}", text)
        if m:
            counts[key] = int(m.group(1))
    return counts, dataset


def main():
    asm_dir, busco_dir, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    rows, datasets = [], set()

    for sample in [f"S{i:02d}" for i in range(1, 15)]:
        path = f"{asm_dir}/{sample}_v2.fasta"
        if not os.path.exists(path):
            print(f"  {sample}: {path} not found, skipped")
            continue
        lengths, strain = read_lengths(path)
        total = sum(lengths)
        n, l = n50(lengths)
        row = [sample, strain, total, len(lengths), total / len(lengths),
               lengths[0], lengths[-1], n, l]

        counts, dataset = read_busco(busco_dir, sample) if os.path.isdir(busco_dir) else (None, None)
        if counts and "n" in counts:
            datasets.add(dataset)
            row += [counts.get("C", ""), counts.get("S", ""), counts.get("D", ""),
                    counts.get("F", ""), counts.get("M", ""), counts.get("n", ""),
                    100.0 * counts["C"] / counts["n"] if counts.get("C") is not None else ""]
        else:
            row += [""] * 7
        rows.append(row)

    if len(datasets) > 1:
        sys.exit(f"  BUSCO summaries use more than one lineage dataset {sorted(datasets)};\n"
                 f"  the percentages would not be comparable. Re-run them against one dataset.")
    if datasets:
        print(f"  BUSCO lineage dataset: {datasets.pop()}")

    tsv = out_path.rsplit(".", 1)[0] + ".tsv"
    with open(tsv, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"  {len(rows)} samples written to {tsv}")

    try:
        import openpyxl
    except ImportError:
        print("  openpyxl not installed; spreadsheet not written (the TSV above has everything)")
        return
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "merged_asm_read_stats"
    ws.append(COLUMNS)
    for r in rows:
        ws.append(r)
    for i, name in enumerate(COLUMNS, 1):
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = max(12, len(name) + 2)
    wb.save(out_path)
    print(f"  spreadsheet written to {out_path}")


main()
