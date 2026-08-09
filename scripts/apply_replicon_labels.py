#!/usr/bin/env python3
"""
Rewrite the delivered assemblies so that every contig carries the label the
evidence supports, and so that contigs which are not Hamiltonella at all are no
longer delivered as though they were.

WHAT THIS CHANGES

Three kinds of correction, all of them driven by results/replicon_table.tsv and
none of them by hand:

  removal      Contigs that align to the aphid genome over their whole length are
               dropped. They are aphid DNA that survived the 2024 contaminant
               filter. Three of them were being delivered with an explicit
               [location=chromosome] tag, which asserts they are Hamiltonella
               chromosome; the rest were merely unplaced. Either way they do not
               belong in a Hamiltonella genome and are removed.

  plasmids     A [plasmid-name=] tag asserts that a separate circular replicon
               exists. That assertion is kept only where the assembler closed the
               contig into a circle and it groups with a plasmid family. Contigs
               that are a fragment of a plasmid, or a mis-joined pair of two
               plasmids, lose the tag and are described as what they are. This
               matters because a fragment carrying its own plasmid name inflates
               the plasmid count and implies a molecule nobody has seen.

  consistency  The same element was labelled a plasmid in one sample and left
               unplaced in another. Labels now come from the family assignment,
               so an element is described the same way wherever it occurs.

WHAT THIS DOES NOT CHANGE

No sequence is altered. Contigs are removed or relabelled; the bases of every
contig that survives are byte-identical to the input. This is checked at the end
rather than assumed, and the check is reported.

USAGE

    python3 scripts/apply_replicon_labels.py results/replicon_table.tsv \\
            assemblies_v2 assemblies_v2

Reading and writing the same directory is intended: the previous version of each
file is kept alongside with a .pre_relabel suffix so the change can be inspected
or undone.

DEPENDENCIES

    python3    standard library only
"""

import collections
import os
import re
import shutil
import sys


def read_fasta_with_headers(path):
    """Return [(header_line, sequence)] preserving the order in the file."""
    entries, header, seq = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq)))
            header, seq = line.rstrip("\n"), []
        else:
            seq.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def main():
    table_path, in_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]

    # Where the 2024 curated assemblies live, so that plasmid names can be read
    # from the original delivery rather than from this script's own output.
    curated_dir = sys.argv[4] if len(sys.argv) > 4 else ""
    if not curated_dir and os.path.exists("data/config.sh"):
        for line in open("data/config.sh"):
            if line.strip().startswith("ASSEMBLY_DIR="):
                curated_dir = line.split("=", 1)[1].strip().strip('"')

    rows = collections.defaultdict(dict)
    with open(table_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[f[col["sample"]]][f[col["contig"]]] = f

    # The topology each phage was found to have when tested against reads. A
    # sample carrying one standalone phage contig takes its answer from here
    # rather than from the assembler.
    phage_topology = {}
    for d in ("results/phage_references", "results/phage_reconstructions"):
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".fasta"):
                continue
            sample = fn.split("_")[0]
            with open(f"{d}/{fn}") as fh:
                header = fh.readline()
            m = re.search(r"\[topology=([^\]]+)\]", header)
            if m and sample not in phage_topology:
                phage_topology[sample] = m.group(1)

    total_removed = total_relabelled = total_kept = 0
    removed_bp = 0

    for sample in sorted(rows):
        path = f"{in_dir}/{sample}_v2.fasta"
        if not os.path.exists(path):
            print(f"  {sample}: {path} not found, skipped")
            continue
        entries = read_fasta_with_headers(path)
        shutil.copyfile(path, f"{path}.pre_relabel")

        strain = ""
        m = re.search(r"\[strain=([^\]]+)\]", entries[0][0]) if entries else None
        if m:
            strain = m.group(1)

        # Plasmid names assigned in 2024 are kept exactly as they are wherever
        # the plasmid itself is confirmed. They are already published in the
        # NCBI submission, so renaming them -- even into a tidier order -- would
        # break the correspondence between this delivery and the public record
        # for no benefit. A plasmid that loses its name here loses it because
        # the replicon turned out not to exist, not because the numbering moved.
        #
        # The names are read from the 2024 curated assembly, never from this
        # script's own output. Reading them from the delivered file would make a
        # demotion permanent: once a contig has lost its name, a later run that
        # re-confirms the plasmid could not give the original name back, and
        # would invent a new one. That happened to S10 contig_10, which was
        # demoted on the assembler's word and then shown by reads to be circular
        # after all.
        existing_name = {}
        curated = f"{curated_dir}/{sample}/{sample}.curated.fasta" if curated_dir else ""
        source_headers = []
        if curated and os.path.exists(curated):
            for line in open(curated):
                if line.startswith(">"):
                    source_headers.append(line.rstrip("\n"))
        else:
            source_headers = [h for h, _ in entries]
        for header in source_headers:
            contig = header[1:].split()[0]
            m = re.search(r"\[plasmid-name=([^\]]+)\]", header)
            if m:
                existing_name[contig] = m.group(1)
                n = re.search(r"(\d+)", contig)
                if n:
                    existing_name[f"contig_{int(n.group(1)):02d}"] = m.group(1)

        # Only a plasmid confirmed here keeps its name; anything newly confirmed
        # that never had one is numbered after the highest existing number, so
        # no name is ever reused for a different molecule.
        used = set(existing_name.values())
        next_number = 1
        plasmid_name = {}
        for contig, f in sorted(rows[sample].items(),
                                key=lambda kv: -int(kv[1][col["length_bp"]])):
            if f[col["label"]] != "plasmid":
                continue
            if contig in existing_name:
                plasmid_name[contig] = existing_name[contig]
            else:
                while f"p{strain or sample}_{next_number}" in used:
                    next_number += 1
                name = f"p{strain or sample}_{next_number}"
                plasmid_name[contig] = name
                used.add(name)

        out, removed, relabelled = [], [], 0
        for header, seq in entries:
            contig = header[1:].split()[0]
            row = rows[sample].get(contig)
            if row is None:
                out.append((header, seq))
                continue
            label = row[col["label"]]

            if label == "aphid":
                removed.append((contig, len(seq)))
                continue

            bits = ["[organism=Hamiltonella defensa]"]
            if strain:
                bits.append(f"[strain={strain}]")

            if label == "APSE phage":
                bits.append("[note=APSE bacteriophage; supplied separately as a "
                            "phage reference]")
            elif label == "chromosome":
                bits.append("[location=chromosome]")
            elif label == "plasmid":
                bits.append(f"[plasmid-name={plasmid_name[contig]}]")
            elif label == "plasmid, not circularised":
                bits.append("[note=plasmid; the assembler did not close this "
                            "contig into a circle, so it is delivered as a "
                            "linear fragment of the plasmid]")
            elif label == "plasmid fragment":
                bits.append("[note=fragment of a plasmid that assembles completely "
                            "in another sample of this panel; not a separate "
                            "replicon]")
            else:
                bits.append("[note=Hamiltonella defensa contig, unplaced]")

            # Topology normally comes from the assembler's own call. A phage
            # contig is the exception: its circularity was tested directly, by
            # rotating the sequence and asking whether reads span the join, and
            # that test overrides the assembler in both directions. Without this
            # the same molecule ends up labelled linear in the genome file and
            # circular in the phage file, which is the annotation fault version
            # 2.0 exists to remove.
            topology = "circular" if row[col["circular"]] == "Y" else "linear"
            if label == "APSE phage" and sample in phage_topology:
                topology = phage_topology[sample]
            bits.append(f"[topology={topology}]")

            # Anything the original header said about an integrated prophage is
            # a result of this project and must survive relabelling.
            for note in re.findall(r"\[note=([^\]]*prophage[^\]]*)\]", header):
                bits.append(f"[note={note}]")

            new_header = f">{contig} " + " ".join(bits)
            if new_header != header:
                relabelled += 1
            out.append((new_header, seq))

        os.makedirs(out_dir, exist_ok=True)
        with open(f"{out_dir}/{sample}_v2.fasta", "w") as fh:
            for header, seq in out:
                fh.write(header + "\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i:i + 60] + "\n")

        # Confirm that nothing but the labels moved.
        before = {h[1:].split()[0]: s for h, s in entries}
        after = {h[1:].split()[0]: s for h, s in out}
        changed = [c for c in after if after[c] != before[c]]
        assert not changed, f"{sample}: sequence changed in {changed}"

        total_removed += len(removed)
        total_relabelled += relabelled
        total_kept += len(out)
        removed_bp += sum(n for _, n in removed)
        note = ""
        if removed:
            note = "  removed: " + ", ".join(f"{c} ({n:,} bp)" for c, n in removed)
        print(f"  {sample}: {len(out):>2} contigs kept, {relabelled:>2} relabelled,"
              f" {len(removed)} removed{note}")

    print(f"\n  {total_kept} contigs delivered, {total_relabelled} relabelled, "
          f"{total_removed} removed ({removed_bp:,} bp of aphid DNA)")
    print("  every surviving contig is byte-identical in sequence to its input")


main()
