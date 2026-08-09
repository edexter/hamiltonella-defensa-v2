#!/usr/bin/env python3
"""
Build the version 2.0 assembly for a sample: correct the headers, and where a
prophage has been spliced in, substitute the rebuilt contig and drop the
now-redundant standalone phage contig.

The header corrections are the Tier A work. Three faults are fixed:

  topology     Every sample delivered in 2024 claims [topology=circular] on its
               chromosome while FLYE called the underlying contig circ=N. Some
               samples label three or four separate contigs as circular
               chromosomes, which no genome can have. Topology is taken from
               FLYE's own call.
  APSE         Seven samples label the phage [plasmid-name=...]. An APSE
               prophage is not a plasmid; it is relabelled.
  descriptors  Contigs delivered with no descriptor at all get one.

Provenance stays out of the headers and in results/v2_provenance.tsv, so the
FASTA remains clean enough to submit.

Usage: python3 scripts/build_v2_assemblies.py <sample> <curated.fa> <assembly_info.txt>
                                              <strain> <out.fa> [replace:CONTIG=file.fa] [drop:CONTIG]
"""

import re
import sys


def read_fa(path):
    d, order, name = {}, [], None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            d[name] = []
            order.append(name)
        else:
            d[name].append(line.strip())
    return order, {k: "".join(v).upper() for k, v in d.items()}


def norm(c):
    m = re.sub(r"[^0-9]", "", c)
    return f"contig_{int(m):02d}" if m else c


def main():
    samp, fa, info, strain, out = sys.argv[1:6]
    replace, drop = {}, set()
    for arg in sys.argv[6:]:
        if arg.startswith("replace:"):
            c, f = arg[8:].split("=")
            replace[c] = list(read_fa(f)[1].values())[0]
        elif arg.startswith("drop:"):
            drop.add(arg[5:])

    circ = {}
    for line in open(info).read().splitlines()[1:]:
        f = line.split("\t")
        circ[norm(f[0])] = f[3]

    order, seqs = read_fa(fa)
    # original descriptors, so plasmid names assigned during curation survive
    desc = {}
    for line in open(fa):
        if line.startswith(">"):
            c = line[1:].split()[0]
            desc[c] = line.rstrip("\n")

    lengths = sorted((len(seqs[c]) for c in order), reverse=True)
    chrom_cut = max(100000, lengths[0] * 0.02)

    with open(out, "w") as fh:
        n_repl = n_drop = 0
        for c in order:
            key = norm(c)
            if key in drop or c in drop:
                n_drop += 1
                continue
            s = replace.get(key) or replace.get(c) or seqs[c]
            if replace.get(key) or replace.get(c):
                n_repl += 1
            d = desc.get(c, "")
            plasmid = re.search(r"\[plasmid-name=([^\]]+)\]", d)
            is_apse = "APSE" in d.upper()
            topo = "circular" if circ.get(key) == "Y" else "linear"
            bits = [f"[organism=Hamiltonella defensa]", f"[strain={strain}]"]
            if is_apse:
                bits.append("[note=APSE bacteriophage; supplied separately as "
                            f"{samp}_APSE_v2.fasta]")
            elif plasmid:
                bits.append(f"[plasmid-name={plasmid.group(1)}]")
            elif len(s) >= chrom_cut:
                bits.append("[location=chromosome]")
            else:
                bits.append("[note=Hamiltonella defensa contig, unplaced]")
            bits.append(f"[topology={topo}]")
            fh.write(f">{key} " + " ".join(bits) + "\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    print(f"  {samp}: {len(order)-n_drop} contigs written "
          f"({n_repl} replaced, {n_drop} dropped)")


main()
