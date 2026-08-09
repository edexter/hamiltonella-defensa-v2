#!/usr/bin/env python3
"""
Splice reconstructed prophages into a bacterial chromosome at att sites that
were determined from split reads.

The att core is a direct repeat. The uninfected chromosome carries one copy
(attB) and the free circular phage carries one copy (attP). Integration
recombines the two, so the prophage ends up flanked by two identical copies
(attL and attR) and the chromosome grows by exactly one phage-genome length.
This script reproduces that geometry:

    chromosome[1 .. attB_end]  +  phage circle opened at attP, core removed
                               +  chromosome[attB_start .. end]

Orientation is read out of the sequence rather than assumed: the chromosomal
core is compared with the phage core in both orientations and the script stops
if neither matches exactly. An exact match is itself a check that the att
coordinates called from the reads are right.

Usage:
    python integrate_prophage.py <chromosome.fa> <phages.fa> <sites.tsv> <out.fa>

The sites file is tab separated with a header line and these columns:
    chr_attB_start  chr_attB_end  phage  phage_attP_start  phage_attP_end
All coordinates 1-based inclusive. Lines beginning with '#' are ignored.
"""
import sys

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def read_fa(path):
    d, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            d[name] = []
        else:
            d[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def read_sites(path):
    sites, cols = [], None
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if cols is None:
            cols = f
            continue
        r = dict(zip(cols, f))
        sites.append((int(r["chr_attB_start"]), int(r["chr_attB_end"]),
                      r["phage"], int(r["phage_attP_start"]),
                      int(r["phage_attP_end"])))
    return sorted(sites)


def main():
    chrfa, phfa, sitesfile, out = sys.argv[1:5]
    chrom = read_fa(chrfa)
    name = list(chrom)[0]
    seq = chrom[name]
    ph = read_fa(phfa)

    pieces, log = [], []
    prev = 0
    for cs, ce, pname, ps, pe in read_sites(sitesfile):
        if pname not in ph:
            sys.exit(f"phage {pname} not found in {phfa}")
        p = ph[pname]
        core_chr = seq[cs - 1:ce]
        core_ph = p[ps - 1:pe]
        if core_chr == core_ph:
            orient = "forward"
        elif core_chr == rc(core_ph):
            orient = "reverse"
        else:
            sys.exit(f"att core mismatch at chromosome {cs}: the core in "
                     f"{pname} does not match in either orientation")
        body = p[pe:] + p[:ps - 1]           # circle opened at the core
        if orient == "reverse":
            body = rc(body)
        pieces.append(seq[prev:ce])          # chromosome up to and including attL
        pieces.append(body)
        prev = cs - 1                        # attR is the same chromosomal core
        log.append((pname, cs, ce, orient, len(body) + (ce - cs + 1)))
    pieces.append(seq[prev:])
    new = "".join(pieces)

    with open(out, "w") as fh:
        fh.write(f">{name}_v2\n")
        for i in range(0, len(new), 60):
            fh.write(new[i:i + 60] + "\n")

    off = 0
    for pname, cs, ce, orient, ins in log:
        print(f"  {pname:<20} attB {cs:>9,}-{ce:<9,} {orient:<7} "
              f"insert {ins:,} bp -> attL {cs+off:,} attR {cs+off+ins:,}")
        off += ins
    print(f"  chromosome {len(seq):,} bp -> {len(new):,} bp")


main()
