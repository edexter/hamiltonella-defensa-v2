#!/usr/bin/env python3
"""
Test whether a standalone APSE contig really exists as a circle in its sample.

An assembler's "circular" flag is a statement about the assembly graph, not about
the molecule. It can be wrong in both directions: a contig flagged circular may
have been closed through a repeat, and one flagged linear may be a genuine circle
the assembler could not close. Sample S06 showed the same distinction from the
other side -- an integrated prophage whose circular form does not exist at all.

The test here is direct. Rotate the contig by half its length, so the point where
its two ends meet now sits in the middle, then map the sample's reads. A read can
only cross that point if the two ends are adjacent in a real molecule.

Also reports whether a known att core is present, which is corroborating: a free
circular phage carries attP, and an att core is what lets it integrate at all.

Usage: python3 scripts/test_phage_circularity.py <name> <contig.fa> <reads.fq> <outdir>
"""

import os
import re
import subprocess
import sys

SPAN = 300          # a read must reach this far past the join on both sides
MIN_READS = 5       # below this we do not claim the circle exists

ATT = {
    "A": "TAACCCCTTGATTTTATTGGTACGCCCTACTGGATTCGAACCAGTGACCTACGGCTTAGAAG",
    "B": "GGGACTTAAAATCCCTTGGCTTATGGCTGTGCGGGTTCAAGTCCCGCCTCGGGTACCA",
}
COMP = str.maketrans("ACGT", "TGCA")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


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


def main():
    name, fa, reads, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    seq = list(read_fa(fa).values())[0]

    half = len(seq) // 2
    rot = seq[half:] + seq[:half]
    rotfa = os.path.join(outdir, f"{name}_rot.fa")
    with open(rotfa, "w") as fh:
        fh.write(f">{name}_rot\n")
        for i in range(0, len(rot), 70):
            fh.write(rot[i:i + 70] + "\n")

    bam = os.path.join(outdir, f"{name}_rot.bam")
    subprocess.run(f"minimap2 -ax map-hifi -t 6 --secondary=no {rotfa} {reads} 2>/dev/null "
                   f"| samtools sort -o {bam} - && samtools index {bam}",
                   shell=True, capture_output=True)

    join = len(seq) - half
    out = subprocess.run(["samtools", "view", bam,
                          f"{name}_rot:{max(1, join-SPAN)}-{join+SPAN}"],
                         capture_output=True, text=True).stdout
    spanning = 0
    for line in out.split("\n"):
        f = line.split("\t")
        if len(f) < 6 or int(f[1]) & 0x904:
            continue
        ref = int(f[3])
        for c, op in CIG.findall(f[5]):
            c = int(c)
            if op in "M=X":
                if ref <= join - SPAN and ref + c - 1 >= join + SPAN:
                    spanning += 1
                ref += c
            elif op in "DN":
                ref += c

    att = [t for t, core in ATT.items() if core in seq or rc(core) in seq]
    dep = subprocess.run(f"samtools depth -a -Q 0 {bam} | awk '{{s+=$3;n++}} END{{if(n)printf \"%.0f\", s/n}}'",
                         shell=True, capture_output=True, text=True).stdout or "0"
    verdict = "CIRCULAR" if spanning >= MIN_READS else "LINEAR_ONLY"
    print(f"{name}\t{len(seq)}\t{verdict}\t{spanning}\t{dep}\t"
          f"{'att-' + ','.join(att) if att else 'no known att core'}")


main()
