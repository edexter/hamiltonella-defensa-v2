#!/usr/bin/env python3
"""
Find chromosome<->phage junctions from split alignments.

Reads SAM on stdin (primary alignments with SA tags). For every read it
reconstructs the full set of alignment blocks (primary + supplementary),
places each on the read in original-read coordinates, and reports every
adjacent pair of blocks that lands on two different reference sequences.

A junction is reported wherever two adjacent blocks of one read sit on two
different reference sequences. When the reference given to the aligner holds both
a chromosome and a phage, those junctions are the points where the prophage joins
the chromosome, and their positions cluster to single bases.

USAGE

    samtools view -F 0x900 <alignments.bam> \
        | python3 scripts/lib/split_read_junctions.py > junctions.tsv

The BAM must have been made with supplementary alignments retained, which is the
default, since those are what carry the second half of a split read. The filter
-F 0x900 drops secondary and supplementary lines themselves while keeping the SA
tags on the primary line, which is where this script reads them from.

OUTPUT

One line per junction, with these columns:

    read   leftRef  leftEnd  leftStrand   rightRef  rightStart  rightStrand   gapOnRead

Positions are given in the coordinates of whichever reference that side of the
junction landed on. Note that they are reference positions, not positions within
the read: a junction seen at different places along different reads is still the
same junction if it maps to the same reference coordinate.
"""
import sys, re, collections

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def blocks(cigar):
    """return (read_start, read_end, ref_span) in the orientation of the alignment"""
    rs = 0
    started = False
    read_len_consumed = 0
    ref = 0
    lead = 0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "SH":
            if not started:
                lead += n
            continue
        started = True
        if op in "M=X":
            read_len_consumed += n
            ref += n
        elif op == "I":
            read_len_consumed += n
        elif op in "DN":
            ref += n
    return lead, lead + read_len_consumed, ref


def total_read_len(cigar):
    n_tot = 0
    for n, op in CIG.findall(cigar):
        if op in "MI=XS":
            n_tot += int(n)
    return n_tot


def main():
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 0x900:            # skip secondary/supplementary lines themselves
            continue
        name, rname, pos, cigar = f[0], f[2], int(f[3]), f[5]
        if rname == "*":
            continue
        sa = None
        for tag in f[11:]:
            if tag.startswith("SA:Z:"):
                sa = tag[5:]
        segs = []
        strand = "-" if flag & 0x10 else "+"
        L = total_read_len(cigar)
        rs, re_, refspan = blocks(cigar)
        if strand == "-":           # convert to original-read coordinates
            rs, re_ = L - re_, L - rs
        segs.append((rs, re_, rname, pos, pos + refspan - 1, strand))
        if sa:
            for rec in sa.rstrip(";").split(";"):
                p = rec.split(",")
                if len(p) < 6:
                    continue
                srn, spos, sstr, scig = p[0], int(p[1]), p[2], p[3]
                sL = total_read_len(scig)
                a, b, sspan = blocks(scig)
                if sstr == "-":
                    a, b = sL - b, sL - a
                segs.append((a, b, srn, spos, spos + sspan - 1, sstr))
        if len(segs) < 2:
            continue
        segs.sort()
        for i in range(len(segs) - 1):
            l, r = segs[i], segs[i + 1]
            if l[2] == r[2]:
                continue
            # breakpoint on each side, in reference coordinates, at the
            # end of the block that faces the junction
            lbp = l[4] if l[5] == "+" else l[3]
            rbp = r[3] if r[5] == "+" else r[4]
            print("\t".join(map(str, [name, l[2], lbp, l[5], r[2], rbp, r[5],
                                      r[0] - l[1]])))


main()
