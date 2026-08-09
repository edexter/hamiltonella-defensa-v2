#!/usr/bin/env python3
"""
Count reads that align partly to one contig and partly to another.

Such a read is direct evidence that the two contigs are adjacent in the real
molecule. This is how a prophage sitting between two chromosome contigs shows
itself, and counting it correctly is what an earlier pass in this project failed
to do.

Alignments starting at position 1 of a contig are discarded: they accumulate at
contig starts as a mapping artefact and are not evidence of adjacency.

Reads SAM on stdin.  Usage: count_contig_links.py <sample>
"""

import collections
import statistics
import sys

MIN_LINKS = 3          # below this a link is not worth reporting
IGNORE_POS = 2         # alignments at position 1 are artefacts


def main():
    sample = sys.argv[1]
    pairs = collections.defaultdict(list)
    for line in sys.stdin:
        if "SA:Z:" not in line:
            continue
        f = line.split("\t")
        here, pos = f[2], int(f[3])
        if pos < IGNORE_POS:
            continue
        sa = line.split("SA:Z:")[1].split(";")[0].split(",")
        there, there_pos = sa[0], int(sa[1])
        if there == here or there_pos < IGNORE_POS:
            continue
        key = tuple(sorted([here, there]))
        pairs[key].append((pos, there_pos) if here == key[0] else (there_pos, pos))

    for (a, b), obs in sorted(pairs.items(), key=lambda kv: -len(kv[1])):
        if len(obs) < MIN_LINKS:
            continue
        print(f"{sample}\t{a}\t{b}\t{len(obs)}\t"
              f"{int(statistics.median(x for x, _ in obs))}\t"
              f"{int(statistics.median(y for _, y in obs))}")


main()
