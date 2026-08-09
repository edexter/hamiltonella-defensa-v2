#!/usr/bin/env python3
"""
Assign reads among N reference lineages by 31-mer containment.

This is the generalisation of kmer_ancestry.py needed when one of the lineages
is itself a recombinant of the others. Diagnostic k-mers -- those unique to one
reference -- barely exist for a mosaic lineage, because almost every one of its
k-mers is shared with one parent or the other. What distinguishes a mosaic is
that it matches ALONG ITS WHOLE LENGTH while each parent matches only part of
it. So instead of asking which reference owns a k-mer, this asks how many of a
read's k-mers each reference contains, and gives the read to whichever contains
most.

A read is placed in every reference within MARGIN of the best score, so reads
from regions where two lineages are identical land in both piles -- correct,
because the sequence there is the same and withholding them fragments both
assemblies.

Usage: python kmer_assign.py <refs.fa> <reads.fq> <outprefix>
Writes <outprefix>.calls.tsv  read, length, per-reference k-mer counts, call
       <outprefix>.<name>.ids one read-name list per reference
"""
import sys, collections

K = 31
MARGIN = 0.002          # within 0.2 % of the best score counts as a tie
MIN_HITS = 50           # a read must match some reference this well to be used

COMP = str.maketrans("ACGT", "TGCA")


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def kmers(seq):
    for i in range(len(seq) - K + 1):
        km = seq[i:i + K]
        if "N" not in km:
            yield canon(km)


def main():
    fa, fq, out = sys.argv[1:4]
    seqs, name = {}, None
    for line in open(fa):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    seqs = {k: "".join(v).upper() for k, v in seqs.items()}
    names = list(seqs)
    sets = {n: set(kmers(seqs[n] + seqs[n][:K - 1])) for n in names}   # circular
    for n in names:
        sys.stderr.write(f"  {n}: {len(seqs[n]):,} bp, {len(sets[n]):,} 31-mers\n")

    fh = open(out + ".calls.tsv", "w")
    fh.write("read\tlen\t" + "\t".join(names) + "\tcall\n")
    piles = {n: open(f"{out}.{n}.ids", "w") for n in names}
    tally = collections.Counter()
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip().upper()
            f.readline(); f.readline()
            ks = list(kmers(seq))
            if not ks:
                continue
            counts = {n: sum(1 for km in ks if km in sets[n]) for n in names}
            best = max(counts.values())
            if best < MIN_HITS:
                continue
            won = [n for n in names if counts[n] >= best * (1 - MARGIN)]
            rid = h[1:].split()[0]
            for n in won:
                piles[n].write(rid + "\n")
            tally["+".join(sorted(won))] += 1
            fh.write(f"{rid}\t{len(seq)}\t"
                     + "\t".join(str(counts[n]) for n in names)
                     + "\t" + "+".join(sorted(won)) + "\n")
    fh.close()
    for p in piles.values():
        p.close()
    tot = sum(tally.values())
    for k, v in tally.most_common():
        sys.stderr.write(f"  {k:>28}: {v:>6}  ({100*v/tot:.1f} %)\n")


main()
