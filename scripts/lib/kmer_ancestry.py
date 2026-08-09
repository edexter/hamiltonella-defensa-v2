#!/usr/bin/env python3
"""
Assembly-independent recombination test.

Build the set of canonical 31-mers that occur in exactly one of the two
reconstructed S07 phages, then walk every read and record the order in which
those parent-diagnostic k-mers appear along the molecule.

A read from a parental phage carries k-mers of one parent only. A read from a
recombinant molecule carries both, in two spatially segregated blocks. Mapping
artefacts cannot produce that pattern: the k-mers are read from the raw
sequence and never touch a reference alignment.

Usage: python kmer_ancestry.py <phages.fa> <reads.fq> <out.tsv>
"""
import sys, collections

K = 31
MIN_PER_PARENT = 15          # k-mers required from each parent to call a read mosaic
MIN_FRACTION = 0.15          # ... and each parent must be this share of the total
MAX_TRANSITIONS = 2          # segregated blocks, not interleaved noise

COMP = str.maketrans("ACGT", "TGCA")


def canon(s):
    r = s.translate(COMP)[::-1]
    return s if s < r else r


def kmers(seq):
    for i in range(len(seq) - K + 1):
        km = seq[i:i + K]
        if "N" in km:
            continue
        yield i, canon(km)


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
    assert len(names) == 2, names

    sets = {}
    for n in names:
        s = seqs[n]
        sets[n] = {km for _, km in kmers(s + s[:K - 1])}   # circular
    diag = {}
    for i, n in enumerate(names):
        other = names[1 - i]
        for km in sets[n] - sets[other]:
            diag[km] = n
    sys.stderr.write(f"diagnostic 31-mers: {names[0]} "
                     f"{sum(1 for v in diag.values() if v==names[0])}, "
                     f"{names[1]} {sum(1 for v in diag.values() if v==names[1])}\n")

    fh = open(out, "w")
    fh.write("read\tlen\tn_%s\tn_%s\ttransitions\tcall\tblocks\n" % tuple(names))
    tally = collections.Counter()
    with open(fq) as f:
        while True:
            h = f.readline()
            if not h:
                break
            seq = f.readline().strip().upper()
            f.readline(); f.readline()
            hits = [(i, diag[km]) for i, km in kmers(seq) if km in diag]
            if not hits:
                continue
            c = collections.Counter(p for _, p in hits)
            n0, n1 = c[names[0]], c[names[1]]
            tot = n0 + n1
            # run-length encode the parent labels along the read
            runs = []
            for pos, p in hits:
                if runs and runs[-1][0] == p:
                    runs[-1][2] = pos
                else:
                    runs.append([p, pos, pos])
            # drop runs supported by fewer than 5 k-mers -- isolated hits are noise
            big = [r for r in runs if sum(1 for pos, p in hits
                                          if r[1] <= pos <= r[2] and p == r[0]) >= 5]
            merged = []
            for r in big:
                if merged and merged[-1][0] == r[0]:
                    merged[-1][2] = r[2]
                else:
                    merged.append(list(r))
            trans = max(len(merged) - 1, 0)
            mn = min(n0, n1)
            if (mn >= MIN_PER_PARENT and mn / tot >= MIN_FRACTION
                    and 1 <= trans <= MAX_TRANSITIONS):
                call = "mosaic"
            elif n0 == 0 or n1 == 0:
                call = names[0] if n1 == 0 else names[1]
            else:
                call = "ambiguous"
            tally[call] += 1
            fh.write(f"{h[1:].split()[0]}\t{len(seq)}\t{n0}\t{n1}\t{trans}\t{call}\t"
                     + ";".join(f"{p}:{a}-{b}" for p, a, b in merged) + "\n")
    fh.close()
    tot = sum(tally.values())
    for k, v in tally.most_common():
        sys.stderr.write(f"  {k:>10}: {v:>6}  ({100*v/tot:.1f} %)\n")


main()
