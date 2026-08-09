#!/usr/bin/env python3
"""
Fast pileup + read-by-site matrix, for one reference sequence.

Pass 1  walk every alignment and tally the base each read contributes to each
        reference position (indels ignored -- they are the dominant HiFi error
        mode in homopolymers and are not used as haplotype markers).
Pass 2  keep positions where a second allele reaches MIN_ALT_COUNT reads and
        MIN_ALT_FRACTION of the depth: these are the informative sites.

Usage:  samtools view <bam> <ref> | python readvars.py <ref.fa> <refname> <outprefix>
Writes  <outprefix>.sites.tsv   pos refBase altBase depth altCount altFrac
        <outprefix>.matrix.tsv  read  then one column per site (R/A/.)
"""
import sys, re, collections

MIN_DEPTH = 20
MIN_ALT_COUNT = 5
MIN_ALT_FRACTION = 0.10

CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def load_ref(path, want):
    seq, name = [], None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif name == want:
            seq.append(line.strip())
    return "".join(seq).upper()


def main():
    ref = load_ref(sys.argv[1], sys.argv[2])
    out = sys.argv[3]
    reads = []                                   # (name, {pos: base})
    depth = collections.Counter()
    alt = collections.defaultdict(collections.Counter)

    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 0x100:                         # drop secondary; keep supplementary
            continue
        if int(f[4]) < 20:                       # MAPQ filter
            continue
        pos, cigar, seq, qual = int(f[3]) - 1, f[5], f[9], f[10]
        calls = {}
        r, q = pos, 0
        for n, op in CIG.findall(cigar):
            n = int(n)
            if op in "M=X":
                for k in range(n):
                    if ord(qual[q + k]) - 33 >= 20:
                        b = seq[q + k]
                        calls[r + k] = b
                        depth[r + k] += 1
                        if b != ref[r + k]:
                            alt[r + k][b] += 1
                r += n
                q += n
            elif op == "I":
                q += n
            elif op in "DN":
                r += n
            elif op == "S":
                q += n
        reads.append((f[0], calls))

    sites = []
    for p, c in sorted(alt.items()):
        d = depth[p]
        if d < MIN_DEPTH:
            continue
        b, n = c.most_common(1)[0]
        if n >= MIN_ALT_COUNT and n / d >= MIN_ALT_FRACTION:
            sites.append((p, ref[p], b, d, n, n / d))

    with open(out + ".sites.tsv", "w") as fh:
        fh.write("pos\tref\talt\tdepth\taltCount\taltFrac\n")
        for p, rb, ab, d, n, fr in sites:
            fh.write(f"{p+1}\t{rb}\t{ab}\t{d}\t{n}\t{fr:.3f}\n")

    with open(out + ".matrix.tsv", "w") as fh:
        fh.write("read\t" + "\t".join(str(p + 1) for p, *_ in sites) + "\n")
        for name, calls in reads:
            row = []
            for p, rb, ab, *_ in sites:
                b = calls.get(p)
                row.append("R" if b == rb else "A" if b == ab else ".")
            if any(x != "." for x in row):
                fh.write(name + "\t" + "\t".join(row) + "\n")

    sys.stderr.write(f"{sys.argv[2]}: {len(sites)} informative sites, "
                     f"{len(reads)} reads, {len(sites)/len(ref)*1000:.2f}/kb\n")


main()
