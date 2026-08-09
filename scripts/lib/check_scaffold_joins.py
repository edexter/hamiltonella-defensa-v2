#!/usr/bin/env python3
"""
Test the joins the 2024 curation made when it scaffolded several assembler
contigs into one chromosome.

The joins are not recorded anywhere, so they are recovered by aligning the
pre-curation assembly to the delivered one: wherever two different assembler
contigs land adjacently on the same delivered contig, a join was made there.
Each join is then tested against the reads.

Usage: check_scaffold_joins.py <sample> <raw_vs_v2.paf> <bam> <span_near> <span_far>
"""

import collections
import re
import subprocess
import sys

# A raw contig must align at least this well to be treated as a real placement
# rather than a stray repeat match.
MIN_ALIGNED_BP = 20000


def main():
    sample, paf, bam, near, far = sys.argv[1:6]
    near, far = int(near), int(far)

    # Where does each assembler contig land in the delivered assembly?
    placed = collections.defaultdict(list)
    for line in open(paf):
        f = line.split("\t")
        if int(f[10]) < MIN_ALIGNED_BP:
            continue
        placed[f[5]].append((int(f[7]), int(f[8]), f[0], int(f[10])))

    for target, hits in placed.items():
        if len(hits) < 2:
            continue
        hits.sort()
        # Consecutive pieces define a join wherever they come from different
        # assembler contigs. The join is placed at the end of the left piece.
        for (s1, e1, n1, _), (s2, e2, n2, _) in zip(hits, hits[1:]):
            if n1 == n2:
                continue
            join = e1 if e1 <= s2 else (e1 + s2) // 2

            reads = []
            out = subprocess.run(
                ["samtools", "view", bam, f"{target}:{max(1, join-20000)}-{join+20000}"],
                capture_output=True, text=True).stdout
            for rec in out.splitlines():
                g = rec.split("\t")
                pos = int(g[3])
                ln = sum(int(x) for x, t in re.findall(r"(\d+)([MDN=X])", g[5]))
                reads.append((pos, pos + ln))

            n_near = sum(1 for a, b in reads if a < join - near and b > join + near)
            n_far = sum(1 for a, b in reads if a < join - far and b > join + far)

            if n_far >= 10:
                verdict = "supported"
            elif n_far >= 3:
                verdict = "weakly supported"
            else:
                verdict = "UNSUPPORTED - inspect"
            print(f"{sample}\t{target}:{join}\t{n1}\t{n2}\t{n_near}\t{n_far}\t{verdict}")


main()
