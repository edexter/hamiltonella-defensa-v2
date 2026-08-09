#!/usr/bin/env python3
"""
Look for a step in read depth inside a contig.

Two plasmids that the assembler fused into one contig are usually present in the
cell at different copy numbers, so the fused contig carries a visible step in
coverage where they meet. A genuine single replicon has flat coverage apart from
noise. This screen finds the interior position that best splits a contig into two
regions of differing depth, and reports the ratio.

The outermost windows are ignored, because coverage always falls at a contig end
where reads cannot extend past it, and that fall is not a step.

Reads `samtools depth -a` on stdin.
Usage: find_depth_steps.py <sample> <min_bp> <max_bp>
"""

import statistics
import sys

WINDOW = 1000          # depth is averaged in windows of this size before testing
EDGE_SKIP = 5          # windows ignored at each end
FLAG_RATIO = 1.5       # a step at or above this is worth a look


def main():
    sample, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    current, depths = None, []

    def emit(name, vals):
        n = len(vals)
        if not (lo <= n <= hi):
            return
        win = [statistics.mean(vals[i:i + WINDOW])
               for i in range(0, n - WINDOW + 1, WINDOW)]
        if len(win) < 2 * EDGE_SKIP + 4:
            return
        core = win[EDGE_SKIP:-EDGE_SKIP]
        best = (1.0, 0, 0, 0)
        for i in range(2, len(core) - 2):
            left, right = statistics.mean(core[:i]), statistics.mean(core[i:])
            if min(left, right) <= 0:
                continue
            ratio = max(left, right) / min(left, right)
            if ratio > best[0]:
                best = (ratio, (EDGE_SKIP + i) * WINDOW, left, right)
        ratio, pos, left, right = best
        verdict = "possible fusion - inspect" if ratio >= FLAG_RATIO else "flat"
        print(f"{sample}\t{name}\t{n}\t{statistics.mean(vals):.1f}\t{pos}\t"
              f"{left:.1f}\t{right:.1f}\t{ratio:.2f}\t{verdict}")

    for line in sys.stdin:
        f = line.split("\t")
        if f[0] != current:
            if current is not None:
                emit(current, depths)
            current, depths = f[0], []
        depths.append(int(f[2]))
    if current is not None:
        emit(current, depths)


main()
