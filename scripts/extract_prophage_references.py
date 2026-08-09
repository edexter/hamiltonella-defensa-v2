#!/usr/bin/env python3
"""
Extract a standalone APSE reference from a sample whose prophage is integrated
in the chromosome, and circularise it only where the reads say the circle exists.

An integrated prophage is flanked by a direct repeat -- the att core, duplicated
when the phage integrated. The repeat is found empirically here rather than
assumed from a known sequence, because the core is not identical in every
sample. Excising between the two copies, keeping one, reconstitutes the phage
genome as it would be when circular.

Whether that circle actually EXISTS in the sample is then a separate question,
and it is answered with reads, not with an assumption. In the chromosome the two
ends of the excised sequence are separated by the whole chromosome, so no read
can span them. In a free circular phage they are adjacent. Reads spanning the
join are therefore direct evidence of the circular form; their absence means the
phage is only ever integrated in this sample, and the reference ships linear.

A good incomplete genome is better than a forced complete one.

Usage: python3 scripts/extract_prophage_references.py <sample> <chrom.fa> <probe.fa> <reads.fq> <outdir>
"""

import os
import subprocess
import sys

MIN_REPEAT = 25          # shortest credible att core
MIN_PHAGE, MAX_PHAGE = 30000, 46000
FLANK = 2500             # window each side of the alignment boundary to search
SPAN = 300               # a read must reach this far past the join to count


def read_fa(path):
    d, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            d[name] = []
        else:
            d[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def longest_shared(left, right):
    """longest exact substring common to both windows, with its offsets"""
    K = MIN_REPEAT
    idx = {}
    for i in range(len(right) - K + 1):
        idx.setdefault(right[i:i + K], []).append(i)
    best = (0, -1, -1)
    for i in range(len(left) - K + 1):
        for j in idx.get(left[i:i + K], []):
            n = K
            while (i + n < len(left) and j + n < len(right)
                   and left[i + n] == right[j + n]):
                n += 1
            if n > best[0]:
                best = (n, i, j)
    return best


def locate(chrom_fa, probe_fa):
    """merged interval on the chromosome covered by same-lineage phage sequence"""
    out = subprocess.run(["minimap2", "-cx", "asm20", "-N", "20", "-p", "0.05",
                          chrom_fa, probe_fa], capture_output=True, text=True).stdout
    blocks = {}
    for line in out.split("\n"):
        f = line.split("\t")
        if len(f) < 12 or int(f[10]) < 2000:
            continue
        blocks.setdefault(f[5], []).append((int(f[7]), int(f[8])))
    best = None
    for c, iv in blocks.items():
        iv.sort()
        cur_s, cur_e = iv[0]
        for s, e in iv[1:]:
            if s > cur_e + 8000:
                if best is None or cur_e - cur_s > best[2] - best[1]:
                    best = (c, cur_s, cur_e)
                cur_s, cur_e = s, e
            else:
                cur_e = max(cur_e, e)
        if best is None or cur_e - cur_s > best[2] - best[1]:
            best = (c, cur_s, cur_e)
    return best


def main():
    samp, chrom_fa, probe_fa, reads, outdir = sys.argv[1:6]
    os.makedirs(outdir, exist_ok=True)
    chrom = read_fa(chrom_fa)

    loc = locate(chrom_fa, probe_fa)
    if loc is None:
        print(f"{samp}\tNO_PHAGE_FOUND")
        return
    cname, a, b = loc
    seq = chrom[cname]
    left = seq[max(0, a - FLANK):a + FLANK]
    right = seq[max(0, b - FLANK):min(len(seq), b + FLANK)]
    n, li, ri = longest_shared(left, right)

    if n < MIN_REPEAT:
        print(f"{samp}\tNO_att_REPEAT\t{cname}:{a}-{b}\t"
              f"aligned {b-a} bp\tno direct repeat found; not excised")
        return
    attL = max(0, a - FLANK) + li
    attR = max(0, b - FLANK) + ri
    phage = seq[attL:attR]                       # one copy of the core retained
    if not (MIN_PHAGE <= len(phage) <= MAX_PHAGE):
        print(f"{samp}\tBAD_LENGTH\t{cname}:{attL}-{attR}\t{len(phage)} bp\t"
              f"att repeat {n} bp; outside {MIN_PHAGE}-{MAX_PHAGE}, not shipped")
        return

    lin = os.path.join(outdir, f"{samp}_APSE_linear.fa")
    with open(lin, "w") as fh:
        fh.write(f">{samp}_APSE\n")
        for i in range(0, len(phage), 70):
            fh.write(phage[i:i + 70] + "\n")

    # does the circular form exist? rotate by half so the join sits mid-sequence,
    # then ask whether any read crosses it
    half = len(phage) // 2
    rot = phage[half:] + phage[:half]
    rotfa = os.path.join(outdir, f"{samp}_rot.fa")
    with open(rotfa, "w") as fh:
        fh.write(f">{samp}_rot\n")
        for i in range(0, len(rot), 70):
            fh.write(rot[i:i + 70] + "\n")
    bam = os.path.join(outdir, f"{samp}_rot.bam")
    subprocess.run(f"minimap2 -ax map-hifi -t 6 --secondary=no {rotfa} {reads} 2>/dev/null "
                   f"| samtools sort -o {bam} - && samtools index {bam}",
                   shell=True, capture_output=True)
    join = len(phage) - half                    # where the original ends meet
    spanning = 0
    out = subprocess.run(["samtools", "view", bam, f"{samp}_rot:{join-SPAN}-{join+SPAN}"],
                         capture_output=True, text=True).stdout
    import re
    CIG = re.compile(r"(\d+)([MIDNSHP=X])")
    for line in out.split("\n"):
        f = line.split("\t")
        if len(f) < 6 or int(f[1]) & 0x904:
            continue
        pos, ref = int(f[3]), int(f[3])
        for c, op in CIG.findall(f[5]):
            c = int(c)
            if op in "M=X":
                if ref <= join - SPAN and ref + c - 1 >= join + SPAN:
                    spanning += 1
                ref += c
            elif op in "DN":
                ref += c
    verdict = "CIRCULAR" if spanning >= 5 else "LINEAR_ONLY"
    print(f"{samp}\t{verdict}\t{cname}:{attL}-{attR}\t{len(phage)} bp\t"
          f"att repeat {n} bp\t{spanning} reads span the join")


main()
