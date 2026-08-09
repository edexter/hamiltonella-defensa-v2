#!/usr/bin/env python3
"""
Compare each disjoint-replicate assembly with the shipped reference.

Circular assemblies are compared rotation-invariantly: the replicate is
doubled head-to-tail before alignment, so a genome that is identical but
started at a different base still aligns as one block. Scoring byte identity
directly would report such a pair as different when it is not.

Differences are split into two classes, because they mean different things:
  homopolymer   inside or adjacent to a run of >=3 identical bases. This is
                HiFi's known error mode and is expected at 150x.
  substantive   everything else. These are the ones that would matter.

USAGE

    cd <the replicate working directory, which holds asm/>
    python3 <repo>/scripts/compare_replicates.py > replicate_comparison.tsv

The directory holding the delivered reference genomes can be given in the
environment variable PHAGE_REFERENCE_DIR; it defaults to the path used by this
project. Run scripts/run_replicate_test.sh first to produce the assemblies.
"""

import glob
import os
import re
import subprocess
import sys

# Where the delivered reference genomes live, relative to the repository root.
REF = os.environ.get("PHAGE_REFERENCE_DIR", "results/phage_reconstructions")
SHIPPED = {
    "S07_RHS":  f"{REF}/S07_APSE_RHS_v2.fasta",
    "S07_CdtB": f"{REF}/S07_APSE_CdtB_v2.fasta",
    "S05_RHS":  f"{REF}/S05_APSE_RHS_v2.fasta",
    "S05_hapA": f"{REF}/S05_APSE_CdtB_hapA_v2.fasta",
    "S05_hapR": f"{REF}/S05_APSE_CdtB_hapR_v2.fasta",
}
CS = re.compile(r"([:*+-])([A-Za-z0-9]+)")


def read_fa(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def longest_contig(d):
    best = max(d.items(), key=lambda kv: len(kv[1]))
    return best[1]


def homopolymer_context(seq, pos, k=3):
    """is position pos inside or beside a run of >=k identical bases?"""
    lo, hi = max(0, pos - 6), min(len(seq), pos + 7)
    w = seq[lo:hi]
    return bool(re.search(r"(A{%d,}|C{%d,}|G{%d,}|T{%d,})" % (k, k, k, k), w))


def compare(refseq, qryseq, tmp):
    """align qry (doubled) against ref; return aligned bp, and difference lists"""
    with open(tmp + ".ref.fa", "w") as fh:
        fh.write(">ref\n" + refseq + "\n")
    with open(tmp + ".qry.fa", "w") as fh:
        fh.write(">qry\n" + qryseq + qryseq + "\n")     # doubled: rotation-proof
    out = subprocess.run(
        ["minimap2", "-cx", "asm5", "--cs", "-N", "5", "-p", "0.1",
         tmp + ".ref.fa", tmp + ".qry.fa"],
        capture_output=True, text=True).stdout

    covered = set()
    subs, indels = [], []
    for line in out.split("\n"):
        f = line.split("\t")
        if len(f) < 12:
            continue
        tstart = int(f[7])
        cs = next((x[5:] for x in f[12:] if x.startswith("cs:Z:")), None)
        if cs is None:
            continue
        p = tstart
        for op, val in CS.findall(cs):
            if op == ":":
                n = int(val)
                covered.update(range(p, p + n))
                p += n
            elif op == "*":
                subs.append(p)
                covered.add(p)
                p += 1
            elif op == "-":                      # deletion from the reference
                indels.append(p)
                covered.update(range(p, p + len(val)))
                p += len(val)
            elif op == "+":                      # insertion relative to reference
                indels.append(p)
    return covered, subs, indels


def main():
    print("lineage\treplicate\tlen\tcircular\tcontigs\tdelta_len\tref_covered_pct"
          "\tsubs\tindels\thomopolymer\tsubstantive")
    for tag, ship in sorted(SHIPPED.items()):
        refseq = longest_contig(read_fa(ship))
        for d in sorted(glob.glob(f"asm/{tag}_rep*")):
            info = os.path.join(d, "assembly_info.txt")
            fa = os.path.join(d, "assembly.fasta")
            rep = os.path.basename(d).split("_rep")[-1]
            if not os.path.exists(info):
                print(f"{tag}\t{rep}\tFAILED\t-\t-\t-\t-\t-\t-\t-\t-")
                continue
            rows = [l.split("\t") for l in open(info).read().splitlines()[1:]]
            rows.sort(key=lambda r: -int(r[1]))
            circ = rows[0][3]
            ncontig = len(rows)
            qryseq = longest_contig(read_fa(fa))
            covered, subs, indels = compare(refseq, qryseq, f"/tmp/cmp_{tag}_{rep}")
            hp = sum(1 for p in subs + indels if homopolymer_context(refseq, p))
            substantive = len(subs) + len(indels) - hp
            print(f"{tag}\t{rep}\t{len(qryseq)}\t{circ}\t{ncontig}\t"
                  f"{len(qryseq)-len(refseq):+d}\t{100*len(covered)/len(refseq):.2f}\t"
                  f"{len(subs)}\t{len(indels)}\t{hp}\t{substantive}")


main()
