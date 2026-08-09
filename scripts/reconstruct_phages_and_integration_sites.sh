#!/bin/bash

################################################################################
# Reconstruction of both APSE bacteriophages in a two-phage sample, and of the
# chromosomal sites into which they are integrated
#
# BACKGROUND
#
# Two samples in this panel, S05 and S07, carry two different APSE phages rather
# than one. The 2024 curation assumed a single phage per sample and deleted the
# contigs of the second one as redundant duplicates, because APSE genomes share
# a highly conserved backbone and a short contig of one phage aligns inside the
# other. Recovering the deleted phage cannot be done by undoing the deletion:
# the discarded contigs were themselves fragments, and the reads of the two
# phages had already been mixed by the assembler.
#
# METHOD
#
# The reads are separated first and assembled afterwards. A read is assigned to
# a phage using 31-mers that occur in exactly one of the two reference
# templates. Those k-mers are read straight out of the raw sequence, so the
# assignment does not depend on an alignment and cannot be biased by a chimeric
# template. Reads carrying no diagnostic k-mer of either phage are placed in
# BOTH piles, because the two phages are identical wherever such a read sits and
# withholding those reads would fragment both assemblies. Reads carrying
# diagnostic k-mers of both phages in two spatially segregated blocks are
# recombinant molecules; they belong to neither parental lineage and are held
# out and counted separately.
#
# When one of the lineages is itself a RECOMBINANT of the others, as in S05,
# unique k-mers barely exist for it -- almost every k-mer it carries is shared
# with one parent or the other. Use scripts/lib/kmer_assign.py instead, which
# asks how many of a read's k-mers each reference contains and gives the read to
# whichever contains most. A mosaic matches along its whole length while each
# parent matches only part of it, so the mosaic wins on its own reads.
#
# THE TEMPLATE IS A READ CLASSIFIER, NEVER A SEQUENCE SOURCE. Each pile is then
# assembled de novo, so no base of any template can reach the output. The
# template only has to SPAN the diversity in the sample; it does not have to be
# free of chimerism. This distinction matters: an earlier attempt used
# chimera-free templates taken from single-phage samples, but both of those
# carry a type-A integrase, so integrase-B reads had no matching template, were
# misassigned, and produced a 56,981 bp assembly carrying both integrases.
#
# The procedure is iterated. The assemblies of round N become the classifier of
# round N+1, and the run stops when the read assignment stops changing. In S07
# it converged after one iteration, with 99.8 % of reads receiving the same call
# from the round-4 and round-5 references.
#
# INTEGRATION SITES
#
# Reads that align partly to the chromosome and partly to a phage are collected
# from the SA tags of the same alignment. Their breakpoints cluster to single
# bases and mark the boundaries of the prophage. The two chromosomal
# breakpoints of one site are the two ends of the att core, the direct repeat
# that integration duplicates; the exact core is recovered as the longest exact
# match between the chromosomal and the phage sequence around the breakpoints.
# Reads that run continuously across the empty att site record the excised state
# and give the lysogenic-to-excised ratio.
#
# ACCEPTANCE
#
# A reconstruction is accepted when the pile it came from produces no
# heterozygous site against it. A clonal element sequenced with HiFi reads is
# monomorphic: sample S12, which genuinely carries one phage, gives 0.00 sites
# per kb. Both S07 phages reach that benchmark, whereas the phage reference
# delivered in 2024 gives 6.94 sites per kb because it is a chimera of the two.
#
# USAGE
#
#   bash scripts/reconstruct_phages_and_integration_sites.sh S07 contig_02 \
#        data/reference/S07_initial_templates.fasta
#
# Arguments: sample ID, name of the chromosome contig in the curated assembly,
# and a two-record FASTA of starting templates. The starting templates are the
# sample's own v1 contigs grouped by integrase type, concatenated per phage.
# Edit data/config.sh first. Run from the root of the repository.
#
# DEPENDENCIES
#
#   minimap2   read and sequence alignment       (tested with 2.31-r1302)
#   samtools   alignment handling                (tested with 1.24)
#   seqtk      read subsetting and subsampling   (tested with 1.5)
#   flye       de novo assembly of each pile     (tested with 2.9.6)
#   python3    the helpers in scripts/lib        (3.8 or later, no modules needed)
################################################################################

set -euo pipefail

SAMPLE="${1:?usage: $0 SAMPLE CHROM_CONTIG TEMPLATES.fasta}"
CHROM_CONTIG="${2:?}"
TEMPLATES="${3:?}"

source data/config.sh
LIB="scripts/lib"
OUT="${WORK_DIR%/*}/reconstruct_${SAMPLE}"
mkdir -p "$OUT"

# Target coverage for assembly. The full phage extraction in this project runs
# to several thousand fold, which crashed hifiasm_meta twice and destabilised
# the machine; Flye is also slower and no more accurate at that depth. 150-fold
# is far above what Flye needs and keeps the run to a few minutes.
TARGET_DEPTH=150

# Maximum number of classify-assemble rounds. Convergence has been reached in
# one round on both two-phage samples; the cap only guards against oscillation.
MAX_ROUNDS=4

READS="${READS_DIR}/${SAMPLE}.fq"
ASM="${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.curated.fasta"

################################################################################
# Pull the chromosome out of the curated assembly. This is the sequence the
# integration coordinates will refer to, so it must be the delivered one.
################################################################################
awk -v c="$CHROM_CONTIG" '/^>/{p=($0 ~ c)} p' "$ASM" \
  | awk -v s="$SAMPLE" 'NR==1{print ">"s"_chr"; next}{print}' > "$OUT/chr.fa"

################################################################################
# Round 0 uses the supplied templates; every later round uses the assemblies of
# the round before it. Each round classifies the reads, assembles each pile and
# compares the classification with the previous one.
################################################################################
cp "$TEMPLATES" "$OUT/ref_0.fa"

for R in $(seq 1 "$MAX_ROUNDS"); do
  PREV=$((R - 1))

  # Reads that touch either phage at all. Everything else is chromosome or
  # plasmid and is irrelevant here.
  minimap2 -t "$THREADS" -cx map-hifi "$OUT/ref_${PREV}.fa" "$READS" 2>/dev/null \
    | cut -f1 | sort -u > "$OUT/phage_reads.ids"
  seqtk subseq "$READS" "$OUT/phage_reads.ids" > "$OUT/phage_reads.fq"

  # Assign every read to a phage from parent-diagnostic 31-mers.
  python3 "$LIB/kmer_ancestry.py" "$OUT/ref_${PREV}.fa" "$OUT/phage_reads.fq" \
          "$OUT/ancestry_${R}.tsv"

  # Build the two piles. A read whose minority parent contributes less than
  # 15 % of its diagnostic k-mers is a parental read with a few stray k-mers
  # from sequencing error; a read with no diagnostic k-mer at all is
  # uninformative and goes to both piles; a mosaic read is held out.
  python3 - "$OUT" "$R" <<'PY'
import sys
out, r = sys.argv[1], sys.argv[2]
seen = set()
with open(f"{out}/pile_A.ids", "w") as a, open(f"{out}/pile_B.ids", "w") as b, \
     open(f"{out}/held_{r}.ids", "w") as h:
    names = None
    for line in open(f"{out}/ancestry_{r}.tsv"):
        f = line.rstrip("\n").split("\t")
        if names is None:
            names = [f[2][2:], f[3][2:]]
            continue
        rid, n0, n1, call = f[0], int(f[2]), int(f[3]), f[5]
        seen.add(rid)
        tot, mn = n0 + n1, min(n0, n1)
        if call == "mosaic" or (tot and mn / tot >= 0.15):
            h.write(rid + "\n")
        else:
            (a if n0 > n1 else b).write(rid + "\n")
    for rid in open(f"{out}/phage_reads.ids"):
        if rid.strip() not in seen:          # no diagnostic k-mer at all
            a.write(rid); b.write(rid)
    open(f"{out}/pile_names.txt", "w").write("\t".join(names) + "\n")
PY
  read -r NAME_A NAME_B < "$OUT/pile_names.txt"

  # Assemble each pile de novo at a workable depth.
  for P in A B; do
    seqtk subseq "$READS" "$OUT/pile_${P}.ids" > "$OUT/pile_${P}.fq"
    BASES=$(awk 'NR%4==2{b+=length($0)} END{print b}' "$OUT/pile_${P}.fq")
    LEN=$(awk '/^>/{next}{n+=length($0)} END{print n}' "$OUT/ref_${PREV}.fa")
    FRAC=$(awk -v b="$BASES" -v l="$LEN" -v d="$TARGET_DEPTH" \
             'BEGIN{f=d*(l/2)/b; print (f>1?1:f)}')
    seqtk sample -s7 "$OUT/pile_${P}.fq" "$FRAC" > "$OUT/pile_${P}.${TARGET_DEPTH}x.fq"
    rm -rf "$OUT/flye_${R}_${P}"
    flye --pacbio-hifi "$OUT/pile_${P}.${TARGET_DEPTH}x.fq" \
         --out-dir "$OUT/flye_${R}_${P}" --threads "$THREADS" --meta \
         > "$OUT/flye_${R}_${P}.log" 2>&1
  done

  # Keep the longest circular contig from each pile; anything else is
  # chromosomal carry-over and is reported but not used.
  python3 - "$OUT" "$R" "$NAME_A" "$NAME_B" <<'PY'
import sys
out, r, na, nb = sys.argv[1:5]
def best(d):
    info = [l.split("\t") for l in open(f"{d}/assembly_info.txt").read().splitlines()[1:]]
    circ = [x for x in info if x[3] == "Y"] or info
    return max(circ, key=lambda x: int(x[1]))[0]
seqs = {}
for pile, name in (("A", na), ("B", nb)):
    d = f"{out}/flye_{r}_{pile}"
    keep, cur = best(d), None
    for line in open(f"{d}/assembly.fasta"):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            if cur == keep: seqs[name] = []
        elif cur == keep: seqs[name].append(line.strip())
with open(f"{out}/ref_{r}.fa", "w") as fh:
    for n, s in seqs.items():
        s = "".join(s)
        fh.write(f">{n}\n" + "\n".join(s[i:i+60] for i in range(0, len(s), 60)) + "\n")
        print(f"  round {r} {n}: {len(s):,} bp")
PY

  # Convergence: has any read changed pile?
  if [ "$R" -gt 1 ]; then
    CHANGED=$(join -t$'\t' \
      <(awk -F'\t' 'NR>1{print $1"\t"$6}' "$OUT/ancestry_$((R-1)).tsv" | sort) \
      <(awk -F'\t' 'NR>1{print $1"\t"$6}' "$OUT/ancestry_${R}.tsv" | sort) \
      | awk -F'\t' '$2!=$3' | wc -l)
    echo "  round ${R}: ${CHANGED} reads changed classification"
    [ "$CHANGED" -lt 50 ] && break
  fi
done
FINAL="$OUT/ref_${R}.fa"

################################################################################
# Integration sites. Map every read against the chromosome and the two finished
# phages together, then read the chromosome-phage junctions off the SA tags.
################################################################################
cat "$OUT/chr.fa" "$FINAL" > "$OUT/combined.fa"
samtools faidx "$OUT/combined.fa"
minimap2 -t "$THREADS" -ax map-hifi --secondary=no -Y "$OUT/combined.fa" "$READS" 2>/dev/null \
  | samtools sort -@2 -o "$OUT/all.bam" -
samtools index "$OUT/all.bam"
samtools view -F 0x900 "$OUT/all.bam" \
  | python3 "$LIB/split_read_junctions.py" > "$OUT/junctions.tsv"

echo "chromosome-side breakpoints (read count, position, phage):"
awk -v c="${SAMPLE}_chr" '$2==c{print $3"\t"$5} $5==c{print $6"\t"$2}' "$OUT/junctions.tsv" \
  | sort | uniq -c | sort -rn | awk '$1>=5{printf "  %6d  %10d  %s\n",$1,$2,$3}'

echo
echo "Junction table written to $OUT/junctions.tsv."
echo "Take the two chromosomal breakpoints of each site as the att core, confirm"
echo "the core is an exact match to the phage, record them in"
echo "results/${SAMPLE}_integration_sites.tsv, then build the integrated"
echo "chromosome with:"
echo
echo "  python3 $LIB/integrate_prophage.py $OUT/chr.fa $FINAL \\"
echo "          results/${SAMPLE}_integration_sites.tsv assemblies_v2/${SAMPLE}_v2_chromosome.fa"
