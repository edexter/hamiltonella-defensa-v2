#!/bin/bash

################################################################################
# Structural audit of the delivered assemblies: are the joins real?
#
# WHY
#
# Three kinds of structural claim in this delivery were made by comparing
# assemblies with each other rather than by looking at reads, and comparison
# between assemblies generates hypotheses well but settles nothing on its own.
# Twice in this project a conclusion reached that way has been overturned once
# reads were consulted. This script consults the reads.
#
# The three claims:
#
#   scaffold joins   Six samples have a chromosome that the 2024 curation built
#                    by joining several assembler contigs on the strength of
#                    BLAST overlaps between their ends. If any of those joins is
#                    wrong, the delivered chromosome contains a rearrangement
#                    that no amount of annotation will reveal. Only S01 has been
#                    checked.
#
#   plasmid fusions  A contig delivered as one plasmid may be two plasmids the
#                    assembler joined. The tell is read depth: two plasmids at
#                    different copy number, fused, leave a step in coverage at
#                    the point where they meet. A genuine single replicon does
#                    not.
#
#   prophage sites   Several samples are documented as having too few reads
#                    crossing a phage-to-chromosome junction to locate the
#                    integration site. That count came from an earlier pass which
#                    is known to have undercounted at least one sample fourfold,
#                    so it is recomputed here from split alignments.
#
# WHAT COUNTS AS SUPPORT
#
# For a join, a read that begins well before it and ends well after it. Reads
# that merely touch the junction prove nothing, because a read ending at a
# junction is what a false join looks like. The script therefore requires the
# read to extend a set distance either side, and reports the count at two
# distances so that a join resting on marginal reads is visible as such.
#
# A depth step is reported as the largest ratio between the mean depth left and
# right of any interior position, ignoring the contig ends where mapping is
# unreliable. A ratio near one means one replicon; a clear step means two.
#
# USAGE
#
#   bash scripts/audit_assembly_structure.sh                # every sample
#   bash scripts/audit_assembly_structure.sh S04 S06        # named samples
#
# Edit data/config.sh first. Run from the root of the repository. Roughly ten to
# twenty minutes per sample, dominated by read mapping. Alignments are deleted as
# it goes because they are large.
#
# OUTPUT
#
#   results/scaffold_join_audit.tsv    one row per scaffold join
#   results/contig_link_audit.tsv      split-read links between contigs
#   results/depth_step_audit.tsv       one row per plasmid-sized contig
#
# DEPENDENCIES
#
#   minimap2, samtools, python3 (standard library only)
################################################################################

set -o pipefail
source data/config.sh

# A read must extend at least this far either side of a join to count as
# spanning it. The larger value is reported alongside so that a join supported
# only by reads that barely reach across is visible.
SPAN_NEAR=150
SPAN_FAR=1000

# Contigs between these sizes are screened for a depth step. Below the floor the
# depth estimate is too noisy; above the ceiling we are looking at chromosome,
# where copy-number steps mean something different.
STEP_MIN_BP=10000
STEP_MAX_BP=200000

WORK="${WORK_DIR%/*}/structure_audit"
mkdir -p "$WORK" results

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
	SAMPLES=(S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 S12 S13 S14)
fi

printf "sample\tjoin_position\tleft_piece\tright_piece\tspanning_reads_150bp\tspanning_reads_1kb\tverdict\n" > results/scaffold_join_audit.tsv
printf "sample\tcontig_a\tcontig_b\tsplit_reads\tmedian_pos_a\tmedian_pos_b\n" > results/contig_link_audit.tsv
printf "sample\tcontig\tlength_bp\tmean_depth\tbest_step_position\tleft_depth\tright_depth\tstep_ratio\tverdict\n" > results/depth_step_audit.tsv

for SAMPLE in "${SAMPLES[@]}"; do
	ASSEMBLY="assemblies_v2/${SAMPLE}_v2.fasta"
	READS="${READS_DIR}/${SAMPLE}.fq"
	# The pre-curation assembly. Four samples name this file "assembly.fasta"
	# rather than "<sample>.raw.fasta", and an earlier version of this script
	# looked only for the latter -- so S07, S09, S10 and S14 were skipped in
	# silence and reported as having no scaffold joins to check. Both names are
	# now tried, and a sample with neither is reported rather than passed over.
	RAW="${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.raw.fasta"
	[ -s "$RAW" ] || RAW="${ASSEMBLY_DIR}/${SAMPLE}/assembly.fasta"
	[ -s "$RAW" ] || echo "  $SAMPLE: no pre-curation assembly found, scaffold joins NOT checked" >&2
	[ -s "$ASSEMBLY" ] && [ -s "$READS" ] || { echo "  $SAMPLE: inputs missing, skipping" >&2; continue; }

	echo "  $SAMPLE: mapping reads"
	minimap2 -ax map-hifi -t "$THREADS" "$ASSEMBLY" "$READS" 2>/dev/null \
		| samtools sort -@ 4 -o "$WORK/${SAMPLE}.bam" - 2>/dev/null
	samtools index "$WORK/${SAMPLE}.bam"

	# ---------------------------------------------------------- scaffold joins
	# Find where each pre-curation assembler contig lands in the delivered
	# chromosome. Where two of them abut, the 2024 curation made a join, and that
	# join is what needs testing.
	if [ -s "$RAW" ]; then
		echo "  $SAMPLE: checking scaffold joins"
		minimap2 -cx asm5 -t "$THREADS" --secondary=no "$ASSEMBLY" "$RAW" 2>/dev/null \
			> "$WORK/${SAMPLE}.raw_vs_v2.paf"
		python3 scripts/lib/check_scaffold_joins.py "$SAMPLE" \
			"$WORK/${SAMPLE}.raw_vs_v2.paf" "$WORK/${SAMPLE}.bam" \
			"$SPAN_NEAR" "$SPAN_FAR" >> results/scaffold_join_audit.tsv
	fi

	# ------------------------------------------------------ split-read linkage
	# Reads aligned partly to one contig and partly to another. This is how a
	# prophage sitting between two contigs announces itself, and it is the count
	# that the earlier pass got wrong.
	echo "  $SAMPLE: counting split-read links between contigs"
	samtools view "$WORK/${SAMPLE}.bam" \
		| python3 scripts/lib/count_contig_links.py "$SAMPLE" \
		>> results/contig_link_audit.tsv

	# ------------------------------------------------------------- depth steps
	echo "  $SAMPLE: screening for depth steps"
	samtools depth -a "$WORK/${SAMPLE}.bam" \
		| python3 scripts/lib/find_depth_steps.py "$SAMPLE" "$STEP_MIN_BP" "$STEP_MAX_BP" \
		>> results/depth_step_audit.tsv

	rm -f "$WORK/${SAMPLE}.bam" "$WORK/${SAMPLE}.bam.bai"
done

echo
echo "Scaffold joins  : results/scaffold_join_audit.tsv"
echo "Contig linkage  : results/contig_link_audit.tsv"
echo "Depth steps     : results/depth_step_audit.tsv"
echo
echo "A join with few or no spanning reads, or a contig with a clear depth step,"
echo "is a candidate error and should be inspected by hand before delivery."
