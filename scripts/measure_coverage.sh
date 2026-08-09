#!/bin/bash

################################################################################
# Read depth for every delivered assembly, reported separately for the
# chromosome, the plasmids, and the APSE prophage
#
# WHY THE PROPHAGE IS SEPARATED OUT
#
# Quoting a single coverage figure for one of these samples is misleading. The
# APSE prophage is actively replicating in most of them, so the reads pile up
# over it at tens to hundreds of times the depth of the surrounding chromosome.
# In S01 the prophage sits at roughly 1300-fold against a chromosome at 33-fold.
# A mean taken across the whole assembly is therefore a number that describes
# neither the chromosome nor the phage: it is inflated by the phage and useless
# as a measure of how well the bacterial genome was sequenced.
#
# This script measures the two separately. The prophage interval is not assumed;
# it is located by aligning the reconstructed phage genome for that sample back
# against the assembly, so that the boundary comes from the sequence itself.
#
# WHY DEPTH IS RECOMPUTED RATHER THAN TAKEN FROM THE ASSEMBLER
#
# The assembler reports a coverage figure per contig, and for these samples it
# happens to agree closely with a direct measurement. It is nonetheless an
# estimate made during assembly, from the reads the assembler chose to use. This
# script maps every read back to the finished assembly and counts, so that the
# figure delivered to the reader describes the assembly actually being delivered.
#
# USAGE
#
#   bash scripts/measure_coverage.sh              # every sample
#   bash scripts/measure_coverage.sh S01 S05      # named samples only
#
# Edit data/config.sh first. Run from the root of the repository. Expect roughly
# five to fifteen minutes per sample depending on depth, and note that the
# intermediate alignment files are large and are deleted as it goes.
#
# OUTPUT
#
#   results/coverage_summary.tsv    one row per sample: chromosome, plasmid and
#                                   prophage depth, and the ratio between them
#   results/coverage_by_contig.tsv  one row per contig
#
# DEPENDENCIES
#
#   minimap2   read alignment                        (tested with 2.31-r1302)
#   samtools   sorting and per-base depth            (tested with 1.24)
#   awk        summarising                           (POSIX awk)
################################################################################

set -o pipefail

source data/config.sh

WORK="${WORK_DIR%/*}/coverage"
BY_CONTIG="results/coverage_by_contig.tsv"
SUMMARY="results/coverage_summary.tsv"
mkdir -p "$WORK" results

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
	SAMPLES=(S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 S12 S13 S14)
fi

printf "sample\tcontig\tlength_bp\tmean_depth\n" > "$BY_CONTIG"
printf "sample\tchromosome_depth\tprophage_depth\tprophage_ratio\tplasmid_depth_range\tprophage_bp\n" > "$SUMMARY"

for SAMPLE in "${SAMPLES[@]}"; do
	ASSEMBLY="assemblies_v2/${SAMPLE}_v2.fasta"
	READS="${READS_DIR}/${SAMPLE}.fq"

	if [ ! -s "$ASSEMBLY" ] || [ ! -s "$READS" ]; then
		echo "  $SAMPLE: assembly or reads not found, skipping" >&2
		continue
	fi

	echo "  $SAMPLE: mapping reads to the delivered assembly"
	minimap2 -ax map-hifi -t "$THREADS" "$ASSEMBLY" "$READS" 2>/dev/null \
		| samtools sort -@ 4 -o "$WORK/${SAMPLE}.bam" - 2>/dev/null
	samtools index "$WORK/${SAMPLE}.bam"
	samtools depth -a "$WORK/${SAMPLE}.bam" > "$WORK/${SAMPLE}.depth"

	# Locate the prophage within the assembly by aligning this sample's own
	# reconstructed phage genome back to it. Several samples carry more than one
	# phage, so every reference belonging to the sample is used and the intervals
	# are pooled.
	: > "$WORK/${SAMPLE}.phage_intervals"
	for PHAGE in results/phage_references/${SAMPLE}_*.fasta \
	             results/phage_reconstructions/${SAMPLE}_*.fasta; do
		[ -s "$PHAGE" ] || continue
		minimap2 -cx asm20 --secondary=no "$ASSEMBLY" "$PHAGE" 2>/dev/null \
			| awk '$11 >= 5000 { print $6 "\t" $8 "\t" $9 }' \
			>> "$WORK/${SAMPLE}.phage_intervals"
	done

	awk -v s="$SAMPLE" '{ sum[$1] += $3; n[$1]++ }
		END { for (c in sum) printf "%s\t%s\t%d\t%.1f\n", s, c, n[c], sum[c]/n[c] }' \
		"$WORK/${SAMPLE}.depth" | sort -k3,3nr >> "$BY_CONTIG"

	# Split the per-base depths into prophage and non-prophage. The chromosome
	# figure excludes the prophage, which is the whole point of the exercise.
	awk -v s="$SAMPLE" '
		FNR == NR { lo[FNR] = $2; hi[FNR] = $3; ctg[FNR] = $1; n = FNR; next }
		{
			isphage = 0
			for (i = 1; i <= n; i++)
				if ($1 == ctg[i] && $2 >= lo[i] && $2 <= hi[i]) { isphage = 1; break }
			if (isphage) { psum += $3; pn++ } else { csum += $3; cn++ }
		}
		END {
			cd = cn ? csum / cn : 0
			pd = pn ? psum / pn : 0
			printf "%s\t%.1f\t%.1f\t%s\t%s\t%d\n", s, cd, pd,
				(cd > 0 && pn) ? sprintf("%.1f", pd / cd) : "n/a",
				"see coverage_by_contig.tsv", pn
		}' "$WORK/${SAMPLE}.phage_intervals" "$WORK/${SAMPLE}.depth" >> "$SUMMARY"

	# The alignments are large and nothing downstream needs them.
	rm -f "$WORK/${SAMPLE}.bam" "$WORK/${SAMPLE}.bam.bai" "$WORK/${SAMPLE}.depth"
done

echo
echo "Per-sample summary written to $SUMMARY"
echo "Per-contig depths written to $BY_CONTIG"
echo "The chromosome figure excludes the prophage; quote the two separately."
