#!/bin/bash

################################################################################
# Audit of the taxonomic read filter applied during the original 2024 assembly
#
# BACKGROUND
#
# The reads used throughout this project are not the raw sequencing reads. They
# are a filtered subset, and the filter was applied in 2024 in the following way.
# Every sample was first assembled as a whole metagenome, containing the aphid,
# the primary symbiont Buchnera, Hamiltonella defensa, and anything else present
# in the insect. The contigs of that assembly were then classified with blobtools
# on the basis of taxonomy, GC content and read coverage, and the contigs that
# did not look like Hamiltonella were discarded. Finally, a read was kept only if
# it aligned to one of the contigs that survived.
#
# Every result in this project inherits that decision, and until now it had never
# been tested. The concern is specific rather than general. Hamiltonella defensa
# is itself an enterobacterium, and the filter excluded contigs whose closest
# taxonomic match fell in several related enterobacterial families. The classes
# of sequence most likely to be excluded on that basis are exactly the ones a
# degraded symbiont genome contains in quantity: horizontally transferred
# regions, insertion sequences and other mobile elements, and short contigs with
# weak or ambiguous database matches.
#
# WHAT THIS SCRIPT DOES
#
# For every sample it compares the assembly as it was before filtering with the
# assembly as it was after, identifies every contig that the filter discarded,
# and asks of each one whether it looks like Hamiltonella or like APSE. A
# discarded contig that aligns convincingly to either is a candidate for having
# been removed in error.
#
# INTERPRETING THE OUTPUT
#
# Two numbers decide whether a flagged contig represents a real loss.
#
# The first is identity. A published Hamiltonella defensa genome is used as the
# yardstick, but it is a different strain from any in this panel, so a perfect
# match should not be expected. Calibrate against the panel itself: the retained
# chromosomes of these samples align to that reference at between 91 and 95 per
# cent. Anything appreciably below that range is not Hamiltonella, however
# alarming a raw alignment might look.
#
# The second is whether the sequence survives anywhere else. A contig can be
# discarded quite correctly if the same sequence is present in a contig that was
# kept, which happens whenever the assembler produced two representations of one
# region. The script therefore aligns each flagged contig back against the
# retained assembly and reports how much of it is already there.
#
# A contig is a genuine loss only if it both looks like Hamiltonella or APSE by
# identity and is absent from the retained assembly.
#
# USAGE
#
#   bash scripts/audit_blobtools_filter.sh                # every sample
#   bash scripts/audit_blobtools_filter.sh S03 S08        # named samples only
#
# Edit data/config.sh first so that the paths point at your copy of the data.
# Run the script from the root of the repository. It takes roughly two minutes
# per sample on a laptop and needs no cluster.
#
# OUTPUT
#
#   results/blobtools_filter_audit.tsv   one row per flagged contig
#
# DEPENDENCIES
#
#   minimap2   sequence alignment                    (tested with 2.31-r1302)
#   samtools   extracting individual contigs         (tested with 1.24)
#   awk        parsing and summarising               (POSIX awk, mawk or gawk)
################################################################################

set -o pipefail

source data/config.sh

# The assembly of each sample as it was BEFORE the blobtools filter was applied.
# This is the whole-metagenome assembly, so it is large: several thousand contigs
# and of the order of a hundred megabases per sample.
PRE_FILTER_DIR="${APHID_DATA}/assemblies"

# The same assemblies AFTER filtering. These hold only the contigs that survived.
POST_FILTER_DIR="${APHID_DATA}/assemblies_filtered"

# A published Hamiltonella defensa genome, used as an external yardstick. Any
# complete H. defensa assembly will serve; the accession is recorded so that the
# comparison can be repeated exactly.
HDEF_REFERENCE="${APHID_DATA}/reference_genomes/GCA_000021705.1_ASM2170v1_genomic.fna"

# The APSE genomes reconstructed by this project, so that a discarded phage
# contig is recognised as a phage rather than as an unknown.
APSE_REFERENCES="results/phage_references results/phage_reconstructions"

# A contig must align over at least this many bases before it is considered at
# all. Shorter matches are common between any two bacterial genomes and say
# nothing about whether a contig belongs to this species.
MIN_ALIGNED_BP=1000

# ... and at least this identity. The value sits deliberately below the 91 to 95
# per cent that the retained chromosomes of this panel show against the published
# reference, so that the audit errs towards flagging too much rather than too
# little. Anything flagged at close to this floor is then judged by hand.
MIN_IDENTITY=0.88

OUT="results/blobtools_filter_audit.tsv"
WORK="${WORK_DIR%/*}/blobtools_audit"
mkdir -p "$WORK" results

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
	SAMPLES=(S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 S12 S13 S14)
fi

################################################################################
# Gather the APSE references into one file, so that a discarded contig can be
# compared against every phage this project has reconstructed in a single pass.
################################################################################
cat $(for d in $APSE_REFERENCES; do ls "$d"/*.fasta 2>/dev/null; done) > "$WORK/apse.fa"

printf "sample\tcontig\tcontig_bp\taligned_bp\tidentity\tpct_present_in_retained\n" > "$OUT"

for SAMPLE in "${SAMPLES[@]}"; do
	PRE="$PRE_FILTER_DIR/$SAMPLE/$SAMPLE.p_ctg.fa"
	POST="$POST_FILTER_DIR/$SAMPLE.p_ctg.filtered.fa"

	if [ ! -s "$PRE" ] || [ ! -s "$POST" ]; then
		echo "  $SAMPLE: assemblies not found, skipping" >&2
		continue
	fi

	# The names of the contigs that survived the filter. Everything in the
	# pre-filter assembly whose name is absent from this list was discarded.
	grep "^>" "$POST" | awk '{print substr($1,2)}' | sort -u > "$WORK/kept_$SAMPLE.txt"

	# Align every contig of the pre-filter assembly against the published
	# Hamiltonella genome and against the APSE references, and keep those with a
	# match long enough and close enough to be worth considering. Alignments are
	# summed per contig because a genuine match is often reported in several
	# blocks rather than one.
	{
		minimap2 -cx asm20 -t "$THREADS" --secondary=no "$HDEF_REFERENCE" "$PRE" 2>/dev/null
		minimap2 -cx asm20 -t "$THREADS" --secondary=no "$WORK/apse.fa"    "$PRE" 2>/dev/null
	} | awk -v minbp="$MIN_ALIGNED_BP" -v minid="$MIN_IDENTITY" '
		$11 >= minbp && $10/$11 >= minid {
			aligned[$1] += $11
			if ($10/$11 > best[$1]) best[$1] = $10/$11
			len[$1] = $2
		}
		END { for (c in aligned) print c "\t" len[c] "\t" aligned[c] "\t" best[c] }
	' > "$WORK/hits_$SAMPLE.tsv"

	# Keep only the hits belonging to contigs the filter discarded.
	awk 'NR==FNR { kept[$1]; next } !($1 in kept)' \
		"$WORK/kept_$SAMPLE.txt" "$WORK/hits_$SAMPLE.tsv" > "$WORK/flagged_$SAMPLE.tsv"

	# For each flagged contig, ask the question that actually decides the matter:
	# is this sequence already present in the assembly that was kept? If it is,
	# the contig was a redundant representation and discarding it lost nothing.
	while IFS=$'\t' read -r CONTIG LENGTH ALIGNED IDENTITY; do
		samtools faidx "$PRE" "$CONTIG" > "$WORK/one.fa" 2>/dev/null || continue
		PRESENT=$(minimap2 -cx asm20 --secondary=no "$POST" "$WORK/one.fa" 2>/dev/null \
			| awk '{ total += $4 - $3 } END { print total + 0 }')
		awk -v s="$SAMPLE" -v c="$CONTIG" -v l="$LENGTH" -v a="$ALIGNED" \
		    -v i="$IDENTITY" -v p="$PRESENT" \
			'BEGIN { printf "%s\t%s\t%d\t%d\t%.3f\t%.0f\n", s, c, l, a, i, (l ? 100*p/l : 0) }' >> "$OUT"
	done < "$WORK/flagged_$SAMPLE.tsv"

	TOTAL=$(grep -c "^>" "$PRE")
	KEPT=$(wc -l < "$WORK/kept_$SAMPLE.txt")
	FLAGGED=$(wc -l < "$WORK/flagged_$SAMPLE.tsv")
	printf "  %-4s %6d contigs assembled, %3d kept, %3d discarded contigs flagged for inspection\n" \
		"$SAMPLE" "$TOTAL" "$KEPT" "$FLAGGED"
done

echo
echo "Results written to $OUT."
echo "A flagged contig represents a real loss only if its identity is within the"
echo "range the retained chromosomes show against the published reference, AND"
echo "the last column shows that little of it survives in the retained assembly."
