#!/bin/bash

################################################################################
# Replicate test: rebuild every phage genome from non-overlapping sets of reads
#
# BACKGROUND
#
# Each phage genome in this project was assembled from a pile of reads that had
# first been sorted by which phage they came from, and that pile was thinned to a
# workable depth before assembly. Two questions follow from that, and a reader is
# entitled to ask both. Does the result depend on which reads happened to be
# drawn? And would a different draw have produced a different genome?
#
# The obvious way to answer this is to repeat the assembly with a different
# random subsample. That is a weak test, because two random subsamples of the
# same pile share most of their reads and are not independent. This script does
# something stronger. It cuts each pile into subsets that do not overlap at all,
# so that no read appears in more than one replicate, and assembles each subset
# on its own. Agreement between replicates then cannot be an artefact of having
# looked at the same reads twice.
#
# WHAT COUNTS AS AGREEMENT
#
# Not byte identity. An assembler chooses an arbitrary starting point when it
# closes a circular genome, so two assemblies of the same molecule routinely
# differ by a rotation and are nonetheless the same sequence. The companion
# script compare_replicates.py handles this by doubling one sequence before
# comparing, and reports differences separately for homopolymer runs, where a
# small number of discrepancies is expected of this sequencing chemistry, and
# everywhere else, where they are not.
#
# WHERE A PILE IS TOO SHALLOW
#
# A pile that cannot supply two non-overlapping subsets at the usual depth is
# given a lower depth rather than being allowed to overlap. Overlapping
# replicates would inflate agreement, which is the one thing this test exists to
# avoid. Where that happens it is recorded in the log.
#
# USAGE
#
#   bash scripts/run_replicate_test.sh
#
# Edit data/config.sh first, and edit the LINEAGES table below so that it names
# the read piles and reference genomes you want to test. Run from the root of the
# repository. Expect roughly three minutes per replicate; the twenty-seven
# replicates reported in this project took a little under two hours.
#
# OUTPUT
#
#   <work>/replicates/asm/<lineage>_rep<n>/   one assembly per replicate
#   <work>/replicates/replicate_log.txt       depths, disjointness check, results
#
# Then run scripts/compare_replicates.py to produce results/replicate_test.tsv.
#
# DEPENDENCIES
#
#   seqtk      subsetting reads                      (tested with 1.5)
#   flye       de novo assembly                      (tested with 2.9.6)
#   awk, shuf  partitioning and reporting            (POSIX awk; shuf from coreutils)
################################################################################

set -o pipefail

source data/config.sh

# Each line names one phage lineage to test, and gives four things: a label, the
# list of read names making up its pile, the FASTQ those reads live in, and how
# many reads to put in each replicate. The read count sets the depth, so it must
# be chosen per lineage: roughly 670 reads of this length give the 150-fold depth
# at which the delivered genomes were assembled.
#
# The last field is how many non-overlapping replicates to cut. Multiplying it by
# the reads per replicate must not exceed the size of the pile, or the replicates
# would have to overlap.
#
#            label       read-name list                    reads FASTQ                 reads/rep  replicates
LINEAGES=(
	"S07_RHS:${WORK_DIR%/*}/s07_v2/pile_rhs.ids:${READS_DIR}/S07.fq:670:10"
	"S07_CdtB:${WORK_DIR%/*}/s07_v2/pile_cdtb.ids:${READS_DIR}/S07.fq:670:7"
	"S05_RHS:${WORK_DIR%/*}/s05_v2/r2_RHS5.ids:${READS_DIR}/S05.fq:670:6"
	"S05_hapA:${WORK_DIR%/*}/s05_v2/r2_hapA5.ids:${READS_DIR}/S05.fq:670:2"
	"S05_hapR:${WORK_DIR%/*}/s05_v2/r2_hapR5.ids:${READS_DIR}/S05.fq:365:2"
)

# The approximate length of an APSE genome, used only to report depth in a form
# that is easy to read. It does not affect the assembly.
PHAGE_LENGTH=39000

WORK="${WORK_DIR%/*}/replicates"
LOG="$WORK/replicate_log.txt"
mkdir -p "$WORK/asm"
cd "$WORK" || exit 1
: > "$LOG"

################################################################################
# Cut each pile into subsets that do not overlap.
#
# The read names are shuffled first, so that a replicate is not systematically
# made of reads from one part of the genome, and the shuffle is given a fixed
# seed so that the whole test can be repeated exactly.
################################################################################
for ENTRY in "${LINEAGES[@]}"; do
	IFS=: read -r LABEL IDS READS PER_REPLICATE N_REPLICATES <<< "$ENTRY"

	if [ ! -s "$IDS" ]; then
		echo "  $LABEL: read list $IDS not found, skipping" | tee -a "$LOG"
		continue
	fi

	echo "### $LABEL: $N_REPLICATES replicates of $PER_REPLICATE reads each" | tee -a "$LOG"
	sort -u "$IDS" | shuf --random-source=<(yes 42) > "${LABEL}_shuffled.ids"

	AVAILABLE=$(wc -l < "${LABEL}_shuffled.ids")
	NEEDED=$(( PER_REPLICATE * N_REPLICATES ))
	if [ "$NEEDED" -gt "$AVAILABLE" ]; then
		echo "  the pile holds $AVAILABLE reads but $NEEDED were requested; reduce the replicate count" | tee -a "$LOG"
		continue
	fi

	for i in $(seq 1 "$N_REPLICATES"); do
		FIRST=$(( (i - 1) * PER_REPLICATE + 1 ))
		LAST=$(( FIRST + PER_REPLICATE - 1 ))
		sed -n "${FIRST},${LAST}p" "${LABEL}_shuffled.ids" > "${LABEL}_rep${i}.ids"
		seqtk subseq "$READS" "${LABEL}_rep${i}.ids" > "${LABEL}_rep${i}.fq"
		awk -v n="$i" -v g="$PHAGE_LENGTH" \
			'NR%4==2 { reads++; bases += length($0) }
			 END { printf "  replicate %d: %d reads, approximately %.0f-fold depth\n", n, reads, bases/g }' \
			"${LABEL}_rep${i}.fq" | tee -a "$LOG"
	done
done

################################################################################
# Confirm that no read ended up in two replicates. This is the assumption the
# whole test rests on, so it is checked rather than assumed, and checked before
# any time is spent on assembly.
################################################################################
echo "### checking that the replicates do not overlap" | tee -a "$LOG"
for ENTRY in "${LINEAGES[@]}"; do
	IFS=: read -r LABEL _ _ _ _ <<< "$ENTRY"
	ls "${LABEL}"_rep*.ids > /dev/null 2>&1 || continue
	SHARED=$(cat "${LABEL}"_rep*.ids | sort | uniq -d | wc -l)
	printf "  %-12s %d reads appear in more than one replicate (this must be zero)\n" \
		"$LABEL" "$SHARED" | tee -a "$LOG"
done

################################################################################
# Assemble each replicate on its own.
################################################################################
echo "### assembling" | tee -a "$LOG"
for ENTRY in "${LINEAGES[@]}"; do
	IFS=: read -r LABEL _ _ _ N_REPLICATES <<< "$ENTRY"
	for i in $(seq 1 "$N_REPLICATES"); do
		[ -s "${LABEL}_rep${i}.fq" ] || continue
		OUTDIR="asm/${LABEL}_rep${i}"
		# Skip anything already assembled, so that an interrupted run can simply
		# be started again without repeating work.
		[ -s "$OUTDIR/assembly_info.txt" ] && continue
		flye --pacbio-hifi "${LABEL}_rep${i}.fq" --out-dir "$OUTDIR" \
			--threads "$THREADS" --meta > "${LABEL}_rep${i}.flye.log" 2>&1
		if [ -s "$OUTDIR/assembly_info.txt" ]; then
			awk -v t="${LABEL} replicate ${i}" \
				'NR>1 { printf "  %-24s %-10s %7d bp, circular = %s\n", t, $1, $2, $4 }' \
				"$OUTDIR/assembly_info.txt" | tee -a "$LOG"
		else
			echo "  ${LABEL} replicate ${i}: assembly failed, see ${LABEL}_rep${i}.flye.log" | tee -a "$LOG"
		fi
	done
done

echo
echo "Assemblies are in $WORK/asm and the log is $LOG."
echo "Run scripts/compare_replicates.py next to compare each replicate with the"
echo "delivered reference and write results/replicate_test.tsv."
