#!/bin/bash

################################################################################
# BUSCO completeness for the version 2.0 assemblies
#
# WHY THE SETTINGS MATTER MORE THAN THE COMPUTE
#
# BUSCO on a two-megabase bacterial genome is cheap: prodigal calls the genes and
# hmmsearch scores them against a few hundred profiles, which takes minutes. The
# part that needs care is comparability. A BUSCO percentage means nothing on its
# own -- it is a fraction of whichever lineage dataset was used, and different
# datasets contain different numbers of markers. The version 1.0 figures were
# produced with BUSCO 5.1.2 against enterobacterales_odb10, which holds 440
# markers. Running version 2.0 against anything else would produce numbers that
# look comparable and are not.
#
# This script therefore pins the version, pins the dataset, and runs offline so
# that BUSCO cannot quietly fetch a newer release of the same lineage.
#
# WHAT TO EXPECT
#
# Most samples should be unchanged. Removing an aphid contig cannot remove a
# bacterial marker, and moving a prophage from its own contig into the chromosome
# does not alter gene content. The two samples that gained sequence are S05 and
# S07, and those are the only ones where a change is expected. A change anywhere
# else is a signal worth chasing rather than a number to report.
#
# USAGE
#
#   bash scripts/run_busco.sh                 # every sample
#   bash scripts/run_busco.sh S05 S07         # named samples only
#
# Edit data/config.sh first. Run from the root of the repository. Roughly three
# to five minutes per sample.
#
# OUTPUT
#
#   <work>/busco/<sample>/short_summary*.txt   BUSCO's own summary per sample
#
# Then run scripts/build_assembly_stats.py, which reads those summaries and fills
# the BUSCO columns of the statistics table.
#
# DEPENDENCIES
#
#   busco      completeness assessment              (5.1.2, to match version 1.0)
#   prodigal   gene calling, invoked by BUSCO       (2.6.3)
#
# BUSCO 5.1.2 must be installed on Python 3.10 or earlier. It opens files in the
# legacy "rU" mode, which was removed in Python 3.11, so on a newer interpreter
# it fails on every sample with an uninformative error about a missing argument.
# The environment used here was built with:
#
#   conda create -n busco512 -c conda-forge -c bioconda python=3.9 busco=5.1.2
#
# Upgrading BUSCO instead would avoid the problem but change the marker set, and
# the version 1.0 percentages would no longer be comparable.
################################################################################

set -o pipefail
source data/config.sh

# The lineage dataset used for version 1.0. It is read from the local copy rather
# than downloaded, so that the marker set cannot change underneath the comparison.
LINEAGE_NAME="enterobacterales_odb10"
LINEAGE_SOURCE="${APHID_DATA}/blobtools/busco/${LINEAGE_NAME}"

WORK="${WORK_DIR%/*}/busco"
mkdir -p "$WORK"

# The dataset lives on an external drive, and BUSCO reads it repeatedly. Copying
# it to local disk once is much faster than letting every sample stream it.
LINEAGE="$WORK/${LINEAGE_NAME}"
if [ ! -d "$LINEAGE" ]; then
	if [ ! -d "$LINEAGE_SOURCE" ]; then
		echo "  lineage dataset not found at $LINEAGE_SOURCE" >&2
		echo "  version 1.0 used ${LINEAGE_NAME}; the comparison needs the same one." >&2
		exit 1
	fi
	echo "  copying $LINEAGE_NAME to local disk"
	cp -r "$LINEAGE_SOURCE" "$LINEAGE"
fi

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
	SAMPLES=(S01 S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 S12 S13 S14)
fi

cd "$WORK" || exit 1

for SAMPLE in "${SAMPLES[@]}"; do
	ASSEMBLY="${OLDPWD}/assemblies_v2/${SAMPLE}_v2.fasta"
	[ -s "$ASSEMBLY" ] || { echo "  $SAMPLE: assembly not found, skipping" >&2; continue; }

	# Skip anything already done, so an interrupted run can simply be restarted.
	if ls "$SAMPLE"/short_summary*.txt >/dev/null 2>&1; then
		echo "  $SAMPLE: already done, skipping"
		continue
	fi

	echo "  $SAMPLE: running BUSCO"
	busco --in "$ASSEMBLY" --out "$SAMPLE" --mode genome \
		--lineage_dataset "$LINEAGE" --offline \
		--cpu "$THREADS" --force > "${SAMPLE}.busco.log" 2>&1

	SUMMARY=$(ls "$SAMPLE"/short_summary*.txt 2>/dev/null | head -1)
	if [ -n "$SUMMARY" ]; then
		grep -o 'C:[0-9.]*%\[S:[0-9.]*%,D:[0-9.]*%\],F:[0-9.]*%,M:[0-9.]*%,n:[0-9]*' "$SUMMARY" \
			| head -1 | awk -v s="$SAMPLE" '{printf "    %s  %s\n", s, $0}'
	else
		echo "    $SAMPLE: BUSCO produced no summary, see ${SAMPLE}.busco.log" >&2
	fi
done

echo
echo "Summaries are in $WORK/<sample>/."
echo "Run scripts/build_assembly_stats.py to fold them into the statistics table."
