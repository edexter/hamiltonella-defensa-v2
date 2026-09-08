#!/bin/bash

################################################################################
# Annotate the Hamiltonella defensa genomes with NCBI PGAP
#
# PGAP is containerised and does not run well on a shared cluster, so it is run
# on a dedicated cloud instance. See ../assembly/original/setup/PGAP_setup.md
# for how that instance is built and how the pipeline is installed.
#
# THREE THINGS THAT WILL STOP THE RUN IF NOT SET CORRECTLY
#
#   the organism   PGAP validates the organism name against its own taxonomy
#                  database and exits before doing any work if the name is not
#                  accepted. The accepted name is "Candidatus Williamhamiltonella
#                  defendens", taxid 138072. PGAP still writes the more familiar
#                  Hamiltonella name into the finished annotation.
#
#   the interpreter  pgap.py uses type syntax introduced in Python 3.10 and
#                  crashes on an older system Python. It is invoked through
#                  python3.11 here.
#
#   the exit code  PGAP can exit zero having produced nothing, so a run counts
#                  as successful below only if annot.gbk and annot.faa were
#                  actually written.
#
# CONCURRENCY
#
# Genomes are run several at a time and the queue refills as each finishes,
# rather than waiting for a whole batch: a fixed batch wastes the tail of every
# round on whichever genome happens to be slowest.
#
# The concurrency figure is a judgement rather than a calculation. A single PGAP
# job was measured driving a load average of 16 on 8 cores, so it oversubscribes
# the processors it is given and jobs contend more than dividing cores between
# them would suggest. Memory is the harder limit: one genome ran comfortably
# inside 30 GB but the peak was never measured, so the value below leaves real
# headroom rather than assuming the observed figure is the maximum.
#
# USAGE
#
#   bash run_pgap_batch.sh              # every sample listed below
#   bash run_pgap_batch.sh S05 S07      # named samples only
#
# Run it from the home directory of the instance, with the genomes already in
# ~/genomes. Launch it under nohup or screen; it runs for hours.
#
# OUTPUT
#
#   ~/<sample>_results/     one directory per genome, as PGAP writes it
#   ~/<sample>.pgap.log     one log per genome, so a single failure can be
#                           diagnosed without reading the others
#   ~/pgap_batch_status.txt a line per genome as it finishes
################################################################################

set -o pipefail

# The current accepted name. Verified against PGAP's own taxonomy database
# rather than assumed: the old name fails, this one is accepted.
ORGANISM="Candidatus Williamhamiltonella defendens"

# How many genomes to have in flight at once. See the note above on why this is
# well below the core count.
CONCURRENCY=7

GENOME_DIR="$HOME/genomes"
STATUS="$HOME/pgap_batch_status.txt"

SAMPLES=("$@")
if [ ${#SAMPLES[@]} -eq 0 ]; then
	SAMPLES=(S02 S03 S04 S05 S06 S07 S08 S09 S10 S11 S12 S13 S14)
fi

: > "$STATUS"
echo "starting $(date '+%H:%M:%S'): ${#SAMPLES[@]} genomes, $CONCURRENCY at a time" | tee -a "$STATUS"

annotate() {
	local SAMPLE="$1"
	local GENOME="$GENOME_DIR/${SAMPLE}_v2.fasta"
	local OUT="$HOME/${SAMPLE}_results"
	local LOG="$HOME/${SAMPLE}.pgap.log"

	if [ ! -s "$GENOME" ]; then
		echo "$SAMPLE SKIPPED - $GENOME not found" >> "$STATUS"
		return
	fi

	# Remove any partial output from an earlier attempt; PGAP refuses to write
	# into a directory that already exists.
	rm -rf "$OUT"

	local START=$(date +%s)
	python3.11 "$HOME/pgap.py" -r -o "$OUT" -g "$GENOME" -s "$ORGANISM" > "$LOG" 2>&1
	local MINUTES=$(( ($(date +%s) - START) / 60 ))

	# A run counts as successful only if the annotation itself was written.
	# An exit code of zero is not sufficient, because PGAP can finish tidily
	# having produced nothing.
	if [ -s "$OUT/annot.gbk" ] && [ -s "$OUT/annot.faa" ]; then
		local PROTEINS=$(grep -c "^>" "$OUT/annot.faa")
		printf "%s OK   %3d min  %5d proteins\n" "$SAMPLE" "$MINUTES" "$PROTEINS" >> "$STATUS"
	else
		printf "%s FAILED after %d min - see %s\n" "$SAMPLE" "$MINUTES" "$LOG" >> "$STATUS"
	fi
}

for SAMPLE in "${SAMPLES[@]}"; do
	# Wait for a slot to free up before starting the next genome.
	while [ "$(jobs -rp | wc -l)" -ge "$CONCURRENCY" ]; do
		sleep 30
	done
	echo "  $(date '+%H:%M:%S') starting $SAMPLE"
	annotate "$SAMPLE" &
done

wait

echo "finished $(date '+%H:%M:%S')" | tee -a "$STATUS"
echo
echo "Per-genome results:"
cat "$STATUS"
