#!/bin/bash

# ---------------------------------------------------------------------------
# Reproduced as it was run, on a SLURM cluster. The absolute paths and the
# `module load` lines are the originals and will not resolve elsewhere: adapt
# the project root, the read paths and the module or conda names to your own
# system before running. The SLURM directives record the resources the step
# actually used and can be ignored if you run it serially.
#
# Note that `module load` pins no version. The versions used are listed in
# scripts/README.md and in the manuscript methods; they were recovered from run
# logs, not from these scripts.
# ---------------------------------------------------------------------------


#SBATCH --job-name=fun_train                  	# Job name
#SBATCH --cpus-per-task=64                     # Number of cores reserved
#SBATCH --mem-per-cpu=2G                       # Memory per core (total: 128GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train_err.log

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate funannotate

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Set path variables
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full"
FASTQ_DIR="RNAseq/fastq"

# Define the number of threads to use
THREADS=64

# Make sure the output directory doesn't already exist. Funnanotate train will throw an error if a previous run exists.
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

# Gather left and right reads as comma-separated lists
LEFT_READS=$(ls "$FASTQ_DIR"/*_1.renamed.fastq.gz | tr '\n' ' ')
RIGHT_READS=$(ls "$FASTQ_DIR"/*_2.renamed.fastq.gz | tr '\n' ' ')

# Join the arrays into space-separated strings
LEFT_READS_STR="${LEFT_READS[@]}"
RIGHT_READS_STR="${RIGHT_READS[@]}"

funannotate train -i "$ASSEMBLY" -o "$OUTDIR" \
	-l "$LEFT_READS_STR" \
	-r "$RIGHT_READS_STR" \
	--species "Aphis fabae" \
    --stranded no \
    --cpus "$THREADS"
