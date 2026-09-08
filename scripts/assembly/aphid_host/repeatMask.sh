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


#SBATCH --job-name=repeatMask                  # Job name
#SBATCH --cpus-per-task=16                     # Number of cores reserved
#SBATCH --mem-per-cpu=4G                       # Memory per core (total: 64GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMask_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMask_err.log

# Define the number of threads to use
THREADS=16

# Define input/output paths
INFILE="/scicore/home/ebertd/dexter0000/aphid/annotation/inputProcessed/aphid.filtered.cleaned.renamed.fa"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output"
REPEAT_LIB="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmodeler_output/consensi.fa.classified"

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Load conda environment
conda activate funannotate

# Confirm environment activation
if [[ "$CONDA_DEFAULT_ENV" != "funannotate" ]]; then
    echo "Conda environment failed to activate" >&2
    exit 1
fi

# Disable usage reporting in BLAST
export BLAST_USAGE_REPORT=false

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Create the output directory if it does not exist
mkdir -p "$OUTDIR"

# Run repeatMasker
RepeatMasker -q -lib "$REPEAT_LIB" -pa "$THREADS" -gff -xsmall -dir "$OUTDIR" "$INFILE" 
