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


#SBATCH --job-name=FLYE2			             # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/FLYE_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/FLYE_err_%A_%a.log
#SBATCH --array=1-14%14                         # Array job specifications

# Define the number of threads to use
THREADS=8

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
module load Flye

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Name index
INDEXNAME=scripts/sampleID.txt

# Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Ensure the output directory exists
mkdir -p FLYE2/"$NAME"

# Assemble genome
flye --pacbio-hifi reads_filtered/"$NAME".fq --out-dir FLYE2/"$NAME" --threads "$THREADS"
