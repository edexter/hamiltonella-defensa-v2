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


#SBATCH --job-name=BLAST			            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 256GB
#SBATCH --time=168:00:00                        # Maximum time the job will run
#SBATCH --qos=1week                             # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/blast_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/blast_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=16

# Load required modules
module load BLAST+

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name index
INDEXNAME=scripts/sampleID.txt

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/FLYE/$NAME/assembly.fasta"

# Navigate to the results directory
cd blobtools

# Set up results directory
mkdir -p "$NAME"

# Perform BLAST against the NT database
blastn -db nt/nt \
       -query "$REF" \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out "$NAME"/"$NAME".ncbi.blastn.out
