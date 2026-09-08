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


#SBATCH --job-name=SAMPLE_BACKMAP	            # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/backmap_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/backmap_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=8

# Load required modules
module load minimap2
module load SAMtools

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#File index
INDEXFILE=scripts/fastq_index.txt

#Name index
INDEXNAME=scripts/sampleID.txt

#File i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

OUTPREFIX=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

# Define the reference genome
REF="assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Map reads to the reference genome
minimap2 -ax map-pb -t "$THREADS" "$REF" "$SAMP" | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o bamsSelfmap/"$NAME".bam

samtools index -c bamsSelfmap/"$NAME".bam
