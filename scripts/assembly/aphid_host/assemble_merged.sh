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


#SBATCH --job-name=HDEF_assemble             	# Job name
#SBATCH --cpus-per-task=16                       # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=168:00:00                         # Maximum time the job will run
#SBATCH --qos=1week                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/meta_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/meta_err_%A_%a.log

# Define the number of threads to use
THREADS=16

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
conda activate eukaryotic_genome_assembly

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Ensure the output directory exists
mkdir -p assembly_merged

# Assemble genome
hifiasm --primary -t "$THREADS" -o assembly_merged/merged \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2177--bc2177.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2178--bc2178.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2179--bc2179.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2180--bc2180.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2181--bc2181.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2182--bc2182.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2183--bc2183.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2184--bc2184.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2185--bc2185.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2186--bc2186.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2187--bc2187.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2188--bc2188.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2189--bc2189.hifi_reads.fastq.gz \
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2190--bc2190.hifi_reads.fastq.gz

#Convert HIFIASM output to fasta format
scripts/gfa2fa assembly_merged/merged.p_ctg.gfa > assembly_merged/merged.p_ctg.fa

#Get some assembly stats
scripts/asmstats assembly_merged/merged.p_ctg.fa > assembly_merged/merged.p_ctg.fa.stats
