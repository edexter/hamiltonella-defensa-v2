# Meta-genome assembly

This script assembles a meta-genome for each of the 14 samples and calculate a few summary statistics. There are multiple species represented in each sample, so the meta-genome assembly is just the first step in the pipeline to assemble the individual genomes. This script runs in less than a day.

* Inputs
  * "fastq_index.txt": An index file of the raw sequencer file names
  * "sampleID.txt": An index file of the sample IDs

* Outputs (the most important ones)
  * "S01.p_ctg.fa": The primary genome assembly (for sample #1 in this case)
  * "S01.p_ctg.fa.stats": Some basic summary stats for the genome assembly (for sample #1 in this case)

````bash
#!/bin/bash

#SBATCH --job-name=HDEF_assemble             	# Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/HDEF_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/HDEF_err_%A_%a.log
#SBATCH --array=1-14%14                         # Array job specifications

# Define the number of threads to use
THREADS=8

# Initialize conda for bash
eval "$(conda shell.bash hook)"

#Load conda environment
conda activate eukaryotic_genome_assembly

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# File index
INDEXFILE=scripts/fastq_index.txt

# Name index
INDEXNAME=scripts/sampleID.txt

# File i from index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

# Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Ensure the output directory exists
mkdir -p assemblies/"$NAME"

# Assemble genome
hifiasm --primary -t "$THREADS" -o assemblies/"$NAME"/"$NAME" $SAMP

#Convert HIFIASM output to fasta format
scripts/gfa2fa assemblies/"$NAME"/"$NAME".p_ctg.gfa > assemblies/"$NAME"/"$NAME".p_ctg.fa

#Get some assembly stats
scripts/asmstats assemblies/"$NAME"/"$NAME".p_ctg.fa > assemblies/"$NAME"/"$NAME".p_ctg.fa.stats
````

