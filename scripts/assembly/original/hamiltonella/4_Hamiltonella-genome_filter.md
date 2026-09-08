# Hamiltonella-genome filtering with BLOBTOOLS 2

Each of the 14 samples represents a genome which should be primarily Hamiltonella, because we pre-filtered out any reads which previously mapped to non-Hamiltonella contigs. BLOBTOOLS 2 will now be used to identify and remove any remaining non-Hamiltonella contigs from each assembly based on coverage, BUSCO completeness, and BLAST against NCBI databases.

# Backmap reads

Backmap the reduced read set the FLYE assemblies. This runs in less than an hour.

````bash
#!/bin/bash

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
REF="FLYE2/$NAME/assembly.fasta"

# Ensure output directory exists
mkdir -p FLYE2/bamsSelfmap

# Map reads to the reference genome
minimap2 -ax map-pb -t "$THREADS" "$REF" "$SAMP" | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o FLYE2/bamsSelfmap/"$NAME".bam

# Index the reference genome (must be .csi format index)
samtools index -c FLYE2/bamsSelfmap/"$NAME".bam
````

# Find BUSCO

Calculate the BUSCO score for each of the FLYE assemblies. This runs in just a few minutes.

````
#!/bin/bash

#SBATCH --job-name=BUSCO			            # Job name
#SBATCH --cpus-per-task=8                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/busco_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/busco_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=8

# Load required modules
module load BUSCO

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name index
INDEXNAME=scripts/sampleID.txt

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/FLYE2/$NAME/assembly.fasta"

# Move into output directory
mkdir -p FLYE2/busco
cd FLYE2/busco

# Calculate BUSCO scores
busco -i "$REF" \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/enterobacterales_odb10 \
-c 8 \
-o $NAME \
-f --offline
````

# BLAST

Perform a BLAST search against the NCBI NT database. This runs in a few hours.

````
#!/bin/bash

#SBATCH --job-name=FL_BLAST			            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 256GB
#SBATCH --time=24:00:00                        # Maximum time the job will run
#SBATCH --qos=1day                             # The job queue (time based)
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
REF="/scicore/home/ebertd/dexter0000/aphid/FLYE2/$NAME/assembly.fasta"

# Navigate to the results directory
cd blobtools

# Set up results directory
mkdir -p FLYE2_"$NAME"

# Perform BLAST against the NT database
blastn -db nt/nt \
       -query "$REF" \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out FLYE2_"$NAME"/"$NAME".ncbi.blastn.out
````

# Run blobtools

Aggregate all of the necessary files to create a new BLOBTOOLS object for each sample. Interactively inspect for any potential assembly issues and filter the genomes as before.

````
# Activate environment
conda activate btk

# Navigate to project folder
cd aphid

mkdir -p FLYE2/blobtools

while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools create with the current sample ID
  blobtools create --replace \
    --fasta FLYE2/$SAMP/assembly.fasta \
    --cov FLYE2/bamsSelfmap/"$SAMP".bam \
    --hits blobtools/FLYE2_"$SAMP"/"$SAMP".ncbi.blastn.out \
    --taxdump blobtools/taxdump \
    --busco FLYE2/busco/"$SAMP"/run_enterobacterales_odb10/full_table.tsv \
    FLYE2/blobtools/"$SAMP"

done < scripts/sampleID.txt
````