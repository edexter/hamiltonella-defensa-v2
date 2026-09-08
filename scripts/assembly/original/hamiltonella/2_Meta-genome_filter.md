# Meta-genome filtering with BLOBTOOLS 2

Each of the 14 samples represents a meta-genome of multiple species, but our goal is to generate cleaned Hamiltonella assemblies. BLOBTOOLS 2 will now be used to identify and remove non-Hamiltonella contigs from each assembly based on coverage, BUSCO completeness, and BLAST against NCBI databases.

# Step 1. Backmap reads to the assembled genomes

BLOBTOOLS needs read coverage for each genome assembly, so we align the raw reads from each sample against each metagenome assembly. This runs in less than 10 minutes.

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
REF="assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Map reads to the reference genome
minimap2 -ax map-pb -t "$THREADS" "$REF" "$SAMP" | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o bamsSelfmap/"$NAME".bam

# Index the reference genome (must be .csi format index)
samtools index -c bamsSelfmap/"$NAME".bam
````



# Step 2: Calculate BUSCO scores for each assembly

BLOBTOOLS will also need BUSCO results for each of our assemblies to help us decide whether we are removing too many contigs. Note that the BUSCO software is really inconsistent about relative file paths, hence the hard coding of file paths in the following script. This runs very quickly.

````bash
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
REF="/scicore/home/ebertd/dexter0000/aphid/assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Move into output directory
cd busco

# Calculate BUSCO scores
busco -i "$REF" \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/enterobacterales_odb10 \
-c 8 \
-o $NAME \
-f --offline
````



# Step 3: BLAST assemblies against NCBI NT database

This requires about one day to run and should be submitted as a batch script.

````
#!/bin/bash

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
REF="/scicore/home/ebertd/dexter0000/aphid/assemblies/"$NAME"/"$NAME".p_ctg.fa"

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

````



# Step 4. BLAST assemblies against Uniprot database

This is requires about one day to run and should be submitted as a batch script.

````
#!/bin/bash

#SBATCH --job-name=DIAMOND			            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 128GB
#SBATCH --time=168:00:00                        # Maximum time the job will run
#SBATCH --qos=1week                             # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/diamond_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/diamond_err_%A_%a.log
#SBATCH --array=1-14%14                        # Array job specifications

# Define the number of threads to use
THREADS=16

# Load required modules
module load DIAMOND 

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name index
INDEXNAME=scripts/sampleID.txt

#Name i from index
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXNAME)

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/assemblies/"$NAME"/"$NAME".p_ctg.fa"

# Navigate to the results directory
cd blobtools

# Set up results directory
mkdir -p "$NAME"

# Perform diamond blast
diamond blastx \
        --query "$REF" \
        --db /scicore/home/ebertd/dexter0000/aphid/blobtools/uniprot/reference_proteomes \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 16 \
        > "$NAME"/"$NAME".diamond.blastx.out
````



# Step 5. Create and populate the blob directory

The minimum requirement to create a new BlobDir is an assembly FASTA file. This command creates a new directory in the location specified by the last argument (in this case “AssemblyName”) that contains a set of files containing values for GC-content (gc.json), length (length.json), number of Ns (ncount.json) and sequence names (identifiers.json) for each sequence in the assembly. A final file (meta.json) contains metadata for the dataset describing the datatypes of the available fields and the ranges of values for each of these fields. We can also add other types of data at the same time. This only requires about 1 minute per sample and can be performed interactively.

````
# Activate environment
conda activate btk

# Navigate to project folder
cd aphid

while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools create with the current sample ID
  blobtools create --replace \
    --fasta assemblies/"$SAMP"/"$SAMP".p_ctg.fa \
    --cov bamsSelfmap/"$SAMP".bam \
    --hits blobtools/"$SAMP"/"$SAMP".ncbi.blastn.out \
    --hits blobtools/"$SAMP"/"$SAMP".diamond.blastx.out \
    --taxdump blobtools/taxdump \
    --busco busco/"$SAMP"/run_enterobacterales_odb10/full_table.tsv \
    blobtools/"$SAMP"

done < scripts/sampleID.txt
````



# Step 6. Interactively filter the assembly

Navigate to the folder where the blob directories are stored and launch the viewer

````bash
blobtools view --remote .
````

Open viewer in a new terminal on your local machine

````
ssh -L 8001:127.0.0.1:8001 -L 8000:127.0.0.1:8000 dexter0000@login.scicore.unibas.ch
````

Navigate to page in a web browser and use the interactive filtering tool. Filter as needed and save the filter parameters using the lists Menu.

````
http://localhost:8001/view/all
````



# Step 7. Create filtered assembly

Upload the JSON files that were created from interactive filtering in BLOBTOOLS and pass them to the filter module.

````
# Request resources for interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Activate environment
conda activate btk

# Navigate to project folder
cd aphid

while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools filter with the current sample ID
  blobtools filter \
     --json blobtools/JSON/"$SAMP".current.json \
     --fasta assemblies/"$SAMP"/"$SAMP".p_ctg.fa \
     --output blobtools/filtered/"$SAMP" \
     blobtools/"$SAMP"

done < scripts/sampleID.txt
````



































