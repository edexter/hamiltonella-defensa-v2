# BLOBTOOLS analysis

The merged sample meta-assembly will be used to create filtered genome assemblies for the Bucherna symbiont and the aphid host with BLOBTOOLS 2. This requires a number of preliminary steps to populate the BLOBTOOLS directory will all of the data that is needed for the filters to operate on



# Step 1: Backmap reads to the meta-genome assembly

Blobtools heavily relies upon a read coverage filter, so we align the raw reads from each sample against the metagenome assembly. This runs in less than 10 minutes.

* Inputs
  * Raw HIFI sequence read fastq files (gz compressed is fine)
  * The meta-genome assembly: "merged.p_ctg.fa"
* Outputs
  * Bam file and index with the base name "backmap.bam"

````bash
#!/bin/bash

#SBATCH --job-name=SAMPLE_BACKMAP	            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 16GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/backmap2_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/backmap2_err_%A_%a.log

# Define the number of threads to use
THREADS=16

# Load required modules
module load minimap2
module load SAMtools

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

# Define the reference genome
REF=assembly_merged/merged.p_ctg.fa

# Map reads to the reference genome
minimap2 -ax map-hifi -t "$THREADS" "$REF" data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2177--bc2177.hifi_reads.fastq.gz \
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
data_round2/fastx/m64141e_240627_145720.hifi_reads.bc2190--bc2190.hifi_reads.fastq.gz \
| samtools view -bS -@ "$THREADS" - | samtools sort -@ "$THREADS" -o assembly_merged/backmap.bam

# Index the reference genome (must be .csi format index)
samtools index -c assembly_merged/backmap.bam
````



# Step 2: Calculate BUSCO scores

Blobtools also needs BUSCO scores, but this needs to be performed for one species at a time, since the reference databases are different. This is fast and can be performed interactively.

* Inputs
  * BUSCO databases for hemiptera (aphid) and enterobacterales (Buchnera)
  * The meta-genome assembly: "merged.p_ctg.fa"
* Outputs
  * Busco results directories

````bash
# Load required module
module load BUSCO

# Navigate to directory
cd assembly_merged

# Calculate BUSCO scores for the Aphid
busco -i merged.p_ctg.fa \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/hemiptera_odb10/ \
-c 16 \
-o merged.busco.out \
-f --offline

# Calculate BUSCO scores for Buchnera
busco -i merged.p_ctg.fa \
-m genome \
-l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/enterobacterales_odb10 \
-c 16 \
-o merged.busco.buchnera.out \
-f --offline
````



# Step 3: Perform BLAST

Blobtools blasts all of the contigs in the assembly against the NCBI database to predict their taxonomic rank. This may run for up to a week. This needs to be performed against a local copy of the database, otherwise it will take even longer.

* Inputs
  * The meta-genome assembly: "merged.p_ctg.fa"
* Outputs
  * BLAST results: "merged.ncbi.blastn.out"

````bash
#!/bin/bash

#SBATCH --job-name=BLAST			            # Job name
#SBATCH --cpus-per-task=16                      # Number of cores reserved
#SBATCH --mem-per-cpu=8G                        # Memory reserved per core
                                                # Total memory reserved: 256GB
#SBATCH --time=168:00:00                        # Maximum time the job will run
#SBATCH --qos=1week                             # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/blast_out_merge.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/blast_err_merge.log

# Define the number of threads to use
THREADS=16

# Load required modules
module load BLAST+

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid

#Name
NAME=merged

# Define the reference genome
REF="/scicore/home/ebertd/dexter0000/aphid/assembly_merged/merged.p_ctg.fa"

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



# Step 4: Populate the BLOBTOOLS directory

Now we have to create a BLOBTOOLS directory that brings all of the inputs together. This only takes a few minutes to to run.

````bash
# Activate environment
conda activate /scicore/home/ebertd/dexter0000/anaconda3/envs/btk

# Navigate to project folder
cd aphid

# Run blobtools create with the current sample ID
blobtools create --replace \
	--fasta assembly_merged/merged.p_ctg.fa \
    --cov assembly_merged/backmap.bam \
    --hits blobtools/merged/merged.ncbi.blastn.out \
    --taxdump blobtools/taxdump \
    --busco assembly_merged/merged.busco.out/run/full_table.tsv \
    viewer/merged
````
