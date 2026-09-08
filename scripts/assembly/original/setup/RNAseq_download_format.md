# Fetch RNA data

My cluster environment won't let me download from the SRA over a compute node, so I just have to log in to a dedicated transfer node and run it from the head with the "nohup" command so that it will continue the download even if my local machine disconnects from the cluster.

### Set up conda environment

We can't use the funannotate conda environment to download the SRA files because it has conflicting dependencies with the sra-tools software. So we'll make a dedicated environment for downloading and preparing RNAseq data. 

````
# Create the new conda environment
conda create -n sra-tools -c bioconda sra-tools=3.0.0
````

### Download and format the SRA data

From this environment we can launch this script to download and pre-process the required files. **SRAdownload.sh**

````bash
#!/bin/bash

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Load conda environment
conda activate sra-tools

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Directory to save downloaded SRA files and output FASTQ files
SRA_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/RNAseq/sra_files"
FASTQ_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/RNAseq/fastq_files"

# Create directories if they do not exist
mkdir -p "$SRA_DIR"
mkdir -p "$FASTQ_DIR"

# List of SRA Run IDs associated with each BioSample ID
SRA_IDS=(
    SAMN10606880 SAMN10606881 SAMN10606882 SAMN10606883
    SAMN10606884 SAMN10606885 SAMN10606886 SAMN10606887
    SAMN10606888 SAMN10606889 SAMN10606890 SAMN10606891
    SAMN10606892 SAMN10606893 SAMN10606894 SAMN10606895
    SAMN10606896 SAMN10606897 SAMN10606898 SAMN10606899
    SAMN10606900 SAMN10606901 SAMN10606902 SAMN10606903
    SAMN10606904 SAMN10606905 SAMN10606906 SAMN10606907
    SAMN10606908 SAMN10606909 SAMN10606910 SAMN10606911
    SAMN10606912 SAMN10606913 SAMN10606914 SAMN10606915
    SAMN10606916 SAMN10606917 SAMN10606918 SAMN10606919
)

# Download each SRA file using prefetch with output file path specified
for SRA_ID in "${SRA_IDS[@]}"; do
    echo "Downloading SRA file for SRA ID: $SRA_ID"
    prefetch "$SRA_ID" --output-file "$SRA_DIR/${SRA_ID}.sra"
done

echo "All SRA files have been downloaded and converted to FASTQ format."

#Make an index of the SRA files that were downloaded
cd /scicore/home/ebertd/dexter0000/aphid/annotation/RNAseq/sra_files
ls *.sra > index.txt
````

My cluster won't allow me to launch this script from a compute node (file downloads are blocked), so I launch it from a dedicated file transfer node that uses a separate login.

````
# Navigate to folder where I stored the SRAdownload.sh script
cd /scicore/home/ebertd/dexter0000/aphid/scripts

# Run the scrip with nohup to avoid interruptions due to disconnects
nohup bash SRAdownload.sh > download_log.out 2>&1 &
````



# Formatting

The headers in the fastq files need to be re-formatted for funannotate with the script **renameRNA.sh**

````bash
#!/bin/bash

#SBATCH --job-name=formatRNA             		# Job name
#SBATCH --cpus-per-task=4                       # Number of cores reserved
#SBATCH --mem-per-cpu=4G                        # Memory reserved per core
                                                # Total memory reserved: 32GB
#SBATCH --time=24:00:00                         # Maximum time the job will run
#SBATCH --qos=1day                              # The job queue (time based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/formatRNA_out_%A_%a.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/formatRNA_err_%A_%a.log
#SBATCH --array=1-40%40                         # Array job specifications

# Define the number of threads to use
THREADS=4

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Load conda environment
conda activate sra-tools

# Navigate to the project folder
cd /scicore/home/ebertd/dexter0000/aphid/annotation/RNAseq

# SRA File index
INDEXFILE=/scicore/home/ebertd/dexter0000/aphid/annotation/RNAseq/sra_files/index.txt

# Get the SRA file name for this job array index
SAMP=$(sed -n ${SLURM_ARRAY_TASK_ID}p $INDEXFILE)

# Extract the base name for each sample from $SAMP
NAME=$(basename "$SAMP" .sra)

# Set output directory
FASTQ_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/RNAseq/fastq"

# Ensure the output directory exists
mkdir -p "$FASTQ_DIR"

# Convert SRA to fastq.gz format with unique headers
echo "Converting $NAME to fastq.gz format"
fastq-dump --split-files --defline-seq '@$sn[_$rn]/$ri' --gzip --outdir "$FASTQ_DIR" sra_files/"$SAMP"
echo "$NAME converted to fastq.gz format"

# Fix headers for Funannotate for left (R1) reads
R1_IN="${FASTQ_DIR}/${NAME}_1.fastq.gz"
R1_OUT="${FASTQ_DIR}/${NAME}_1.renamed.fastq.gz"

# Add sample prefix and replace quality score header with "+" only
zcat "$R1_IN" | \
    sed -e "s/^@/@${NAME}_/" -e "s/^+.*$/+/" | \
    gzip > "$R1_OUT"

# Fix headers for Funannotate for right (R2) reads
R2_IN="${FASTQ_DIR}/${NAME}_2.fastq.gz"
R2_OUT="${FASTQ_DIR}/${NAME}_2.renamed.fastq.gz"

# Add sample prefix and replace quality score header with "+" only
zcat "$R2_IN" | \
    sed -e "s/^@/@${NAME}_/" -e "s/^+.*$/+/" | \
    gzip > "$R2_OUT"

echo "$NAME headers have been renamed and saved to $FASTQ_DIR."
````



