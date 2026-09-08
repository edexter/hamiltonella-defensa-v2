# Overview of annotation pipelineprocess

**Step 1: Remove duplicate contigs**

* Inputs: 
  * "aphid.filtered.fa": Aphid genome with non-target organisms filtered
* Outputs: 
  * "aphid.filtered.cleaned.fa":De-duplicated genome 
  * "cleaninglog.txt": Log file

**Step 2: Format FASTA headers**

* Inputs: 
  * "aphid.filtered.cleaned.fa": De-duplicated genome
* Outputs: 
  * "aphid.filtered.cleaned.renamed.fa": Analysis-ready genome

**Step 3A: Create repeat model**

* Inputs: 
  * "aphid.filtered.cleaned.renamed.fa": Analysis-ready genome
* Outputs: 
  * "consensi.fa.classified": Repeat library
  * "rmod.log": Log file

**Step 3B: Repeat mask genome**

* Inputs: 
  * "aphid.filtered.cleaned.renamed.fa": Analysis-ready genome
  * "consensi.fa.classified": Repeat library
* Outputs: 
  * "aphid.filtered.cleaned.renamed.fa.masked": Repeat masked genome
  * "aphid.filtered.cleaned.renamed.fa.tbl": Logfile

**Step 4: FUNANNOTATE train:**

* Inputs:
  * "aphid.filtered.cleaned.renamed.fa.masked": Repeat masked genome
  * "[...]_1.renamed.fastq.gz" and "[...]2.renamed.fastq.gz" Formatted RNAseq left and right reads.
* Outputs:
  * The initial FUNANNOTATE project directory

**Step 5: FUNANNOTATE predict**

* Inputs:
  * "aphid.filtered.cleaned.renamed.fa.masked": Repeat masked genome
  * The initial FUNANNOTATE project directory
* Outputs:
  * The updated FUNANNOTATE project directory

**Step 6: Interproscan**

* Inputs:
  * "Aphis_fabae.proteins.fa": Predicted proteins from FUNANNOTATE predict module
* Outputs:
  * "Aphis_fabae.proteins.fa.xml": INTERPROSCAN results file

**Step 6: FUNANNOTATE annotate**

* Inputs:

  * The FUNANNOTATE project directory
  * "Aphis_fabae.proteins.fa.xml": INTERPROSCAN results file

* Outputs:

  * The completed FUNANNOTATE project directory

    

# 1. Remove duplicate contigs

We want to identify and remove repetitive contigs that are contained in other scaffolds of the assembly, as those are likely assembly errors. If the repeats are indeed unique, then we want to keep them in the assembly. Funannotate has a module to clean up repetitive contigs  This is done using a “leave one out” approach using minimap2 or mummer (nucmer), where the the shortest contigs/scaffolds are aligned to the rest of the assembly to determine if it is repetitive. The script loops through the contigs starting with the shortest and workings its way to the N50 of the assembly, dropping contigs/scaffolds that are greater than the percent coverage of overlap (`--cov`) and the percent identity of overlap (`--pident`). This scripts sorts contigs by size, starting with shortest contigs it uses minimap2 to find contigs duplicated elsewhere, and then removes duplicated contigs. This step can be performed interactively because it only takes a few minutes.

````bash
./funannotate-docker clean \
--input inputRaw/aphid.filtered.fa \
--out inputProcessed/aphid.filtered.cleaned.fa
````

# 2. Sort/Rename FASTA Headers

Augustus has problems with long contig/scaffold names. Funannotate provides a script to sort by size and rename the headers, but we don't need to sort and rename them for this assembly, since we want to be able to link it back to the original meta-genome assembly. Instead we'll just remove the metadata from the fasta headers.

````bash
# Remove all metadata from contig headers
sed '/^>/ s/ .*//' inputProcessed/aphid.filtered.cleaned.fa > inputProcessed/aphid.filtered.cleaned.renamed.fa
````

# 3. Repeat Modeling/Masking

Annotation requires that repetitive regions are soft-masked (lower-case letters) in the genome assembly. Otherwise it will produce a lot of spurious gene predictions. The funannotate pipeline includes the option to pre-process the genome using repeatModeler and repeatMasker, but it doesn't use the new and improved versions of those tools, so we're going to call them separately. This step requires at least several days to one week of processing and so it submitted as a script on a cluster.

### RepeatMod.sh

Repeat modeling and repeat masking are closely related tasks, but we'll run them separately to make the pipeline modular. First we generate a de novo repeat library for the species. If we start with a good genome assembly, then we can reuse this repeat library for other genome assemblies for the same species, or other assemblies of this species. The may take up to one week to run, but only has to be performed once.

````bash
#!/bin/bash

#SBATCH --job-name=repeatModel                 # Job name
#SBATCH --cpus-per-task=16                     # Number of cores reserved
#SBATCH --mem-per-cpu=4G                       # Memory per core (total: 64GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMod_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMod_err.log

# Define the number of threads to use
THREADS=16

# Define input/output paths
INFILE="/scicore/home/ebertd/dexter0000/aphid/annotation/inputProcessed/aphid.filtered.cleaned.renamed.fa"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmodeler_output"
REPEAT_LIB="$OUTDIR/consensi.fa.classified"

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

# Step 1: Build the database
BuildDatabase -name "$OUTDIR/aphid_repeat_db" "$INFILE"

# Step 2: Run RepeatModeler with LTRStruct, directing output to a fixed location
RepeatModeler -LTRStruct -database "$OUTDIR/aphid_repeat_db" -threads "$THREADS" -dir "$OUTDIR"

# Verify RepeatModeler output
if [[ -f "$REPEAT_LIB" ]]; then
    echo "RepeatModeler completed successfully. Repeat library saved at $REPEAT_LIB"
else
    echo "Error: RepeatModeler did not generate consensi.fa.classified" >&2
    exit 1
fi
````

### RepeatMask.sh

We can use the repeat model that we developed in the previous step to soft-mask (lowercase letters) repeat elements. This is mandatory for gene prediction and annotation, and it take a very long time to run (up to a week).

````bash
#!/bin/bash

#SBATCH --job-name=repeatMask                  # Job name
#SBATCH --cpus-per-task=32                     # Number of cores reserved
#SBATCH --mem-per-cpu=2G                       # Memory per core (total: 64GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMask_quick_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/repeatMask_quick_err.log

# Define the number of threads to use
THREADS=32

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
RepeatMasker -qq -lib "$REPEAT_LIB" -pa "$THREADS" -gff -xsmall -dir "$OUTDIR" "$INFILE" -noisy
````

# 4. FUNANNOTATE train

````bash
#!/bin/bash

#SBATCH --job-name=fun_train                  	# Job name
#SBATCH --cpus-per-task=64                     # Number of cores reserved
#SBATCH --mem-per-cpu=2G                       # Memory per core (total: 128GB)
#SBATCH --time=168:00:00                       # Max runtime
#SBATCH --qos=1week                            # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_train_err.log

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate conda environment
conda activate funannotate

# Set working directory
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Set path variables
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full"
FASTQ_DIR="RNAseq/fastq"

# Define the number of threads to use
THREADS=64

# Create a custom temporary directory in your project folder
SCRATCH_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/tmp"
mkdir -p "$SCRATCH_DIR"
export TMPDIR="$SCRATCH_DIR"

# Make sure the output directory doesn't already exist. Funnanotate train will throw an error if a previous run exists.
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

# Gather left and right reads as individual variables
LEFT_READS=( "$FASTQ_DIR"/*_1.renamed.fastq.gz )
RIGHT_READS=( "$FASTQ_DIR"/*_2.renamed.fastq.gz )

funannotate train -i "$ASSEMBLY" -o "$OUTDIR" \
    --left "${LEFT_READS[@]}" \
    --right "${RIGHT_READS[@]}" \
    --species "Aphis fabae" \
    --stranded no \
    --cpus "$THREADS"
    
# Clean up the temporary directory after the job completes
rm -rf "$SCRATCH_DIR"
````



# 5. FUNNANOTATE predict

This runs in ~8 hours.

````bash
./funannotate-docker predict -i /home/ec2-user/aphid.filtered.cleaned.renamed.fa.masked \
    -o "funannotate_train_full" \
    -s "Aphis fabae" \
    --cpus 16 \
    --organism other \
    --repeats2evm \
    --ploidy 2 \
    --optimize_augustus \
    --max_intronlen 50000
````

# 6. INTERPROSCAN

This takes about 8-12 hours with 12 CPUs. Run it in a screen to prevent disconnect issues. DO NOT USE all available CPUs or it will crash the server. DO NOT launch using the INTERPROSCAN wrapper function in FUNANNOTATE, run INTERPROSCAN DIRECTLY INSTEAD. The input file "[genome].proteins.fa" is located in the "predict_results" folder that INTERPROSCAN predict previously generated.

````
# Start a new screen
screen -S interpro

# Launch interproscan (with output log from screen)
# Detach from the screen with Ctrl + A, then D
interproscan.sh -i funannotate_train_full/predict_results/Aphis_fabae.proteins.fa -f XML -goterms -pa --cpu 12 > interpro.log 2>&1

# Reattach to the screen (if needed)
screen -r interpro
````

# 7. FUNNOTATE annotate

Will need to pass the output files from IPRSCAN to this function.

````bash
# Start a new screen
screen -S annotate

./funannotate-docker annotate \
    --input "funannotate_train_full" \
    --cpus 14 \
    --iprscan Aphis_fabae.proteins.fa.xml \
    > annotate.log 2>&1

# Reattach to the screen (if needed)
screen -r annotate
````

