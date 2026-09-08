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


#SBATCH --job-name=fun_predict                    # Job name
#SBATCH --cpus-per-task=32                        # Number of cores
#SBATCH --mem-per-cpu=2G                          # Memory per CPU (total 64GB)
#SBATCH --time=168:00:00                          # Max runtime
#SBATCH --qos=1week                               # Queue (time-based)
#SBATCH --output=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred_out.log
#SBATCH --error=/scicore/home/ebertd/dexter0000/aphid/logs/fun_pred_err.log

# Load environment with Funannotate
eval "$(conda shell.bash hook)"
conda activate funannotate

# Set paths to your data files
ASSEMBLY="/scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked"
TRAIN_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_short/training"
OUTDIR="/scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_predict"
SCRATCH_DIR="/scicore/home/ebertd/dexter0000/aphid/annotation/scratch"

# Define the number of threads to use
THREADS=32

# Navigate to project folder
cd /scicore/home/ebertd/dexter0000/aphid/annotation

# Create scratch directory if it doesn't exist
if [ ! -d "$SCRATCH_DIR" ]; then
    mkdir -p "$SCRATCH_DIR"
fi

# Make sure the output directory doesn't already exist. Funannotate will throw an error if a previous run exists.
if [ -d "$OUTDIR" ]; then
    rm -r "$OUTDIR"
fi

# Run funannotate predict with updated parameters, using scratch directory as tmpdir
funannotate predict -i "$ASSEMBLY" \
    -o "$OUTDIR" \
    -s "Aphis fabae" \
    --cpus "$THREADS" \
    --augustus_species "insect" \
    --transcript_evidence "$TRAIN_DIR/trinity.fasta.clean"  \
    --genemark_mode ET \
    --other_gff "$TRAIN_DIR/pasa.step1.gff3:10" \
    --max_intronlen 3000 \
    --organism other \
    --repeats2evm \
    --busco_db insecta \
    --tmpdir "$SCRATCH_DIR"

# Optional: clean up scratch directory if desired after run completion
rm -r "$SCRATCH_DIR"
