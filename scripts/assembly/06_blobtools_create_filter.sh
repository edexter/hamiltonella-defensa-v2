#!/bin/bash

# ---------------------------------------------------------------------------
# Unlike the other scripts in this directory, this one was not preserved as a
# job script: the step was carried out interactively. The commands below are
# the ones that were run, taken from the protocol document in original/, with
# the interactive contig selection replaced by the saved JSON filters so that
# it can be repeated without the viewer. Paths are the original cluster paths.
# ---------------------------------------------------------------------------


################################################################################
# Build a BlobToolKit dataset per sample and apply the saved contig filter.
#
# This is the step that separates Hamiltonella contigs from everything else in
# the metagenome assembly. It combines four kinds of evidence produced by the
# preceding scripts: read coverage (02), gene-space completeness (03),
# nucleotide-level taxonomy (04) and protein-level taxonomy (05).
#
# Contig selection was originally made by hand in the BlobToolKit viewer, and
# the criteria were saved as one JSON file per sample. Those files are
# distributed with this repository in data/blobtools_JSON/, so the selection is
# replayed here non-interactively and does not have to be repeated by eye.
#
# The recorded criteria are: GC between 0.35 and 0.50; superkingdom not
# Eukaryota; family not Erwiniaceae, Yersiniaceae or Morganellaceae; genus not
# Providencia; phylum Uroviricota retained, so that tailed bacteriophage
# sequence is kept deliberately rather than discarded as non-target. Minimum
# length and minimum coverage were set per sample and are in the JSON.
#
# To inspect or re-derive the filter interactively instead:
#     blobtools view --remote .
# then forward the ports and open http://localhost:8001/view/all
################################################################################

set -euo pipefail

conda activate btk

while IFS= read -r SAMP; do
	echo "Processing sample ID: $SAMP"

	blobtools create --replace \
		--fasta assemblies/"$SAMP"/"$SAMP".p_ctg.fa \
		--cov bamsSelfmap/"$SAMP".bam \
		--hits blobtools/"$SAMP"/"$SAMP".ncbi.blastn.out \
		--hits blobtools/"$SAMP"/"$SAMP".diamond.blastx.out \
		--taxdump blobtools/taxdump \
		--busco busco/"$SAMP"/run_enterobacterales_odb10/full_table.tsv \
		blobtools/"$SAMP"

	blobtools filter \
		--json blobtools/JSON/"$SAMP".current.json \
		--fasta assemblies/"$SAMP"/"$SAMP".p_ctg.fa \
		--output blobtools/filtered/"$SAMP" \
		blobtools/"$SAMP"

done < scripts/sampleID.txt
