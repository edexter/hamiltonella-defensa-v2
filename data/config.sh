#!/bin/bash

################################################################################
# Configuration for the APSE / Hamiltonella multiplicity screen
#
# This is the ONLY file that needs to be edited to run the analysis on a
# different machine. Everything else in this repository is path independent.
#
# The sequencing reads and genome assemblies are NOT distributed with this
# repository because of their size (~9 GB of reads, ~8 GB of alignments). They
# are available from the authors and from the accessions listed in the
# publication. Set the two paths below to point at your local copy.
################################################################################

# Root directory holding the original project data. Expected to contain the
# sub-directories "reads_filtered" and "FLYE2" described below.
APHID_DATA="/mnt/d/aphid"

# Directory of Hamiltonella-filtered HiFi reads, one FASTQ per sample, named
# after the sample ID (for example "S07.fq"). These are the same reads that
# were supplied to the FLYE assembler in the original project.
READS_DIR="${APHID_DATA}/reads_filtered"

# Directory of raw (pre-curation) genome assemblies, one sub-directory per
# sample. Two file names are in use and the script accepts either:
# "SXX.raw.fasta" or "assembly.fasta".
#
# These are the assemblies that were delivered and to which the manual curation
# was applied, so they are the authoritative record and are what this screen
# examines. Note that the FLYE working directories on the analysis drive are NOT
# a safe substitute: for three samples they hold the output of a later or failed
# re-run rather than the delivered assembly. For S05 the working directory
# contains only 109 kb of a 2.27 Mb genome. The script verifies contig count and
# total length for every sample and records them in the log.
ASSEMBLY_DIR="/mnt/c/Users/ericd/Dropbox/Eric Work/DPA/APHIDS/data"

# Scratch directory for intermediate alignments. This must NOT be inside the
# repository: the script writes one sorted BAM per sample and these total
# several gigabytes. Only the small summary tables are kept in "results".
WORK_DIR="${HOME}/aphid_work/multiplicity_screen"

# Number of CPU threads to use for read mapping and sorting.
THREADS=8
