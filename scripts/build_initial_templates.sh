#!/bin/bash

################################################################################
# Build the starting templates for read partitioning.
#
# The templates are the sample's own pre-curation contigs, grouped by integrase
# type and toxin family, one FASTA record per phage. They are used only to sort
# reads: no base of a template reaches the assembled output, so a template need
# not be free of misassembly. What it must do is span the diversity present in
# the sample. A lineage with no matching template does not merely get polished
# badly, it has its reads misassigned to another lineage, which is the one
# failure this step cannot recover from.
#
# Grouping a contig is a two-part question with a cheap answer in both parts:
# which toxin does it carry (align data/reference/toxins_CDS.fasta against it),
# and which integrase (align data/reference/integrases.fasta against it, and see
# scripts/detect_multiple_APSE_and_strains.sh, which reports both per contig).
#
# The grouping used for S07:
#
#     RHS phage    contig_3 contig_4
#     CdtB phage   contig_5 contig_1
#
# contig_6 is a true duplicate of part of contig_5 and is left out; including it
# would have been harmless.
#
# S05 carries three lineages, one of which is a recombinant of the other two.
# A recombinant has almost no unique k-mers, so it cannot be separated by the
# same route; reads are assigned there by containment instead, using
# scripts/lib/kmer_assign.py. See the header of
# reconstruct_phages_and_integration_sites.sh.
#
# Usage:
#   bash scripts/build_initial_templates.sh S07 out.fa RHS:contig_3,contig_4 CdtB:contig_5,contig_1
#
# Arguments: sample ID, output FASTA, then one NAME:contig,contig,... group per
# phage. Contigs are taken from the pre-curation assembly named in
# data/config.sh.
################################################################################

set -euo pipefail

. data/config.sh

SAMPLE="$1"
OUT="$2"
shift 2

ASM="${ASSEMBLY_DIR}/${SAMPLE}/${SAMPLE}.raw.fasta"
[ -s "$ASM" ] || ASM="${ASSEMBLY_DIR}/${SAMPLE}/assembly.fasta"
[ -s "$ASM" ] || { echo "no pre-curation assembly for $SAMPLE under $ASSEMBLY_DIR" >&2; exit 1; }

: > "$OUT"

for group in "$@"; do
	name="${group%%:*}"
	contigs="${group#*:}"
	echo ">${SAMPLE}_${name}" >> "$OUT"
	IFS=',' read -ra cs <<< "$contigs"
	for c in "${cs[@]}"; do
		# Concatenating the contigs of one phage is deliberate: the template is
		# a k-mer source, not a genome, so contig order and junctions between
		# them carry no meaning and are never assembled.
		python3 - "$ASM" "$c" >> "$OUT" <<'PY'
import sys
asm, want = sys.argv[1], sys.argv[2]
keep, seq = False, []
for line in open(asm):
    if line.startswith(">"):
        keep = line[1:].split()[0] == want
    elif keep:
        seq.append(line.strip())
if not seq:
    sys.exit(f"contig {want} not found in {asm}")
sys.stdout.write("".join(seq) + "\n")
PY
	done
done

printf "%s: %d templates written to %s\n" "$SAMPLE" "$#" "$OUT"
grep -c "^>" "$OUT" > /dev/null
