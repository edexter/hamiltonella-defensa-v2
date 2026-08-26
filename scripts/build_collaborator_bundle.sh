#!/bin/bash

################################################################################
# Assemble a working bundle of the version 2.0 results for a collaborator
#
# WHAT THIS IS FOR
#
# A colleague who wants to start analysing these genomes needs the sequences and
# the annotation, and very little else. This builds that: assemblies, phage
# genomes, PGAP annotation, and the few tables needed to know what each contig
# is. It is deliberately not the full delivery -- the scripts, the method
# documentation and the audit trail are a separate handover.
#
# WHAT IS LEFT OUT, AND WHY
#
# The PGAP output is sent complete, including the assembly that was fed in. That
# duplicates what is in assemblies/, but it means each annotation directory
# carries the exact input it was produced from, so there is never any question
# about which version of an assembly a given annotation describes.
#
# One file is dropped: cwltool.log, at roughly 20 MB per sample. It is a record
# of how the pipeline ran rather than a result, and across fourteen samples it
# would account for about 280 MB of the bundle on its own. If a run ever needs
# debugging, the log is still on the instance and in the full handover.
#
# USAGE
#
#   bash scripts/build_collaborator_bundle.sh <user@host> <ssh-key> [outdir]
#
# The first two arguments say where the PGAP results are. The annotation is
# trimmed on the remote machine before transfer, so the 20 MB logs are never
# copied across the network.
#
# OUTPUT
#
#   <outdir>/hamiltonella_v2_working/       the bundle
#   <outdir>/hamiltonella_v2_working.zip    the same, compressed
#
# DEPENDENCIES
#
#   ssh, rsync, zip
################################################################################

set -o pipefail

REMOTE="$1"
KEY="$2"
OUTDIR="${3:-$HOME/aphid_work/bundle}"

if [ -z "$REMOTE" ] || [ -z "$KEY" ]; then
	echo "usage: bash scripts/build_collaborator_bundle.sh <user@host> <ssh-key> [outdir]" >&2
	exit 1
fi

# Everything PGAP produced is sent except the pipeline log, which is large and
# is a record of the run rather than a result.
DROP="cwltool.log"

BUNDLE="$OUTDIR/hamiltonella_v2_working"
rm -rf "$BUNDLE" "$BUNDLE.zip"
mkdir -p "$BUNDLE"/{assemblies,phages,annotation,docs}

################################################################################
# Trim the annotation on the remote machine, so that the large files are never
# transferred at all.
################################################################################
echo "  trimming annotation on the remote machine"
ssh -i "$KEY" -o ConnectTimeout=25 "$REMOTE" "
	rm -rf ~/bundle_staging && mkdir -p ~/bundle_staging
	for d in ~/S*_results; do
		[ -d \"\$d\" ] || continue
		s=\$(basename \"\$d\" _results)
		[ -s \"\$d/annot.gbk\" ] || { echo \"    \$s: no annotation, skipping\"; continue; }
		cp -r \"\$d\" ~/bundle_staging/\$s
		for f in $DROP; do rm -f ~/bundle_staging/\$s/\$f; done
	done
	du -sh ~/bundle_staging | awk '{print \"    staged: \"\$1}'
	ls ~/bundle_staging | wc -l | awk '{print \"    samples: \"\$1}'
"

echo "  transferring"
rsync -a --info=progress2 -e "ssh -i $KEY -o ConnectTimeout=25" \
	"$REMOTE:~/bundle_staging/" "$BUNDLE/annotation/" 2>/dev/null | tail -1

################################################################################
# The local parts: assemblies, phage genomes, and the tables that say what each
# contig is.
################################################################################
echo "  adding assemblies and phage genomes"
cp assemblies_v2/S*_v2.fasta "$BUNDLE/assemblies/"
cp results/phage_references/*.fasta results/phage_reconstructions/*.fasta \
	"$BUNDLE/phages/" 2>/dev/null
rm -f "$BUNDLE/phages/"*.fai

cp results/contig_inventory.txt "$BUNDLE/docs/" 2>/dev/null
cp results/phage_lineage_names.tsv "$BUNDLE/docs/" 2>/dev/null
cp docs/curated_assembly_stats_v2.xlsx "$BUNDLE/docs/" 2>/dev/null
cp "$HOME/aphid_work/docs/Assembly summaries.docx" "$BUNDLE/docs/" 2>/dev/null

################################################################################
# A note saying what the bundle is, so it does not have to be explained twice.
################################################################################
cat > "$BUNDLE/README.txt" <<'NOTE'
Hamiltonella defensa, version 2.0 -- working results
====================================================

This is an interim bundle so that analysis can start. It is not the final
handover: the scripts, the method documentation and the audit trail follow
separately once they have been tidied.

  assemblies/   One genome per sample, fourteen in all. Chromosome, plasmids,
                and the APSE prophage in place where the evidence supports it.

  phages/       The APSE genomes as standalone sequences, one file per phage per
                sample. S05 carries three lineages plus a second allele of one of
                them; S07 carries two; the rest carry one each.

  annotation/   PGAP annotation, one directory per sample. The NCBI submission
                file and the pipeline logs have been left out; everything used
                for analysis is here.

  docs/         Assembly summaries.docx      what changed per sample, and why
                curated_assembly_stats_v2.xlsx  assembly and BUSCO statistics
                contig_inventory.txt         one line per contig saying what it is
                phage_lineage_names.tsv      which lineage each phage belongs to

Two things worth knowing before using these
-------------------------------------------

Topology in the FASTA headers is a claim about the molecule, tested against
reads, not a copy of what the assembler reported. The assembler's own flag was
wrong in both directions in this panel. No chromosome is circular: none could be
closed, and that is a limitation of the data rather than an oversight.

A shared toxin gene does not mean a shared phage lineage. The panel contains
three distinct CdtB-carrying lineages and three distinct RHS-carrying ones, and
two phages sharing both toxin and integrase type can share as little as 59% of
their 31-mers. The numeral in the lineage name is what distinguishes them.

Known limitations are listed in the summaries document, per sample, rather than
being left for the reader to discover.
NOTE

################################################################################
# Compress, and report what went in.
################################################################################
echo "  compressing"
( cd "$OUTDIR" && zip -qr hamiltonella_v2_working.zip hamiltonella_v2_working )

echo
echo "  bundle contents:"
for d in assemblies phages annotation docs; do
	printf "    %-12s %3d items  %8s\n" "$d" \
		"$(find "$BUNDLE/$d" -maxdepth 1 -mindepth 1 | wc -l)" \
		"$(du -sh "$BUNDLE/$d" | cut -f1)"
done
echo
du -sh "$BUNDLE" | awk '{print "    uncompressed: "$1}'
du -sh "$BUNDLE.zip" | awk '{print "    zipped:       "$1}'
echo
echo "  $BUNDLE.zip"
