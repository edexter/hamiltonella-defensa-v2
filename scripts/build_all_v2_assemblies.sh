#!/bin/bash

################################################################################
# Build all fourteen genome assemblies, in order, from the curated assemblies
# and the phage genomes.
#
# WHY THIS SCRIPT EXISTS
#
# The per-sample arguments to build_v2_assemblies.py -- which contig is
# substituted, which contigs are dropped, and the strain name -- cannot be
# derived from the inputs. They are recorded in data/assembly_manifest.tsv and
# this script applies them. Without the two together the assemblies cannot be
# regenerated at all.
#
# WHAT IT DOES, IN ORDER
#
#   1. For each sample that has an integration-site table, splice the prophage
#      into the chromosome at the att core (lib/integrate_prophage.py).
#   2. Build the assembly, substituting spliced contigs, dropping the contigs
#      the manifest names, and correcting the headers (build_v2_assemblies.py).
#   3. Classify every contig from taxonomy, cross-panel comparison, topology and
#      depth (lib/classify_replicons.py).
#   4. Apply those labels to the assemblies (apply_replicon_labels.py). This is
#      the last step to touch the sequence files.
#
# ONE STEP IS NOT AUTOMATED. S10 needs two scaffold gaps filled from the
# assembly graph before it is built. That was done once, by hand, and is
# described in results/v2_provenance.tsv. Point S10_GAPFILLED at the resulting
# chromosome, or the script will build S10 with the gaps still in place and say
# so.
#
# VERIFYING THE RESULT
#
# Compare the rebuilt assemblies against results/v2_checksums.txt. Every file
# should match. The samples that need no substitution and no drop are the
# control: if one of those changes, something is wrong with the pipeline rather
# than with that sample.
#
# Usage:  bash scripts/build_all_v2_assemblies.sh [outdir]
################################################################################

set -euo pipefail

. data/config.sh

OUTDIR="${1:-assemblies_v2}"
MANIFEST="data/assembly_manifest.tsv"
PHAGES="results/phage_reconstructions"
REFS="results/phage_references"
WORK="${WORK_DIR%/*}/build_v2"
S10_GAPFILLED="${S10_GAPFILLED:-}"

mkdir -p "$OUTDIR" "$WORK"

# The phage genomes that get spliced into chromosomes, gathered into one file
# per sample so integrate_prophage.py can look them up by name.
gather_phages () {
	local samp="$1" out="$WORK/${samp}_phages.fa"
	cat "$PHAGES/${samp}_"*.fasta "$REFS/${samp}_"*.fasta 2>/dev/null > "$out" || true
	[ -s "$out" ] || { echo "  no phage genome found for $samp" >&2; return 1; }
	echo "$out"
}

built=0
skipped=0

while IFS=$'\t' read -r samp strain replace drop note; do
	case "$samp" in \#*|sample|"") continue ;; esac

	curated="${ASSEMBLY_DIR}/${samp}/${samp}.curated.fasta"
	info="${ASSEMBLY_DIR}/${samp}/assembly_info.txt"
	if [ ! -s "$curated" ]; then
		echo "$samp SKIPPED: $curated not found"
		skipped=$((skipped+1))
		continue
	fi

	args=()

	if [ "$replace" != "-" ]; then
		# Each entry is contig=source, where source is "splice" or "gapfill".
		IFS=';' read -ra items <<< "$replace"
		for item in "${items[@]}"; do
			contig="${item%%=*}"
			source="${item##*=}"

			if [ "$source" = "splice" ]; then
				sites="results/${samp}_integration_sites.tsv"
				[ -s "$sites" ] || { echo "$samp: $sites missing" >&2; exit 1; }
				phages=$(gather_phages "$samp")
				spliced="$WORK/${samp}_${contig}.fa"
				# integrate_prophage.py verifies the att core against both
				# sequences and stops rather than guessing if it does not match.
				python3 scripts/lib/integrate_prophage.py \
					"$curated" "$phages" "$sites" "$spliced"
				args+=("replace:${contig}=${spliced}")

			elif [ "$source" = "gapfill" ]; then
				if [ -n "$S10_GAPFILLED" ] && [ -s "$S10_GAPFILLED" ]; then
					args+=("replace:${contig}=${S10_GAPFILLED}")
				else
					echo "$samp: building WITHOUT the scaffold gap fill;" \
					     "set S10_GAPFILLED to include it"
				fi
			fi
		done
	fi

	if [ "$drop" != "-" ]; then
		IFS=';' read -ra dc <<< "$drop"
		for c in "${dc[@]}"; do args+=("drop:${c}"); done
	fi

	python3 scripts/build_v2_assemblies.py \
		"$samp" "$curated" "$info" "$strain" "$OUTDIR/${samp}_v2.fasta" "${args[@]}"

	printf "%s built (%d contigs)\n" "$samp" "$(grep -c '^>' "$OUTDIR/${samp}_v2.fasta")"
	built=$((built+1))
done < "$MANIFEST"

echo
echo "$built assemblies built, $skipped skipped"
[ "$skipped" -eq 0 ] || echo "SOME SAMPLES WERE SKIPPED - the run is not complete"

echo
echo "Classifying contigs and applying labels"
python3 scripts/lib/classify_replicons.py data/config.sh results/replicon_table.tsv
python3 scripts/apply_replicon_labels.py results/replicon_table.tsv "$OUTDIR" "$OUTDIR"

echo
echo "Now compare against results/v2_checksums.txt:"
echo "    md5sum -c results/v2_checksums.txt"
