#!/bin/bash

################################################################################
# Detection of multiple APSE bacteriophages and multiple Hamiltonella defensa
# strains within a single sequenced sample
#
# BACKGROUND
#
# Genome assemblies of the facultative aphid symbiont Hamiltonella defensa are
# routinely curated by removing redundant contigs: where one contig aligns
# entirely inside another, the shorter one is deleted as a duplicate. That rule
# is correct when a sample contains a single APSE prophage, which is commonly
# assembled twice (once integrated in the chromosome and once as a circular
# excised element). It is NOT correct when a sample contains two different APSE
# phages, because APSE genomes share a highly conserved backbone and differ
# mainly in a single toxin-encoding module. Two distinct phages therefore look
# like one phage plus a redundant duplicate, and the second phage is discarded.
#
# This script screens every sample for that situation using two independent
# tests, because neither is sufficient on its own:
#
#   TEST 1 (allelic discordance) detects a second phage only if the assembler
#          collapsed it into the first. Two phages divergent enough to assemble
#          separately produce no discordance and are missed.
#
#   TEST 2 (toxin content) detects a second phage only if it carries a
#          different, previously catalogued toxin. Two phages carrying the same
#          toxin family are missed.
#
# The same discordance measurement applied to the bacterial chromosome answers a
# separate question: whether the sample contains one Hamiltonella strain or
# several. This matters because two strains, each carrying a different phage,
# would produce the same phage-level signal as one strain carrying two phages.
# The chromosome is therefore both a second result and the control for the first.
#
# INTERPRETATION
#
# A clonal element sequenced with PacBio HiFi reads should be essentially
# monomorphic. In the validation sample S12, which carries a single APSE, the
# phage yields zero discordant sites at 1382-fold coverage. The assay therefore
# has no measurable false-positive rate at phage-level depth, and any non-zero
# discordance on an APSE contig is meaningful. Discordance that is dispersed
# across a contig indicates two collapsed sequences; discordance confined to one
# short block indicates a repeat, and every Hamiltonella chromosome in this
# panel contains one such block of roughly 10 kb.
#
# USAGE
#
#   bash scripts/detect_multiple_APSE_and_strains.sh
#
# Edit data/config.sh first to point at your copy of the reads and assemblies.
# Run from the root of the repository. Runtime is approximately one minute per
# sample on a laptop.
#
# DEPENDENCIES
#
#   minimap2   read and sequence alignment          (tested with 2.31-r1302)
#   samtools   alignment handling and pileup        (tested with 1.24)
#   awk        all parsing and summarisation        (POSIX awk, mawk or gawk)
#
# NOTE ON CONVENTION
#
# The other scripts in this project are SLURM array jobs written for the sciCore
# cluster. This one is deliberately a portable loop with no scheduler directives
# and no hard-coded absolute paths, so that it can be re-run by reviewers and
# collaborators on an ordinary workstation.
################################################################################

# Fail if any stage of a pipeline fails, rather than only the last one. A single
# sample that fails is reported and skipped; it does not abort the whole run.
set -o pipefail

################################################################################
# Configuration
################################################################################

# Load the machine-specific paths (APHID_DATA, READS_DIR, ASSEMBLY_DIR,
# WORK_DIR, THREADS)
source data/config.sh

# Index of sample IDs, one per line
INDEXNAME=data/sampleID.txt

# Reference panel of complete APSE genomes, used to decide which contigs in each
# assembly are phage derived
APSE_PANEL=data/reference/APSE_panel.fasta

# Coding sequences of every toxin gene catalogued across the 14 samples, used to
# determine which toxin families a sample carries
TOXIN_CDS=data/reference/toxins_CDS.fasta

# Directory for the small summary tables that are kept with the repository
RESULTS=results

################################################################################
# Analysis thresholds
#
# Every threshold is defined here with the reason for its value, so that the
# sensitivity of the screen can be assessed and adjusted by a reader.
################################################################################

# Minimum read depth at a site before its allele fractions are considered.
# Below roughly ten reads an apparent minor allele cannot be distinguished from
# sampling noise.
MIN_DEPTH=10

# Minimum number of reads carrying the alternative allele. Requiring three
# excludes isolated sequencing errors, which in HiFi data are overwhelmingly
# single-read events.
MIN_ALT_COUNT=3

# Minimum fraction of reads carrying the alternative allele. HiFi base accuracy
# is above 99.5 percent, so a genuine minor allele at 15 percent is far above
# the error floor.
MIN_ALT_FRACTION=0.15

# Minimum mapping quality. Reads that align equally well to two near-identical
# copies of the same sequence are assigned a mapping quality of zero and are
# excluded here, so that a phage present twice in one assembly does not create
# spurious discordance.
MIN_MAPQ=20

# Minimum base quality at a site
MIN_BASEQ=20

# Two discordant sites separated by more than this distance are treated as
# belonging to separate blocks. Used to distinguish a localised repeat from
# variation dispersed along a contig.
CLUSTER_GAP=20000

# A contig whose discordant sites span more than this fraction of its length is
# called dispersed rather than clustered. Used for the per-contig report.
DISPERSAL_FRACTION=0.50

# Fraction of the one-kilobase windows in a region that must contain at least
# one discordant site before the region is called dispersed. This is the measure
# the verdicts use, because it behaves correctly when a region is split across
# several contigs and when a small phage sits inside a large chromosome. A
# repeat block of ten kilobases inside a two megabase chromosome occupies well
# under one per cent of its windows; two collapsed phages typically occupy more
# than half.
MIN_OCCUPANCY=0.10

# Fraction of a contig that must align to the APSE reference panel before the
# contig is classified as phage derived
APSE_MIN_COVERAGE=0.30

# Two APSE contigs sharing more than this fraction of their sequence are treated
# as the same phage in two states, not as two different phages
SAME_PHAGE_IDENTITY=0.95

# Fraction of a toxin gene that must be covered before the toxin is called
# present, whether in an assembly or in the reads
TOXIN_MIN_BREADTH=80

# Minimum contig length to be considered part of the bacterial chromosome rather
# than a plasmid or phage
MIN_CHROM_LENGTH=100000

# Mean depth below which the discordance test cannot support a negative result.
# Samples below this report insufficient power rather than a clean negative.
MIN_POWER_DEPTH=30

################################################################################
# Preparation
################################################################################

# Confirm that the required software is available before doing any work
for TOOL in minimap2 samtools awk; do
	if ! command -v "$TOOL" > /dev/null 2>&1; then
		echo "ERROR: required program '$TOOL' was not found in PATH" >&2
		exit 1
	fi
done

# Confirm that the reference files distributed with this repository are present
for FILE in "$INDEXNAME" "$APSE_PANEL" "$TOXIN_CDS"; do
	if [ ! -s "$FILE" ]; then
		echo "ERROR: required input '$FILE' is missing or empty" >&2
		exit 1
	fi
done

# Confirm that the externally stored data has been located correctly
if [ ! -d "$READS_DIR" ] || [ ! -d "$ASSEMBLY_DIR" ]; then
	echo "ERROR: could not find the read or assembly directories." >&2
	echo "       Edit data/config.sh so that APHID_DATA points at your copy." >&2
	exit 1
fi

# Ensure the output directories exist
mkdir -p "$WORK_DIR"
mkdir -p "$RESULTS"/logs

# Write the column headings of each output table
echo -e "sample\tcontig\tlength\tclass\tmeanDepth\tsitesTested\tdiscordantSites\tper_kb\tmedianMAF\tsiteSpan\tdispersal\tblocks\tpattern" \
	> "$RESULTS"/per_contig_discordance.tsv

echo -e "sample\tregion\tmeanDepth\tsitesTested\tdiscordantSites\tper_kb\toccupancy\tpattern" \
	> "$RESULTS"/per_region_discordance.tsv

echo -e "sample\ttoxin\tfamily\tin_assembly_pct\tin_reads_breadth_pct\tin_reads_meanDepth\tcall" \
	> "$RESULTS"/toxin_screen.tsv

echo -e "sample\tcontigA\tcontigB\tsplitReads\tmodalPositionA" \
	> "$RESULTS"/split_read_links.tsv

echo -e "sample\tapseRegions\tapseMeanDepth\tsamePhageDuplicate\ttoxinFamilies\ttoxinsInReadsOnly\tapse_per_kb\tapsePattern\tchromMeanDepth\tchrom_per_kb\tchromPattern\tAPSE_verdict\tstrain_verdict\tpower" \
	> "$RESULTS"/multiplicity_summary.tsv

################################################################################
# Per sample analysis
################################################################################

while IFS= read -r NAME; do

	# Skip blank lines and comments in the sample index
	case "$NAME" in ""|\#*) continue ;; esac

	echo "Processing sample ID: $NAME"

	# Locate this sample's reads and assembly. Two assembly file names are in
	# use across the project, so both are accepted.
	READS="$READS_DIR/$NAME.fq"
	LOG="$RESULTS/logs/$NAME.log"
	ASSEMBLY=""
	for CANDIDATE in "$ASSEMBLY_DIR/$NAME/$NAME.raw.fasta" "$ASSEMBLY_DIR/$NAME/assembly.fasta"; do
		if [ -s "$CANDIDATE" ]; then
			ASSEMBLY="$CANDIDATE"
			break
		fi
	done

	# Skip the sample rather than aborting if either input is unavailable
	if [ ! -s "$READS" ] || [ -z "$ASSEMBLY" ]; then
		echo "  WARNING: missing reads or assembly for $NAME, skipping" | tee "$LOG"
		continue
	fi

	# Record the run parameters so that each result can be traced to the
	# settings that produced it
	{
		echo "sample: $NAME"
		echo "reads: $READS"
		echo "assembly: $ASSEMBLY"
		echo "minimap2: $(minimap2 --version)"
		echo "samtools: $(samtools --version | head -1)"
		echo "thresholds: depth>=$MIN_DEPTH alt>=$MIN_ALT_COUNT maf>=$MIN_ALT_FRACTION mapq>=$MIN_MAPQ"
	} > "$LOG"

	############################################################################
	# Step 1. Align the reads to the sample's own assembly
	#
	# Mapping a sample's reads back to its own assembly is what exposes
	# collapsed sequence: where two similar elements were merged into one
	# contig, reads from both pile up on that contig and disagree with each
	# other at every position where the two elements differ.
	############################################################################

	BAM="$WORK_DIR/$NAME.bam"

	if [ ! -s "$BAM" ]; then
		echo "  aligning reads to assembly"
		minimap2 -ax map-hifi -t "$THREADS" "$ASSEMBLY" "$READS" 2>> "$LOG" \
			| samtools sort -@ 4 -m 512M -o "$BAM" - 2>> "$LOG"
		samtools index "$BAM" 2>> "$LOG"
	else
		echo "  reusing existing alignment"
	fi

	############################################################################
	# Step 2. Identify which parts of the assembly are APSE derived
	#
	# Each assembly is aligned against a panel of complete APSE genomes and the
	# aligned intervals are recorded. Overlapping alignments are merged so that
	# repeated hits to the same region are not counted twice.
	#
	# Intervals are used rather than whole contigs because an APSE prophage is
	# often integrated in the chromosome rather than assembled as a separate
	# circular element. In those samples the phage occupies two per cent of a
	# two megabase contig, so any test applied to the whole contig would dilute
	# the phage signal to nothing. Restricting the analysis to the aligned
	# intervals keeps integrated and excised phages equally detectable.
	############################################################################

	echo "  locating APSE regions"

	# Record the length of every contig in the assembly
	samtools faidx "$ASSEMBLY" 2>> "$LOG"
	cut -f1,2 "$ASSEMBLY".fai > "$WORK_DIR/$NAME.contig_lengths.tsv"

	# Align the reference phage panel to the assembly and merge the aligned
	# intervals for each contig. Writes one line per merged interval.
	minimap2 -x asm10 -t "$THREADS" "$ASSEMBLY" "$APSE_PANEL" 2>> "$LOG" \
		| awk -F'\t' '{print $6"\t"$8"\t"$9}' \
		| sort -k1,1 -k2,2n \
		| awk -F'\t' '
			# Emit each merged interval as it is completed
			{
				if ($1 != contig || $2 > end) {
					if (contig != "") print contig"\t"start"\t"end
					contig = $1; start = $2; end = $3
				} else if ($3 > end) {
					end = $3
				}
			}
			END { if (contig != "") print contig"\t"start"\t"end }' \
		> "$WORK_DIR/$NAME.apse_intervals.tsv"

	# Total the phage-derived bases on each contig, used only for classifying
	# contigs in the per-contig report
	awk -F'\t' '{ covered[$1] += $3 - $2 }
		END { for (c in covered) print c"\t"covered[c] }' \
		"$WORK_DIR/$NAME.apse_intervals.tsv" > "$WORK_DIR/$NAME.apse_covered.tsv"

	# Classify every contig as APSE, chromosome or other
	awk -F'\t' -v mincov="$APSE_MIN_COVERAGE" -v minchrom="$MIN_CHROM_LENGTH" '
		FNR == NR { covered[$1] = $2; next }
		{
			frac = ($1 in covered) ? covered[$1] / $2 : 0
			if (frac >= mincov)          class = "APSE"
			else if ($2 >= minchrom)     class = "chromosome"
			else                         class = "other"
			printf "%s\t%s\t%s\t%.3f\n", $1, $2, class, frac
		}' "$WORK_DIR/$NAME.apse_covered.tsv" "$WORK_DIR/$NAME.contig_lengths.tsv" \
		> "$WORK_DIR/$NAME.contig_class.tsv"

	############################################################################
	# Step 3. Distinguish one phage in two states from two different phages
	#
	# A single APSE is frequently assembled twice, as an integrated copy and as
	# a circular excised copy. Those two contigs are nearly identical along
	# their whole length. Two different phages instead share only their backbone
	# and diverge across the toxin module. Aligning the APSE contigs against one
	# another separates these cases.
	############################################################################

	# Extract every phage region as its own sequence. Regions are used instead
	# of whole contigs so that an integrated prophage can be compared against an
	# excised one, which is exactly the case this test needs to recognise. Only
	# regions of appreciable length are considered, so that short scattered
	# alignments to conserved backbone fragments are not mistaken for a phage.
	awk -F'\t' '$3 - $2 >= 5000 {printf "%s:%d-%d\n", $1, $2 + 1, $3}' \
		"$WORK_DIR/$NAME.apse_intervals.tsv" > "$WORK_DIR/$NAME.apse_regions.txt"

	APSE_COUNT=$(wc -l < "$WORK_DIR/$NAME.apse_regions.txt")
	SAME_PHAGE="no"

	if [ "$APSE_COUNT" -gt 1 ]; then

		samtools faidx "$ASSEMBLY" $(tr '\n' ' ' < "$WORK_DIR/$NAME.apse_regions.txt") \
			> "$WORK_DIR/$NAME.apse_regions.fasta" 2>> "$LOG"

		# Align the phage regions against themselves. A pair whose alignment
		# covers most of the shorter region is the same phage twice, once
		# integrated and once excised, rather than two different phages.
		SAME_PHAGE=$(minimap2 -x asm10 -t "$THREADS" \
				"$WORK_DIR/$NAME.apse_regions.fasta" \
				"$WORK_DIR/$NAME.apse_regions.fasta" 2>> "$LOG" \
			| awk -F'\t' -v cut="$SAME_PHAGE_IDENTITY" '
				# Column 1 is the query name, 2 its length, 11 the alignment
				# block length, 6 the target name. Self hits are ignored.
				$1 != $6 && $11 / $2 >= cut { found = 1 }
				END { print (found ? "yes" : "no") }')
	fi

	############################################################################
	# Step 4. Determine which toxin families are present in the assembly
	#
	# The toxin module is the most variable part of an APSE genome, so a
	# permissive alignment preset is used. A toxin is called present when most
	# of its coding sequence is covered.
	############################################################################

	echo "  screening toxin genes"

	minimap2 -x asm20 -t "$THREADS" "$ASSEMBLY" "$TOXIN_CDS" 2>> "$LOG" \
		| awk -F'\t' '{ if ($11 > best[$1]) best[$1] = $11; len[$1] = $2 }
			END { for (t in best) printf "%s\t%.1f\n", t, 100 * best[t] / len[t] }' \
		| sort > "$WORK_DIR/$NAME.toxin_in_assembly.tsv"

	############################################################################
	# Step 5. Determine which toxin families are present in the reads
	#
	# This is the test that detects a phage which is present in the sample but
	# absent from the assembly. Reads are aligned directly to the toxin coding
	# sequences, so the result does not depend on any assembly or curation
	# decision. Breadth of coverage is reported alongside depth because a
	# partially covered toxin usually reflects cross-mapping from a related
	# family rather than genuine presence.
	############################################################################

	minimap2 -ax map-hifi -t "$THREADS" "$TOXIN_CDS" "$READS" 2>> "$LOG" \
		| samtools sort -@ 2 -m 256M -o "$WORK_DIR/$NAME.toxins.bam" - 2>> "$LOG"
	samtools index "$WORK_DIR/$NAME.toxins.bam" 2>> "$LOG"

	samtools coverage "$WORK_DIR/$NAME.toxins.bam" 2>> "$LOG" \
		| awk 'NR > 1 {print $1"\t"$6"\t"$7}' \
		| sort > "$WORK_DIR/$NAME.toxin_in_reads.tsv"

	# Combine the assembly and read evidence for each toxin and record the call
	join -t $'\t' -a 2 -e "0.0" -o 0,1.2,2.2,2.3 \
			"$WORK_DIR/$NAME.toxin_in_assembly.tsv" "$WORK_DIR/$NAME.toxin_in_reads.tsv" \
		| awk -F'\t' -v s="$NAME" -v br="$TOXIN_MIN_BREADTH" -v md="$MIN_DEPTH" '
			{
				# The toxin name encodes the family after the sample prefix.
				# The name is normalised so that variants of one family are
				# counted once: the sample prefix is removed, the suffixes
				# marking a toxin split across two contigs are removed, and the
				# result is upper-cased because the catalogue spells the same
				# family both as "CdtB" and as "cdtB".
				family = $1
				sub(/^S[0-9]+_/, "", family)
				sub(/_(First|Second)Part$/, "", family)
				family = toupper(family)
				inreads = ($3 + 0 >= br && $4 + 0 >= md)
				inasm   = ($2 + 0 >= br)
				if      (inasm && inreads) call = "present_both"
				else if (inreads)          call = "READS_ONLY"
				else if (inasm)            call = "assembly_only"
				else                       call = "absent"
				printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", s, $1, family, $2, $3, $4, call
			}' >> "$RESULTS"/toxin_screen.tsv

	############################################################################
	# Step 6. Count discordant sites along every contig
	#
	# A discordant site is a position where a minority of reads consistently
	# disagree with the assembled sequence. Isolated sequencing errors are
	# excluded by requiring several supporting reads and a minimum allele
	# fraction. Every qualifying site is written out individually so that the
	# summary statistics below can be audited and recomputed.
	############################################################################

	echo "  counting discordant sites"

	SITES="$WORK_DIR/$NAME.discordant_sites.tsv"

	samtools mpileup -f "$ASSEMBLY" -q "$MIN_MAPQ" -Q "$MIN_BASEQ" \
			-d 5000 --no-BAQ "$BAM" 2>> "$LOG" \
		| awk -F'\t' -v mind="$MIN_DEPTH" -v minalt="$MIN_ALT_COUNT" \
			-v minfrac="$MIN_ALT_FRACTION" -v minchrom="$MIN_CHROM_LENGTH" \
			-v depthfile="$WORK_DIR/$NAME.contig_depth.tsv" \
			-v regionfile="$WORK_DIR/$NAME.region_depth.tsv" '
			# The first input file lists the phage intervals, the second the
			# contig lengths; the pileup itself arrives on standard input.
			FNR == NR && FILENAME ~ /apse_intervals/ {
				n = ++nint[$1]; istart[$1, n] = $2; iend[$1, n] = $3; next
			}
			FILENAME ~ /contig_lengths/ { clen[$1] = $2; next }
			{
				contig = $1; pos = $2; depth = $4; bases = $5

				# Assign this position to a region. Phage intervals take
				# precedence, so a prophage integrated in the chromosome is
				# measured as phage rather than as chromosome.
				region = "other"
				if (clen[contig] >= minchrom) region = "chromosome"
				for (k = 1; k <= nint[contig]; k++) {
					if (pos > istart[contig, k] && pos <= iend[contig, k]) {
						region = "APSE"
						break
					}
				}

				# Accumulate depth across every position so that mean coverage
				# can be reported per contig and per region
				sum[contig] += depth
				covered[contig]++
				rsum[region] += depth
				rcovered[region]++

				if (depth < mind) next

				# Count only the positions that carry enough reads to be
				# assessed. Discordance rates are expressed per kilobase of
				# this analysable sequence, not per kilobase of contig, so that
				# poorly covered regions do not dilute the rate.
				tested[contig]++
				rtested[region]++

				# Walk the pileup string. A caret introduces a read start and is
				# followed by a mapping quality character; a dollar marks a read
				# end; a plus or minus introduces an indel whose length is given
				# by the following number.
				n = length(bases); i = 1; matches = 0
				delete count
				while (i <= n) {
					c = substr(bases, i, 1)
					if (c == "^") { i += 2; continue }
					if (c == "$") { i += 1; continue }
					if (c == "+" || c == "-") {
						j = i + 1; num = ""
						while (j <= n && substr(bases, j, 1) ~ /[0-9]/) {
							num = num substr(bases, j, 1); j++
						}
						i = j + num + 0
						continue
					}
					if (c == "." || c == ",") {
						matches++
					} else {
						u = toupper(c)
						if (u == "A" || u == "C" || u == "G" || u == "T") count[u]++
					}
					i++
				}

				# Identify the most frequent alternative allele
				total = matches; best = 0; allele = "."
				for (b in count) {
					total += count[b]
					if (count[b] > best) { best = count[b]; allele = b }
				}
				if (total < mind || best < minalt) next

				frac = best / total
				if (frac < minfrac) next

				rdiscord[region]++
				printf "%s\t%d\t%s\t%d\t%d\t%.4f\t%s\n",
					contig, pos, allele, best, total, frac, region
			}
			END {
				# Emit mean depth and the number of analysable positions per
				# contig
				for (c in covered)
					printf "%s\t%.1f\t%d\n", c, sum[c] / covered[c],
						(c in tested ? tested[c] : 0) > depthfile

				# Emit the same figures per region, which is what the verdicts
				# are based on
				for (r in rcovered)
					printf "%s\t%.1f\t%d\t%d\n", r, rsum[r] / rcovered[r],
						(r in rtested ? rtested[r] : 0),
						(r in rdiscord ? rdiscord[r] : 0) > regionfile
			}' "$WORK_DIR/$NAME.apse_intervals.tsv" "$WORK_DIR/$NAME.contig_lengths.tsv" - \
		> "$SITES"

	############################################################################
	# Step 7. Summarise each contig and describe the spatial pattern
	#
	# The distinction that matters is between discordance dispersed along a
	# contig, which indicates two collapsed sequences, and discordance confined
	# to one short block, which indicates a repeat. Every Hamiltonella
	# chromosome examined here contains one such block of roughly 10 kb, and
	# without this step it would be reported as evidence of multiple strains in
	# every sample.
	############################################################################

	sort -k1,1 -k2,2n "$SITES" > "$SITES".sorted

	awk -F'\t' -v s="$NAME" -v gap="$CLUSTER_GAP" -v dispfrac="$DISPERSAL_FRACTION" '
		# First file: contig classification and length
		FNR == NR && FILENAME ~ /contig_class/ { len[$1] = $2; class[$1] = $3; next }
		# Second file: mean depth per contig
		FILENAME ~ /contig_depth/ { depth[$1] = $2; tested[$1] = $3; next }
		# Third file: the discordant sites themselves
		{
			c = $1
			n[c]++
			if (!(c in first)) { first[c] = $2; blocks[c] = 1 }
			else if ($2 - last[c] > gap) blocks[c]++
			last[c] = $2
			mafs[c] = mafs[c] " " $6
		}
		END {
			for (c in len) {
				count = (c in n) ? n[c] : 0
				span  = (count > 0) ? last[c] - first[c] : 0
				disp  = (len[c] > 0) ? span / len[c] : 0
				perkb = (tested[c] > 0) ? 1000 * count / tested[c] : 0

				# Median minor allele fraction, computed by sorting the values
				med = 0
				if (count > 0) {
					k = split(mafs[c], v, " ")
					for (a = 2; a <= k; a++) {
						tmp = v[a]; b = a - 1
						while (b >= 1 && v[b] + 0 > tmp + 0) { v[b+1] = v[b]; b-- }
						v[b+1] = tmp
					}
					med = v[int((k + 1) / 2)] + 0
				}

				if (count == 0)                    pattern = "clean"
				else if (disp >= dispfrac)         pattern = "DISPERSED"
				else if (blocks[c] >= 3)           pattern = "DISPERSED"
				else                               pattern = "clustered"

				printf "%s\t%s\t%d\t%s\t%.1f\t%d\t%d\t%.2f\t%.2f\t%d\t%.3f\t%d\t%s\n",
					s, c, len[c], class[c], depth[c], tested[c], count, perkb, med,
					span, disp, (count > 0 ? blocks[c] : 0), pattern
			}
		}' "$WORK_DIR/$NAME.contig_class.tsv" "$WORK_DIR/$NAME.contig_depth.tsv" "$SITES".sorted \
		| sort -k4,4 -k3,3nr >> "$RESULTS"/per_contig_discordance.tsv

	############################################################################
	# Step 7b. Summarise each region
	#
	# The verdicts are based on regions rather than contigs, because a phage may
	# be integrated inside a chromosome contig and a chromosome may be split
	# across several contigs. Occupancy, the fraction of one-kilobase windows
	# holding at least one discordant site, is used in place of the per-contig
	# span because it behaves correctly in both situations.
	############################################################################

	awk -F'\t' -v s="$NAME" -v occ="$MIN_OCCUPANCY" '
		# First file: mean depth, analysable positions and site count per region
		FNR == NR { depth[$1] = $2; tested[$1] = $3; sites[$1] = $4; next }
		# Second file: the discordant sites, tagged with their region
		{ bins[$7, $1 "_" int($2 / 1000)] = 1 }
		END {
			for (r in depth) {
				# Count the distinct kilobase windows containing a site
				hit = 0
				for (k in bins) {
					split(k, f, SUBSEP)
					if (f[1] == r) hit++
				}
				kb = tested[r] / 1000
				occupancy = (kb > 0) ? hit / kb : 0
				perkb = (tested[r] > 0) ? 1000 * sites[r] / tested[r] : 0
				pattern = (sites[r] == 0) ? "clean" \
					: ((occupancy >= occ) ? "DISPERSED" : "clustered")
				printf "%s\t%s\t%.1f\t%d\t%d\t%.2f\t%.4f\t%s\n",
					s, r, depth[r], tested[r], sites[r], perkb, occupancy, pattern
			}
		}' "$WORK_DIR/$NAME.region_depth.tsv" "$SITES" \
		| sort -k2,2 >> "$RESULTS"/per_region_discordance.tsv

	############################################################################
	# Step 8. Record reads that span two contigs
	#
	# A read whose alignment is split between two contigs physically connects
	# them. Splits between a phage contig and the chromosome locate the site at
	# which the prophage is integrated; splits between two phage contigs mark
	# the boundary at which they diverge.
	############################################################################

	samtools view -q 1 "$BAM" 2>> "$LOG" \
		| awk -F'\t' -v s="$NAME" '
			{
				sa = ""
				for (i = 12; i <= NF; i++) if ($i ~ /^SA:Z:/) sa = $i
				if (sa == "") next
				split(sa, a, ":"); n = split(a[3], parts, ";")
				for (j = 1; j <= n; j++) {
					if (parts[j] == "") continue
					split(parts[j], p, ",")
					if (p[1] != $3) {
						key = $3 SUBSEP p[1]
						count[key]++
						# Record the modal alignment start, rounded to 1 kb
						posbin[key SUBSEP int($4 / 1000)]++
					}
				}
			}
			END {
				for (k in count) {
					split(k, f, SUBSEP)
					best = 0; bestpos = 0
					for (pk in posbin) {
						split(pk, g, SUBSEP)
						if (g[1] == f[1] && g[2] == f[2] && posbin[pk] > best) {
							best = posbin[pk]; bestpos = g[3]
						}
					}
					if (count[k] >= 5)
						printf "%s\t%s\t%s\t%d\t%d\n", s, f[1], f[2], count[k], bestpos * 1000
				}
			}' \
		| sort -k4,4nr >> "$RESULTS"/split_read_links.tsv

	############################################################################
	# Step 9. Reach a verdict for this sample
	#
	# The two verdicts are deliberately conservative. A sample is only called
	# negative when there was enough sequencing depth for a negative result to
	# be meaningful; otherwise it is reported as having insufficient power, so
	# that absence of evidence is not recorded as evidence of absence.
	############################################################################

	# Read the phage and chromosome figures from the region summary
	read -r APSE_DEPTH APSE_PERKB APSE_PATTERN <<< "$(awk -F'\t' -v s="$NAME" '
		$1 == s && $2 == "APSE" { printf "%.1f %.2f %s", $3, $6, $8; found = 1 }
		END { if (!found) printf "0.0 0.00 absent" }' \
		"$RESULTS"/per_region_discordance.tsv)"

	read -r CHROM_DEPTH CHROM_PERKB CHROM_PATTERN <<< "$(awk -F'\t' -v s="$NAME" '
		$1 == s && $2 == "chromosome" { printf "%.1f %.2f %s", $3, $6, $8; found = 1 }
		END { if (!found) printf "0.0 0.00 absent" }' \
		"$RESULTS"/per_region_discordance.tsv)"

	# Count the distinct toxin families found. A family counts as present when
	# either the assembly or the reads support it: a toxin detected only in the
	# assembly is still a toxin the sample carries, and in the most weakly
	# covered samples the assembly is the only evidence available. Counting is
	# done at the level of the family rather than the individual gene, so that
	# near-identical variants of one toxin are not mistaken for two phages.
	TOXIN_FAMILIES=$(awk -F'\t' -v s="$NAME" '$1 == s && $7 != "absent" {print $3}' \
		"$RESULTS"/toxin_screen.tsv | sort -u | paste -sd, -)
	TOXIN_FAMILY_COUNT=$(awk -F'\t' -v s="$NAME" '$1 == s && $7 != "absent" {print $3}' \
		"$RESULTS"/toxin_screen.tsv | sort -u | wc -l)
	TOXIN_READS_ONLY=$(awk -F'\t' -v s="$NAME" '$1 == s && $7 == "READS_ONLY" {print $3}' \
		"$RESULTS"/toxin_screen.tsv | sort -u | paste -sd, -)

	# Decide whether the sample carries more than one phage
	if [ "$TOXIN_FAMILY_COUNT" -ge 2 ]; then
		APSE_VERDICT="MULTIPLE_APSE_distinct_toxins"
	elif [ -n "$TOXIN_READS_ONLY" ]; then
		APSE_VERDICT="MULTIPLE_APSE_toxin_missing_from_assembly"
	elif [ "$APSE_PATTERN" = "DISPERSED" ] && [ "$SAME_PHAGE" = "no" ]; then
		APSE_VERDICT="MULTIPLE_APSE_collapsed"
	elif [ "$APSE_PATTERN" = "DISPERSED" ] && [ "$SAME_PHAGE" = "yes" ]; then
		APSE_VERDICT="review_same_phage_duplicated"
	elif awk -v d="$APSE_DEPTH" -v m="$MIN_POWER_DEPTH" 'BEGIN{exit !(d < m)}'; then
		APSE_VERDICT="insufficient_power"
	else
		APSE_VERDICT="single_APSE"
	fi

	# Decide whether the sample contains more than one Hamiltonella strain.
	# Dispersed chromosomal discordance is the signature of a strain mixture;
	# discordance confined to one block is the repeat present in every sample.
	if [ "$CHROM_PATTERN" = "DISPERSED" ]; then
		STRAIN_VERDICT="MULTIPLE_STRAINS_possible"
	elif awk -v d="$CHROM_DEPTH" -v m="$MIN_POWER_DEPTH" 'BEGIN{exit !(d < m)}'; then
		STRAIN_VERDICT="insufficient_power"
	else
		STRAIN_VERDICT="single_strain"
	fi

	# Record whether the depth was adequate for these conclusions
	POWER=$(awk -v a="$APSE_DEPTH" -v c="$CHROM_DEPTH" -v m="$MIN_POWER_DEPTH" \
		'BEGIN { printf "%s", (a >= m && c >= m) ? "adequate" : "LIMITED" }')

	# Append this sample's row to the summary table
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
		"$NAME" "$APSE_COUNT" "$APSE_DEPTH" "$SAME_PHAGE" \
		"${TOXIN_FAMILIES:-none}" "${TOXIN_READS_ONLY:-none}" \
		"$APSE_PERKB" "$APSE_PATTERN" \
		"$CHROM_DEPTH" "$CHROM_PERKB" "$CHROM_PATTERN" \
		"$APSE_VERDICT" "$STRAIN_VERDICT" "$POWER" \
		>> "$RESULTS"/multiplicity_summary.tsv

	echo "  APSE: $APSE_VERDICT   strains: $STRAIN_VERDICT   power: $POWER"

done < "$INDEXNAME"

################################################################################
# Report
################################################################################

echo
echo "Screen complete. Summary tables written to $RESULTS:"
echo "  multiplicity_summary.tsv     one row per sample, with both verdicts"
echo "  per_contig_discordance.tsv   one row per contig"
echo "  toxin_screen.tsv             toxin presence in assembly and in reads"
echo "  split_read_links.tsv         reads spanning two contigs"
echo
echo "Intermediate alignments were kept in $WORK_DIR and can be deleted."
