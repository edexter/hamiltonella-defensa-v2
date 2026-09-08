#!/bin/bash

################################################################################
# Report the version of every tool the analysis uses, and say plainly which ones
# are missing or differ from the versions the work was run with.
#
# Run this before anything else. A version mismatch is not necessarily a
# problem, but discovering one after a six-hour assembly is worse than
# discovering it now.
################################################################################

EXPECTED_minimap2="2.31"
EXPECTED_samtools="1.24"
EXPECTED_seqtk="1.5"
EXPECTED_flye="2.9.6"
EXPECTED_busco="5.1.2"

status=0

check () {
	local name="$1" cmd="$2" want="$3" got
	if ! command -v "$cmd" >/dev/null 2>&1; then
		printf "  %-10s MISSING   (expected %s)\n" "$name" "$want"
		status=1
		return
	fi
	got=$("${@:4}" 2>&1 | head -5 | grep -oE "[0-9]+\.[0-9]+(\.[0-9]+)?(-r[0-9]+)?" | head -1)
	if [ -z "$got" ]; then
		printf "  %-10s present, version not determined (expected %s)\n" "$name" "$want"
	elif [[ "$got" == "$want"* ]]; then
		printf "  %-10s %s\n" "$name" "$got"
	else
		printf "  %-10s %s  <- differs from %s used here\n" "$name" "$got" "$want"
	fi
}

echo "Tools used by the analysis scripts:"
check minimap2 minimap2 "$EXPECTED_minimap2" minimap2 --version
check samtools samtools "$EXPECTED_samtools" samtools --version
check seqtk    seqtk    "$EXPECTED_seqtk"    seqtk
check flye     flye     "$EXPECTED_flye"     flye --version

echo
echo "Tool used for completeness, which lives in its own environment:"
check busco busco "$EXPECTED_busco" busco --version

echo
echo "Python and the one non-standard-library import:"
python3 --version 2>&1 | sed 's/^/  /'
python3 -c "import matplotlib; print('  matplotlib', matplotlib.__version__)" 2>/dev/null \
	|| echo "  matplotlib MISSING (needed only for the figure scripts)"

echo
if [ -f data/config.sh ]; then
	# shellcheck disable=SC1091
	. data/config.sh
	echo "Data paths from data/config.sh:"
	for v in READS_DIR ASSEMBLY_DIR WORK_DIR; do
		p="${!v}"
		if [ -d "$p" ]; then printf "  %-12s %s\n" "$v" "$p"
		else printf "  %-12s %s  <- does not exist\n" "$v" "$p"; status=1; fi
	done
else
	echo "data/config.sh not found; copy it and set the two paths at the top."
	status=1
fi

exit $status
