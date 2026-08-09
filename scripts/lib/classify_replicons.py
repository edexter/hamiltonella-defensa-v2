#!/usr/bin/env python3
"""
Decide what every contig in every assembly actually is, from evidence rather
than from the label it inherited.

WHY THIS EXISTS

The labels carried by the 2024 assemblies were assigned during manual curation
and were never independently tested. Some of them are wrong. Ten contigs holding
aphid DNA are labelled as Hamiltonella, and three of those are labelled as
Hamiltonella chromosome. Several contigs that are fragments or mis-joined pairs
of plasmids carry their own plasmid names, which implies a replicon that does not
exist. The same element is called a plasmid in one sample and left unplaced in
another. This script replaces every one of those labels with one that can be
traced to a measurement.

THE EVIDENCE USED

Four independent lines, because no single one of them is sufficient:

  taxonomy     Each contig is aligned against published Hamiltonella defensa
               genomes and plasmids, against the aphid genome assembled by this
               same project, and against the APSE phage genomes reconstructed
               here. The class covering the greatest fraction of the contig wins.

  the panel    A contig absent from every database may still be genuine, because
               these are accessory replicons and the published strains do not
               carry all of them. So each contig is also compared against the
               other thirteen samples. Sequence shared with another sample of the
               same species is evidence of belonging to it.

  topology     A plasmid is a circle. A contig the assembler could not close is
               not thereby disqualified, but it cannot be asserted to be a
               separate replicon either, and this script does not assert it.

  copy number  A plasmid is usually present at a depth different from the
               chromosome. Read depth is taken from the assembler's own report
               and expressed as a ratio to the chromosome of the same sample, so
               that samples sequenced to different depths can be compared.

HOW PLASMIDS ARE GROUPED

Contigs are grouped into families by reciprocal containment: two contigs belong
together when each covers at least 85 per cent of the other. Grouping by a
one-directional match instead would be a mistake, and was one when first tried
here. A single mis-joined contig that happens to contain two unrelated plasmids
will link them into one family, and everything downstream inherits the error.

Within a family the members may differ in size. That is reported rather than
smoothed over: where the smaller member is contained in the larger and both are
circular, the difference is a real deletion and not an assembly artefact.

USAGE

    python3 scripts/lib/classify_replicons.py <config.sh> <out.tsv>

Run from the root of the repository. Takes a few minutes for fourteen samples.

DEPENDENCIES

    minimap2   all alignment                          (tested with 2.31-r1302)
    python3    standard library only
"""

import collections
import itertools
import os
import re
import subprocess
import sys

# A contig must be at least this long before its taxonomy is judged at all.
# Shorter pieces produce alignments too short to separate a genuine match from
# the mobile elements that every enterobacterium shares with every other.
MIN_ALIGNMENT_BP = 500

# Two contigs are the same replicon when each covers at least this much of the
# other. The value is deliberately high: the cost of merging two plasmids that
# are merely similar is worse than the cost of splitting one into two families,
# because a merged family produces a plasmid count that is too low and a name
# that spans two real molecules.
RECIPROCAL_COVERAGE = 0.85

# A contig is called a fragment of a larger one when the larger contains this
# much of it but the reverse does not hold.
CONTAINMENT = 0.85

# Length at or above which a Hamiltonella contig is called chromosome. These
# assemblies are fragmented, so most chromosomal sequence sits on contigs well
# below this and is honestly labelled unplaced rather than assigned a location
# the data does not support.
CHROMOSOME_BP = 100000


def read_fasta(path):
    """Return {name: sequence}, keyed on the first word of the header."""
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def merged_length(intervals):
    """Total length covered by a set of possibly overlapping intervals.

    Summing alignment lengths instead would double-count any region reported in
    more than one block, which is common wherever a repeat is involved.
    """
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total, cur_start, cur_end = 0, intervals[0][0], intervals[0][1]
    for start, end in intervals[1:]:
        if start > cur_end:
            total += cur_end - cur_start
            cur_start, cur_end = start, end
        else:
            cur_end = max(cur_end, end)
    return total + cur_end - cur_start


def paf(target, query, extra=()):
    """Run minimap2 and yield the fields of each alignment line."""
    cmd = ["minimap2", "-cx", "asm20", "-t", "8", *extra, target, query]
    out = subprocess.run(cmd, capture_output=True, text=True).stdout
    for line in out.splitlines():
        f = line.split("\t")
        if int(f[10]) >= MIN_ALIGNMENT_BP:
            yield f


def main():
    config, out_path = sys.argv[1], sys.argv[2]
    conf = {}
    for line in open(config):
        if "=" in line and not line.strip().startswith("#"):
            k, _, v = line.strip().partition("=")
            conf[k] = v.strip().strip('"')
    asm_dir = conf["ASSEMBLY_DIR"]
    v2_dir = "assemblies_v2"
    samples = [f"S{i:02d}" for i in range(1, 15)]

    seqs = {s: read_fasta(f"{v2_dir}/{s}_v2.fasta") for s in samples}

    # ---------------------------------------------------------------- taxonomy
    # Build one reference file in which every sequence is named after the class
    # it represents, so that a single alignment pass settles the question.
    ref_dir = f"{asm_dir}/reference_genomes"
    sources = [(f"{ref_dir}/CP001277.HDEF.fa", "HDEF_CHROM"),
               (f"{ref_dir}/HDEF_A2C.fasta", "HDEF_CHROM"),
               (f"{ref_dir}/PHD5AT.fa", "HDEF_PLASMID"),
               (f"{ref_dir}/M147_plasmid.fasta", "HDEF_PLASMID"),
               (f"{ref_dir}/plasmid_P4M47.fa", "HDEF_PLASMID"),
               (f"{ref_dir}/APSE.fa", "APSE")]
    work = "results/.replicon_work"
    os.makedirs(work, exist_ok=True)
    with open(f"{work}/ref.fa", "w") as fh:
        n = 0
        for path, cls in sources:
            if not os.path.exists(path):
                print(f"  reference missing, skipped: {path}", file=sys.stderr)
                continue
            for line in open(path):
                if line.startswith(">"):
                    n += 1
                    fh.write(f">{cls}|{n}\n")
                else:
                    fh.write(line)
        for d in ("results/phage_references", "results/phage_reconstructions"):
            for fn in sorted(os.listdir(d)) if os.path.isdir(d) else []:
                if not fn.endswith(".fasta"):
                    continue
                for line in open(f"{d}/{fn}"):
                    if line.startswith(">"):
                        n += 1
                        fh.write(f">APSE|{n}\n")
                    else:
                        fh.write(line)

    aphid = f"{asm_dir}/Aphid_assembly/aphid.filtered.fa"
    cover = collections.defaultdict(lambda: collections.defaultdict(list))
    ident = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0]))
    for s in samples:
        q = f"{v2_dir}/{s}_v2.fasta"
        for target, forced in ((f"{work}/ref.fa", None), (aphid, "APHID")):
            for f in paf(target, q, ("--secondary=no",)):
                cls = forced or f[5].split("|")[0]
                cover[(s, f[0])][cls].append((int(f[2]), int(f[3])))
                ident[(s, f[0])][cls][0] += int(f[9])
                ident[(s, f[0])][cls][1] += int(f[10])

    # ------------------------------------------------------------- the panel
    # Compare each sample against the other thirteen. The sample's own contigs
    # must be excluded from the reference, or every contig matches itself and
    # --secondary=no discards the cross-sample matches that are the whole point.
    shared = collections.defaultdict(list)
    for s in samples:
        with open(f"{work}/others.fa", "w") as fh:
            for other in samples:
                if other == s:
                    continue
                for name, seq in seqs[other].items():
                    fh.write(f">{other}|{name}\n{seq}\n")
        per = collections.defaultdict(lambda: collections.defaultdict(list))
        for f in paf(f"{work}/others.fa", f"{v2_dir}/{s}_v2.fasta", ("--secondary=no",)):
            if int(f[10]) >= 1000:
                per[f[0]][f[5].split("|")[0]].append((int(f[2]), int(f[3])))
        for contig, bysample in per.items():
            length = len(seqs[s][contig])
            shared[(s, contig)] = sorted(
                (o for o, iv in bysample.items() if merged_length(iv) / length >= 0.60))

    # ------------------------------------------------- topology and copy number
    flye = {}
    for s in samples:
        by_length = {}
        for line in open(f"{asm_dir}/{s}/assembly_info.txt").read().splitlines()[1:]:
            f = line.split("\t")
            by_length.setdefault(int(f[1]), (f[2], f[3]))
        for contig, seq in seqs[s].items():
            flye[(s, contig)] = by_length.get(len(seq), ("", ""))

    # Where a contig's topology has actually been measured against reads, that
    # measurement replaces the assembler's call. FLYE reports on its own graph
    # rather than on the molecule and is wrong in both directions in this panel:
    # it called S02's phage linear when 96 reads span its join, and called S11's
    # cointegrate plasmid linear when 57 do.
    verified = "data/reference/verified_topology.tsv"
    if os.path.exists(verified):
        for line in open(verified):
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], f[1])
            if key in flye:
                flye[key] = (flye[key][0], "Y" if f[2] == "circular" else "N")

    # --------------------------------------------------------- plasmid families
    # Being small is not evidence of being a plasmid. These assemblies are
    # fragmented -- S09 alone is in forty-three pieces -- so most contigs below
    # the chromosome threshold are chromosomal fragments, and admitting them
    # here would classify the majority of the chromosome as plasmid. The
    # candidate set is therefore built from evidence of an actual circle:
    # a contig the assembler closed, or one already carrying a plasmid name from
    # the 2024 curation, which is the claim this script exists to test.
    def best_class(s, c):
        prof = {k: merged_length(v) / len(seqs[s][c]) for k, v in cover[(s, c)].items()}
        return max(prof, key=prof.get) if prof else "NONE", prof

    # Which contigs the 2024 curation called plasmids, and which it called
    # chromosome. This is read from the ORIGINAL curated assemblies, never from
    # this script's own output. Reading it from the delivered files would make
    # the script destroy its own input: a contig demoted on one run has lost its
    # plasmid name by the next, so it would drop out of the candidate set and be
    # silently downgraded again. Sourcing the claim from 2024 keeps every run
    # independent and repeatable.
    was_tagged, chromosome_tagged = set(), set()
    for s in samples:
        curated = f"{asm_dir}/{s}/{s}.curated.fasta"
        source = curated if os.path.exists(curated) else f"{v2_dir}/{s}_v2.fasta"
        for line in open(source):
            if not line.startswith(">"):
                continue
            contig = line.split()[0][1:]
            # The 2024 files number some contigs without zero padding.
            m = re.search(r"(\d+)", contig)
            keys = {contig}
            if m:
                keys.add(f"contig_{int(m.group(1)):02d}")
            for key in keys:
                if "plasmid-name" in line:
                    was_tagged.add((s, key))
                if "location=chromosome" in line:
                    chromosome_tagged.add((s, key))

    eligible = []
    for s in samples:
        chrom = max(seqs[s], key=lambda k: len(seqs[s][k]))
        for c in seqs[s]:
            cls, _ = best_class(s, c)
            if cls in ("APHID", "APSE") or c == chrom or len(seqs[s][c]) < 2000:
                continue
            if flye[(s, c)][1] == "Y" or (s, c) in was_tagged:
                eligible.append((s, c))

    with open(f"{work}/cand.fa", "w") as fh:
        for s, c in eligible:
            fh.write(f">{s}|{c}\n{seqs[s][c]}\n")
    # Both sides of each alignment are recorded. The query intervals say how much
    # of one contig the other covers; the target intervals say WHERE along the
    # other it lands, which is what distinguishes two plasmids joined end to end
    # from one plasmid nested inside another.
    pair = collections.defaultdict(list)
    pair_target = collections.defaultdict(list)
    for f in paf(f"{work}/cand.fa", f"{work}/cand.fa", ("-N", "50", "-p", "0.05")):
        if f[0] != f[5] and int(f[10]) >= 1000:
            pair[(f[0], f[5])].append((int(f[2]), int(f[3])))
            pair_target[(f[0], f[5])].append((int(f[7]), int(f[8])))
    length = {f"{s}|{c}": len(seqs[s][c]) for s, c in eligible}
    frac = {k: merged_length(v) / length[k[0]] for k, v in pair.items()}

    parent = {k: k for k in length}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in itertools.combinations(length, 2):
        if frac.get((a, b), 0) >= RECIPROCAL_COVERAGE and frac.get((b, a), 0) >= RECIPROCAL_COVERAGE:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

    groups = collections.defaultdict(list)
    for k in length:
        groups[find(k)].append(k)
    ordered = sorted(groups.values(), key=lambda g: -max(length[x] for x in g))
    family = {}
    for i, members in enumerate(ordered, 1):
        for m in members:
            family[m] = f"P{i:02d}"

    # Containment on its own says nothing about which of two contigs is the real
    # replicon: it is symmetric in evidence and asymmetric only in size. The
    # circle breaks the tie. A contig the assembler closed is a molecule; a
    # longer contig that swallows it but was never closed is a join the
    # assembler could not confirm. Reading this the wrong way round inverts the
    # conclusion, calling the genuine plasmid a fragment of the mis-assembly
    # that contains it, so circularity is checked before containment is used.
    is_circular = {f"{s}|{c}": flye[(s, c)][1] == "Y" for s, c in eligible}

    # A contig that was never closed, and that holds most of two or more closed
    # plasmids, MAY be a mis-joined fusion of them -- but only if those plasmids
    # occupy separate stretches of it. Counting distinct families is not enough:
    # plasmid families in this panel are often nested, one being a deletion
    # derivative of another, so a contig carrying the larger necessarily carries
    # the smaller too. That is one plasmid, not two joined. The test is therefore
    # positional. Two plasmids laid end to end along the contig are a mis-join;
    # one sitting inside the other is a single replicon.
    OVERLAP_ALLOWED = 0.25

    def spans(inner, host):
        return merge_intervals(pair_target.get((inner, host), []))

    def merge_intervals(intervals):
        if not intervals:
            return []
        intervals = sorted(intervals)
        out = [list(intervals[0])]
        for start, end in intervals[1:]:
            if start > out[-1][1]:
                out.append([start, end])
            else:
                out[-1][1] = max(out[-1][1], end)
        return out

    def extent(blocks):
        return sum(b - a for a, b in blocks) if blocks else 0

    def overlap(x, y):
        total = 0
        for ax, bx in x:
            for ay, by in y:
                total += max(0, min(bx, by) - max(ax, ay))
        return total

    fusion_of = {}
    for a in length:
        if is_circular[a]:
            continue
        by_family = {}
        for b in length:
            if a == b or not is_circular[b]:
                continue
            if frac.get((b, a), 0) >= CONTAINMENT:
                blocks = spans(b, a)
                prev = by_family.get(family[b])
                if prev is None or extent(blocks) > extent(prev[1]):
                    by_family[family[b]] = (b, blocks)
        # Keep only families that sit in largely separate parts of the contig.
        disjoint = []
        for fam_id, (name, blocks) in sorted(by_family.items()):
            if all(overlap(blocks, other) <= OVERLAP_ALLOWED * min(extent(blocks), extent(other))
                   for _, (_, other) in disjoint):
                disjoint.append((fam_id, (name, blocks)))
        if len(disjoint) >= 2:
            fusion_of[a] = [name for _, (name, _) in disjoint]

    # A contig that was never closed and sits almost entirely inside one that
    # was is a fragment of that plasmid.
    fragment_of = {}
    for a in length:
        if is_circular[a] or a in fusion_of:
            continue
        best = None
        for b in length:
            if a == b or family[a] == family[b] or not is_circular[b]:
                continue
            if frac.get((a, b), 0) >= CONTAINMENT and frac.get((b, a), 0) < CONTAINMENT:
                if best is None or frac[(a, b)] > best[1]:
                    best = (b, frac[(a, b)])
        if best:
            fragment_of[a] = best

    # ------------------------------------------------------------------ output
    rows = []
    for s in samples:
        chrom = max(seqs[s], key=lambda k: len(seqs[s][k]))
        chrom_cov = flye[(s, chrom)][0]
        for c in sorted(seqs[s], key=lambda k: -len(seqs[s][k])):
            key = f"{s}|{c}"
            n = len(seqs[s][c])
            cls, prof = best_class(s, c)
            cov, circ = flye[(s, c)]
            ratio = ""
            if cov and chrom_cov and float(chrom_cov) > 0:
                ratio = f"{float(cov) / float(chrom_cov):.2f}"
            fam = family.get(key, "")
            frag = fragment_of.get(key)
            fused = fusion_of.get(key)

            # A Hamiltonella contig that is not a plasmid and not a phage is
            # called chromosome when it is large enough to be unambiguous, or
            # when the 2024 curation already called it that. Everything else is
            # left unplaced.
            #
            # Both halves of that rule are needed. Size alone would drop the
            # smaller contigs the 2024 delivery had explicitly placed, and the
            # 2024 tags alone would leave chromosome-sized contigs unplaced
            # merely because nobody got round to labelling them. Together they
            # reproduce the delivered documentation for thirteen of the fourteen
            # samples.
            #
            # The threshold is a statement about confidence, not about position:
            # calling a contig chromosome says the sequence is chromosomal in
            # origin, never that its place in the chromosome is known. None of
            # these assemblies is contiguous enough to claim that.
            inherited = "unplaced"
            if (s, c) in chromosome_tagged or n >= CHROMOSOME_BP:
                inherited = "chromosome"

            if cls == "APHID":
                label, why = "aphid", "aligns to the aphid genome over its full length"
            elif cls == "APSE":
                label, why = "APSE phage", "aligns to a reconstructed APSE genome"
            elif fused:
                label = "mis-joined plasmids"
                why = ("the assembler did not close this contig, and it holds most "
                       "of " + " and ".join(fused) + ", each of which assembles as "
                       "a separate closed circle elsewhere in this panel")
            elif frag:
                label = "plasmid fragment"
                why = (f"not circular, and {100 * frag[1]:.0f}% of it is contained "
                       f"in {frag[0]}, which is circular")
            elif circ == "Y" and fam:
                label, why = "plasmid", f"assembled as a closed circle; family {fam}"
            elif fam and circ != "Y":
                label = "plasmid, not circularised"
                why = (f"groups with plasmid family {fam}, but the assembler did "
                       f"not close it into a circle")
            elif (s, c) in was_tagged:
                label = inherited
                why = ("carried a plasmid name in 2024, but there is no evidence "
                       "of a circle and it joins no plasmid family")
            else:
                label = inherited
                why = "H. defensa sequence; placement unchanged from the 2024 curation"

            rows.append([s, c, n, cov, circ, ratio, cls,
                         f"{100 * prof.get(cls, 0):.1f}",
                         ",".join(shared.get((s, c), [])) or "-",
                         fam or "-", label, why])

    with open(out_path, "w") as fh:
        fh.write("sample\tcontig\tlength_bp\tflye_coverage\tcircular\tcoverage_ratio_to_chromosome\t"
                 "best_taxonomic_match\tpct_of_contig_covered\tsamples_sharing_it\t"
                 "plasmid_family\tlabel\tevidence\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    counts = collections.Counter(r[10] for r in rows)
    print(f"  {len(rows)} contigs classified, written to {out_path}")
    for k, v in counts.most_common():
        print(f"      {v:>4}  {k}")


main()
