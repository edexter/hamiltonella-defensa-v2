# Scripts

Everything here can be re-run. Nothing depends on a path inside our working
environment: the only file that needs editing is `data/config.sh`, which says
where your copy of the reads and assemblies lives. All scripts are run from the
root of the repository, not from inside this directory.

The sequencing reads and the original assemblies are far too large to distribute
with a repository and are not included. They are available from the authors and
from the accessions listed in the publication.

---

## Before anything else

Open `data/config.sh` and set the two paths at the top of it. Then check that the
tools below are on your `PATH`. Every version given is the one this project was
run with; nearby versions should behave the same way.

| Tool | Used for | Version tested |
|---|---|---|
| `minimap2` | all sequence and read alignment | 2.31-r1302 |
| `samtools` | alignment handling and pileups | 1.24 |
| `seqtk` | subsetting and subsampling reads | 1.5 |
| `flye` | de novo assembly | 2.9.6 |
| `python3` | the analysis and figure scripts | 3.8 or later |
| `matplotlib` | figures only | 3.10.6 |

The Python scripts use nothing outside the standard library except `matplotlib`,
and that only for drawing figures. If you do not want the figures you do not need
it.

---

## The order things run in

The pipeline divides into four stages. Each depends only on the ones above it, so
you can stop after any of them and still have something coherent.

### Stage 1 — Ask how many phages each sample carries

```
bash scripts/detect_multiple_APSE_and_strains.sh
```

Screens all fourteen samples for evidence that a sample carries more than one
APSE phage, and separately for evidence that it carries more than one strain of
*Hamiltonella*. Those are different questions and the script keeps them apart.
Read the header of the script before running it; it explains why two independent
tests are used and what each of them is blind to.

Output: `results/multiplicity_summary.tsv` and the per-contig tables beside it.

### Stage 2 — Reconstruct the phages of a sample that carries two

```
bash scripts/reconstruct_phages_and_integration_sites.sh S07 contig_02 <templates.fasta>
```

Only needed for a sample that stage 1 flags. This is the substantial piece of
work and its header explains the method in full. In outline: the reads are sorted
into one pile per phage before anything is assembled, each pile is assembled from
scratch, and the result is tested against the reads that produced it. The
reference used for sorting contributes no sequence to the output.

The script calls four helpers in `scripts/lib/`, which are documented
individually and can be used on their own:

| Helper | What it does |
|---|---|
| `kmer_ancestry.py` | sorts reads by k-mers unique to one phage, and identifies recombinant reads |
| `kmer_assign.py` | the same job where one lineage is itself a recombinant of the others, which defeats the unique-k-mer approach |
| `split_read_junctions.py` | finds the points where a read crosses from chromosome into phage |
| `read_variant_matrix.py` | counts positions where reads disagree with an assembly |
| `integrate_prophage.py` | splices a phage into a chromosome at an *att* site |

### Stage 3 — Build the delivered genomes

```
python3 scripts/extract_prophage_references.py <sample> <chromosome.fa> <probe.fa> <reads.fq> <outdir>
python3 scripts/test_phage_circularity.py      <name> <contig.fa> <reads.fq> <outdir>
python3 scripts/build_v2_assemblies.py         <sample> <curated.fa> <assembly_info.txt> <strain> <out.fa>
```

The first cuts a prophage out of a chromosome, finding the boundary repeat from
the sequence itself rather than assuming a known one. The second decides whether
a phage genome really exists as a circle, by asking whether any read crosses the
point where its two ends meet. The third writes the delivered assembly with
corrected headers.

Use the second one on anything you intend to label circular. An assembler's
circularity flag describes its own graph, not the molecule, and in this panel it
was wrong in both directions.

### Stage 4 — Check the result

```
bash   scripts/run_replicate_test.sh
python3 scripts/compare_replicates.py > replicate_comparison.tsv
bash   scripts/audit_blobtools_filter.sh
python3 scripts/compare_phages_across_samples.py <phages.fa> <modules.tsv>
```

The replicate test rebuilds each genome from sets of reads that do not overlap,
which is the strongest check available without new sequencing. The filter audit
asks whether the read filtering done in 2024, upstream of everything else, threw
away anything belonging to *Hamiltonella* or APSE. The comparison script places
every phage in the panel against every other, so that two samples are not called
the same lineage merely because they share a toxin gene.

---

## Figures

Each of these regenerates one figure from measurements already in `results/`, and
each is named after what it draws. They take no arguments and are safe to re-run.

```
python3 scripts/plot_S07_method.py            the method, and what each round changed
python3 scripts/plot_S07_rounds.py            one page per round, with its diagnostics
python3 scripts/plot_S07_one_round.py         a single round in detail
python3 scripts/plot_S05_architecture.py      the three phage lineages of S05
python3 scripts/plot_S07_architecture.py      the two phages of S07 and their att sites
python3 scripts/plot_replicate_test.py        the replicate test
python3 scripts/plot_S07_contig_coverage.py   read depth against the coverage Flye reported
python3 scripts/plot_APSE_dotplot.py          two phage genomes aligned to each other
python3 scripts/plot_method_iteration.py      the iterative loop, and when to stop
```

---

## A note on modifying these

Thresholds are declared as named variables at the top of each script, with a
comment saying why that value and not another. If you change one, the comment is
the thing to read first: several of the values are set by a property of this
dataset rather than by convention, and the reasoning matters more than the
number.

Two in particular are worth understanding before altering them. The identity
floor in the filter audit is set relative to how similar this panel's own
chromosomes are to a published genome, so it must move if you apply the script to
a different species. And the depth targets in the assembly scripts exist because
the untrimmed read piles here run to several thousand fold, which is far more
than the assembler needs and enough to destabilise it.
