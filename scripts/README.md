# Scripts

Everything needed to reproduce the genomes in this repository, from the raw
sequencing reads to the annotated assemblies.

The reads and the intermediate alignments are far too large to distribute here
and are not included; they are available from the accessions given in the
publication. The only file that needs editing is `data/config.sh`, which says
where your copy of the reads and assemblies lives. All scripts are run from the
root of the repository, not from inside this directory.

---

## Before anything else

Open `data/config.sh` and set the paths at the top. Then create the environment
and check it:

```
conda env create -f environment.yml
conda activate aphid
bash scripts/check_environment.sh
```

`check_environment.sh` reports the version of every tool and says which are
missing or differ from the versions used here. Read its output before starting a
long step.

### Versions

Two stages of this work ran under different toolchains, and the same tool
appears at two versions where that is the case. This is recorded rather than
smoothed over, because the phage assemblies were accepted against the output of
one Flye version and the genome assemblies produced by another.

| Tool | Used for | Version |
|---|---|---|
| `hifiasm` | metagenome assembly | not recorded, see below |
| `blobtools` (BlobToolKit) | contig classification and filtering | not recorded, see below |
| `minimap2` | alignment, assembly stage | not recorded, see below |
| `minimap2` | alignment, phage reconstruction and validation | 2.31-r1302 |
| `samtools` | alignment handling, assembly stage | not recorded, see below |
| `samtools` | alignment handling, phage reconstruction | 1.24 |
| `seqtk` | read subsetting, assembly stage | not recorded, see below |
| `seqtk` | read subsetting, phage reconstruction | 1.5 |
| `diamond` | protein-level taxonomy | not recorded, see below |
| `blastn` (BLAST+) | nucleotide-level taxonomy | 2.16.0 |
| `flye` | genome assembly, thirteen samples | 2.9.2 |
| `flye` | genome assembly, S07 | 2.9.3 |
| `flye` | phage assembly from partitioned reads | 2.9.6 |
| `busco` | completeness, with prodigal 2.6.3 | 5.1.2 |
| `pgap.py` | annotation | 2026-06-18.build8602 |
| `python3` | the analysis scripts | 3.8 or later |
| `matplotlib` | figures only | 3.10.6 |

The assembly-stage scripts call `module load` without pinning a version, so
those scripts record which tool was used and not which build. The versions that
could be recovered came from run logs and from the BUSCO summary headers, not
from the scripts. The ones marked "not recorded" have not been established; they
are stated as unknown rather than guessed.

Reference data: NCBI `nt` and `new_taxdump`, UniProt reference proteomes, and
the BUSCO lineage datasets `enterobacterales_odb10` (2020-03-06, 440 markers)
and `hemiptera_odb10` (2020-08-05), all retrieved in 2024.

The Python scripts use nothing outside the standard library except `matplotlib`
for figures and `openpyxl` for one spreadsheet output. If you want neither, you
do not need either.

---

## The pipeline, in order

### Stage 0 — Assembly (`scripts/assembly/`)

Reads to curated genome. These steps ran on a SLURM cluster and are reproduced
as they were run, with the original absolute paths and `module load` lines
intact; adapt them to your system before running. Each has a header saying so.

```
01_metagenome_assembly.sh        hifiasm --primary, per sample
02_backmap_reads.sh              minimap2 + samtools, for coverage
03_busco.sh                      BUSCO against enterobacterales_odb10
04_blastn_nt.sh                  blastn against NCBI nt
05_diamond_uniprot.sh            diamond blastx against UniProt
06_blobtools_create_filter.sh    build the BlobDir and apply the saved filter
07_extract_filtered_reads.sh     reads mapping to retained contigs, via seqtk
08_flye_assembly.sh              flye --pacbio-hifi, the genome assembly
09_backmap_flye.sh               coverage for the second classification
10_busco_flye.sh                 completeness for the second classification
11_blastn_flye.sh                taxonomy for the second classification
```

Contig selection in step 06 was originally made by eye in the BlobToolKit
viewer. The criteria were saved, and the fourteen saved filters are in
`data/blobtools_JSON/`, so the step replays without the viewer.

`original/` holds the protocol documents these scripts came from, for the
*Hamiltonella*, aphid and *Buchnera* tracks and for the software setup. They are
the record of what was done, including the manual curation in
`original/hamiltonella/5_Curate_genome_assemblies.md`, which was carried out by
hand and has no script. `aphid_host/` holds the host and primary-symbiont
scripts, which are included for completeness; that track is summarised in the
methods rather than reproduced here. They are `assemble_merged.sh` (the pooled
hifiasm assembly of all fourteen read sets), `repeatMod.sh` and `repeatMask.sh`
(RepeatModeler and RepeatMasker), and `funTrain.sh` and `funPred.sh`
(funannotate training on RNA-seq and gene prediction).

### Stage 1 — How many phages, and how many strains

```
bash scripts/detect_multiple_APSE_and_strains.sh
```

Screens all fourteen samples for evidence of more than one APSE phage, and
separately for more than one strain of *Hamiltonella*. Those are different
questions and the script keeps them apart. Read its header before running it: it
explains why independent tests are used and what each is blind to.

Output: `results/multiplicity_summary.tsv` and the per-contig tables beside it.

### Stage 2 — Separate and reassemble the phages of a sample that carries two

```
bash scripts/build_initial_templates.sh S07 S07_templates.fa RHS:contig_3,contig_4 CdtB:contig_5,contig_1
bash scripts/reconstruct_phages_and_integration_sites.sh S07 contig_02 S07_templates.fa
```

The first command assembles the starting templates from the sample's own
pre-curation contigs, grouped by integrase type and toxin family. Its header
explains how a contig is assigned to a group, and gives the grouping used here.

Only needed for a sample stage 1 flags, which here is S05 and S07. This is the
substantial piece of work and its header explains the method in full. In
outline: reads are sorted into one set per phage before anything is assembled,
each set is assembled from scratch, and the result is tested against the reads
that produced it. The reference used for sorting contributes no sequence to the
output.

It calls four helpers in `lib/`, each usable on its own:

| Helper | What it does |
|---|---|
| `kmer_ancestry.py` | sorts reads by k-mers unique to one phage, and identifies recombinant reads |
| `kmer_assign.py` | the same job where one lineage is itself a recombinant of the others, which defeats the unique-k-mer approach |
| `split_read_junctions.py` | finds the points where a read crosses from chromosome into phage |
| `integrate_prophage.py` | splices a phage into a chromosome at an *att* site |

### Stage 3 — Phage reference genomes

```
python3 scripts/extract_prophage_references.py <sample> <chromosome.fa> <probe.fa> <reads.fq> <outdir>
python3 scripts/test_phage_circularity.py      <name> <contig.fa> <reads.fq> <outdir>
```

The first cuts a prophage out of a chromosome, finding the boundary repeat from
the sequence itself rather than assuming a known one. The second decides whether
a phage genome really exists as a circle, by asking whether any read crosses the
point where its two ends meet.

Use the second on anything you intend to label circular. An assembler's
circularity flag describes its own graph, not the molecule, and in this panel it
was wrong in both directions.

### Stage 4 — Build the assemblies

```
bash scripts/build_all_v2_assemblies.sh
```

Splices the prophages in, drops the contigs that should not be there, corrects
the headers, classifies every contig and applies the labels. It reads
`data/assembly_manifest.tsv`, which records the per-sample arguments; those
arguments cannot be derived from the inputs, so the manifest and the driver are
needed together.

The driver calls, in order, `lib/integrate_prophage.py`, `build_v2_assemblies.py`,
`lib/classify_replicons.py` and `apply_replicon_labels.py`. The last of these is
the final step to touch a sequence file.

### Stage 5 — Annotate

```
bash scripts/annotate/run_pgap_batch.sh
```

Runs NCBI PGAP over the finished genomes on a cloud instance;
`assembly/original/setup/PGAP_setup.md` describes how that instance is built.
Note the organism name: PGAP validates it against its own taxonomy and refuses
to start if it is not the currently accepted one.

### Stage 6 — Check the result

```
bash    scripts/run_replicate_test.sh
python3 scripts/compare_replicates.py > replicate_comparison.tsv
bash    scripts/audit_blobtools_filter.sh
bash    scripts/audit_assembly_structure.sh
bash    scripts/measure_coverage.sh
bash    scripts/run_busco.sh
python3 scripts/build_assembly_stats.py assemblies_v2 <busco_dir> deliverables/assembly_stats.xlsx
python3 scripts/compare_phages_across_samples.py <phages.fa> <modules.tsv>
python3 scripts/build_contig_inventory.py results/replicon_table.tsv assemblies_v2 results/contig_inventory.txt
```

The replicate test rebuilds each phage genome from sets of reads that do not
overlap, which is the strongest check available without new sequencing. The
filter audit asks whether the taxonomic filtering applied before assembly threw
away anything belonging to *Hamiltonella* or APSE. The structure audit tests the
scaffold joins, contig adjacency and within-contig coverage steps against reads
rather than against other assemblies. The comparison script places every phage
in the panel against every other, so that two samples are not called the same
lineage merely because they share a toxin gene.

The structure audit calls three helpers in `lib/`, each of which answers one
question against reads rather than against another assembly:

| Helper | What it tests |
|---|---|
| `check_scaffold_joins.py` | whether reads span each join made during curation, at two distances |
| `count_contig_links.py` | how many reads align partly to one contig and partly to another |
| `find_depth_steps.py` | whether coverage steps within a contig suggest two elements at different copy number |

### Figures (`scripts/figures/`)

```
python3 scripts/figures/plot_S05_architecture.py     the three phage lineages of S05
python3 scripts/figures/plot_S07_architecture.py     the two phages of S07 and their att sites
python3 scripts/figures/plot_replicate_test.py       the replicate test
python3 scripts/figures/plot_APSE_dotplot.py         two phage genomes aligned to each other
```

Most of these carry their measured values in the source rather than reading them
from `results/`, so they redraw a figure rather than recompute it. If you change
an underlying number, change it here too.

---

## Which script produces which delivered file

A row with no script is a file that is maintained by hand. There are two, and
both are documentation of what was done rather than derived output.

| File | Produced by |
|---|---|
| `assemblies_v2/S*_v2.fasta` | `build_all_v2_assemblies.sh` |
| `results/phage_references/*.fasta` | `extract_prophage_references.py`, topology from `test_phage_circularity.py` |
| `results/phage_reconstructions/*.fasta` | `reconstruct_phages_and_integration_sites.sh` |
| `results/multiplicity_summary.tsv` and the per-contig tables | `detect_multiple_APSE_and_strains.sh` |
| `results/S0{3,5,7}_integration_sites.tsv` | `lib/split_read_junctions.py`, transcribed by hand (see below) |
| `results/replicon_table.tsv` | `lib/classify_replicons.py` |
| `results/contig_inventory.txt` | `build_contig_inventory.py` |
| `results/cross_sample_containment.tsv` | `compare_phages_across_samples.py` |
| `results/replicate_test.tsv` | `run_replicate_test.sh` then `compare_replicates.py` |
| `results/blobtools_filter_audit.tsv` | `audit_blobtools_filter.sh` |
| `results/scaffold_join_audit.tsv`, `contig_link_audit.tsv`, `depth_step_audit.tsv` | `audit_assembly_structure.sh` |
| `results/coverage_summary.tsv`, `coverage_by_contig.tsv` | `measure_coverage.sh` |
| `deliverables/assembly_stats.tsv` | `run_busco.sh` then `build_assembly_stats.py` |
| `results/v2_checksums.txt` | `md5sum` over the delivered files |
| `results/v2_provenance.tsv` | **hand-maintained** |
| `results/phage_lineage_names.tsv` | **hand-maintained**, from `cross_sample_containment.tsv` |

---

## Steps that were done by hand

These are recorded rather than automated, and each has a file that records what
was decided.

- **Contig selection during taxonomic filtering.** Judged in the BlobToolKit
  viewer; the criteria are saved in `data/blobtools_JSON/` and replay
  non-interactively.
- **Curation.** Scaffolding and contig removal, described in
  `assembly/original/hamiltonella/5_Curate_genome_assemblies.md`. The outcome is
  captured in `data/assembly_manifest.tsv` and in the curated assemblies
  themselves.
- **Choosing the starting references for read partitioning.** The criterion is
  that they span the diversity in the sample; the choice is a list of contig
  names and is recorded in the reconstruction script's header.
- **Transcribing *att* coordinates.** `lib/split_read_junctions.py` reports
  junctions; clustering them to a site and reading off the core is a judgement,
  and the results are written into `results/S0*_integration_sites.tsv`.
- **The S10 scaffold gap fill.** Done once from the assembly graph and described
  in `results/v2_provenance.tsv`. `build_all_v2_assemblies.sh` takes the result
  through `S10_GAPFILLED` and says so if it is not set.
- **Choosing which iteration round to keep for the recombinant S05 lineage.**
  The stop rule and the length series are in the methods.

---

## Thresholds

Thresholds are declared as named variables at the top of each script, with a
comment saying why that value and not another. If you change one, the comment is
the thing to read first: several of the values are set by a property of this
dataset rather than by convention, and the reasoning matters more than the
number.

Two are worth understanding before altering them. The identity floor in the
filter audit is set relative to how similar this panel's own chromosomes are to
a published genome, so it must move if the script is applied to another species.
And the depth targets in the assembly scripts exist because the untrimmed read
sets here run to several thousand fold, which is far more than the assembler
needs and enough to destabilise it.

---

## Verifying a run

```
md5sum -c results/v2_checksums.txt
```

Every one of the delivered files should match. The samples that need no
substitution and no contig dropped are the control: S01, S04, S06, S11, S12,
S13 and S14 must come out unchanged, and if one of them differs the fault is in
the pipeline rather than in that sample.

Then re-run stage 1 against the finished assemblies and confirm that each phage
reference now reports a single phage, and every sample a single strain.

---

## Attribution

`data/reference/toxins_CDS.fasta` was compiled by **Giacomo Moggioli** and is
included so the analysis is reproducible from this repository alone.

`APSE-1` (`NC_000935.1`) is from RefSeq. The remaining sequences in the APSE
panel were generated in this study.

The assembly-stage scripts and protocol documents in `scripts/assembly/` come
from https://github.com/edexter/H_defensa_assembly, which is their canonical
home.
