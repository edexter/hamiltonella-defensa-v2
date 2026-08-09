# Detection of multiple APSE phages and multiple *Hamiltonella* strains

This script screens each sequenced sample for two things that a conventional genome assembly
pipeline will silently hide: the presence of **more than one APSE bacteriophage**, and the
presence of **more than one *Hamiltonella defensa* strain**.

It exists because the standard curation rule used to clean these assemblies — *where one contig
aligns entirely inside another, delete the shorter one as a duplicate* — is correct for a single
APSE assembled twice (integrated and excised), but destroys the evidence when a sample carries
two different APSE phages. APSE genomes share a highly conserved backbone and differ mainly in a
single toxin-encoding module, so a second phage presents as a redundant fragment of the first.

* **Inputs**
  * `data/sampleID.txt`: index of sample IDs, one per line
  * `data/config.sh`: paths to the reads and assemblies on the local machine — **the only file
    that needs editing**
  * `data/reference/APSE_panel.fasta`: four complete APSE genomes used to decide which contigs
    are phage derived — APSE-1 (`NC_000935.1`) plus the H76, H101 and H244 phages from this study
  * `data/reference/toxins_CDS.fasta`: coding sequences of all 15 catalogued toxin genes across
    the 14 samples, spanning five families (RHS, CdtB, MAC/perforin, Shiga-like, leucine-rich
    repeat). Supplied by G. Moggioli
  * `SXX.fq`: *Hamiltonella*-filtered HiFi reads per sample, held externally
  * `FLYE2/SXX/assembly.fasta`: FLYE genome assembly per sample, held externally

* **Outputs**
  * `results/multiplicity_summary.tsv`: one row per sample carrying both verdicts, the evidence
    behind them, and an explicit statement of whether sequencing depth was adequate
  * `results/per_contig_discordance.tsv`: one row per contig — length, class, depth, discordant
    site count and rate, median minor allele fraction, and the spatial pattern of the sites
  * `results/toxin_screen.tsv`: presence of each toxin gene in the assembly and in the reads
  * `results/split_read_links.tsv`: reads whose alignment is split between two contigs, which
    locate prophage integration sites
  * `results/logs/SXX.log`: per-sample run log recording software versions and thresholds

* **Dependencies**
  * `minimap2`: read and sequence alignment (tested with 2.31-r1302)
  * `samtools`: alignment handling and pileup (tested with 1.24)
  * `awk`: all parsing and summarisation (POSIX awk; mawk and gawk both work)

## How to run

```bash
# Edit data/config.sh so that APHID_DATA points at your copy of the data, then:
bash scripts/detect_multiple_APSE_and_strains.sh
```

Runtime is roughly two minutes per sample on an ordinary laptop, so about half an hour for the
full panel of 14. Intermediate alignments are written to the scratch directory named in
`config.sh` and can be deleted afterwards; only the small summary tables are kept.

## Why two independent tests are used

Neither test alone is sufficient, and a sample can be positive by one and negative by the other.

**Test 1 — allelic discordance.** Reads are mapped back to the sample's own assembly and
positions are counted where a consistent minority of reads disagree with the assembled sequence.
Where two similar elements were collapsed into one contig, reads from both pile up on it and
disagree at every position where the two differ. This test only detects a second phage **if the
assembler collapsed it**. Two phages divergent enough to assemble separately produce no
discordance at all.

**Test 2 — toxin content.** Reads are mapped directly to a catalogue of APSE toxin coding
sequences. Two distinct toxin families in one sample means two phages, regardless of what the
assembly shows. Crucially this test is independent of any assembly or curation decision, so it
detects a phage that is present in the sample but **absent from the assembly**. It only detects a
second phage **if that phage carries a different, previously catalogued toxin**; two phages
bearing the same family are invisible to it.

## Interpreting the discordance results

A clonal element sequenced with PacBio HiFi reads should be essentially monomorphic. In the
validation sample **S12, which carries a single APSE, the phage yields zero discordant sites at
1394-fold coverage**. The assay therefore has no measurable false-positive rate at phage-level
depth, and any non-zero discordance on an APSE contig is meaningful. This is why the thresholds
below can be simple rather than statistical.

The **spatial pattern** of the discordant sites carries as much information as their number:

| Pattern | Meaning |
|---|---|
| `DISPERSED` — spread across the contig | two collapsed sequences |
| `clustered` — confined to one short block | a repeat or paralogous region |
| `clean` — no qualifying sites | a single, clonal element |

This distinction is essential rather than cosmetic. **Every *Hamiltonella* chromosome in this
panel contains one repeat block of roughly 10 kb** that generates hundreds of discordant sites.
Without the clustering step it would be reported as evidence of multiple strains in all 14
samples. In S12 that block spans 9.9 kb of a 987 kb contig, giving a dispersal fraction of 0.010
and a correct `clustered` call.

## Thresholds

Every threshold is a named variable at the top of the script, each with the reason for its value
recorded in a comment.

| Variable | Value | Reason |
|---|---|---|
| `MIN_DEPTH` | 10 | below ten reads a minor allele cannot be separated from sampling noise |
| `MIN_ALT_COUNT` | 3 | excludes isolated sequencing errors, which in HiFi data are almost always single-read events |
| `MIN_ALT_FRACTION` | 0.15 | HiFi base accuracy exceeds 99.5 %, so 15 % is far above the error floor |
| `MIN_MAPQ` | 20 | excludes reads that map equally well to two near-identical copies of the same sequence |
| `MIN_BASEQ` | 20 | standard |
| `CLUSTER_GAP` | 20000 | separates a localised repeat from variation dispersed along a contig |
| `DISPERSAL_FRACTION` | 0.50 | a contig whose sites span more than half its length is called dispersed |
| `APSE_MIN_COVERAGE` | 0.30 | fraction of a contig aligning to the phage panel before it is called phage derived |
| `SAME_PHAGE_IDENTITY` | 0.95 | two phage contigs sharing more than this are one phage in two states |
| `TOXIN_MIN_BREADTH` | 80 | percentage of a toxin gene that must be covered before it is called present |
| `MIN_POWER_DEPTH` | 30 | below this a negative result is reported as insufficient power |

## Known limitations

These are deliberate and are stated so that results can be read correctly.

**Toxin variants within one family cross-map.** The H101 and H244 CdtB genes differ by a single
nucleotide, so reads from one align acceptably to the other. The script normalises toxin names to
their family before counting, which prevents this from producing a false positive — but it also
means the toxin test cannot distinguish two phages that carry different variants of the *same*
family. Discordance (Test 1) covers that case.

**One phage assembled twice is not two phages.** Seven samples in this panel contain APSE in both
an integrated and a circular excised form. These are near-identical along their full length and
are detected by mutual alignment of the phage contigs; the sample is flagged
`review_same_phage_duplicated` rather than being called positive. This is the principal
false-positive route and is handled explicitly.

**Depth varies by two orders of magnitude across the panel.** APSE coverage ranges from 14-fold to
2159-fold. Samples below `MIN_POWER_DEPTH` report `insufficient_power` rather than a negative
result, so that absence of evidence is never recorded as evidence of absence. The `power` column
should be read before any negative verdict is relied upon.

**Discordance cannot by itself distinguish two phages in one genome from two strains each
carrying one phage.** That is precisely why the same measurement is applied to the chromosome.
Dispersed chromosomal discordance indicates a strain mixture; discordance confined to one block
is the repeat present in every sample.

## Validation

The screen was calibrated against three samples whose composition was established independently
before it was applied to the rest of the panel.

| Sample | Strain | Expected | Basis |
|---|---|---|---|
| **S12** | H244 | negative — one APSE | only CdtB detectable in its reads |
| **S07** | H101 | strongly positive | RHS present in reads, absent from the assembly |
| **S05** | H76 | positive | CdtB present as assembly fragments alongside a complete RHS phage |

The chromosome test must return a single strain for both S07 and S12, with the repeat block
correctly classified as `clustered` rather than as strain-level variation.

## Attribution

`data/reference/toxins_CDS.fasta` was compiled by **Giacomo Moggioli** and is included here so
that the analysis is reproducible from the repository alone. This repository is private; the
decision on whether to release these sequences publicly rests with him.

`APSE-1` (`NC_000935.1`) is from RefSeq. The three remaining sequences in the APSE panel were
generated in this study.
