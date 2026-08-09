# *Hamiltonella defensa* genome assemblies, version 2.0

Fourteen samples. For each one: a genome assembly with the APSE prophage in place wherever the
evidence supports it, and a separate reference genome for every APSE phage the sample carries.

Version 2.0 supersedes the 2024 delivery. It is an improved version of the same assemblies, not a
new assembly: no sample was resequenced and no sample was reassembled from scratch except where
stated.

---

## What is here

```
assemblies_v2/            one file per sample: chromosome, plasmids, and prophages in place
results/phage_references/ one file per APSE genome, per sample
results/v2_provenance.tsv how every sequence was produced
results/phage_lineage_names.tsv   which lineage each sample's phage belongs to
```

Eighteen phage reference files cover the fourteen samples. S05 carries four and S07 two; the rest
carry one each.

---

## Naming

Files keep the **sample** identifier, even where two samples carry an identical phage genome.
S03, S04 and S13 all carry the same 35,455 bp phage, and each ships its own copy under its own
sample name. Retaining the sample information was preferred to deduplicating it.

The **lineage** is recorded in the FASTA header and in `results/phage_lineage_names.tsv`. Lineages
are named `APSE-<toxin>-<numeral>`, and **the numeral matters**. A shared toxin does not mean a
shared lineage: this panel contains three distinct CdtB-carrying lineages and three distinct
RHS-carrying ones. Two phages sharing both toxin and integrase type can share as little as 59 % of
their 31-mers.

---

## How to read the headers

**`[topology=...]`** is a statement about the molecule, tested against reads, not a copy of the
assembler's flag. For each phage the sequence was rotated so its two ends met in the middle and the
sample's reads were mapped; only a genuine circle produces reads spanning that join. The flag was
wrong in both directions:

- **S02 and S08** were flagged linear by the assembler and are circular — 96 and 110 spanning reads.
- **S06** has a complete phage sequence whose circular form does not exist. Zero reads span the
  join, and its phage sits at only twice the chromosome depth, while S01 carries the same lineage
  at forty times. Delivered linear.

Every chromosome is delivered `[topology=linear]`. The 2024 delivery marked chromosomes circular in
all fourteen samples; FLYE had called every one of those contigs `circ=N`, and S05 and S02 labelled
three or four separate contigs as circular chromosomes each. **Chromosome circularisation is out of
scope** — no read links the chromosome ends, and closing them needs a different sequencing
technology.

---

## Prophage placement

The prophage is in place in **eleven of fourteen** samples.

| | Samples |
|---|---|
| Prophage in place | S01, S03, S04, S05, S06, S07, S10, S11, S12, S13 |
| **Not placed** | **S02, S08, S09, S14** |

For the four not placed, the phage is delivered as a separate reference and the chromosome is
unchanged. The reasons differ and are worth stating:

- **S02, S09, S14** — too few reads cross the junction to locate an integration site (3, 0 and 7).
  Limited by coverage, not by biology.
- **S08** — 65 reads place the phage against the end of one chromosome contig, but the two flanking
  contigs share no sequence at all, so there is no *att* repeat to splice at. The phage probably
  sits between two contigs; asserting that would mean claiming a scaffold join on one-sided
  evidence. Left unplaced deliberately.

---

## Sample-specific notes

**S05 carries four phage files.** Three lineages plus one allele:

| File | What it is |
|---|---|
| `S05_APSE_RHS_v2` | complete, circular |
| `S05_APSE_CdtB_hapR_v2` | complete, circular |
| `S05_APSE_CdtB_hapA_v2` | a **recombinant** of the other two — real, but **partial and not circular** |
| `S05_APSE_RHS_alt_v2` | an **alternative allele** of the RHS phage, not a separate phage |

hapA cannot be closed, and this is structural rather than a shortcoming of effort: the only thing
that distinguishes it from its two parents is a single recombination junction, so reads far from
that junction are indistinguishable from a parent and cannot be assigned to it.

`S05_APSE_RHS_alt` differs from `S05_APSE_RHS` at one 436 bp cassette, where the population carries
two unrelated sequences at roughly **58:42**. Both are supplied because mapping reads against only
the majority version loses 42 % of them at that locus.

**S10** is a draft. Only 46 reads match its phage, at about tenfold coverage. The assembly reached
39,246 bp and agrees with the previously delivered sequence exactly over their overlap, but the
terminal 5,151 bp rested on two- to fourfold coverage and was removed. What ships is the 34,095 bp
that is well supported.

**S14** is a 19,800 bp fragment at thirteenfold coverage. It is too short to assign to a lineage,
because the sequence it does contain is shared between two candidate lineages.

---

## What changed from the 2024 delivery

**Ten samples are byte-identical in sequence** — S01, S02, S04, S06, S08, S09, S10, S11, S12, S13.
Only their headers changed. That was the control: a pipeline that altered them would have been
wrong.

Four samples changed, and only where intended:

| Sample | Change |
|---|---|
| S03 | prophage moved into the chromosome; net **0 bp** |
| S05 | +15,853 bp — two prophages in, three APSE fragments out |
| S07 | +39,657 bp — two prophages in, one chimeric APSE contig out |

The single largest correction is **S07**, where a complete APSE phage had been deleted during
curation as a redundant duplicate. It was not redundant: the sample carries two different phages
whose backbones are 97–99 % identical, so one looked like a copy of the other. The reference
delivered in 2024 was a chimera of both, which is why mapping against it produced apparent
polymorphism — 6.94 heterozygous sites per kb — that exists in neither phage. Both reconstructions
now give 0.00.

---

## How far to trust these

Every phage genome was rebuilt from **non-overlapping** subsets of reads and reassembled
independently: 27 replicate assemblies in total, no read shared between any two. **Twenty-one
reproduced the delivered sequence exactly**, at 100 % coverage with zero differences.

The exceptions are informative rather than contradictory: two replicates recovered the alternative
S05 allele described above, and two are hapA, which necessarily shortens when its pile is halved.

Three further checks support the method:

- `S05_APSE_RHS_v2` is **identical to the contig delivered in 2024** over 39,286 bp, reached by a
  completely different route.
- S04 and S13, excised from their chromosomes, are **identical to S03's independently assembled
  circular contig** — 35,455 bp, three samples, two methods.
- The read filter applied in 2024, upstream of everything here, was audited across all 68,863
  pre-filter contigs. Nothing belonging to *Hamiltonella* or APSE was lost.

---

## Known limitations

1. No chromosome is closed. Out of scope; needs long-range data.
2. Four samples have their phage unplaced (above).
3. S10 and S14 have incomplete phage genomes, limited by coverage.
4. S05's hapA is partial and cannot be closed.
5. S14 cannot be assigned to a lineage.
6. Occupancy of S05's shared integration site is not determined: the two lineages that compete for
   it give 5 and 7 junction reads, and an earlier count on more reads disagreed in direction.

None of these is fixable with the present data. Each is stated rather than papered over.
