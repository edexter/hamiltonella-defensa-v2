# Version 2.0 — per-sample summaries

Draft text for the version 2.0 half of each sample entry. Coverage and BUSCO are
left as placeholders and will be filled from `results/coverage_summary.tsv` once
the batch measurement finishes; every other number here is measured.

Coverage is quoted with the prophage separated from the chromosome. A single
figure would be meaningless for these samples: the APSE is actively replicating
in most of them and sits at tens to hundreds of times the chromosomal depth, so
an assembly-wide mean describes neither.

---

## S01

**Coverage:** 32X (*Hamiltonella*), 1313X (APSE)
**BUSCO:** 93.9% (unchanged)
**APSE:** 1 genotype (complete + circularized), 39,092 bp, lineage APSE-RHS-I

The version 2.0 assembly is identical in sequence to version 1.0 — the same seven
contigs and the same 2,290,146 bp, base for base. S01 was one of the samples used
as a control on the reprocessing: it was already correct, and any change to it
would have indicated a fault in the new pipeline rather than a correction to the
old assembly. What changed is the annotation, and the amount of evidence standing
behind it.

Three claims that version 1.0 asserted are now tested rather than assumed. The
chromosome was labelled circular; FLYE had in fact called that contig `circ=N`,
no read links the two ends, and it is now labelled linear. Version 1.0 assumed a
single APSE per sample; all fourteen samples have now been screened for a second
phage genotype and a second *Hamiltonella* strain, and S01 returns one of each,
with zero heterozygous sites per kb at coverage sufficient to have detected a
second genotype had one been present. The circularity of the phage is now
established by mapping reads across the join rather than taken from the
assembler's flag, which proved wrong in both directions elsewhere in this panel.

The six plasmids are confirmed as independent closed circles and keep their
original names pH15_1 to pH15_6.

---

## S02

**Coverage:** 15X (*Hamiltonella*), 168X (APSE)
**BUSCO:** 81.6%
**APSE:** 1 genotype (complete + circularized), lineage APSE-LeucineRich-I

Thirteen contigs, one fewer than version 1.0. The contig that was removed,
`contig_06`, is 12,126 bp of aphid DNA that survived the 2024 contaminant filter.
It aligns to the aphid genome assembled by this same project over its entire
length at 98% identity, and to no *Hamiltonella* sequence at all. Apart from that
removal the assembly is unchanged in sequence.

The APSE genome was delivered in version 1.0 as a linear contig. It is in fact
circular: when the sequence is rotated so that its two ends meet in the middle,
96 reads span the join. A linear phage genome and a circular one are different
molecules, and the assembler's flag was simply wrong here.

The prophage has not been placed in the chromosome. Only three reads cross the
junction between phage and chromosome, which is too few to locate an integration
site. This is a limit of coverage rather than of biology, and the phage is
delivered as a separate reference instead.

Coverage is not sufficient to test whether more than one *Hamiltonella* strain is
present, so that question is open rather than answered in the negative.

---

## S03

**Coverage:** 8X (*Hamiltonella*), 379X (APSE)
**BUSCO:** 75.2%
**APSE:** 1 genotype (complete + circularized), 35,455 bp, lineage APSE-CdtB-I

Seventeen contigs, down from twenty-two. Four of the removed contigs — 43,204 bp
in total — are aphid DNA rather than *Hamiltonella*, confirmed both by alignment
to the aphid genome and, independently, by tetranucleotide composition, which
separates them from every *Hamiltonella* contig in the sample without reference
to any database. The fifth change is the prophage, which has been moved into the
chromosome at its integration site; because the phage was already present in the
assembly as a separate contig, this is a rearrangement rather than an addition
and the total assembly length is unchanged.

The APSE genome is circular, with 217 reads spanning the join. It is identical to
the phage genomes carried by S04 and S13 — 35,455 bp in all three — which were
reached by a different route, being excised from those chromosomes rather than
assembled standalone. Three samples and two methods agreeing exactly is a useful
check on both.

Coverage of the chromosome is low at roughly 8-fold, so this sample cannot be
tested for the presence of a second *Hamiltonella* strain.

---

## S04

**Coverage:** 63X (*Hamiltonella*), 2162X (APSE)
**BUSCO:** 93.4%
**APSE:** 1 genotype (complete + circularized), 35,455 bp, lineage APSE-CdtB-I

Unchanged in sequence from version 1.0: two contigs, 2,151,406 bp. Like S01, this
sample served as a control on the reprocessing.

The APSE genome was excised from the chromosome at an integration site located
from the sequence itself rather than from an assumed boundary, and it is
identical to the phage of S03 and S13. The chromosome topology annotation is
corrected from circular to linear for the reason given throughout.

The screen finds a single phage genotype and a single *Hamiltonella* strain, at
coverage adequate to have detected a second of either.

---

## S05

**Coverage:** 15X (*Hamiltonella*), 577X (APSE)
**BUSCO:** 90.2%
**APSE:** 3 genotypes — two complete and circularized, one partial

This is the most substantially changed assembly in the panel, and the one where
version 1.0 was most incomplete. Twelve contigs, down from nineteen: three
contigs totalling 25,957 bp were aphid DNA and have been removed, and three
fragmentary APSE contigs have been replaced by two complete phage genomes spliced
into the chromosome at their respective integration sites. The assembly gains
15,853 bp of *Hamiltonella* sequence overall.

Version 1.0 recorded one APSE. The sample in fact carries three distinct phage
lineages, which the original assembly had collapsed together:

- **APSE-RHS-II**, 39,289 bp, complete and circular.
- **APSE-CdtB-II (hapR)**, 38,980 bp, complete and circular.
- **hapA**, 36,248 bp, a recombinant of the other two — genuine, but partial and
  not circularized.

The recombinant cannot be closed, and the reason is structural rather than a
shortcoming of effort. The only feature distinguishing it from its two parents is
a single recombination junction, so reads drawn from anywhere else in its genome
are indistinguishable from parental reads and cannot be assigned to it. Supplying
it partial was preferred to omitting it, since it is a real molecule.

A fourth file is supplied for completeness: `S05_APSE_RHS_alt`. It is not a
fourth phage but a second allele of the RHS phage, differing at one 436 bp
cassette where the population carries two unrelated sequences in roughly 58:42
proportion. Both are provided because mapping reads against only the majority
version discards 42% of them at that locus.

---

## S06

**Coverage:** 97X (*Hamiltonella*), 213X (APSE)
**BUSCO:** 93.9%
**APSE:** 1 genotype (complete, **not** circularized)

Unchanged in sequence from version 1.0: seven contigs, 2,294,161 bp.

S06 and S01 are the same strain of *Hamiltonella defensa*. Their chromosomes are
99.988% identical over 1.97 Mb, they share five of six plasmids, and they carry
the same APSE lineage. This was not known in version 1.0 and is worth recording,
since the two samples are not independent observations of the symbiont.

The phage is nonetheless delivered differently from S01's. Its sequence is
complete, but the circular form of the molecule does not appear to exist in this
sample: no read spans the join when the sequence is rotated, whereas 96 do in
S02 and 110 in S08. The phage also sits at only about twice chromosomal depth
here, against roughly forty times in S01, which carries the same lineage. Taken
together this indicates the prophage is integrated and quiescent rather than
replicating, so it is delivered linear. Version 1.0 would have labelled it
circular on the strength of the assembler's flag alone.

---

## S07

**Coverage:** 38X (*Hamiltonella*), 1785X (APSE)
**BUSCO:** 93.6%
**APSE:** 2 genotypes (both complete + circularized)

Four contigs, and 39,657 bp longer than version 1.0. This is the single largest
correction in the panel.

S07 carries two different APSE phages, **APSE-RHS-II** and **APSE-CdtB-II**,
whose backbones are 97–99% identical. During the 2024 curation one of them was
identified as a redundant duplicate of the other and deleted. It was not
redundant. The consequence was that the reference delivered in version 1.0 was a
chimera of the two, assembled from reads belonging to both, and matching neither.

The symptom was visible in the original data without being recognised: mapping
the sample's reads against that reference produced 6.94 heterozygous sites per
kb, spread across the whole element rather than clustered. Apparent polymorphism
of that kind, dispersed rather than localised, is the signature of two similar
sequences being forced onto one template. Both reconstructed genomes now give
0.00 heterozygous sites per kb against their own reads.

Both phages are complete and circular, and both integration sites have been
located to base resolution, each in a tRNA gene. The two phages carry different
integrase types, which is what directs them to different chromosomal sites — the
two integrases share no 21-mer sequence at all.

---

## S08

**Coverage:** 6X (*Hamiltonella*), 174X (APSE)
**BUSCO:** 85.0%
**APSE:** 1 genotype (complete + circularized), lineage APSE-MAC-I

Thirty contigs, one fewer than version 1.0. The removed contig, `contig_21`, is
11,121 bp of aphid DNA. Otherwise unchanged in sequence.

The APSE was delivered linear in version 1.0 and is circular: 110 reads span the
join.

The prophage has deliberately **not** been placed in the chromosome, and the
reason differs from the other unplaced cases. Here there is no shortage of
evidence that the phage is integrated — 65 reads place it against the end of one
chromosomal contig. But the two contigs flanking it share no sequence at all, so
there is no *att* repeat to splice at. The phage almost certainly sits between
two contigs that the assembly has not joined, and asserting the join would mean
claiming a scaffold on one-sided evidence. It is left unplaced on purpose.

---

## S09

**Coverage:** 4X (*Hamiltonella*), 168X (APSE)
**BUSCO:** 64.1%
**APSE:** 1 genotype (complete + circularized), lineage APSE-Shiga-I

Forty-two contigs, one fewer than version 1.0 after removal of `contig_36`, 8,850
bp of aphid DNA. Otherwise unchanged in sequence.

This is the most fragmented assembly in the panel and the shallowest, at roughly
4-fold chromosomal coverage. That limits what can be concluded rather than what
was assembled: the phage is complete and circular, with 96 reads spanning the
join, but no read at all crosses the junction between phage and chromosome, so
the integration site cannot be located and the prophage is left unplaced.

Coverage is far too low to test for a second *Hamiltonella* strain.

The phage carries a Shiga-like toxin, as does the published APSE-1 genome, but
the two share only 62% of their 31-mers. They are not the same lineage, and this
is the clearest illustration in the panel of why a shared toxin gene is not
evidence of a shared phage lineage.

---

## S10

**Coverage:** 15X (*Hamiltonella*), 19X (APSE)
**BUSCO:** 94.1%
**APSE:** 1 genotype (**draft**, incomplete), 34,095 bp, lineage APSE-MAC-II

Unchanged in sequence from version 1.0: five contigs.

The phage genome is a draft and is labelled as such. Only 46 reads match it, at
about tenfold coverage. The assembly reached 39,246 bp and agrees exactly with
the sequence delivered in version 1.0 over the region they share, but its
terminal 5,151 bp rested on two- to fourfold coverage and has been removed. What
is delivered is the 34,095 bp that is well supported. A shorter genome that is
right was preferred to a longer one that is partly guesswork.

Coverage is insufficient to establish whether a second phage genotype or a second
*Hamiltonella* strain is present; both questions are open.

One contig previously delivered as a plasmid, `contig_10`, is a plasmid that the
assembler did not close into a circle. The same plasmid closes cleanly in S14, so
its existence is not in doubt, but it is now described as an unclosed fragment
rather than as a complete circular replicon.

---

## S11

**Coverage:** 53X (*Hamiltonella*), 51X (APSE)
**BUSCO:** 93.4%
**APSE:** 1 genotype (complete + circularized), 40,511 bp, lineage APSE-RHS-III

Unchanged in sequence from version 1.0: six contigs. The corrections here are
entirely in what the contigs are said to be, and they are substantial.

Two of the five contigs delivered as separate plasmids in version 1.0 are not
separate plasmids, and a third is not what it first appeared to be:

- **`contig_02`** (136,360 bp) carries sequence that assembles as two separate
  closed circles in S01 and S06, meeting at position ~78,178 with no overlap
  between them, which initially suggested that the assembler had joined two
  plasmids in error. Reads say otherwise. The internal junction is spanned by
  107 reads at 150 bp either side, 80 of them by more than a kilobase, with only
  4 of 24 nearby reads soft-clipped and no depth step across it (100.9x, 114.6x
  and 105.1x before, across and after). The contig itself closes into a circle,
  with 57 reads spanning the join. It is therefore a single circular molecule in
  this strain, a **cointegrate** of two plasmids that are separate elsewhere in
  the panel, and it is delivered as one plasmid rather than as an assembly
  error. The recorded topology is in `data/reference/verified_topology.tsv`.
- **`contig_04`** (28,430 bp) and **`contig_01`** (2,960 bp) are fragments of a
  third plasmid, 86% and 100% contained respectively within a contig that closes
  as a circle in S01 and S06.

`contig_04` and `contig_01` were not closed by the assembler and do not close
against reads, so delivering them under their own plasmid names asserted two
replicons that have not been observed. They are now described as fragments.

The `contig_02` case is worth stating plainly because it was nearly got wrong in
the other direction. Comparing assemblies to each other made it look like an
assembly error, and only the reads settled it. Containment and position can
establish that two elements sit in one contig; they cannot distinguish a
mis-join from a genuine cointegrate, and a cointegrated plasmid pair is a
described phenomenon in this system.

A separate note on `contig_01`: it matches *Arsenophonus* at about 92% identity,
which at first appeared to indicate a second symbiont. It does not. It is the
portion of a *Hamiltonella* plasmid that resembles an *Arsenophonus* mobile
element — a single element inside a plasmid, not a second organism. No
*Arsenophonus* genome is present in this sample.

---

## S12

**Coverage:** 64X (*Hamiltonella*), 1922X (APSE)
**BUSCO:** 93.6%
**APSE:** 1 genotype (complete + circularized), lineage APSE-CdtB-III

Unchanged in sequence from version 1.0: four contigs.

This sample served as the benchmark against which the other assemblies were
judged. Its phage returns 0.00 heterozygous sites per kb at 1,382-fold coverage,
which establishes what a single, correctly assembled phage genome looks like when
its own reads are mapped back to it. That figure is what made the 6.94 per kb
seen in S07's version 1.0 reference interpretable as a defect rather than as
biological variation.

Both the phage and the *Hamiltonella* strain are single, at coverage adequate to
have detected a second of either.

---

## S13

**Coverage:** 22X (*Hamiltonella*), 1445X (APSE)
**BUSCO:** 93.4%
**APSE:** 1 genotype (complete + circularized), 35,455 bp, lineage APSE-CdtB-I

Unchanged in sequence from version 1.0: two contigs.

The phage is identical to those of S03 and S04. All three are 35,455 bp, and the
agreement is exact despite S03's having been assembled as a standalone contig
while S04's and S13's were excised from their chromosomes.

Chromosomal coverage is not sufficient to test for a second *Hamiltonella*
strain, so that question is left open.

---

## S14

**Coverage:** 19X (*Hamiltonella*), 13X (APSE)
**BUSCO:** 86.4%
**APSE:** 1 genotype (**partial**, 19,800 bp, lineage unresolved)

Unchanged in sequence from version 1.0: eight contigs.

The phage genome is a 19,800 bp fragment recovered at about 15-fold coverage, and
it cannot be assigned to a lineage. This is not simply a matter of having too
little sequence. The sequence that was recovered lies in the region shared
between two candidate lineages — it matches APSE-MAC-I at 0.979 and APSE-MAC-II
at 0.998 — so the fragment is nearly identical to both. A partial genome cannot
be assigned to a lineage when the candidates share the backbone it happens to
cover, and naming it would imply a discrimination the data does not support.

The prophage is not placed in the chromosome; seven reads cross the junction,
which is too few to locate an integration site.

One contig previously delivered as a plasmid, `contig_04`, was not closed into a
circle by the assembler and is described as an unclosed plasmid rather than a
complete one.
