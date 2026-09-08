# Genome curation

The filtered genomes can still benefit from a degree of manual curation for the purposes of:

* Identifying duplicate contigs
* Formatting FASTA headers
* Identifying plasmids

While this process is more difficult to annotate than programmatic approaches, it be described by the following algorithms

* Genomes are self-aligned using miniMap and possible duplicate contigs are identified. These contigs are then BLASTed against each other. In cases when 1 contig is completely nested inside of another, the shorter contig is removed. This frequently occurred in the case of the APSE phage, which could exist as both the linear lysogenic and circular lytic form within a single sample
* Circular contigs were BLASTed against the NCBI database to ascertain if they were plasmids.
* The assembly graph was assembled for each sample to determine if genome misassembles required contigs to be split or scaffolded together
* Contig headers were formatted according to the NCBI guidelines at: https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/

The curated genomes were then saved under the name scheme of **S01.curated.fasta**

