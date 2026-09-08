# *Hamiltonella defensa* genome assemblies

Fourteen *Candidatus* Hamiltonella defensa genomes, one per aphid line, each with the APSE
prophage placed where the evidence supports it, plus a standalone reference genome for every APSE
phage in the panel and the annotation of all fourteen.

**Start with [`deliverables/README.md`](deliverables/README.md).** It says what was produced, how
to read the files, and what the limitations are.

| Where | What |
|---|---|
| [`deliverables/`](deliverables/) | the written documents: the delivery note, the methods, the per-sample summaries and the assembly statistics |
| [`assemblies_v2/`](assemblies_v2/) | the fourteen genome assemblies |
| [`results/`](results/) | the phage genomes and every table produced by the analysis |
| [`scripts/`](scripts/) | everything needed to reproduce the above from the raw reads; see [`scripts/README.md`](scripts/README.md) |
| [`data/`](data/) | reference sequences, the sample and barcode maps, the saved contamination filters, and the one file you must edit to run anything |
| [`docs/`](docs/) | documentation for the screening tool, and files received from collaborators |

The sequencing reads and the intermediate alignments are far too large to distribute here and are
not included. Set the paths in `data/config.sh` to point at your own copy.
