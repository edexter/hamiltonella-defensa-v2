# Prepare HIFI sequencing reads

The raw reads for each of the 14 samples come from the sequencing center in the form of a tarball. This needs to be decompressed. Store the absolute paths of the decompressed FASTQ files in an index file named **fastq_index.txt** and the sample names in an index file named **sampleID.txt**. These index files will be need for many different steps of the pipeline

```bash
# Request computing resources
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Navidate to data directory
cd data_round2

#Decompress the data tarball from the sequencing center
tar -xvf p30633_o35228_multiplexed_1_A01.dmxData.tar.gz
```

























