# Species-specific genome filtering

Now we use BLOBTOOLS 2 to interactively filter out any contigs which we can identify as non-aphid in origin (for the aphid assembly), or non-Buchnera in origin (for the Buchnera assembly).

# Step 1. Interactively filter the assemblies

````bash
# Navigate to the folder when the BLOBTOOLS directories are stored
cd /scicore/home/ebertd/dexter0000/aphid/viewer

#Launch the remote viewer
blobtools view --remote .
````

Open local viewer from a new terminal on your local machine

````bash
ssh -L 8001:127.0.0.1:8001 -L 8000:127.0.0.1:8000 dexter0000@login.scicore.unibas.ch
````

Navigate to page in a web browser and use the interactive filtering tool. Filter as needed and save the filter parameters as JSON files using the lists Menu.

````bash
http://localhost:8001/view/all
````



# Step 2: Create filtered assemblies

This only takes a moment and can be performed interactively

````
# Request resources for interactive node
srun --nodes=1 --cpus-per-task=8 --mem=32G --pty bash

# Activate environment
conda activate btk

# Navigate to project folder
cd aphid

# Run blobtools filter to extract Buchnera contigs
blobtools filter \
     --json blobtools/JSON/buchnera_blobtools_filter.json \
     --fasta assembly_merged/merged.p_ctg.fa \
     --output blobtools/filtered/buchnera \
     viewer/merged

# Run blobtools filter to extract Aphid contigs
blobtools filter \
     --json blobtools/JSON/aphid_blobtools_filter.json \
     --fasta assembly_merged/merged.p_ctg.fa \
     --output blobtools/filtered/aphid \
     viewer/merged
````



# Step 3: Append species names to FASTA headers

For the NCBI submission we need to append some metadata to the FASTA headers. We will do this manually for the Buchnera genome because it contains few contigs, but we need to do it programmatically for the Aphid genome.

````
sed '/^>/ s/$/ [organism=Aphae fabae] [strain=407]/' aphid.filtered.fa > aphid.filtered_with_info.fa
````

