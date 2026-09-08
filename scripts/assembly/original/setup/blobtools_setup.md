# Set up conda environment and install blobtools

Create and activate the conda environment

````bash
conda create -y -n btk -c conda-forge python=3.9
conda activate btk
````

Install blobtoolkit

````
pip install "blobtoolkit[full]"
````

Install personal copy of BLAST+ because the server installation is broken

````
# Download BLAST+
wget ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz

# Decompress BLAST+
tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz
````



# Prepare NCBI databases for BLAST

First install the NT Database. This uses FTP so sciCore won't allow the transfer on an interactive node. Need to use to use the transfer node.

````
# Need to load Perl separately
module load Perl

# Need to use a local instal of Blast+ because the cluster install is broken.
/scicore/home/ebertd/dexter0000/aphid/blobtools/blast/ncbi-blast-2.16.0+/bin/update_blastdb.pl --decompress nt [*]
````

Next install and prepare the UniProt database.

````
mkdir -p uniprot
wget -q -O uniprot/reference_proteomes.tar.gz \
 ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$(curl \
     -vs ftp.ebi.ac.uk/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/ 2>&1 | \
     awk '/tar.gz/ {print $9}')
cd uniprot
tar xf reference_proteomes.tar.gz

touch reference_proteomes.fasta.gz
find . -mindepth 2 | grep "fasta.gz" | grep -v 'DNA' | grep -v 'additional' | xargs cat >> reference_proteomes.fasta.gz

echo -e "accession\taccession.version\ttaxid\tgi" > reference_proteomes.taxid_map
zcat */*/*.idmapping.gz | grep "NCBI_TaxID" | awk '{print $1 "\t" $1 "\t" $3 "\t" 0}' >> reference_proteomes.taxid_map

diamond makedb -p 16 --in reference_proteomes.fasta.gz --taxonmap reference_proteomes.taxid_map --taxonnodes ../taxdump/nodes.dmp -d reference_proteomes.dmnd
cd -
````

Finally, fetch the NCBI Taxdump.

````

mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
````



# Add BUSCO databases

````bash
mkdir -p busco

# Hemiptera database
wget -q -O hemiptera_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/hemiptera_odb10.2020-08-05.tar.gz" \
        && tar xf hemiptera_odb10.gz -C busco

# Most specific Hamiltonela and Buchneria database
wget -q -O enterobacterales_odb10.gz "https://busco-data.ezlab.org/v4/data/lineages/enterobacterales_odb10.2020-03-06.tar.gz" \
        && tar xf enterobacterales_odb10.gz -C busco

````








