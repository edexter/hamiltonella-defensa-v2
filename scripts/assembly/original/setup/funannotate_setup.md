# FUNNANOTATE setup

FUNANNOTATE is an extremely difficult pipeline to run because it has many dependencies, some of which have unique and troublesome installation requirements. A docker container is available for FUNANNOTATE which eases this process somewhat, but most cluster environments will not allow the necessary permissions to run docker, so life is much easier if an EC2 instance is created to run the more difficult parts of the pipeline. Preliminary steps in the pipeline such as repeat masking can still be run on a cluster, which tends to be much cheaper.

# Setting up the EC2 instance

An EC instance with 16 CPUs, 64 GB of memory, and 300 GB of storage works well. **THE IP ADDRESS WILL CHANGE EVERY TIME THE INSTANCE IS RESTARTED**. 

````
# Set permission level for key
chmod 400 funannotate-server-key.pem

# Login with keypair
ssh -i funannotate-server-key.pem ec2-user@34.207.122.120

# Start with the usual update 
sudo yum update -y

# Install docker
sudo yum install -y docker

# Start docker
sudo service docker start

# Check docker version
docker --version

# Make sure docker always starts when you start this instance
sudo systemctl enable docker

# Give the ec2-user permission to use docker (requires log-out to take effect)
sudo usermod -a -G docker ec2-user
exit

# Log back in
ssh -i funannotate-server-key.pem ec2-user@18.206.254.237

# Check docker version
docker --version

# Check docker access
docker ps
````

# Install difficult external dependencies

Genemark is difficult to install because of licensing issues. The instructions are described here: http://topaz.gatech.edu/GeneMark/license_download.cgi. Don't forget to separately download, extract, and store the license key in the home directory (not the installation directory)!

````
# Send genemark key to server from local
scp -i ~/funannotate-server-key.pem  \
-r /mnt/c/Users/ericd/Downloads/gm_key.gz \
ec2-user@18.206.254.237:/home/ec2-user/

# Send genemark program to server from local
scp -i ~/funannotate-server-key.pem  \
-r /mnt/c/Users/ericd/Downloads/gmes_linux_64.tar.gz \
ec2-user@18.206.254.237:/home/ec2-user/

# Send signalp program to server from local
scp -i ~/funannotate-server-key.pem  \
-r /mnt/c/Users/ericd/Downloads/signalp-4.1g.Linux.tar.gz \
ec2-user@18.206.254.237:/home/ec2-user/

# Send repbase to local server
scp -i ~/funannotate-server-key.pem  \
-r /mnt/c/Users/ericd/Downloads/RepBaseRepeatMaskerEdition-20170127.tar.gz \
ec2-user@18.206.254.237:/home/ec2-user/

# Log back in
ssh -i funannotate-server-key.pem ec2-user@18.206.254.237

# Rename the genemark key to the name docker is expecting
mv gm_key.gz gm_key_64.gz

# Rename the signalP file to the name docker is expecting
mv signalp-4.1g.Linux.tar.gz signalp-4.1f.Linux.tar.gz

#Optional
tar -xzf gmes_linux_64_4.tar.gz
````

# Set up FUNANNOTATE docker environment

````bash
# Get newest docker image
docker pull nextgenusfs/funannotate

# Get all-in-one docker script
wget -O funannotate-docker https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker

# Make the script executable
chmod +x funannotate-docker

# Run on a test file
./funannotate-docker test -t predict --cpus 12

# Add the insecta busco database
sudo ./funannotate-docker setup -b insecta

# Import data from scicore cluster (where earlier steps of annotation were performed for a better price)
rsync -avz -e "ssh -i funannotate-server-key.pem" /scicore/home/ebertd/dexter0000/aphid/annotation/output/funannotate_train_full ec2-user@18.206.254.237:/home/ec2-user

# Move the input genome to EC2
rsync -avvz -e "ssh -i funannotate-server-key.pem" /scicore/home/ebertd/dexter0000/aphid/annotation/output/repeatmasker_output/aphid.filtered.cleaned.renamed.fa.masked ec2-user@18.206.254.237:/home/ec2-user
````

### Add the insecta BUSCO database to the docker instance

````bash
docker run -it --name funannotate_container nextgenusfs/funannotate /bin/bash
funannotate setup -b insecta
docker commit funannotate_container funannotate-with-insecta
docker ps

# Test it
docker exec -it funannotate_container funannotate test -t predict --cpus 12 --busco_db insecta
````



# Download and setup INTERPROSCAN

````bash
# Make directory for program files
mkdir my_interproscan

# Navigate into program folder
cd my_interproscan

# Download the program
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz

# Download the checksum file
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.71-102.0/interproscan-5.71-102.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful
# Must return *interproscan-5.71-102.0-64-bit.tar.gz: OK*
# If not - try downloading the file again as it may be a corrupted copy.
md5sum -c interproscan-5.71-102.0-64-bit.tar.gz.md5

# Extract the tarball
tar -pxvzf interproscan-5.71-102.0-*-bit.tar.gz

# Add path to bash profile
nano ~/.bashrc
export PATH=$PATH:/home/ec2-user/my_interproscan/interproscan-5.71-102.0
source ~/.bashrc

# Index hmm models (takes a few minutes - be patient)
python setup.py -f interproscan.properties

# Install various dependencies
sudo yum install -y java-11-amazon-corretto
sudo yum install -y perl-FindBin
sudo yum install -y perl-lib
sudo yum install -y perl-Data-Dumper

# Perform two test runs
interproscan.sh -i test_all_appl.fasta -f tsv -dp
interproscan.sh -i test_all_appl.fasta -f tsv
````

# Set up conda environment for REPEATMASKER and REPEATMODELER

````
# Create conda environment
conda create -n repeatMaskModel

conda install -c bioconda repeatmodeler
conda install -c bioconda repeatmasker
````

# Upload input files to server

This is a placeholder as a reminder that if you're not using a cluster at all and performing the entire annotation on the EC2 instance, then you need to upload some additional files to the EC2 instance.

````bash
# The genome assembly

# The RNA seq files

# The NCBI metadata file (optional)
````

# Install Eggnog mapper

This dependency is required for the FUNANNOTATE annotate module

````
conda activate eggnog-mapper

# Install some dependencies that the Amazon linx version is missing
conda install -c bioconda -c conda-forge eggnog-mapper

# Add to path
export PATH=~/miniconda3/bin:$PATH
echo 'export PATH=~/miniconda3/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

# Make database directory
mkdir -p ~/eggnog-mapper-data

# Set database directory
export EGGNOG_DATA_DIR=~/eggnog-mapper-data
echo 'export EGGNOG_DATA_DIR=~/eggnog-mapper-data' >> ~/.bashrc
source ~/.bashrc
````

