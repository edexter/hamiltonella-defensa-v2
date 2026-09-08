# PGAP setup overview

Annotation of the Hamiltonella the Buchnera genomes will be performed using the NCBI Prokaryote Genome Annotation Pipeline (PGAP). NCBI provides this as an optional service when submitting genomes, but it's much easier and faster to run the pipeline locally because NCBI is extremely slow. The pipeline is based on DOCKER which doesn't work well on shared computer cluster environments. So the first step will be to set up a dedicated Amazon EC2 instance to PGAP.



## Step 1: Set up Amazon EC2 instance

First, set up the EC2 instance from the AWS developer console: 8 CPUs / 32 GB memory / 100 GB storage works well (for example: m5.2xlarge). Create the access key file and place it into the local home directory of the machine that you will use to connect to the instance. Note that all of the commands shown here are for "Amazon 2023 Linux", and are quite different from the syntax that one finds on the internet for the older "Amazon 2 Linux".



## Step 2: Establish connection to EC2 instance

On your local machine, assuming the key is named "pgap-key.pem" execute the following commands. **Note that the IP address  will change every time that the instance is restarted!**

````bash
# Navigate to directory where key is stored
cd ~

# Set permission level for key
chmod 400 pgap-key.pem

# Login with keypair
ssh -i pgap-key.pem ec2-user@3.80.203.237
````



## Step 3: Set up PGAP pipeline on EC2 instance

SSH into the EC2 instance and then:

````bash
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
ssh -i pgap-key.pem ec2-user@3.80.203.237

# Check docker version
docker --version

# Check docker access
docker ps

# Download the PGAP docker image
docker pull ncbi/pgap-utils:2024-07-18.build7555

# Confirm image is available
docker images

# Install curl (Amazon 2023 OS specific)
sudo dnf swap curl-minimal curl-full

# Get the rest of the PGAP install that is not included with the docker image
curl -OL https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py

# Make the script executable
chmod +x pgap.py

# Run the update script to get the rest of the files (takes ~5 minutes)
./pgap.py --update

# Run the test data that comes with PGAP
./pgap.py -r -o mg37_results -g $HOME/.pgap/test_genomes/MG37/ASM2732v1.annotation.nucleotide.1.fasta -s 'Mycoplasmoides genitalium'
````



