# PGAP Pipeline execution

This requires setup on a dedicated server using Amazon EC2 - see the associated PGAP setup script. Once the PGAP pipeline has been installed and validated on the dedicated EC2 server, the genome assemblies can be annotated. **Don't forget that the EC2 server IP address will change every time that the instance is restarted!**



# Step 1: Upload genomes

Upload the genome assembly FASTA files to the EC2 instance from wherever they happen to be stored.

````bash
scp -i ~/pgap-key.pem \
-r /mnt/c/Users/ericd/Downloads/genomes \
ec2-user@3.80.203.237:/home/ec2-user/
````



# Step 2: Annotate genomes

The pipeline takes approximately 1 hour per genome. You could write a script to loop through all the genomes, but writing it out for each genome like shown makes it simpler to re-run individual genomes. Don't forget to launch jobs using SCREEN in the linux terminal so an accidental disconnect doesn't interrupt the job.

````bash
# Buchnera
./pgap.py -r \
-o buchnera_results \
-g $HOME/genomes/buchnera.curated.fa \
-s 'Buchnera aphidicola'

#S01
./pgap.py -r \
-o S01_results \
-g $HOME/genomes/S01.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S02
./pgap.py -r \
-o S02_results \
-g $HOME/genomes/S02.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S03
./pgap.py -r \
-o S03_results \
-g $HOME/genomes/S03.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S04
./pgap.py -r \
-o S04_results \
-g $HOME/genomes/S04.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S05
./pgap.py -r \
-o S05_results \
-g $HOME/genomes/S05.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S06
./pgap.py -r \
-o S06_results \
-g $HOME/genomes/S06.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S07
./pgap.py -r \
-o S07_results \
-g $HOME/genomes/S07.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S08
./pgap.py -r \
-o S08_results \
-g $HOME/genomes/S08.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S09
./pgap.py -r \
-o S09_results \
-g $HOME/genomes/S09.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S10
./pgap.py -r \
-o S10_results \
-g $HOME/genomes/S10.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S11
./pgap.py -r \
-o S11_results \
-g $HOME/genomes/S11.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S12
./pgap.py -r \
-o S12_results \
-g $HOME/genomes/S12.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S13
./pgap.py -r \
-o S13_results \
-g $HOME/genomes/S13.curated.fasta \
-s 'Candidatus Hamiltonella defensa'

#S14
./pgap.py -r \
-o S14_results \
-g $HOME/genomes/S14.curated.fasta \
-s 'Candidatus Hamiltonella defensa'
````



# Step 3: Download results

The EC2 server isn't meant to be used for permanent storage, so download the results after the job is finished.

````bash
# Upload all of the Hamiltonella assemblies
rsync -av -e "ssh -i ~/pgap-key.pem" ec2-user@3.80.203.237:/home/ec2-user/S* /mnt/c/Users/ericd/Downloads/results/

rsync -av -e "ssh -i ~/pgap-key.pem" ec2-user@3.80.203.237:/home/ec2-user/buchnera* /mnt/c/Users/ericd/Downloads/results/
````

