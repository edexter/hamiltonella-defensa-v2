# Curated assembly summary statistics

Now that the second round of genome assembly is completed (first with HIFASM and then with FLYE) and the genomes have been filtered a second time using BLOBTOOLS, some basic summary statistics can be calculated on the curated genome assemblies



### Calculate ASM assembly stats 

````
# Activate environment
conda activate eukaryotic_genome_assembly

# Navigate to project folder
cd aphids

# Calculate ASM stats across all filtered assemblies
while IFS= read -r SAMP; do
  echo "Processing sample ID: $SAMP"

  # Run blobtools filter with the current sample ID
  scripts/asmstats curated/"$SAMP".curated.fasta > curated/"$SAMP".stats

done < scripts/sampleID.txt
````



### Aggregate ASM assembly stats

````
RESULTS_DIR="curated"
OUTPUT_FILE="curated/merged_asm_read_stats.csv"

# Find all files ending in "asm_read_stats" and save them to a temporary file
find "$RESULTS_DIR" -name "*.stats" -print > temp_results.txt

# Add a header to the output file
echo "Filename,sum,n,mean,largest,smallest,N50,L50" > "$OUTPUT_FILE"

# Process each file to extract the required fields and append to the output file
while read -r file; do
  sum=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $2); print $2}' "$file")
  n=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $4); print $4}' "$file")
  mean=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $6); print $6}' "$file")
  largest=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $8); print $8}' "$file")
  smallest=$(awk -F '[=,]' '/sum =/ {gsub(/ /, "", $10); print $10}' "$file")
  N50=$(awk -F '[=,]' '/N50 =/ {gsub(/ /, "", $2); print $2}' "$file")
  L50=$(awk -F '[=,]' '/N50 =/ {gsub(/ /, "", $4); print $4}' "$file")

  echo "$file,$sum,$n,$mean,$largest,$smallest,$N50,$L50" >> "$OUTPUT_FILE"
done < temp_results.txt

# Clean up temporary file
rm temp_results.txt

echo "File processing complete. Output saved to $OUTPUT_FILE"
````



## BUSCO

Calculate the BUSCO scores. Only takes a minute to run.

````
# Load required module
module load BUSCO

# Navigate to directory
cd curated

# Calculate BUSCO scores
while read -r SAMP; do
	
	#Print start message
	echo "Processing sample ID: $SAMP"
	
	#Run BUSCO
    busco -i "$SAMP.curated.fasta" \
    -m genome \
    -l /scicore/home/ebertd/dexter0000/aphid/blobtools/busco/enterobacterales_odb10 \
    -c 8 \
    -o $SAMP.busco.out \
    -f --offline

done < /scicore/home/ebertd/dexter0000/aphid/scripts/sampleID.txt
````



## Aggregate BUSCO results

Aggregate the BUSCO results together. This script needs to be saved and launched.

````
#!/bin/bash

# Directory to search for "busco.out.txt" files
search_dir=$1

# Output CSV file
output_file="busco_summary.csv"

# Write the header for the CSV file
echo "Filename,Complete BUSCOs (C),Single-copy BUSCOs (S),Duplicated BUSCOs (D),Fragmented BUSCOs (F),Missing BUSCOs (M),Total BUSCO groups searched" > $output_file

# Function to extract information from the file
extract_info() {
    file=$1

    # Extract the filename (the last part of the path)
    filename=$(grep "Summarized benchmarking in BUSCO notation for file" "$file" | awk '{print $NF}' | xargs basename)

    # Extract the numbers
    complete=$(grep -P "\d+\s+Complete BUSCOs \(C\)" "$file" | awk '{print $1}')
    single=$(grep -P "\d+\s+Complete and single-copy BUSCOs \(S\)" "$file" | awk '{print $1}')
    duplicated=$(grep -P "\d+\s+Complete and duplicated BUSCOs \(D\)" "$file" | awk '{print $1}')
    fragmented=$(grep -P "\d+\s+Fragmented BUSCOs \(F\)" "$file" | awk '{print $1}')
    missing=$(grep -P "\d+\s+Missing BUSCOs \(M\)" "$file" | awk '{print $1}')
    total=$(grep -P "\d+\s+Total BUSCO groups searched" "$file" | awk '{print $1}')

    # Append the data to the CSV file
    echo "$filename,$complete,$single,$duplicated,$fragmented,$missing,$total" >> $output_file
}

# Find all "busco.out.txt" files and process each one
find "$search_dir" -type f -name "*busco.out.txt" | while read -r file; do
    extract_info "$file"
done

echo "Results saved to $output_file"
````



This is how the BUSCO aggregate script is launched.

````
chmod +x scripts/busco_merge.sh
scripts/busco_merge.sh curated
````

