#!/bin/bash

# Check if the input argument is provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_bam_file> <output_tsv_folder>"
    exit 1
fi

# Assign the input file and unblocked flag to variables
input_bam="$1"
output_tsv_folder="$2"

base=${input_bam##*/}

# Check if the input file exists
if [ ! -f "$input_bam" ]; then
    echo "Error: Input file '$input_bam' not found."
    exit 1
fi

# Remove 
# List of chromosomes including MT, X, and Y
chromosomes=( $(seq -f "chr%.0f" 1 22) "chrMT" "chrX" "chrY" )

for chr in "${chromosomes[@]}"; do
  echo "Extracting ${chr} to table..."
  # Define the output file name
  output_tsv="$output_tsv_folder/${base%.bam}_$chr.tsv"
  ./modkit extract --threads 16 --force $input_bam $output_tsv --region $chr --mapped --ignore 'h'
done
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract table from BAM file."
    exit 1
fi

echo "Extraction completed successfully."