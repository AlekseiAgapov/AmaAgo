#!/usr/bin/bash

# This script preprocesses and aligns sequencing reads to a reference genome.
# It performs the following steps:
# 1. Takes two optional arguments:
#    -p: number of threads to use (default: number of cores on the machine)
#    -d: path to the working directory containing a single fq.gz file with reads
# 2. Validates input arguments and checks for the existence of the working directory and input file
# 3. Generates a FastQC report for the raw data
# 4. Trims adapters and filters out reads shorter than 15 nt using Cutadapt
# 5. Generates a FastQC report for the processed data
# 6. Aligns reads to the reference genome using Bowtie
# 7. Converts the SAM alignment file to a BAM file and removes the SAM file
# 8. Extracts relevant columns from the aligned BAM file and saves them as a TSV file

# Set default values for arguments
num_threads=$(nproc)
working_directory="non_existent_default_directory"

# Parse input arguments
while getopts p:d: flag
do
    case "${flag}" in
        p) num_threads=${OPTARG};;
        d) working_directory=${OPTARG};;
    esac
done

# Validate the number of threads
if [ $num_threads -gt 0 ]; then
	echo "$num_threads threads"
else 
	echo "Indicate the number of threads! -p argument should be > 0."
	exit
fi

# Validate the working directory
if [ -d $working_directory ]; then
	echo "Working directory exists, starting the calculations."
else 
	echo "Working directory does not exist! Check -d argument."
	exit
fi

# Change to the working directory
cd $working_directory

file_name=trimmed.fastq.gz

# Check if the output file already exists
if [ -f $file_name ]; then
	echo "trimmed.fastq file already exists. The script must have already been run. Exiting."
	exit
else 
	echo "Starting the FASTQ file preprocessing."
fi

# Generate FastQC report for the raw data
mkdir raw_fastqc_report
fastqc -o raw_fastqc_report *fastq.gz

# Trim adapters and filter out reads shorter than 15 nt
cutadapt -a AGATCGGAAGAG -m 15 -o trimmed.fastq.gz *fastq.gz

# Generate FastQC report for the processed data
mkdir proc_fastqc_report
fastqc -o proc_fastqc_report trimmed.fastq.gz

echo 'The FASTQ file is ready for alignment. Check FastQC report.'

# Align reads to the genome (specify the path to the Bowtie index)
bowtie -k 1 -p $num_threads ../../BL21_genome/BL21index ./trimmed.fastq.gz -S ./aligned.sam

# Convert SAM to BAM and remove the SAM file
samtools view -@ $num_threads -b aligned.sam > aligned.bam
rm aligned.sam

# Extract relevant columns from the aligned BAM file and save as TSV
samtools view -@ $num_threads aligned.bam | cut -f 1,2,4,10 > aligned.tsv