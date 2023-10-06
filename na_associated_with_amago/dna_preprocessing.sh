#!/usr/bin/bash

# Set default values for arguments
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

# Trim adapters and filter out reads shorter than 14 nt
cutadapt -a TGGAATTCTCGGGTGCCAAGG -m 14 -o trimmed.fastq.gz *fastq.gz

# Generate FastQC report for the processed data
mkdir proc_fastqc_report
fastqc -o proc_fastqc_report trimmed.fastq.gz

echo 'The FASTQ file is ready for alignment. Check FastQC report.'
