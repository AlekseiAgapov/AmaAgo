#!/usr/bin/bash

# DESCRIPTION

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

echo "Starting the job"

# Change to the working directory
cd $working_directory

# Create a directory for alignment results
mkdir alignment

# Map the reads unaligned to the plasmid to the genome
bowtie -k 1 -n 2 --best -p $num_threads ./genome_reference/ref ./*.fastq.gz -S ./alignment/aligned2genome.sam

echo "Made the alignment"

# Change to the alignment directory
cd ./alignment

# Convert SAM to BAM, sort and make index files
samtools view -b aligned2genome.sam | samtools sort -@ $num_threads > sorted_aligned2genome.bam
samtools index sorted_aligned2genome.bam

echo "Converted it to a BAM file"

# Intersect alignment with genes and report only those intersections that are on the same strand
bedtools intersect -a ../genes.gtf -b ./sorted_aligned2genome.bam -c -s > genes_coverage.tsv

echo "Intersected with genes"

# Calculate number of reads mapped to the genome
samtools view -F 4 ./sorted_aligned2genome.bam | cut -f 1 | sort | uniq | wc -l > n_reads_genome.txt

echo "Done!"
