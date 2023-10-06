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

# Map the reads to the plasmid DNA
bowtie -k 1 -n 1 --best -p $num_threads ./plasmid_reference/ref ./trimmed.fastq.gz -S ./alignment/aligned2plasmid.sam --un ./no_plasmid_reads.fastq

# Map the reads unaligned to the plasmid to the genome
bowtie -k 1 -n 1 --best -p $num_threads ./genome_reference/ref ./no_plasmid_reads.fastq -S ./alignment/aligned2genome.sam
rm no_plasmid_reads.fastq

echo "Made the alignments"

# Change to the alignment directory
cd ./alignment

# Convert SAM to BAM, sort and make index files
samtools view -b aligned2plasmid.sam | samtools sort -@ $num_threads > sorted_aligned2plasmid.bam
samtools view -b aligned2genome.sam | samtools sort -@ $num_threads > sorted_aligned2genome.bam
samtools index sorted_aligned2plasmid.bam
samtools index sorted_aligned2genome.bam

echo "Converted them to BAM files"

# Intersect alignment with genes and report only those intersections that are on the same strand
bedtools intersect -a ../genes.gtf -b ./sorted_aligned2genome.bam -c -s > nt_genes_coverage.tsv
bedtools intersect -a ../genes.gtf -b ./sorted_aligned2genome.bam -c -S > t_genes_coverage.tsv

echo "Intersected with genes"

# Calculate coverage of each position in the plasmid and the genome for positive and negative strands
# Create BAM file with reads aligned to positive strand (or unaligned at all)
samtools view -b -F 16 -@ $num_threads ./sorted_aligned2plasmid.bam > pos_sorted_aligned2plasmid.bam
samtools index pos_sorted_aligned2plasmid.bam

samtools view -b -F 16 -@ $num_threads ./sorted_aligned2genome.bam > pos_sorted_aligned2genome.bam
samtools index pos_sorted_aligned2genome.bam

# Create BAM file only with reads aligned to negative strand
samtools view -b -f 16 -@ $num_threads ./sorted_aligned2plasmid.bam > neg_sorted_aligned2plasmid.bam
samtools index neg_sorted_aligned2plasmid.bam

samtools view -b -f 16 -@ $num_threads ./sorted_aligned2genome.bam > neg_sorted_aligned2genome.bam
samtools index neg_sorted_aligned2genome.bam

# Calculate coverage
samtools depth -a ./pos_sorted_aligned2plasmid.bam > pos_strand_plasmid_cov.tsv
samtools depth -a ./pos_sorted_aligned2genome.bam > pos_strand_genome_cov.tsv
samtools depth -a ./neg_sorted_aligned2plasmid.bam > neg_strand_plasmid_cov.tsv
samtools depth -a ./neg_sorted_aligned2genome.bam > neg_strand_genome_cov.tsv

echo "Calculated coverage of each position in genome and chromosome"

# Calculate number of reads mapped to the plasmid and the genome
# Plasmid
samtools view -F 4 ./sorted_aligned2plasmid.bam | cut -f 1 | sort | uniq | wc -l > n_reads_plasmid.txt
# Genome
samtools view -F 4 ./sorted_aligned2genome.bam | cut -f 1 | sort | uniq | wc -l > n_reads_genome.txt

echo "Done!"
