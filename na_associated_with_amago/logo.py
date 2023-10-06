#!/usr/bin/env python
# coding: utf-8


# Import libraries

import pysam
import logomaker
import numpy as np
import matplotlib.pyplot as plt


# Specify how many nucleotides to show upstream and downstream from the 5'-end of the mapped read
# Upstream
left_from_5_prime_end = 5
# Downstream
right_from_5_prime_end = 9


# Function to convert a sequence to reverse and complement
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))


# Function to fetch a sequence from a circular DNA molecule
def fetch_circular_sequence(fasta_file, contig, start, end, length):
    if end > length:
        wrap_around_length = end - length
        return fasta_file.fetch(contig, start, length) + fasta_file.fetch(contig, 0, wrap_around_length)
    elif start < 0:
        wrap_around_start = length + start
        return fasta_file.fetch(contig, wrap_around_start, length + 1) + fasta_file.fetch(contig, 0, end)
    else:
        return fasta_file.fetch(contig, start, end)


# Function to generate a list of sequences for Logo plot
def fetch_seqs_for_logo(bamfile, reference):
    logo_seqs = []
    for read in bamfile:
        # Ignore unmapped reads
        if read.is_unmapped:
            continue
        # Fetch sequences surrounding 5'-ends of the reads mapped to the positive strand
        elif not read.is_reverse:
            five_prime_end = read.reference_start
            extracted_sequence = fetch_circular_sequence(reference, read.reference_name, five_prime_end - left_from_5_prime_end, five_prime_end + 1 + right_from_5_prime_end, reference.get_reference_length(read.reference_name)).upper()
            logo_seqs.append(extracted_sequence)
        # Fetch sequences surrounding 5'-ends of the reads mapped to the negative strand
        elif read.is_reverse:
            five_prime_end = read.reference_end - 1
            extracted_sequence = reverse_complement(fetch_circular_sequence(reference, read.reference_name, five_prime_end - right_from_5_prime_end, five_prime_end + 1 + left_from_5_prime_end, reference.get_reference_length(read.reference_name)).upper())
            logo_seqs.append(extracted_sequence)
    return logo_seqs


# Function to create and save a sequence logo from a list of sequences
def create_sequence_logo_prob(sequences, output_filename):
    # Create a count matrix from the list of sequences
    count_matrix = logomaker.alignment_to_matrix(sequences)
    # Normalize the count matrix to create a probability matrix
    prob_matrix = count_matrix.divide(count_matrix.sum(axis=1), axis=0)

    # Create a sequence logo using the weighted probability matrix
    logo = logomaker.Logo(prob_matrix, figsize=(12, 4), color_scheme='colorblind_safe')

    # Set font size for x-axis and y-axis labels
    logo.ax.set_xlabel('Position', fontsize=30)
    logo.ax.set_ylabel('Probability', fontsize=30)

    # Set font size for axis tick labels
    logo.ax.tick_params(axis='both', labelsize=26, width=3)

    # Set the y-axis limits
    logo.ax.set_ylim(0, 1)

    # Create new labels excluding zero
    new_labels = [i for i in range(left_from_5_prime_end * -1, right_from_5_prime_end + 2) if i != 0]

    # Set new tick positions and labels
    logo.ax.set_xticks(range(len(new_labels)))
    logo.ax.set_xticklabels(new_labels)

    # Remove frame
    logo.ax.spines['right'].set_visible(False)
    logo.ax.spines['top'].set_visible(False)
    logo.ax.spines['left'].set_linewidth(3)  # Set left spine linewidth to 1.5
    logo.ax.spines['bottom'].set_linewidth(3)  # Set bottom spine linewidth to 1.5
    
    # Save the plot as a PNG file
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')


# Load the BAM files
plasmid_alignment = pysam.AlignmentFile("../dago_agap_rna/alignment/sorted_aligned2plasmid.bam", "rb") # Specify a path to the data
genome_alignment = pysam.AlignmentFile("../dago_agap_rna/alignment/sorted_aligned2genome.bam", "rb") # Specify a path to the data


# Load the FASTA files
genome = pysam.FastaFile("../dago_agap_rna/genome_reference/genome.fa") # Specify a path to the data
plasmid = pysam.FastaFile("../dago_agap_rna/plasmid_reference/plasmid.fa") # Specify a path to the data


# Generate list of sequences for Logo plot
plasmid_logo_seqs = fetch_seqs_for_logo(plasmid_alignment, plasmid)
genome_logo_seqs = fetch_seqs_for_logo(genome_alignment, genome)

# Merge them into one list
logo_seqs = plasmid_logo_seqs + genome_logo_seqs

# Change 'T' for 'U'
logo_seqs = [seq.replace('T', 'U') for seq in logo_seqs]


# Make a probability plot and save it into PNG file
create_sequence_logo_prob(logo_seqs, "sequence_prob_plot.png")