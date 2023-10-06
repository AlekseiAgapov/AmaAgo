#!/usr/bin/env python
# coding: utf-8

# This script processes alignment data from AgAP-treated RNA and control samples to create
# sequence logos for AgAP cleavage sites. The script loads data, processes it, identifies
# cleavage sites, and generates sequence logos for visualization.


# Specify threshold for comparing data from AgAP treated RNA and control samples
threshold = 100


# 1. Importing libraries

import pandas as pd
import logomaker
import numpy as np
import matplotlib.pyplot as plt


# 2. Reading data

# Function to load alignment data from the specified directory.
def load_alignment(path_to_directory):
    # Construct the file path by concatenating the directory path and the file name
    file_path = path_to_directory + 'aligned.tsv'
    
    # Read the TSV file into a pandas DataFrame
    df = pd.read_csv(file_path, sep='\t', header=None)
    
    # Assign column names to the DataFrame
    df.columns = ['read_name', 'flag', 'position', 'sequence']
    
    # Calculate the length of each sequence and store it in a new column called 'length'
    df['length'] = df['sequence'].apply(lambda x: len(x))
    
    # Drop the 'sequence' column as it's not needed for further analysis
    df = df.drop(columns='sequence')
    
    # Return the processed DataFrame
    return df

# Load alignment data from 6 libraries into DataFrames
control_1 = load_alignment('./reads/control_1/') # Specify a path to the data
control_2 = load_alignment('./reads/control_2/') # Specify a path to the data
agap_1 = load_alignment('./reads/AgAP_1/') # Specify a path to the data
agap_2 = load_alignment('./reads/AgAP_2/') # Specify a path to the data
agap_3 = load_alignment('./reads/AgAP_3/') # Specify a path to the data
agap_ago_1 = load_alignment('./reads/AgAP_Ago_1/') # Specify a path to the data


# Read the genome file as a pandas DataFrame
genome_df = pd.read_csv('./BL21_genome/BL21.fa', sep='\t') # Specify a path to the data

# Rename the column to 'fasta'
genome_df.columns = ['fasta']

# Concatenate all the sequences in the 'fasta' column into a single string
genome_sequence = ''.join(genome_df['fasta'])

# Calculate the length of the genome sequence
genome_length = len(genome_sequence)


# 3. Processing data to get a DataFrame with coverage of reads starting positions expressed in RPM for each library

def process_data(df):
    """
    This function processes the input DataFrame containing alignment data by separating reads based on their
    orientation (plus or minus) and calculating the coverage at each start position.

    Args:
        df (pd.DataFrame): A DataFrame containing alignment data with columns: read_name, flag, position, and length.

    Returns:
        plus_pivot (pd.DataFrame): A DataFrame containing coverage data for plus strand reads, with columns: start_position and coverage.
        minus_pivot (pd.DataFrame): A DataFrame containing coverage data for minus strand reads, with columns: start_position and coverage.
    """
    # Filter the DataFrame for records with flag equal to 0 and create a copy
    plus_df = df.query('flag == 0').copy()
    
    # Add a new column 'start_position' with values equal to the 'position' column
    plus_df.loc[:, 'start_position'] = plus_df['position']
    
    # Create a pivot table to count the number of reads per start position, and reset the index
    plus_pivot = plus_df.pivot_table(index='start_position', values='read_name', aggfunc='count').reset_index()
    
    # Rename the columns of the resulting DataFrame
    plus_pivot.columns = ['start_position', 'coverage']
    
    # Filter the DataFrame for records with flag equal to 16 and create a copy
    minus_df = df.query('flag == 16').copy()
    
    # Add a new column 'start_position' with values equal to the 'position' plus 'length' minus 1
    minus_df.loc[:, 'start_position'] = minus_df['position'] + minus_df['length'] - 1
    
    # Create a pivot table to count the number of reads per start position, and reset the index
    minus_pivot = minus_df.pivot_table(index='start_position', values='read_name', aggfunc='count').reset_index()
    
    # Rename the columns of the resulting DataFrame
    minus_pivot.columns = ['start_position', 'coverage']
    
    # Return the processed plus and minus DataFrames
    return plus_pivot, minus_pivot


# Separating reads from each library to two datasets based on their orientation and calculating the coverage at each start position.
control_1_plus_pivot, control_1_minus_pivot = process_data(control_1)
control_2_plus_pivot, control_2_minus_pivot = process_data(control_2)
agap_1_plus_pivot, agap_1_minus_pivot = process_data(agap_1)
agap_2_plus_pivot, agap_2_minus_pivot = process_data(agap_2)
agap_3_plus_pivot, agap_3_minus_pivot = process_data(agap_3)
agap_ago_1_plus_pivot, agap_ago_1_minus_pivot = process_data(agap_ago_1)


def merge_coverages(pivot_dataframes, genome_length):
    """
    This function creates a DataFrame with coverage data from multiple samples for plus and minus strands.
    
    Args:
        pivot_dataframes (list): A list of DataFrames containing the coverage data for each sample.
        genome_length (int): The length of the genome.

    Returns:
        coverage_df (pd.DataFrame): A DataFrame containing merged coverage data for all samples.
    """
    # Create a DataFrame with a range of positions from 1 to genome_length + 1
    positions = list(range(1, genome_length + 1))
    coverage_df = pd.DataFrame({'start_position': positions})

    # Iterate over the pivot_dataframes and merge them with the coverage_df
    for pivot_df in pivot_dataframes:
        coverage_df = coverage_df.merge(pivot_df, how='left', on='start_position')

    # Fill NaN values with 0
    coverage_df.fillna(0, inplace=True)

    return coverage_df


# List of DataFrames for plus and minus strands
plus_pivot_dataframes = [control_1_plus_pivot, control_2_plus_pivot, agap_1_plus_pivot, agap_2_plus_pivot, agap_3_plus_pivot, agap_ago_1_plus_pivot]
minus_pivot_dataframes = [control_1_minus_pivot, control_2_minus_pivot, agap_1_minus_pivot, agap_2_minus_pivot, agap_3_minus_pivot, agap_ago_1_minus_pivot]

# Merge coverage data for plus and minus strands
plus_coverage = merge_coverages(plus_pivot_dataframes, genome_length)
minus_coverage = merge_coverages(minus_pivot_dataframes, genome_length)

# Set column names for merged DataFrames
column_names = ['position', 'control_1', 'control_2', 'agap_1', 'agap_2', 'agap_3', 'agap_ago_1']
plus_coverage.columns = column_names
minus_coverage.columns = column_names


# Normalize read counts to reads per million (RPM)
for column in ['control_1', 'control_2', 'agap_1', 'agap_2', 'agap_3', 'agap_ago_1']:
    new_column_name = column + '_rpm'
    reads_number = plus_coverage[column].sum() + minus_coverage[column].sum()
    plus_coverage[column] = plus_coverage[column].fillna(1)
    plus_coverage[new_column_name] = plus_coverage[column] / reads_number * 1000000
    minus_coverage[column] = minus_coverage[column].fillna(1)
    minus_coverage[new_column_name] = minus_coverage[column] / reads_number * 1000000

# Calculate average RPM values for control, agap, and agap_ago for plus and minus strands
plus_coverage['control'] = (plus_coverage['control_1_rpm'] + plus_coverage['control_2_rpm']) / 2
plus_coverage['agap'] = (plus_coverage['agap_1_rpm'] + plus_coverage['agap_2_rpm'] + plus_coverage['agap_3_rpm']) / 3
plus_coverage['agap_ago'] = plus_coverage['agap_ago_1_rpm']
plus_coverage.replace(0, 1, inplace=True)

minus_coverage['control'] = (minus_coverage['control_1_rpm'] + minus_coverage['control_2_rpm']) / 2
minus_coverage['agap'] = (minus_coverage['agap_1_rpm'] + minus_coverage['agap_2_rpm'] + minus_coverage['agap_3_rpm']) / 3
minus_coverage['agap_ago'] = minus_coverage['agap_ago_1_rpm']
minus_coverage.replace(0, 1, inplace=True)


# 4. Retrieving the sequences of AgAP cleavage sites (in the absence or presence of AmaAgo protein)

# Identify positions where agap and agap_ago RPM values are greater than or equal to the control RPM values multiplied by the threshold
agap_plus_coordinates = plus_coverage.query('agap >= control * @threshold')['position']
agap_minus_coordinates = minus_coverage.query('agap >= control * @threshold')['position']

agap_ago_plus_coordinates = plus_coverage.query('agap_ago >= control * @threshold')['position']
agap_ago_minus_coordinates = minus_coverage.query('agap_ago >= control * @threshold')['position']

# Get sequences for the identified positions
agap_plus_sites = [genome_sequence[position-6 : position+4] for position in agap_plus_coordinates]
agap_ago_plus_sites = [genome_sequence[position-6 : position+4] for position in agap_ago_plus_coordinates]

# Function to compute the reverse complement of a given DNA sequence
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# Get reverse complement sequences for the identified positions on the minus strand
agap_minus_sites = [reverse_complement(genome_sequence[position-5 : position+5]) for position in agap_minus_coordinates]
agap_ago_minus_sites = [reverse_complement(genome_sequence[position-5 : position+5]) for position in agap_ago_minus_coordinates]

# Unite sequences from positive and negative strands for each experimental setup
agap_sites = agap_plus_sites + agap_minus_sites
agap_ago_sites = agap_ago_plus_sites + agap_ago_minus_sites


# 5. Making sequence logos of the cleavage site

# Function to create and save a sequence logo from a list of sequences
def create_sequence_logo(sequences, output_filename):
    # Create a dataframe from the list of sequences
    df_sites = pd.DataFrame({'seqs': sequences})
    # Convert DNA into RNA by replacing Ts to Us
    df_sites['seqs'] = df_sites['seqs'].apply(lambda x : x.replace('T', 'U'))
    # Create a count matrix from the dataframe of sequences
    count_matrix = logomaker.alignment_to_matrix(df_sites['seqs'])
    # Normalize the count matrix to create a probability matrix
    prob_matrix = count_matrix.divide(count_matrix.sum(axis=1), axis=0)
    # Compute the entropy for each position
    entropy = -(prob_matrix * np.log2(prob_matrix)).sum(axis=1).fillna(0)
    # Compute the information content matrix by subtracting the entropy from the maximum possible entropy (2 bits for DNA)
    info_matrix = 2 - entropy
    # Expand the information content matrix into a weighted probability matrix
    weighted_matrix = prob_matrix.multiply(info_matrix, axis=0)

    # Create a sequence logo using the weighted probability matrix
    logo = logomaker.Logo(weighted_matrix, figsize=(12, 4), color_scheme='colorblind_safe')

    # Set font size for x-axis and y-axis labels
    logo.ax.set_xlabel('Position', fontsize=22)
    logo.ax.set_ylabel('Bits', fontsize=22)

    # Set font size for axis tick labels
    logo.ax.tick_params(axis='both', labelsize=18, width=3)

    # Set the y-axis limits
    logo.ax.set_ylim(0, 2)

    # Update the x-axis tick labels to start with 1
    logo.ax.set_xticks(range(weighted_matrix.shape[0]))
    logo.ax.set_xticklabels(range(1, weighted_matrix.shape[0] + 1))

    # Remove frame
    logo.ax.spines['right'].set_visible(False)
    logo.ax.spines['top'].set_visible(False)
    logo.ax.spines['left'].set_linewidth(3)  # Set left spine linewidth to 1.5
    logo.ax.spines['bottom'].set_linewidth(3)  # Set bottom spine linewidth to 1.5

    # Save the plot as a PNG file
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')

# Create and save sequence logos for AgAP and AgAP_Ago
create_sequence_logo(agap_sites, 'sequence_logo_AgAP.png')
create_sequence_logo(agap_ago_sites, 'sequence_logo_AgAP_Ago.png')