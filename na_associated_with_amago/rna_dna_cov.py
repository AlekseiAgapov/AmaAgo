#!/usr/bin/env python
# coding: utf-8


# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

# Load the data

rna_pos = pd.read_csv('../dago_agap_rna/alignment/pos_strand_genome_cov.tsv', sep='\t', header=None) # Specify a path to the data
rna_pos.columns = ['chr', 'position', 'rna_pos']

rna_neg = pd.read_csv('../dago_agap_rna/alignment/neg_strand_genome_cov.tsv', sep='\t', header=None) # Specify a path to the data
rna_neg.columns = ['chr', 'position', 'rna_neg']

dna_pos = pd.read_csv('../dago_agap_dna/alignment/pos_strand_genome_cov.tsv', sep='\t', header=None) # Specify a path to the data
dna_pos.columns = ['chr', 'position', 'dna_pos']

dna_neg = pd.read_csv('../dago_agap_dna/alignment/neg_strand_genome_cov.tsv', sep='\t', header=None) # Specify a path to the data
dna_neg.columns = ['chr', 'position', 'dna_neg']

df = rna_pos[['position', 'rna_pos']]
df['rna_neg'] = rna_neg['rna_neg']
df['dna_pos'] = dna_pos['dna_pos']
df['dna_neg'] = dna_neg['dna_neg']
df['rna'] = df['rna_pos'] + df['rna_neg']

# Choose the threshold
threshold = 1

# Choose the length threshold
len_thr = 100

positions = list(df.query('rna > @threshold')['position'])

# Initialize variables
intervals = []
start = positions[0]
end = positions[0]

# Loop through the list to find intervals
for i in range(1, len(positions)):
    if positions[i] == end + 1:
        end = positions[i]
    else:
        intervals.append([start, end])
        start = positions[i]
        end = positions[i]
        
# Don't forget the last interval
intervals.append([start, end])

# Convert to DataFrame
intervals_df = pd.DataFrame(intervals, columns=['start', 'end'])
intervals_df = intervals_df.query('end - start > @len_thr')


# Initialize an empty list to hold the results
results = []

# Loop through each interval and calculate coverage
for i, row in intervals_df.iterrows():
    start, end = row['start'], row['end']
    
    # Filter the coverage DataFrame to only include positions within the current interval
    filtered_coverage = df[(df['position'] >= start) & (df['position'] <= end)]
    
    # Calculate the total and average coverage for the current interval
    rna_pos = filtered_coverage['rna_pos'].sum()
    rna_neg = filtered_coverage['rna_neg'].sum()
    dna_pos = filtered_coverage['dna_pos'].sum()
    dna_neg = filtered_coverage['dna_neg'].sum()
    
    # Append results to the list
    results.append({
        'start': start,
        'end': end,
        'rna_pos': rna_pos,
        'rna_neg': rna_neg,
        'dna_pos': dna_pos,
        'dna_neg': dna_neg,
    })

# Convert the results list to a DataFrame
intervals_df = pd.DataFrame(results)



intervals_df.replace(0, 1, inplace=True)
intervals_df['rna_ratio'] = intervals_df['rna_pos'] / intervals_df['rna_neg']
intervals_df['dna_ratio'] = intervals_df['dna_neg'] / intervals_df['dna_pos']


genes = pd.read_csv('../dago_agap_rna/genes.gtf', sep='\t', header=None)
genes.columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
genes['gene_name'] = genes['attribute'].str.extract('Name=([^;]+)')
genes['gene_biotype'] = genes['attribute'].str.extract('gene_biotype=([^;]+)')
genes = genes[['gene_name', 'gene_biotype', 'start', 'end']]


biotypes = []

for i, row in intervals_df.iterrows():
    bg = ''
    for j, row_genes in genes.iterrows():
        if (row_genes['start'] < row['start'] < row_genes['end']) or (row_genes['start'] < row['end'] < row_genes['end']):
            gb = row_genes['gene_biotype']
            break
    biotypes.append(gb)

intervals_df['gene_biotype'] = biotypes


colours = intervals_df['gene_biotype'].map({'protein_coding': '#84a59d', 'pseudogene': '#84a59d', 'rRNA': 'red', 'tRNA': 'black', 'ncRNA': 'black'})

plt.figure(figsize=(8, 8))

# Create scatter plot
plt.scatter(intervals_df['rna_ratio'], intervals_df['dna_ratio'], s=50, edgecolors=colours, facecolors='none', marker='o', alpha=0.8, linewidths=2)

ax = plt.gca()  # get current axis

# Set font size for x-axis and y-axis labels
plt.xlabel('RNA ratio +/-', size=30)
plt.ylabel('DNA ratio -/+', size=30)

ax.tick_params(axis='both', which='major', labelsize=22, length=8, width=3)
ax.tick_params(axis='both', which='minor', length=5, width=1.5)

# Set axis to log scale
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')

# Add lines
plt.axvline(x=1, color='grey', linestyle='--')
plt.axhline(y=1, color='grey', linestyle='--')

# Remove frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(3)  # Set left spine linewidth to 1.5
ax.spines['bottom'].set_linewidth(3)  # Set bottom spine linewidth to 1.5

# Save the plot as a PNG file
plt.savefig('ratio.png', dpi=600, bbox_inches='tight')



spearman_corr, spearman_p = spearmanr(intervals_df['rna_ratio'], intervals_df['dna_ratio'])
print('spearman_corr =' spearman_corr)
print('spearman_p =' spearman_p)