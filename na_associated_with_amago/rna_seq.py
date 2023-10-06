#!/usr/bin/env python
# coding: utf-8


# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# Read TSV files into a pandas DataFrame
rna_nt = pd.read_csv('../dago_agap_rna/alignment/nt_genes_coverage.tsv', sep='\t', header=None) # Specify a path to the data
rna_nt.columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute', 'rna_nt']
rna_nt['gene_name'] = rna_nt['attribute'].str.extract('Name=([^;]+)')
rna_nt['gene_biotype'] = rna_nt['attribute'].str.extract('gene_biotype=([^;]+)')

rna_seq = pd.read_csv('../rna_seq/alignment/genes_coverage.tsv', sep='\t', header=None) # Specify a path to the data
rna_seq.columns = ['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute', 'rna_seq']


df = rna_nt[['gene_name', 'gene_biotype', 'start', 'end', 'rna_nt']]
df['rna_seq'] = rna_seq['rna_seq']


# Load the number of reads
with open('../dago_agap_rna/alignment/n_reads_genome.txt', 'r') as f: # Specify a path to the data
    n_reads_rna = int(f.read().strip())

with open('../rna_seq/alignment/n_reads_genome.txt', 'r') as f: # Specify a path to the data
    n_reads_rna_seq = int(f.read().strip())


df_plot_rna_seq = df.query('gene_biotype != "rRNA"')

# Create a color map based on the 'condition' column
colours = df_plot_rna_seq['gene_biotype'].map({'protein_coding': '#84a59d', 'pseudogene': '#84a59d', 'tRNA': 'black', 'ncRNA': '#eb6424'})


plt.figure(figsize=(8, 8))

# Create scatter plot
plt.scatter(df_plot_rna_seq['rna_seq_rpkm'], df_plot_rna_seq['rna_nt_rpkm'], s=50, edgecolors=colours, facecolors='none', marker='o', alpha=0.8, linewidths=2)

ax = plt.gca()  # get current axis

# Set font size for x-axis and y-axis labels
plt.xlabel('RNA-seq, RPKM', size=30)
plt.ylabel('RNA guides, RPKM', size=30)

ax.tick_params(axis='both', which='major', labelsize=22, length=8, width=3)
ax.tick_params(axis='both', which='minor', length=5, width=1.5)

# Set axis to log scale
plt.gca().set_xscale('log')
plt.gca().set_yscale('log')

# Remove frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(3)  # Set left spine linewidth to 1.5
ax.spines['bottom'].set_linewidth(3)  # Set bottom spine linewidth to 1.5

# Save the plot as a PNG file
plt.savefig('rna_seq.png', dpi=600, bbox_inches='tight')


spearman_corr, spearman_p = spearmanr(df_plot_rna_seq['rna_seq_rpkm'], df_plot_rna_seq['rna_nt_rpkm'])
print('spearman_corr =' spearman_corr)
print('spearman_p =' spearman_p)