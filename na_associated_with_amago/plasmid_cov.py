#!/usr/bin/env python
# coding: utf-8


# Import libraries
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data
rna_pos = pd.read_csv('../dago_agap_rna/alignment/pos_strand_plasmid_cov.tsv', sep='\t', header=None) # Specify a path to the data
rna_pos.columns = ['chr', 'position', 'rna_pos']

rna_neg = pd.read_csv('../dago_agap_rna/alignment/neg_strand_plasmid_cov.tsv', sep='\t', header=None) # Specify a path to the data
rna_neg.columns = ['chr', 'position', 'rna_neg']

dna_pos = pd.read_csv('../dago_agap_dna/alignment/pos_strand_plasmid_cov.tsv', sep='\t', header=None) # Specify a path to the data
dna_pos.columns = ['chr', 'position', 'dna_pos']

dna_neg = pd.read_csv('../dago_agap_dna/alignment/neg_strand_plasmid_cov.tsv', sep='\t', header=None) # Specify a path to the data
dna_neg.columns = ['chr', 'position', 'dna_neg']

df = rna_pos[['position', 'rna_pos']]
df['rna_neg'] = rna_neg['rna_neg']
df['dna_pos'] = dna_pos['dna_pos']
df['dna_neg'] = dna_neg['dna_neg']


df.replace(0, 1, inplace=True)
df['rna_pos_n'] = df['rna_pos'] / (df['rna_pos'].sum() + df['rna_neg'].sum()) * 1
df['rna_neg_n'] = df['rna_neg'] / (df['rna_pos'].sum() + df['rna_neg'].sum()) * -1
df['dna_pos_n'] = df['dna_pos'] / (df['dna_pos'].sum() + df['dna_neg'].sum()) * 1
df['dna_neg_n'] = df['dna_neg'] / (df['dna_pos'].sum() + df['dna_neg'].sum()) * -1

df['rna_dif'] = np.log2(df['rna_pos_n'] / df['rna_neg_n'] * -1)
df['dna_dif'] = np.log2(df['dna_pos_n'] / df['dna_neg_n'] * -1)


# Parameters
window_size = 200  # Width of the sliding window
step_size = 20  # Step size for each slide

# Calculate rolling mean with a sliding window
df['dna_pos_n_rolling'] = df['dna_pos_n'].rolling(window=window_size, min_periods=1).mean()
df['dna_neg_n_rolling'] = df['dna_neg_n'].rolling(window=window_size, min_periods=1).mean()
df['rna_pos_n_rolling'] = df['rna_pos_n'].rolling(window=window_size, min_periods=1).mean()
df['rna_neg_n_rolling'] = df['rna_neg_n'].rolling(window=window_size, min_periods=1).mean()

df['rna_dif_rolling'] = df['rna_dif'].rolling(window=window_size, min_periods=1).mean()
df['dna_dif_rolling'] = df['dna_dif'].rolling(window=window_size, min_periods=1).mean()

# Take only every Nth row based on step_size
df_plot = df.iloc[::step_size, :]


plt.figure(figsize=(10, 6))

plt.fill_between(df_plot['position'], df_plot['rna_dif_rolling'], color='skyblue', alpha=0.4)
plt.plot(df_plot['position'], df_plot['rna_dif_rolling'], color='Slateblue', alpha=0.6)

plt.fill_between(df_plot['position'], df_plot['dna_dif_rolling'], color='pink', alpha=0.4)
plt.plot(df_plot['position'], df_plot['dna_dif_rolling'], color='red', alpha=0.6)

plt.xlabel('Plasmid coordinate, kb', size=30)
plt.ylabel('+ strand / - strand', size=30)

plt.axhline(y=0, color='black', linewidth=2)

ax = plt.gca()  # get current axis

ax.tick_params(axis='both', which='major', labelsize=22, length=8, width=3)

# Set the ticks on the y-axis and label them
plt.yticks([-4, -2, 0, 2, 4, 6], ['1/16', '1/4', '1', '4', '16', '64'])

# Set the ticks on the x-axis and label them
plt.xticks([0 ,1000, 2000, 3000, 4000, 5000, 6000, 7000], ['0', '1', '2', '3', '4', '5', '6', '7'])

# Add vertical lines
plt.axvline(x=560, color='black', linestyle='--')
plt.axvline(x=1375, color='black', linestyle='--')
plt.axvline(x=4600, color='black', linestyle='--')
plt.axvline(x=7323, color='black', linestyle='--')

# Remove frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(3)  # Set left spine linewidth to 1.5
ax.spines['bottom'].set_linewidth(3)  # Set bottom spine linewidth to 1.5

# Save the plot as a PNG file
plt.savefig('plasmid_coverage.png', dpi=600, bbox_inches='tight')