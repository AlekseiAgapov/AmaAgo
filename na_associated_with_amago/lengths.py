#!/usr/bin/env python
# coding: utf-8


# Import libraries

import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt


lengths = []

with pysam.AlignmentFile("../dago_agap_rna/alignment/sorted_aligned2plasmid.bam", "rb") as bamfile: # Specify a path to the data
    # Loop through each read in the BAM file
    for read in bamfile:
        # Check if the read is mapped
        if not read.is_unmapped:
            # Get the read length
            lengths.append(read.query_length)

with pysam.AlignmentFile("../dago_agap_rna/alignment/sorted_aligned2genome.bam", "rb") as bamfile: # Specify a path to the data
    # Loop through each read in the BAM file
    for read in bamfile:
        # Check if the read is mapped
        if not read.is_unmapped:
            # Get the read length
            lengths.append(read.query_length)


lengths_df = pd.DataFrame(lengths, columns=['lengths'])
df = lengths_df['lengths'].value_counts().reset_index()
df.columns = ['length', 'counts']
df = df.sort_values('length').reset_index(drop=True)



# Create bar plot
plt.figure(figsize=(12, 4))  # Set the figure size
plt.bar(df['length'], df['counts'], color='#629c99')

ax = plt.gca()  # get current axis

# Set font size for x-axis and y-axis labels
plt.xlabel('Read length', fontsize=30)
plt.ylabel('Number of reads', fontsize=30)

# Set major ticks
major_ticks = [15, 20, 25, 30, 35, 40, 45, 50]
ax.set_xticks(major_ticks)

# Set minor ticks
minor_ticks = [x for x in np.arange(12, 51, 1) if x not in major_ticks]
ax.set_xticks(minor_ticks, minor=True)

# Customize tick labels
major_tick_labels = [str(x) if x != 50 else "≥50" for x in major_ticks]
ax.set_xticklabels(major_tick_labels)

# Set ticks settings
ax.tick_params(axis='both', which='major', labelsize=26, length=8, width=3)
ax.tick_params(which='minor', length=4, width=2, direction='out')

# Remove frame
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(3)  # Set left spine linewidth to 1.5
ax.spines['bottom'].set_linewidth(3)  # Set bottom spine linewidth to 1.5

# Save the plot
plt.savefig('lengths.png', dpi=600, bbox_inches='tight')