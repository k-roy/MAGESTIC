"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy
"""

import pandas as pd
from sys import argv
import matplotlib.pyplot as plt

infilename = argv[1]
outfilename = argv[2]

# infilename = '/path/to/processed_data/by_project/NNS/20240411_Twist_200mer_subpool_1-16_step_1_libraries/guide_donor_bc0/V536_matched_guide_donor_bc0_counts.tsv'
# outfilename = '/path/to/processed_data/by_project/NNS/20240411_Twist_200mer_subpool_1-16_step_1_libraries/guide_donor_bc0_purity_filtered/V536_matched_guide_donor_bc0_purity_filtered.tsv'
df = pd.read_csv(infilename, sep='\t')

# If the guide and bc0 are the same, sum the counts.
summed_guide_bc0_df = df.groupby(['guide', 'bc0'])['counts'].sum().reset_index()
summed_guide_bc0_df.rename(columns={'counts': 'grouped_guide_bc0_counts'}, inplace=True)

# Retain only the top guide-bc0 combination from the initial table.
idx = df.groupby(['guide', 'bc0'])['counts'].idxmax()
top_guide_bc0_df = df.loc[idx]

# Merge
merged_guide_bc0_df = top_guide_bc0_df.merge(summed_guide_bc0_df)

# Assess guide purity by taking the guide-bc0 total counts and dividing by the top guide-donor-bc0 counts.
merged_guide_bc0_df.columns
merged_guide_bc0_df['guide_bc0_purity'] = merged_guide_bc0_df['counts'] / merged_guide_bc0_df['grouped_guide_bc0_counts']

# filter for 3 read minimum
filtered_guide_bc0_df = merged_guide_bc0_df.query('counts > 2')
# filtering brings the file size from 600k to 300k, so perhaps it is unnecessary to filter.
merged_guide_bc0_df.columns
merged_guide_bc0_df.to_csv(outfilename, sep='\t', index=False)

'''
# Create the density plot
merged_guide_bc0_df['guide_bc0_purity'].plot.density()
# Add labels and title
plt.xlabel('guide_bc0_purity')
plt.ylabel('Density')
plt.title('Density Plot of guide_purity')
# Save the plot
PLOTS_DIR = '/path/to/scripts_and_keyfiles/by_project/NNS/guide_donor_bc0/plots/'
plt.savefig(PLOTS_DIR + 'density_plot.png')

# Clear the current figure and plot histogram
plt.clf()  
df_merged['guide_purity'].plot.hist(bins=100)
# Add labels and title
plt.xlabel('guide_purity')
plt.ylabel('Counts')
plt.title('Histogram of guide_purity')
# Save the plot
plt.savefig(PLOTS_DIR + 'hist_plot.png')

# Clear the current figure and plot purity by read counts
plt.clf()  
df_merged.plot.hexbin(x='guide_purity', y='grouped_guide_bc0_counts')
plt.show()
# Add labels and title
plt.xlabel('guide_purity')
plt.ylabel('grouped_guide_bc0_counts')
plt.title('hexbin plot of counts vs guide_purity')
# Save the plot
plt.savefig(PLOTS_DIR + 'hexbin_plot.png')

# Filter the guide table for guide_purity > 10%. Expect PCR recombination to result in many guide combinations with <10%.
'''