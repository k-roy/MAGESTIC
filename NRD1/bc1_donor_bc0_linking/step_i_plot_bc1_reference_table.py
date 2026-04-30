"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


BASE_DIR = "/path/to/processed_data/by_project/"

PLOTS_DIR = os.path.join(BASE_DIR, "NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/plots")
if not os.path.exists(PLOTS_DIR):
    os.makedirs(PLOTS_DIR)

BC1_DONOR_BC0_DIR = os.path.join(
    BASE_DIR, "NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/data_tables"
)

BC1_REFERENCE_TABLE_PATH = os.path.join(
    BC1_DONOR_BC0_DIR, "20250820_NRD1_bc1_reference_table.tsv"
)

# Plot the number of unique bc1 per aa position
bc1_reference_df = pd.read_csv(BC1_REFERENCE_TABLE_PATH, sep="\t")

mapped_oligo_names_df = bc1_reference_df[bc1_reference_df["aa_nums_str"].notnull()]
mapped_oligo_names_df.columns

# If alt_aas_str == *, then edit_type is "PTC"
# If ref_aas_str == alt_aas_str, then edit_type is "synonymous"
# If ref_aas_str != alt_aas_str, then edit_type is "missense"
mapped_oligo_names_df["edit_type"] = mapped_oligo_names_df.apply(
    lambda row: (
        "PTC"
        if row["alt_aas_str"] == "*"
        else ("synonymous" if row["ref_aas_str"] == row["alt_aas_str"] else "missense")
    ),
    axis=1,
)
# Ensure all AA positions from 1 to the length of NRD1 protein are present
nrd1_length = mapped_oligo_names_df["aa_nums_str"].astype(int).max()
aa_positions = pd.Series(range(1, nrd1_length + 1), name="aa_nums_str")
# Merge to ensure all positions are present
plot_df = aa_positions.to_frame().merge(
    mapped_oligo_names_df, on="aa_nums_str", how="left"
)
edit_types = ["PTC", "synonymous", "missense"]

fig, axes = plt.subplots(
    len(edit_types), 1, figsize=(12, 6 * len(edit_types)), sharex=True
)

for i, edit_type in enumerate(edit_types):
    subset = plot_df[plot_df["edit_type"] == edit_type]
    # Count unique bc1 per position, fill missing positions with zero
    counts = (
        subset.groupby("aa_nums_str")["bc1"]
        .nunique()
        .reindex(range(1, nrd1_length + 1), fill_value=0)
    )
    axes[i].scatter(
        counts.index, counts.values, color=sns.color_palette("viridis")[i], s=30
    )
    axes[i].set_title(f"Number of Unique bc1 per AA Position ({edit_type})")
    axes[i].set_xlabel("AA Position")
    axes[i].set_ylabel("Count")
    xticks = list(range(1, nrd1_length + 1, 50))
    axes[i].set_xticks(xticks)
    axes[i].set_xticklabels([str(x) for x in xticks], rotation=45)

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "bc1_per_aa_position_by_edit_type_scatter.png"))
plt.close()

# Scatterplot for total_bc1_donor_bc0_counts per edit type in subplots
fig, axes = plt.subplots(
    len(edit_types), 1, figsize=(12, 6 * len(edit_types)), sharex=True
)

for i, edit_type in enumerate(edit_types):
    subset = plot_df[plot_df["edit_type"] == edit_type]
    # Sum total_bc1_donor_bc0_counts per position, fill missing positions with zero
    total_counts = (
        subset.groupby("aa_nums_str")["total_bc1_donor_bc0_counts"]
        .sum()
        .reindex(range(1, nrd1_length + 1), fill_value=0)
    )
    axes[i].scatter(
        total_counts.index,
        total_counts.values,
        color=sns.color_palette("viridis")[i],
        s=30,
    )
    axes[i].set_title(f"Total bc1 Donor bc0 Counts at Each AA Position ({edit_type})")
    axes[i].set_xlabel("AA Position")
    axes[i].set_ylabel("Total Counts")
    axes[i].set_xticks(xticks)
    axes[i].set_xticklabels([str(x) for x in xticks], rotation=45)

plt.tight_layout()
plt.savefig(
    os.path.join(
        PLOTS_DIR, "total_bc1_donor_bc0_counts_per_aa_position_by_edit_type_scatter.png"
    )
)
plt.close()


# Make a plot for how many aa positions are represented by how many unique bc1 for each edit type
unique_bc1_counts = (
    mapped_oligo_names_df.groupby("aa_nums_str")["bc1"]
    .nunique()
    .reindex(range(1, nrd1_length + 1), fill_value=0)
)
edit_type_counts = (
    mapped_oligo_names_df.groupby("aa_nums_str")["edit_type"]
    .first()
    .reindex(range(1, nrd1_length + 1), fill_value="unknown")
)
unique_bc1_per_edit_type = pd.DataFrame(
    {"unique_bc1_count": unique_bc1_counts, "edit_type": edit_type_counts}
).reset_index()
# Count how many aa positions have each unique_bc1_count
unique_bc1_distribution = (
    unique_bc1_per_edit_type.groupby("unique_bc1_count")["edit_type"]
    .count()
    .reset_index()
)
# Plot the distribution of unique bc1 counts per AA position, separated by edit_type
unique_bc1_per_edit_type["edit_type"] = unique_bc1_per_edit_type["edit_type"].replace(
    "unknown", pd.NA
)
edit_types = ["PTC", "synonymous", "missense"]

fig, axes = plt.subplots(
    len(edit_types), 1, figsize=(10, 6 * len(edit_types)), sharex=True
)
# Define bins: 0, 1, 2, 3, 4, 5, 6-10, 11-20, 21-30, 31-40, 41-50, 51-100, 101-200, 201-500, 501+
bins = [0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 200, 500, np.inf]
bin_labels = [
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6-10",
    "11-20",
    "21-30",
    "31-40",
    "41-50",
    "51-100",
    "101-200",
    "201-500",
]

for i, edit_type in enumerate(edit_types):
    subset = unique_bc1_per_edit_type[
        unique_bc1_per_edit_type["edit_type"] == edit_type
    ].copy()
    # Bin the unique_bc1_count values
    subset["bin"] = pd.cut(
        subset["unique_bc1_count"],
        bins=bins,
        labels=bin_labels,
        right=True,
        include_lowest=True,
    )
    bin_counts = (
        subset.groupby("bin")["aa_nums_str"]
        .count()
        .reset_index()
        .rename(columns={"aa_nums_str": "aa_position_count"})
    )
    sns.barplot(
        x="bin",
        y="aa_position_count",
        data=bin_counts,
        palette="viridis",
        ax=axes[i],
    )
    axes[i].set_title(
        f"Distribution of Unique bc1 Counts per AA Position ({edit_type})"
    )
    axes[i].set_xlabel("Unique bc1 Count Bin")
    axes[i].set_ylabel("Number of AA Positions")
    axes[i].set_xticklabels(bin_labels, rotation=45)

plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "unique_bc1_count_distribution_by_edit_type.png"))
plt.close()
