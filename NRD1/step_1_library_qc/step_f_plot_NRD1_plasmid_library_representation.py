"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy
"""

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

PROCESSED_DATA_DIR = (
    "/path/to/processed_data/by_project/"
)
DIR = PROCESSED_DATA_DIR + "NNS/20250628_repeat_step_1_library_cloning/"
PLOTS_DIR = DIR + "plots/"

nrd1_df = pd.read_csv(DIR + "NRD1_libraries_step_1_oligo_counts.tsv", sep="\t")


# Plot all libraries on separate subplots (facets) using seaborn.relplot with col wrapping (max 4 columns)
import scipy.stats as stats
import numpy as np
import seaborn as sns
import numpy as np

sorted(nrd1_df["library"].unique())
# Let's set an order for comparison:
# ['V536', 'V622', 'V623', 'V624', 'V625',

# Plot in a single row
library_order = [
    "V536",
    "V622",
    "V623",
    "V624",
    "V625",
]

nrd1_df.columns
# ['library', 'oligo_name', 'grouped_oligo_name_counts', 'unique_bc0',
#    'oligo_seq', 'designed_donor', 'gc_content', 'oligo_seq_designed',
#    'designed_donor_designed', 'gc_content_designed', 'SPS']

# Oligo names look like this:
# SpCas9_YNL251C_NRD1_3_Q_A_CAG_GCA_2_US_syn_changes_1_DS_syn_changes_-_chrXIV_+_174318_TGG_6_174308_174310
# In this example, the aa num is 3, the ref_aa is Q, the edit_aa is A, and the edit_type is missense.

# Extract the relevant information from the oligo_name column.
nrd1_df["aa_num"] = nrd1_df["oligo_name"].str.split("_").str[3].astype(int)
nrd1_df["ref_aa"] = nrd1_df["oligo_name"].str.split("_").str[4]
nrd1_df["edit_aa"] = nrd1_df["oligo_name"].str.split("_").str[5]
nrd1_df["edit_type"] = np.where(
    nrd1_df["ref_aa"] == nrd1_df["edit_aa"], "synonymous", "missense"
)
nrd1_df["edit_type"] = np.where(nrd1_df["edit_aa"] == "*", "PTC", nrd1_df["edit_type"])

# Write the annotated nrd1_df to a new file for record keeping.
nrd1_df.to_csv(
    DIR + "NRD1_libraries_step_1_oligo_counts_annotated.tsv", sep="\t", index=False
)

# Plot the aa_num on the x-axis and the grouped_oligo_name_counts on the y-axis.
# Facet columns by library and rows by edit_type.
# Use a linear scale for both axes.
# Color the points by edit_type (for legend clarity, but each row is a type).
g = sns.relplot(
    data=nrd1_df,
    x="aa_num",
    y="grouped_oligo_name_counts",
    col="library",
    row="edit_type",
    kind="scatter",
    height=4,
    aspect=1.5,
    palette={"missense": "blue", "synonymous": "green", "PTC": "red"},
    facet_kws={"sharey": False},
)

g.set_axis_labels("Amino Acid Position", "Oligo Counts")
g.set(yscale="linear")
g.set_titles(row_template="{row_name}", col_template="{col_name} Library")
plt.subplots_adjust(top=0.9)
g.fig.suptitle("Oligo Counts by Amino Acid Position and Edit Type in NRD1 Libraries")
# Add GC_content as faint black points with its own axes on the right for each subplot
for ax, (library, edit_type) in zip(
    g.axes.flat,
    [(lib, et) for et in nrd1_df["edit_type"].unique() for lib in library_order],
):
    sub_df = nrd1_df[
        (nrd1_df["library"] == library) & (nrd1_df["edit_type"] == edit_type)
    ]
    ax2 = ax.twinx()
    ax2.scatter(
        sub_df["aa_num"],
        sub_df["gc_content"],
        color="black",
        marker="o",
        label="GC Content",
        alpha=0.05,  # Make points much fainter
    )
    ax2.set_ylabel("GC Content")
    ax2.set_ylim(0, 1)
    ax2.tick_params(axis="y", colors="black")
    # Optionally, add legend only to first subplot
    if ax == g.axes.flat[0]:
        ax2.legend(loc="upper right")

plt.savefig(
    PLOTS_DIR + "NRD1_libraries_oligo_counts_by_aa_position_facet_edit_type.png",
    dpi=300,
)
plt.savefig(
    PLOTS_DIR + "NRD1_libraries_oligo_counts_by_aa_position_facet_edit_type.pdf"
)
plt.close()


# Do another version where we only plot V625, and omit the GC content.
v625_df = nrd1_df[nrd1_df["library"] == "V625"]

# Plot all nonzero points as usual, then overlay zero points in red
g = sns.relplot(
    data=v625_df[v625_df["grouped_oligo_name_counts"] != 0],
    x="aa_num",
    y="grouped_oligo_name_counts",
    col="edit_type",
    kind="scatter",
    col_wrap=3,
    height=4,
    aspect=1.5,
    palette={"missense": "blue", "synonymous": "green", "PTC": "red"},
    facet_kws={"sharey": False},
)

# Overlay zero points in red
for ax, edit_type in zip(g.axes.flat, v625_df["edit_type"].unique()):
    zero_df = v625_df[
        (v625_df["edit_type"] == edit_type)
        & (v625_df["grouped_oligo_name_counts"] == 0)
    ]
    ax.scatter(
        zero_df["aa_num"],
        zero_df["grouped_oligo_name_counts"],
        color="red",
        marker="o",
        label="Zero Count",
    )

g.set_axis_labels("Amino Acid Position", "Oligo Counts")
g.set(yscale="linear")
g.set_titles(col_template="{col_name} edits in V625 Library")
plt.subplots_adjust(top=0.85)
g.fig.suptitle("Oligo Counts by Amino Acid Position and Edit Type in V625 Library")
plt.savefig(
    PLOTS_DIR + "V625_library_oligo_counts_by_aa_position_facet_edit_type.png", dpi=300
)
plt.savefig(PLOTS_DIR + "V625_library_oligo_counts_by_aa_position_facet_edit_type.pdf")
plt.close()
