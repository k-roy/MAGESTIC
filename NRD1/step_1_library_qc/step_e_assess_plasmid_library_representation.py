"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy
"""

import glob
import os
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

PROCESSED_DATA_DIR = (
    "/path/to/processed_data/by_project/"
)
PROJECT_PATH = (
    "NNS/202250628_repeat_step_1_library_cloning/guide_donor_bc0_purity_filtered/"
)
DIR = PROCESSED_DATA_DIR + PROJECT_PATH
OUT_DIR = PROCESSED_DATA_DIR + "NNS/202250628_repeat_step_1_library_cloning/"
suffix = "_guide_donor_bc0_purity_filtered.tsv"

# Find all files in the directory that end with the suffix and read them into a single DataFrame.
# Put the prefix as a column in the DataFrame.

files = glob.glob(DIR + "*" + suffix)
df_list = []
for file in files:
    # Extract the prefix from the file name
    prefix = os.path.basename(file).replace(suffix, "")
    print(f"Reading file: {file} with prefix: {prefix}")
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep="\t")
    # Add the prefix as a column in the DataFrame
    df["prefix"] = prefix
    # Append the DataFrame to the list
    df_list.append(df)

# For comparison, we should also plot the previous V536 and V545 data.
# Only load these libraries if they exist.
PREVIOUS_DATA_DIR = (
    PROCESSED_DATA_DIR
    + "NNS/20240411_Twist_200mer_subpool_1-16_and_42-65_step_1_libraries/guide_donor_bc0_purity_filtered/"
)
previous_files = glob.glob(PREVIOUS_DATA_DIR + "*_guide_donor_bc0_purity_filtered.tsv")
previous_files = [f for f in previous_files if "V536" in f or "V545" in f]
if previous_files:
    for file in previous_files:
        prefix = os.path.basename(file).replace(
            "_guide_donor_bc0_purity_filtered.tsv", ""
        )
        print(f"Reading previous file: {file} with prefix: {prefix}")
        df = pd.read_csv(file, sep="\t")
        df["prefix"] = prefix
        df_list.append(df)

# Concatenate all DataFrames in the list into a single DataFrame.
combined_df = pd.concat(df_list, ignore_index=True)

# Rename prefix as library.
combined_df.rename(columns={"prefix": "library"}, inplace=True)

combined_df.columns


# Get GC content for each donor sequence.
def gc_content(seq):
    """Calculate the GC content of a sequence."""
    seq_str = str(seq)
    g_count = seq_str.count("G")
    c_count = seq_str.count("C")
    return (g_count + c_count) / len(seq_str) if len(seq_str) > 0 else 0


combined_df["gc_content"] = combined_df["donor"].apply(gc_content)

# How many donors are not string objects?
non_string_donors = combined_df[
    ~combined_df["donor"].apply(lambda x: isinstance(x, str))
]
print(f"Number of non-string donors: {len(non_string_donors)}")
if not non_string_donors.empty:
    print(non_string_donors)

# Write these non-string donors to a file.
non_string_donors.to_csv(OUT_DIR + "non_string_donors.tsv", sep="\t", index=False)

# Read in the designed oligos.
# /path/to/scripts_and_keyfiles/by_project/NNS/20240411_Twist_200mer_subpool_1-16_and_42-65_step_1_libraries_guide_donor_bc0/keyfiles
KEYFILE_DIR = "/path/to/scripts_and_keyfiles/by_project/NNS/20240411_Twist_200mer_subpool_1-16_and_42-65_step_1_libraries_guide_donor_bc0/keyfiles/"
designed_oligos_filename = KEYFILE_DIR + "20240411_Twist_200mer_oligo_array_order.tsv"
designed_oligos = pd.read_csv(designed_oligos_filename, sep="\t")

# Calculate GC content for designed oligos. The donor begins at pos 51 and is 129 bases long.
# ctgcgattggcaggcgcgccGTTATGGTTTTCATCTGACAgtttgaagagcATTGACTTACGTTTTTCCCGGACTGTCCTTTCATAATATAATAACCATCTGCAAGCCtaaagcgacgagaatCATAACAGTGATGTTCAAGATATTCCTTCACCTGAACTATCCGTCGATAGTAACTCTgaactacgtgtctcgtcagg

designed_oligos["designed_donor"] = (
    designed_oligos["oligo_seq"].str[51:180].str.upper()
)  # Extract the donor sequence from the designed oligos
designed_oligos["SPS"] = (
    designed_oligos["oligo_seq"].str[180:].str.upper()
)  # Extract the donor sequence from the designed oligos

designed_oligos["gc_content"] = designed_oligos["designed_donor"].apply(gc_content)

# Group by library, oligo_name, sum the grouped_guide_bc0_counts, and count unique bc0s.
oligo_name_grouped_df = combined_df.groupby(
    ["library", "oligo_name"], as_index=False
).agg(
    grouped_guide_bc0_counts=("grouped_guide_bc0_counts", "sum"),
    unique_bc0=("bc0", "nunique"),
)

# NRD1_mut	2308	cgatagacgactggacagca
# SpG NGNG Bloom QTL subpool	10430	cttgactcgacagtgagagc
# ['V536', 'V622', 'V623', 'V624', 'V625' --> merge with SPS == cgatagacgactggacagca
# 'V545',  'V626', 'V627', 'V628', 'V629' --> merge with SPS == cttgactcgacagtgagagc

# Do a merge the designed oligos with the grouped DataFrame on the oligo_name column.
oligo_name_grouped_df = oligo_name_grouped_df.merge(
    designed_oligos[["oligo_name", "oligo_seq", "designed_donor", "gc_content"]],
    on="oligo_name",
    how="left",
    suffixes=("", "_designed"),
)

# For each library_ID, check if any designed_oligos were not matched.
# Add these and assign with 0 counts.
# Define the mapping of libraries to their expected SPS (use uppercase for consistency)
library_sps_map = {
    "V536": "CGATAGACGACTGGACAGCA",
    "V622": "CGATAGACGACTGGACAGCA",
    "V623": "CGATAGACGACTGGACAGCA",
    "V624": "CGATAGACGACTGGACAGCA",
    "V625": "CGATAGACGACTGGACAGCA",
    "V545": "CTTGACTCGACAGTGAGAGC",
    "V626": "CTTGACTCGACAGTGAGAGC",
    "V627": "CTTGACTCGACAGTGAGAGC",
    "V628": "CTTGACTCGACAGTGAGAGC",
    "V629": "CTTGACTCGACAGTGAGAGC",
}
# Ensure SPS and oligo_seq columns are uppercase for consistency
designed_oligos["SPS"] = designed_oligos["SPS"].str.upper()
designed_oligos["oligo_seq"] = designed_oligos["oligo_seq"].str.upper()
if "SPS" in oligo_name_grouped_df.columns:
    oligo_name_grouped_df["SPS"] = oligo_name_grouped_df["SPS"].str.upper()
if "oligo_seq" in oligo_name_grouped_df.columns:
    oligo_name_grouped_df["oligo_seq"] = oligo_name_grouped_df["oligo_seq"].str.upper()

# For each library, find missing oligos with the correct SPS and add them with 0 counts
for library, sps in library_sps_map.items():
    # Oligos with correct SPS for this library
    oligos_for_library = designed_oligos[designed_oligos["SPS"] == sps]
    # Oligos not already present for this library
    existing_oligos = oligo_name_grouped_df[
        oligo_name_grouped_df["library"] == library
    ]["oligo_name"]
    missing_oligos = oligos_for_library[
        ~oligos_for_library["oligo_name"].isin(existing_oligos)
    ]
    if not missing_oligos.empty:
        # Check if any missing oligos have the same oligo_seq as an oligo already present for this library
        existing_oligo_seqs = oligo_name_grouped_df[
            oligo_name_grouped_df["library"] == library
        ]["oligo_seq"]
        overlap = missing_oligos[missing_oligos["oligo_seq"].isin(existing_oligo_seqs)]
        if not overlap.empty:
            print(
                f"Warning: For library {library}, {len(overlap)} missing oligos have the same oligo_seq as an existing oligo."
            )
            print(overlap[["oligo_name", "oligo_seq"]])
        else:
            print(
                f"Library {library}: Adding {len(missing_oligos)} missing oligos with SPS {sps}"
            )
            # Add missing oligos with 0 counts
            to_append = missing_oligos.assign(
                library=library, grouped_guide_bc0_counts=0
            )[
                [
                    "library",
                    "oligo_name",
                    "grouped_guide_bc0_counts",
                    "oligo_seq",
                    "designed_donor",
                    "gc_content",
                    "SPS",
                ]
            ]
            oligo_name_grouped_df = pd.concat(
                [oligo_name_grouped_df, to_append], ignore_index=True
            )

oligo_name_grouped_df.query("grouped_guide_bc0_counts == 0").to_csv(
    OUT_DIR + "missing_oligos_with_zero_counts.tsv", sep="\t", index=False
)

# Rename grouped_guide_bc0_counts to grouped_oligo_name_counts
oligo_name_grouped_df.rename(
    columns={"grouped_guide_bc0_counts": "grouped_oligo_name_counts"}, inplace=True
)

# Extract the last 20 bases of the oligo_seq as SPS (Subpool Sequence) for each oligo.
oligo_name_grouped_df["SPS"] = oligo_name_grouped_df["oligo_seq"].str[-20:]

# Groupby library and SPS, and sum the grouped_oligo_name_counts and unique_bc0 counts.
sps_grouped_df = oligo_name_grouped_df.groupby(["library", "SPS"], as_index=False).agg(
    grouped_oligo_name_counts=("grouped_oligo_name_counts", "sum"),
    unique_bc0=("unique_bc0", "sum"),
)

# Only take the top SPS for each library based on the maximum grouped_oligo_name_counts.
sps_grouped_df = sps_grouped_df.sort_values(
    by=["library", "grouped_oligo_name_counts"], ascending=[True, False]
)
sps_grouped_df = sps_grouped_df.drop_duplicates(subset=["library"])

# Filter oligo_name_grouped_df to only include the top SPS for each library.
oligo_name_filtered_df = oligo_name_grouped_df.merge(
    sps_grouped_df[["library", "SPS"]], on=["library", "SPS"], how="inner"
)

# Write the oligo_name_filtered_df to a file.
oligo_name_filtered_df.to_csv(
    OUT_DIR + "all_libraries_oligo_name_filtered.tsv", sep="\t", index=False
)

# Also write only V536 and V625 libraries to a separate file.
v536_v625_df = oligo_name_filtered_df[
    oligo_name_filtered_df["library"].isin(["V536", "V625"])
]

v536_v625_df.to_csv(
    OUT_DIR + "V536_V625_libraries_oligo_name_filtered.tsv", sep="\t", index=False
)

# Also write all the NRD1 libraries to a separate file.
nrd1_libraries = ["V536", "V622", "V623", "V624", "V625"]
nrd1_df = oligo_name_filtered_df[oligo_name_filtered_df["library"].isin(nrd1_libraries)]
nrd1_df.to_csv(
    OUT_DIR + "NRD1_libraries_step_1_oligo_counts.tsv", sep="\t", index=False
)
