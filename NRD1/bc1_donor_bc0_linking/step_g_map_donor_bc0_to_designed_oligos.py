"""
Created on Fri Jun 19 02:23:35 2020

@author: kevinroy

This script maps donor barcode fragments to designed oligo names using perfect
and partial sequence matches.
"""

import pandas as pd
import os
from tqdm import tqdm

tqdm.pandas()

# Constants
BASE_DIR = "/path/to/processed_data/by_project/"
BC1_DONOR_BC0_DIR = os.path.join(
    BASE_DIR, "NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/data_tables"
)
KEYFILE_DIR = (
    "/path/to/scripts_and_keyfiles/"
    "by_project/NNS/20250628_repeat_step_1_library_cloning/keyfiles/"
)
OLIGO_POOL_FILENAME = os.path.join(
    KEYFILE_DIR, "20240411_Twist_200mer_oligo_array_order.tsv"
)
SpCas9_CLONING_SITE = "GTTTGAAGAGC"
MIDDLE_DONOR_START_IDX = 30
MIDDLE_DONOR_END_IDX = 90
MIDDLE_DONOR_LENGTH = MIDDLE_DONOR_END_IDX - MIDDLE_DONOR_START_IDX

# Load data
donor_bc0_df = pd.read_csv(
    os.path.join(BC1_DONOR_BC0_DIR, "combined_donor_bc0_counts.tsv"), sep="\t"
)
oligo_pool_df = pd.read_csv(OLIGO_POOL_FILENAME, sep="\t")

# Extract guide and donor sequences
oligo_pool_df["guide"] = oligo_pool_df["oligo_seq"].str.slice(20, 40)
oligo_pool_df["donor"] = oligo_pool_df["oligo_seq"].str.slice(51, 180).str.upper()
oligo_pool_df["middle_donor"] = (
    oligo_pool_df["oligo_seq"].str.slice(81, 141).str.upper()
)

# Build lookup dictionaries
guide_to_oligo_names = {}
perfect_donor_to_oligo_names = {}
partial_donor_to_oligo_names = {}
guide_donor_to_oligo_names = {}

for _, row in oligo_pool_df.iterrows():
    guide = row["guide"]
    donor = row["donor"]
    oligo_name = row["oligo_name"]
    guide_donor = guide + SpCas9_CLONING_SITE + donor

    guide_donor_to_oligo_names[guide_donor] = oligo_name
    guide_to_oligo_names.setdefault(guide, []).append(oligo_name)
    perfect_donor_to_oligo_names.setdefault(donor, []).append(oligo_name)

    for idx in range(0, len(donor) - MIDDLE_DONOR_LENGTH):
        partial_donor = donor[idx : idx + MIDDLE_DONOR_LENGTH]
        partial_donor_to_oligo_names.setdefault(partial_donor, []).append(oligo_name)


def map_oligo(row):
    donor = row["donor"]
    donor_bc0_fragment = row["donor_bc0_fragment"]
    counts = row["counts"]

    perfect_donor = donor in perfect_donor_to_oligo_names
    middle_donor = donor[MIDDLE_DONOR_START_IDX:MIDDLE_DONOR_END_IDX]
    perfect_middle_donor_region = middle_donor in partial_donor_to_oligo_names

    oligo_name = "NA"
    if perfect_donor:
        oligo_names = perfect_donor_to_oligo_names[donor]
        if len(oligo_names) == 1:
            oligo_name = oligo_names[0]
    elif perfect_middle_donor_region:
        oligo_names = partial_donor_to_oligo_names[middle_donor]
        if len(oligo_names) == 1:
            oligo_name = oligo_names[0]

    return pd.Series(
        {
            "donor_bc0_fragment": donor_bc0_fragment,
            "db_counts": counts,
            "db_perfect_donor": perfect_donor,
            "db_perfect_middle_donor_region": perfect_middle_donor_region,
            "db_donor": donor,
            "db_oligo_name": oligo_name,
        }
    )


# Let's first remove duplicate donor_bc0_fragment, favoring the one with the highest counts
# In the case of ties, we will keep all occurrences and then use db_perfect_donor or
# db_perfect_middle_donor_region to determine the best match later.

# Sort donor_bc0_df by counts in descending order in place
donor_bc0_df.sort_values("counts", ascending=False, inplace=True)

# Deduplicate donor_bc0_fragment, keeping all rows tied for the highest counts (faster version)
idx = (
    donor_bc0_df.groupby("donor_bc0_fragment")["counts"].transform("max")
    == donor_bc0_df["counts"]
)
deduplicated_donor_bc0_df = donor_bc0_df[idx].reset_index(drop=True)

# Vectorized mapping for speed

# Precompute sets for fast lookup
perfect_donor_set = set(perfect_donor_to_oligo_names.keys())
partial_donor_set = set(partial_donor_to_oligo_names.keys())

donors = deduplicated_donor_bc0_df["donor"].values
donor_bc0_fragments = deduplicated_donor_bc0_df["donor_bc0_fragment"].values
counts = deduplicated_donor_bc0_df["counts"].values

middle_donors = [d[MIDDLE_DONOR_START_IDX:MIDDLE_DONOR_END_IDX] for d in donors]

db_perfect_donor = [d in perfect_donor_set for d in donors]
db_perfect_middle_donor_region = [md in partial_donor_set for md in middle_donors]

"CAACCGTATGGTTATGCCCCCAACCAACCATTACCTTCTCAGGGACCTGCCGCAGCGGCTCCACCAGTCCCTCAGCAACAATTTGACCCCACTGCTCAATTGAATTCTTTGATGAATATGCTTAACCAA" in perfect_donor_set

db_oligo_name = []
unambiguous_oligo_name_lst = []
for i in range(len(donors)):
    oligo_name = "NA"
    unambiguous_oligo_name = False
    if db_perfect_donor[i]:
        oligo_names = perfect_donor_to_oligo_names[donors[i]]
        if len(oligo_names) == 1:
            unambiguous_oligo_name = True
        oligo_name = oligo_names[0]
    elif db_perfect_middle_donor_region[i]:
        oligo_names = partial_donor_to_oligo_names[middle_donors[i]]
        if len(oligo_names) == 1:
            unambiguous_oligo_name = True
        oligo_name = oligo_names[0]
    db_oligo_name.append(oligo_name)
    unambiguous_oligo_name_lst.append(unambiguous_oligo_name)

result_df = pd.DataFrame(
    {
        "donor_bc0_fragment": donor_bc0_fragments,
        "db_counts": counts,
        "db_perfect_donor": db_perfect_donor,
        "db_perfect_middle_donor_region": db_perfect_middle_donor_region,
        "db_donor": donors,
        "db_oligo_name": db_oligo_name,
        "unambiguous_oligo_name": unambiguous_oligo_name_lst,
    }
)

# Remove duplicate rows based on 'donor_bc0_fragment',
# In case of ties based on total counts, favor those with 'db_perfect_donor' first
# and then 'db_perfect_middle_donor_region' second.
result_df = (
    result_df.sort_values(
        by=[
            "donor_bc0_fragment",
            "db_counts",
            "db_perfect_donor",
            "db_perfect_middle_donor_region",
        ],
        ascending=[True, False, False, False],
    )
    .drop_duplicates(subset=["donor_bc0_fragment"])
    .reset_index(drop=True)
)

# Sort by 'db_counts'
result_df = result_df.sort_values(by=["db_counts"], ascending=False)

# Check number of unique donor_bc0_fragments
num_unique_donor_bc0_fragments = result_df["donor_bc0_fragment"].nunique()
print(f"Number of unique donor_bc0_fragments: {num_unique_donor_bc0_fragments}")

output_path = os.path.join(BC1_DONOR_BC0_DIR, "donor_bc0_to_oligo_mapping.tsv")
result_df.to_csv(output_path, sep="\t", index=False)
