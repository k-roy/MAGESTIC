import os
import pandas as pd
from tqdm import tqdm

# Define directories
BASE_DIR = "/path/to/processed_data/by_project/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1"
COMPLETE_DONOR_BC0_DIR = os.path.join(BASE_DIR, "complete_donor_bc0_counts")
DATA_TABLES_DIR = os.path.join(BASE_DIR, "data_tables")

KEYFILE_DIR = os.path.join(
    "/path/to/scripts_and_keyfiles/by_project/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1",
    "keyfiles",
)
COMBINED_KEY_FILE = os.path.join(KEYFILE_DIR, "bc1_donor_bc0_combined_key.tsv")

# Load and filter keyfile
combined_key = pd.read_csv(COMBINED_KEY_FILE, sep="\t")
combined_key = combined_key[combined_key["description"] == "donor-bc0"]

# Collect sample dataframes
dfs = []
for _, row in tqdm(
    combined_key.iterrows(),
    total=combined_key.shape[0],
    desc="Combining donor bc0 data",
):
    sample_number = row["sample_number"]
    filename = os.path.join(
        COMPLETE_DONOR_BC0_DIR,
        f"20250815_plate_1C_{sample_number}_complete_donor_bc0_counts.tsv",
    )
    if os.path.exists(filename):
        df = pd.read_csv(filename, sep="\t")
        df["sample_number"] = sample_number
        dfs.append(df)

# Concatenate all sample data
if dfs:
    combined_df = pd.concat(dfs, ignore_index=True)
else:
    combined_df = pd.DataFrame()

# Merge with keyfile
combined_df = combined_df.merge(combined_key, how="left")

# Drop unnecessary columns
drop_cols = ["sequencing_date", "inner_product_size_bp", "outer_product_size_bp"]
combined_df = combined_df.drop(columns=drop_cols, errors="ignore")

# Group and sum counts
grouping_cols = [
    "RE_site_found",
    "RE_site",
    "donor",
    "SPS",
    "bc0",
    "donor_bc0_fragment",
]
grouped_df = combined_df.groupby(grouping_cols, as_index=False).sum(numeric_only=True)

# Save result
output_file = os.path.join(DATA_TABLES_DIR, "combined_donor_bc0_counts.tsv")
grouped_df.to_csv(output_file, sep="\t", index=False)
