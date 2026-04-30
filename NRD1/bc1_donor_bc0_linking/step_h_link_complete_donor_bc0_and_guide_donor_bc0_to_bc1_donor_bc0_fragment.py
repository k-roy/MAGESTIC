import os
import pandas as pd

# =========================
# Define filepaths at top
# =========================

BASE_DIR = "/path/to/processed_data/by_project/"

BC1_DONOR_BC0_DIR = os.path.join(
    BASE_DIR, "NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/data_tables"
)
GUIDE_DONOR_BC0_DIR = os.path.join(
    BASE_DIR,
    "NNS/20250628_repeat_step_1_library_cloning/guide_donor_bc0_purity_filtered",
)
KEYFILE_DIR = os.path.join(
    "/path/to/scripts_and_keyfiles/"
    "by_project/NNS/20250628_repeat_step_1_library_cloning/keyfiles"
)

OLIGO_POOL_PATH = os.path.join(
    KEYFILE_DIR, "20240411_Twist_200mer_oligo_array_order.tsv"
)
NRD1_OLIGO_INFO_PATH = os.path.join(
    KEYFILE_DIR,
    "20240409_NRD1_mut_SpCas9_NGG_codon_change_guide_donor_200_bp_oligos_info.txt",
)
BC1_DONOR_BC0_PATH = os.path.join(BC1_DONOR_BC0_DIR, "top_bc1_donor_bc0_fragments.tsv")
DONOR_BC0_PATH = os.path.join(BC1_DONOR_BC0_DIR, "donor_bc0_to_oligo_mapping.tsv")
GUIDE_DONOR_BC0_PATH = os.path.join(
    GUIDE_DONOR_BC0_DIR, "V625_guide_donor_bc0_purity_filtered.tsv"
)
BC1_REFERENCE_TABLE_PATH = os.path.join(
    BC1_DONOR_BC0_DIR, "20250820_NRD1_bc1_reference_table.tsv"
)

# =========================
# Data loading
# =========================

oligo_pool_df = pd.read_csv(OLIGO_POOL_PATH, sep="\t")
NRD1_oligo_info_df = pd.read_csv(NRD1_OLIGO_INFO_PATH, sep="\t")
bc1_donor_bc0_df = pd.read_csv(BC1_DONOR_BC0_PATH, sep="\t")
donor_bc0_df = pd.read_csv(DONOR_BC0_PATH, sep="\t")
guide_donor_bc0_df = pd.read_csv(GUIDE_DONOR_BC0_PATH, sep="\t")

# =========================
# Data filtering
# =========================

filtered_bc1_donor_bc0_df = bc1_donor_bc0_df[
    bc1_donor_bc0_df["bc1_purity"] > 0.05
].copy()

filtered_bc1_donor_bc0_df
# Do value counts of value counts on bc1 to check for duplicates
bc1_counts = filtered_bc1_donor_bc0_df["bc1"].value_counts()
bc1_counts.value_counts()


# =========================
# Prefix columns for the guide_donor_bc0_df
# =========================

guide_donor_bc0_df = guide_donor_bc0_df.rename(
    columns=lambda x: "gdb_" + x if x != "donor_bc0_fragment" else x
)

# =========================
# Merge steps
# =========================

merged_df = pd.merge(
    filtered_bc1_donor_bc0_df,
    donor_bc0_df,
    how="left",
)

result_df = pd.merge(
    merged_df,
    guide_donor_bc0_df,
    how="left",
    on="donor_bc0_fragment",
)

# =========================
# Matching checks
# =========================

result_df["donor_bc0_df_matched"] = result_df["db_donor"].notnull()
result_df["guide_donor_bc0_df_matched"] = result_df["gdb_donor"].notnull()

# Merge result_df with NRD1 oligo info, first use db_oligo_name, then gdb_oligo_name if not found
result_df["db_oligo_name"] = result_df["db_oligo_name"].fillna(
    result_df["gdb_oligo_name"]
)
result_df["gdb_oligo_name"] = result_df["gdb_oligo_name"].fillna(
    result_df["db_oligo_name"]
)
result_df = pd.merge(
    result_df,
    NRD1_oligo_info_df,
    how="left",
    left_on="db_oligo_name",
    right_on="oligo_name",
)
# Now only add gdb_oligo_name if db_oligo_name is not found
result_df["gdb_oligo_name"] = result_df.apply(
    lambda row: (
        row["gdb_oligo_name"] if pd.isnull(row["oligo_name"]) else row["gdb_oligo_name"]
    ),
    axis=1,
)


# Make a column indicating if the bc1 is unique in the table or not
result_df["multiple_donor_bc0_fragments_passing_filter"] = ~result_df["bc1"].duplicated(
    keep=False
)

result_df.query("~multiple_donor_bc0_fragments_passing_filter")
result_df.query('bc1 == "GCGTCGTCAAAACATCTTAT"')

# Summarize multiple_donor_bc0_fragments_passing_filter
num_multiple_donor_bc0_fragments_passing_filter = result_df[
    "multiple_donor_bc0_fragments_passing_filter"
].sum()
print(
    f"Number of unique bc1 in the table: {num_multiple_donor_bc0_fragments_passing_filter}"
)

result_df
# Do value counts of value counts on bc1 to check for duplicates
bc1_counts = result_df["bc1"].value_counts()
bc1_counts.value_counts()
result_df["multiple_donor_bc0_fragments_passing_filter"].value_counts()

result_df.to_csv(BC1_REFERENCE_TABLE_PATH, sep="\t", index=False)

result_df.columns

# =========================
# Summary statistics
# =========================
num_bc1_with_donor_match = result_df["donor_bc0_df_matched"].sum()
print(
    f"Number of bc1 with a matching donor_bc0_fragment in donor_bc0_df: {num_bc1_with_donor_match}"
)

num_bc1_with_guide_match = result_df["guide_donor_bc0_df_matched"].sum()
print(
    f"Number of bc1 with a matching donor_bc0_fragment in guide_donor_bc0_df: {num_bc1_with_guide_match}"
)

num_bc1_with_both_matches = result_df[
    result_df["donor_bc0_df_matched"] & result_df["guide_donor_bc0_df_matched"]
].shape[0]
print(
    f"Number of bc1 with a matching donor_bc0_fragment in both donor_bc0_df and guide_donor_bc0_df: {num_bc1_with_both_matches}"
)

num_bc1_with_oligo_name = result_df["db_oligo_name"].notnull().sum()
print(f"Number of bc1 with a db_oligo_name: {num_bc1_with_oligo_name}")

num_bc1_with_guide_oligo_name = result_df["gdb_oligo_name"].notnull().sum()
print(f"Number of bc1 with a gdb_oligo_name: {num_bc1_with_guide_oligo_name}")

num_bc1_with_any_oligo_name = result_df[
    result_df["db_oligo_name"].notnull() | result_df["gdb_oligo_name"].notnull()
].shape[0]
print(
    f"Number of bc1 with either a db_oligo_name or gdb_oligo_name: {num_bc1_with_any_oligo_name}"
)

num_bc1_with_perfect_donor = result_df[
    result_df["db_perfect_donor"] | result_df["gdb_perfect_donor"]
].shape[0]
print(
    f"Number of bc1 with either db_perfect_donor or gdb_perfect_donor: {num_bc1_with_perfect_donor}"
)

num_bc1_with_perfect_middle_donor_region = result_df[
    result_df["db_perfect_middle_donor_region"]
    | result_df["gdb_perfect_middle_donor_region"]
].shape[0]
print(
    f"Number of bc1 with either db_perfect_middle_donor_region or gdb_perfect_middle_donor_region: {num_bc1_with_perfect_middle_donor_region}"
)
