import os
import pandas as pd
from tqdm import tqdm
from Levenshtein import distance
import seaborn as sns
import multiprocessing
import matplotlib.pyplot as plt

DATA_DIR = (
    "/path/to/processed_data/by_project/"
    "NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/"
)
DATA_TABLES_DIR = os.path.join(DATA_DIR, "data_tables")
PLOT_DIR = os.path.join(DATA_DIR, "plots")
os.makedirs(DATA_TABLES_DIR, exist_ok=True)
os.makedirs(PLOT_DIR, exist_ok=True)

FRAGMENT_COUNTS_DIR = os.path.join(DATA_DIR, "bc1_donor_bc0_fragment_counts")
KEYFILE_DIR = os.path.join(
    "/path/to/scripts_and_keyfiles/"
    "by_project/NNS/20250813_NRD1_MAGESTIC_bc0_donor_bc1/keyfiles"
)
COMBINED_KEY_FILE = os.path.join(KEYFILE_DIR, "bc1_donor_bc0_combined_key.tsv")
EDIT_DISTANCE_THRESHOLD = 3


def process_group(args):
    """
    Collapse similar sequences within a group by edit distance threshold.

    Args:
        args (tuple): Contains group DataFrame, group_cols, seq_col, count_col, edit_distance_threshold.

    Returns:
        list: List of dicts representing collapsed rows.
    """
    group, group_cols, seq_col, count_col, edit_distance_threshold = args
    group = group.copy()
    result_rows = []
    while not group.empty:
        top_idx = (
            group[count_col].sum(axis=1).idxmax()
            if isinstance(count_col, list)
            else group[count_col].idxmax()
        )
        top_seq = group.loc[top_idx, seq_col]
        group["edit_distance"] = group[seq_col].apply(lambda x: distance(x, top_seq))
        collapse_mask = group["edit_distance"] <= edit_distance_threshold
        if isinstance(count_col, list):
            collapsed_counts = {
                col: group.loc[collapse_mask, col].sum() for col in count_col
            }
        else:
            collapsed_counts = {count_col: group.loc[collapse_mask, count_col].sum()}
        result_rows.append(
            {
                **{col: group.iloc[0][col] for col in group_cols},
                seq_col: top_seq,
                **collapsed_counts,
            }
        )
        group = group.loc[~collapse_mask].drop(columns=["edit_distance"])
    return result_rows


def collapse_sequences(
    df,
    group_cols,
    seq_col,
    count_col,
    edit_distance_threshold,
    show_progress=True,
    n_jobs=1,
):
    """
    Collapse similar sequences within each group by edit distance threshold.

    Args:
        df (pd.DataFrame): Input DataFrame.
        group_cols (list): Columns to group by.
        seq_col (str): Sequence column name.
        count_col (str or list): Count column(s).
        edit_distance_threshold (int): Edit distance threshold.
        show_progress (bool): Show progress bar.
        n_jobs (int): Number of parallel jobs.

    Returns:
        pd.DataFrame: Collapsed DataFrame.
    """
    groups = [group for _, group in df.groupby(group_cols)]
    groups_iter = tqdm(groups, desc="Collapsing sequences") if show_progress else groups

    if n_jobs > 1:
        args_list = [
            (group, group_cols, seq_col, count_col, edit_distance_threshold)
            for group in groups
        ]
        result_rows = []
        if show_progress:
            with multiprocessing.Pool(n_jobs) as pool, tqdm(
                total=len(args_list), desc="Collapsing sequences"
            ) as pbar:
                for group_result in pool.imap_unordered(process_group, args_list):
                    result_rows.extend(group_result)
                    pbar.update(1)
        else:
            with multiprocessing.Pool(n_jobs) as pool:
                for group_result in pool.imap_unordered(process_group, args_list):
                    result_rows.extend(group_result)
    else:
        result_rows = []
        for group in groups_iter:
            result_rows.extend(
                process_group(
                    (group, group_cols, seq_col, count_col, edit_distance_threshold)
                )
            )

    result = pd.DataFrame(result_rows)
    return result


def plot_edit_distance_distribution(df, plot_dir):
    """
    Plot distribution of edit distances by polymerase.

    Args:
        df (pd.DataFrame): DataFrame with 'edit_distance' and 'polymerase'.
        plot_dir (str): Directory to save plot.
    """
    sns.set(style="whitegrid")
    plt.figure(figsize=(10, 6))
    sns.histplot(
        data=df,
        x="edit_distance",
        hue="polymerase",
        multiple="stack",
        bins=30,
        kde=True,
    )
    plt.title("Distribution of Edit Distances by Polymerase")
    plt.xlabel("Edit Distance")
    plt.ylabel("Frequency")
    plt.legend(title="Polymerase")
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "edit_distance_distribution.png"))
    plt.close()


def plot_bc1_purity_density(df, plot_dir):
    """
    Plot histogram of BC1 purity by polymerase, faceted by sample number.

    Args:
        df (pd.DataFrame): DataFrame with 'bc1_purity', 'polymerase', 'sample_number', 'counts'.
        plot_dir (str): Directory to save plot.
    """
    sns.set(style="whitegrid")
    bins = [i / 100 for i in range(1, 101)]  # 1% to 100% bins
    g = sns.FacetGrid(
        df.query("counts >= 10"),
        col="sample_number",
        hue="polymerase",
        col_wrap=4,
        sharex=True,
        sharey=True,
        height=4,
    )
    g.map(
        sns.histplot,
        "bc1_purity",
        bins=bins,
        stat="density",
        element="step",
        common_norm=False,
        alpha=0.7,
    )
    g.add_legend(title="Polymerase")
    g.set_axis_labels("BC1 Purity", "Density")
    g.fig.suptitle(
        "Histogram of BC1 Purity by Polymerase (Faceted by Sample Number)", y=1.02
    )
    plt.tight_layout()
    g.savefig(os.path.join(plot_dir, "bc1_purity_density_facet_plot.png"))
    plt.close()


def main():
    """
    Main analysis pipeline for collapsing and summarizing bc1-donor-bc0 fragment data.
    """
    combined_key = pd.read_csv(COMBINED_KEY_FILE, sep="\t")
    combined_key = combined_key[combined_key["description"] == "bc1-donor-bc0"]

    dfs = []
    for _, row in tqdm(
        combined_key.iterrows(),
        total=combined_key.shape[0],
        desc="Combining bc1 donor bc0 data",
    ):
        fname = (
            f"{row['sequencing_date']}_{row['outer_primers']}_{row['sample_number']}_"
            "bc1_donor_bc0_fragment_counts.tsv"
        )
        filepath = os.path.join(FRAGMENT_COUNTS_DIR, fname)
        if os.path.exists(filepath):
            df = pd.read_csv(filepath, sep="\t")
            df["sample_number"] = row["sample_number"]
            dfs.append(df)

    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df = combined_df.merge(combined_key, how="left")

    # Group and sum counts
    combined_df = (
        combined_df.groupby(
            ["sample_number", "polymerase", "bc1", "donor_bc0_fragment"]
        )["counts"]
        .sum()
        .reset_index()
    )

    num_unique_bc1s = combined_df["bc1"].nunique()
    print(f"Number of unique bc1s: {num_unique_bc1s}")

    # Collapse similar donor_bc0_fragment sequences within each sample_number, bc1
    collapsed_df = collapse_sequences(
        combined_df,
        group_cols=["sample_number", "bc1"],
        seq_col="donor_bc0_fragment",
        count_col="counts",
        edit_distance_threshold=EDIT_DISTANCE_THRESHOLD,
        show_progress=True,
        n_jobs=24,
    )
    collapsed_df["total_bc1_counts"] = collapsed_df.groupby(["sample_number", "bc1"])[
        "counts"
    ].transform("sum")
    collapsed_df["bc1_purity"] = (
        collapsed_df["counts"] / collapsed_df["total_bc1_counts"]
    )
    collapsed_df.sort_values(by=["counts"], ascending=False, inplace=True)
    collapsed_df = collapsed_df.merge(combined_key, how="left")

    # Top donor fragment per bc1 per sample for plotting purity
    top_bc1_donor_fragment_df = collapsed_df.drop_duplicates(
        subset=["sample_number", "bc1"], keep="first"
    )

    plot_bc1_purity_density(top_bc1_donor_fragment_df, PLOT_DIR)

    # Collapse across all samples for each bc1
    all_samples_df = (
        collapsed_df.groupby(["bc1", "donor_bc0_fragment"])["counts"]
        .sum()
        .reset_index()
    )
    all_samples_df["num_PCR_replicates"] = (
        collapsed_df.groupby(["bc1", "donor_bc0_fragment"])["sample_number"]
        .nunique()
        .reset_index(drop=True)
    )

    all_samples_df.sort_values(by=["num_PCR_replicates"], ascending=False, inplace=True)

    collapsed_donor_bc0_fragments = collapse_sequences(
        all_samples_df,
        group_cols=["bc1"],
        seq_col="donor_bc0_fragment",
        count_col=["counts", "num_PCR_replicates"],
        edit_distance_threshold=EDIT_DISTANCE_THRESHOLD,
        show_progress=True,
        n_jobs=24,
    )

    collapsed_donor_bc0_fragments.sort_values(
        by=["counts"], ascending=False, inplace=True
    )

    collapsed_donor_bc0_fragments.rename(
        columns={"counts": "total_bc1_donor_bc0_counts"}, inplace=True
    )

    top_donor_bc0_fragment_per_bc1 = collapsed_donor_bc0_fragments.drop_duplicates(
        subset=["bc1"], keep="first"
    ).reset_index(drop=True)

    top_donor_bc0_fragment_per_bc1.rename(
        columns={
            "donor_bc0_fragment": "top_bc1_donor_fragment",
            "total_bc1_donor_bc0_counts": "top_total_bc1_donor_bc0_counts",
            "num_PCR_replicates": "top_num_PCR_replicates",
        },
        inplace=True,
    )

    collapsed_donor_bc0_fragments = collapsed_donor_bc0_fragments.merge(
        top_donor_bc0_fragment_per_bc1, how="left"
    )

    num_PCR_replicates_counts = (
        collapsed_donor_bc0_fragments["num_PCR_replicates"].value_counts().reset_index()
    )
    num_PCR_replicates_counts.columns = ["num_PCR_replicates", "count"]
    print(num_PCR_replicates_counts)

    filtered_df = collapsed_donor_bc0_fragments[
        collapsed_donor_bc0_fragments["num_PCR_replicates"] >= 2
    ]

    filtered_df["edit_distance"] = filtered_df.apply(
        lambda row: distance(row["donor_bc0_fragment"], row["top_bc1_donor_fragment"]),
        axis=1,
    )

    num_bc1s_with_multiple_donor_fragments = (
        filtered_df.groupby("bc1")["donor_bc0_fragment"]
        .nunique()
        .reset_index()
        .query("donor_bc0_fragment > 1")
        .shape[0]
    )
    print(
        f"Number of bc1s with more than one donor fragment: {num_bc1s_with_multiple_donor_fragments}"
    )

    num_unique_bc1s = filtered_df["bc1"].nunique()
    print(f"Number of unique bc1s: {num_unique_bc1s}")

    num_unique_bc1s_before_filtering = collapsed_donor_bc0_fragments["bc1"].nunique()
    print(
        f"Number of unique bc1s before filtering on num_PCR_replicates: {num_unique_bc1s_before_filtering}"
    )

    bc1_counts = (
        filtered_df.groupby("bc1")["total_bc1_donor_bc0_counts"].sum().reset_index()
    )
    bc1_counts.sort_values(
        by="total_bc1_donor_bc0_counts", ascending=False, inplace=True
    )

    total_bc1_counts = bc1_counts["total_bc1_donor_bc0_counts"].sum()
    print(f"Total counts for all bc1s: {total_bc1_counts}")

    filtered_df["num_distinct_donor_fragments"] = filtered_df.groupby("bc1")[
        "donor_bc0_fragment"
    ].transform("nunique")

    filtered_df["total_bc1_counts"] = filtered_df.groupby("bc1")[
        "total_bc1_donor_bc0_counts"
    ].transform("sum")

    filtered_df["bc1_purity"] = (
        filtered_df["total_bc1_donor_bc0_counts"] / filtered_df["total_bc1_counts"]
    )

    filtered_df.to_csv(
        os.path.join(DATA_TABLES_DIR, "top_bc1_donor_bc0_fragments.tsv"),
        sep="\t",
        index=False,
    )


if __name__ == "__main__":
    main()
