#!/usr/bin/env python3
"""
Integrate bc1-donor-bc0 data from multiple sequencing runs.

This script combines data from 2x150 bp and 2x300 bp runs to maximize bc1 recovery:
- 2x300 bp runs capture full donor in merged reads
- 2x150 bp runs capture partial donor; full donor recovered via donor-bc0 linkage

Integration strategy:
1. Load collapsed data from each run
2. For 2x150 data, link bc1 to full donor via donor_bc0_fragment
3. Combine datasets
4. Collapse similar sequences across all runs
5. Apply quality filters

Usage:
    python integrate_multiple_runs.py \\
        --input-files <file1.tsv> <file2.tsv> ... \\
        --output <integrated.tsv> \\
        [options]

Author: Kevin Roy Lab
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from Levenshtein import distance
import multiprocessing


def process_group(args):
    """Collapse similar sequences within a group."""
    group, group_cols, seq_col, count_cols, edit_distance_threshold = args
    group = group.copy()
    result_rows = []

    while not group.empty:
        if isinstance(count_cols, list):
            top_idx = group[count_cols].sum(axis=1).idxmax()
        else:
            top_idx = group[count_cols].idxmax()

        top_seq = group.loc[top_idx, seq_col]
        group["edit_distance"] = group[seq_col].apply(
            lambda x: distance(x, top_seq) if pd.notna(x) else 999
        )
        collapse_mask = group["edit_distance"] <= edit_distance_threshold

        if isinstance(count_cols, list):
            collapsed_counts = {col: group.loc[collapse_mask, col].sum() for col in count_cols}
        else:
            collapsed_counts = {count_cols: group.loc[collapse_mask, count_cols].sum()}

        result_rows.append({
            **{col: group.iloc[0][col] for col in group_cols},
            seq_col: top_seq,
            **collapsed_counts,
        })

        group = group.loc[~collapse_mask].drop(columns=["edit_distance"])

    return result_rows


def collapse_sequences(df, group_cols, seq_col, count_cols, edit_distance_threshold, n_jobs=1):
    """Collapse similar sequences within each group."""
    grouped = df.groupby(group_cols)
    n_groups = grouped.ngroups

    def generate_args():
        for _, group in grouped:
            yield (group, group_cols, seq_col, count_cols, edit_distance_threshold)

    with multiprocessing.Pool(n_jobs) as pool:
        chunksize = max(1, n_groups // (n_jobs * 4))
        iterator = pool.imap_unordered(process_group, generate_args(), chunksize=chunksize)
        iterator = tqdm(iterator, total=n_groups, desc=f"Collapsing {seq_col}")
        result_rows = [row for group_result in iterator for row in group_result]

    return pd.DataFrame(result_rows)


def main():
    parser = argparse.ArgumentParser(description="Integrate bc1 data from multiple runs")
    parser.add_argument("--input-files", nargs="+", required=True, help="Input collapsed TSV files")
    parser.add_argument("--output", required=True, help="Output integrated TSV")
    parser.add_argument("--donor-edit-distance", type=int, default=5)
    parser.add_argument("--bc1-edit-distance", type=int, default=1)
    parser.add_argument("--min-pcr-replicates", type=int, default=2)
    parser.add_argument("--n-jobs", type=int, default=24)

    args = parser.parse_args()

    print("=" * 60)
    print("BC1-DONOR-BC0 INTEGRATION PIPELINE")
    print("=" * 60)

    # Load all input files
    all_dfs = []
    for f in args.input_files:
        print(f"Loading: {f}")
        df = pd.read_csv(f, sep="\t")
        df["source_file"] = Path(f).name
        all_dfs.append(df)
        print(f"  {len(df)} rows, {df['bc1'].nunique()} unique bc1s")

    combined = pd.concat(all_dfs, ignore_index=True)
    print(f"\nCombined: {len(combined)} total rows, {combined['bc1'].nunique()} unique bc1s")

    # Standardize columns
    if "donor" not in combined.columns and "top_donor" in combined.columns:
        combined["donor"] = combined["top_donor"]
    if "bc0" not in combined.columns and "top_bc0" in combined.columns:
        combined["bc0"] = combined["top_bc0"]

    # Filter for entries with full donor
    combined = combined[combined["donor"].notna()].copy()
    print(f"After filtering for full donor: {combined['bc1'].nunique()} unique bc1s")

    # Aggregate
    agg_df = combined.groupby(["bc1", "donor", "bc0"]).agg({
        "counts": "sum",
        "num_PCR_replicates": "sum",
        "source_file": lambda x: ",".join(sorted(set(x)))
    }).reset_index()

    print(f"\nAfter aggregation: {len(agg_df)} rows")

    # Collapse similar sequences
    print("\nCollapsing similar donors...")
    collapsed_donors = collapse_sequences(
        agg_df,
        group_cols=["bc1", "bc0"],
        seq_col="donor",
        count_cols=["counts", "num_PCR_replicates"],
        edit_distance_threshold=args.donor_edit_distance,
        n_jobs=args.n_jobs
    )

    print("\nCollapsing similar bc1s...")
    collapsed_bc1 = collapse_sequences(
        collapsed_donors,
        group_cols=["bc0", "donor"],
        seq_col="bc1",
        count_cols=["counts", "num_PCR_replicates"],
        edit_distance_threshold=args.bc1_edit_distance,
        n_jobs=args.n_jobs
    )

    print("\nCollapsing similar bc0s...")
    collapsed_bc0 = collapse_sequences(
        collapsed_bc1,
        group_cols=["bc1", "donor"],
        seq_col="bc0",
        count_cols=["counts", "num_PCR_replicates"],
        edit_distance_threshold=1,
        n_jobs=args.n_jobs
    )

    # Sort and get top per bc1
    collapsed_bc0 = collapsed_bc0.sort_values("counts", ascending=False)

    top_per_bc1 = collapsed_bc0.drop_duplicates(subset=["bc1"], keep="first").copy()
    top_per_bc1 = top_per_bc1.rename(columns={
        "donor": "top_donor",
        "bc0": "top_bc0",
        "counts": "top_counts",
        "num_PCR_replicates": "top_num_PCR_replicates"
    })

    result = collapsed_bc0.merge(
        top_per_bc1[["bc1", "top_donor", "top_bc0", "top_counts", "top_num_PCR_replicates"]],
        on="bc1",
        how="left"
    )

    # Filter
    filtered = result[result["num_PCR_replicates"] >= args.min_pcr_replicates].copy()
    print(f"\nAfter filtering for >= {args.min_pcr_replicates} PCR replicates:")
    print(f"  {filtered['bc1'].nunique()} unique bc1s")

    # Add metrics
    filtered["total_bc1_counts"] = filtered.groupby("bc1")["counts"].transform("sum")
    filtered["bc1_purity"] = filtered["counts"] / filtered["total_bc1_counts"]
    filtered["num_distinct_donors"] = filtered.groupby("bc1")["donor"].transform("nunique")

    # Save
    filtered.to_csv(args.output, sep="\t", index=False)
    print(f"\nSaved integrated data to: {args.output}")

    # Summary
    print("\n" + "=" * 60)
    print("INTEGRATION SUMMARY")
    print("=" * 60)
    print(f"Total unique bc1s: {filtered['bc1'].nunique()}")
    print(f"Total rows: {len(filtered)}")
    print(f"Total counts: {filtered['counts'].sum()}")


if __name__ == "__main__":
    main()
