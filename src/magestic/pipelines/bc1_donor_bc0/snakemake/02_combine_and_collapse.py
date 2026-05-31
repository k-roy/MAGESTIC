#!/usr/bin/env python3
"""
Combine and collapse bc1-donor-bc0 data across samples for a single sequencing run.

This script:
1. Loads parsed bc1-donor-bc0 data from multiple samples
2. Joins trafo_ID from mapping keyfile (if provided) for sublibrary tracking
3. Aggregates counts across PCR replicates
4. Collapses similar sequences using edit distance thresholds
5. Calculates bc1 purity metrics
6. Filters for reproducible bc1s (found in multiple PCR replicates)

Usage:
    python 02_combine_and_collapse.py \\
        --parsed-dir <dir> \\
        --sample-key <file> \\
        --plate-key <file> \\
        --run-date <date> \\
        --output-collapsed <file> \\
        --output-purity <file> \\
        [--trafo-key <file>] \\
        [options]

trafo_ID Propagation (CRITICAL for yL437/yL442 and yL406):
    For libraries with multiple trafos pooled together, V634 bc1s cannot be
    distinguished between trafo_3 (SpG_NGG) and trafo_4 (SpCas9_NGG) from SPS alone.

    Solution: Use --trafo-key to provide sample→trafo mapping from individual trafo sequencing.
    The trafo_key file should be TSV with columns: sample_number, trafo_ID
    (and optionally: V_library, nuclease_PAM, yeast_sublibrary)

    This enables bc1s to be assigned to their correct yeast_sublibrary for DESeq2 analysis.

Author: Kevin Roy Lab
"""

import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from Levenshtein import distance
import multiprocessing


def process_group(args):
    """
    Collapse similar sequences within a group by edit distance threshold.

    Args:
        args: tuple of (group_df, group_cols, seq_col, count_cols, edit_distance_threshold, preserve_cols)

    Returns:
        list of dicts representing collapsed rows
    """
    # Unpack args - preserve_cols is optional for backward compatibility
    if len(args) == 6:
        group, group_cols, seq_col, count_cols, edit_distance_threshold, preserve_cols = args
    else:
        group, group_cols, seq_col, count_cols, edit_distance_threshold = args
        preserve_cols = []

    group = group.copy()
    result_rows = []

    while not group.empty:
        # Find sequence with highest count
        if isinstance(count_cols, list):
            top_idx = group[count_cols].sum(axis=1).idxmax()
        else:
            top_idx = group[count_cols].idxmax()

        top_seq = group.loc[top_idx, seq_col]

        # Calculate edit distance to top sequence
        group["edit_distance"] = group[seq_col].apply(
            lambda x: distance(x, top_seq) if pd.notna(x) else 999
        )
        collapse_mask = group["edit_distance"] <= edit_distance_threshold

        # Aggregate counts for similar sequences
        if isinstance(count_cols, list):
            collapsed_counts = {col: group.loc[collapse_mask, col].sum() for col in count_cols}
        else:
            collapsed_counts = {count_cols: group.loc[collapse_mask, count_cols].sum()}

        # Preserve additional columns (take from top-count row)
        preserved = {}
        for col in preserve_cols:
            if col in group.columns:
                preserved[col] = group.loc[top_idx, col]

        result_rows.append({
            **{col: group.iloc[0][col] for col in group_cols},
            seq_col: top_seq,
            **collapsed_counts,
            **preserved,
        })

        group = group.loc[~collapse_mask].drop(columns=["edit_distance"])

    return result_rows


def collapse_sequences(df, group_cols, seq_col, count_cols, edit_distance_threshold,
                       show_progress=True, n_jobs=1, preserve_cols=None):
    """
    Collapse similar sequences within each group based on edit distance threshold.

    Args:
        df: Input DataFrame
        group_cols: Columns to group by
        seq_col: Sequence column to collapse
        count_cols: Count column(s) to aggregate
        edit_distance_threshold: Maximum edit distance for collapsing
        show_progress: Show progress bar
        n_jobs: Number of parallel jobs
        preserve_cols: Additional columns to preserve (take value from top-count row)

    Returns:
        DataFrame with collapsed sequences
    """
    if preserve_cols is None:
        preserve_cols = []

    grouped = df.groupby(group_cols)
    n_groups = grouped.ngroups

    if n_jobs > 1:
        def generate_args():
            for _, group in grouped:
                yield (group, group_cols, seq_col, count_cols, edit_distance_threshold, preserve_cols)

        with multiprocessing.Pool(n_jobs) as pool:
            chunksize = max(1, n_groups // (n_jobs * 4))
            iterator = pool.imap_unordered(process_group, generate_args(), chunksize=chunksize)
            if show_progress:
                iterator = tqdm(iterator, total=n_groups, desc=f"Collapsing {seq_col}")
            result_rows = [row for group_result in iterator for row in group_result]
    else:
        groups = [group for _, group in grouped]
        groups_iter = tqdm(groups, desc=f"Collapsing {seq_col}") if show_progress else groups
        result_rows = [
            row
            for group in groups_iter
            for row in process_group(
                (group, group_cols, seq_col, count_cols, edit_distance_threshold, preserve_cols)
            )
        ]

    return pd.DataFrame(result_rows)


def main():
    parser = argparse.ArgumentParser(
        description="Combine and collapse bc1-donor-bc0 data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
trafo_ID Propagation:
    For libraries like yL437/yL442 where V634 is shared between trafo_3 (SpG_NGG) and
    trafo_4 (SpCas9_NGG), use --trafo-key to provide sample→trafo mapping.

    trafo_key format (TSV):
        sample_number   trafo_ID    V_library   nuclease_PAM
        1               trafo_1     V629        SpG_NGNG
        2               trafo_2     V631        SpG_NGNH
        ...

    This enables accurate yeast_sublibrary assignment for DESeq2 analysis.
        """
    )
    parser.add_argument("--parsed-dir", required=True, help="Directory with parsed TSV files")
    parser.add_argument("--sample-key", required=True, help="Sample key TSV file")
    parser.add_argument("--plate-key", required=True, help="Plate key TSV file")
    parser.add_argument("--run-date", required=True, help="Sequencing run date")
    parser.add_argument("--output-collapsed", required=True, help="Output collapsed TSV")
    parser.add_argument("--output-purity", required=True, help="Output purity summary TSV")
    parser.add_argument("--trafo-key", help="Optional trafo mapping TSV (sample_number → trafo_ID)")
    parser.add_argument("--donor-edit-distance", type=int, default=5)
    parser.add_argument("--bc1-edit-distance", type=int, default=1)
    parser.add_argument("--bc0-edit-distance", type=int, default=1)
    parser.add_argument("--min-pcr-replicates", type=int, default=2)
    parser.add_argument("--n-jobs", type=int, default=24)

    args = parser.parse_args()

    print(f"Processing run: {args.run_date}")

    # Load keyfiles
    plate_key = pd.read_csv(args.plate_key, sep="\t")
    sample_key = pd.read_csv(args.sample_key, sep="\t")

    # Check for trafo_ID in various sources:
    # 1. From sample_key (via combined_key merge) - preferred
    # 2. From separate --trafo-key file - for supplementation/override
    trafo_key = None
    has_trafo_in_sample_key = "trafo_ID" in combined_key.columns or "trafo_id" in combined_key.columns
    trafo_col_name = "trafo_ID" if "trafo_ID" in combined_key.columns else ("trafo_id" if "trafo_id" in combined_key.columns else None)

    if has_trafo_in_sample_key:
        print(f"trafo_ID found in sample_key (column: {trafo_col_name})")
        print(f"  trafo_ID distribution: {combined_key[trafo_col_name].value_counts().to_dict()}")

    # Load additional trafo_key if provided (for supplementation or when sample_key lacks trafo_ID)
    if args.trafo_key:
        trafo_key = pd.read_csv(args.trafo_key, sep="\t")
        # Ensure sample_number is the join key
        if "sample_number" not in trafo_key.columns:
            print(f"WARNING: trafo_key missing 'sample_number' column. Available: {list(trafo_key.columns)}")
            trafo_key = None
        else:
            print(f"Loaded supplementary trafo_key with {len(trafo_key)} entries")
            if "trafo_ID" in trafo_key.columns:
                print(f"  trafo_ID distribution: {trafo_key['trafo_ID'].value_counts().to_dict()}")

    # Filter for bc1-donor-bc0 amplicons
    plate_key = plate_key[plate_key["description"] == "bc1-donor-bc0"]

    # Create combined key
    combined_key = plate_key.merge(sample_key, how="cross")
    combined_key["sequencing_date"] = args.run_date

    print(f"Loading data from {len(combined_key)} sample files...")

    # Load all parsed files
    dfs = []
    for _, row in tqdm(combined_key.iterrows(), total=len(combined_key), desc="Loading samples"):
        fname = f"{args.run_date}_{row['outer_primers']}_{row['sample_number']}_parsed_bc1_donor_bc0.tsv"
        filepath = Path(args.parsed_dir) / fname

        if filepath.exists():
            df = pd.read_csv(filepath, sep="\t")
            df["sample_number"] = row["sample_number"]
            df["outer_primers"] = row["outer_primers"]

            # Join trafo_ID from sample_key (via combined_key) if available
            trafo_cols_to_join = ["trafo_ID", "trafo_id", "V_library", "nuclease_PAM", "yeast_sublibrary"]
            for col in trafo_cols_to_join:
                if col in row.index and pd.notna(row[col]):
                    # Normalize trafo_ID column name
                    target_col = "trafo_ID" if col == "trafo_id" else col
                    df[target_col] = row[col]

            # Supplement/override from separate trafo_key if provided
            if trafo_key is not None:
                sample_trafo = trafo_key[trafo_key["sample_number"] == row["sample_number"]]
                if len(sample_trafo) == 1:
                    for col in ["trafo_ID", "V_library", "nuclease_PAM", "yeast_sublibrary"]:
                        if col in sample_trafo.columns:
                            df[col] = sample_trafo[col].iloc[0]
                elif len(sample_trafo) > 1:
                    print(f"WARNING: Multiple trafo entries for sample {row['sample_number']}")

            dfs.append(df)

    if not dfs:
        print("No data files found!")
        return

    combined_df = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(combined_df)} total rows")

    # Report trafo_ID coverage if applicable
    if has_trafo and "trafo_ID" in combined_df.columns:
        trafo_coverage = combined_df["trafo_ID"].notna().sum() / len(combined_df) * 100
        print(f"trafo_ID coverage: {trafo_coverage:.1f}%")
        print(f"trafo_ID distribution: {combined_df['trafo_ID'].value_counts().to_dict()}")

    # Determine columns based on read type (2x300 has 'donor', 2x150 has 'donor_bc0_fragment')
    has_full_donor = "donor" in combined_df.columns
    seq_col = "donor" if has_full_donor else "donor_bc0_fragment"

    # Identify trafo columns to preserve through collapsing
    trafo_cols = ["trafo_ID", "V_library", "nuclease_PAM", "yeast_sublibrary"]
    trafo_cols_present = [col for col in trafo_cols if col in combined_df.columns]

    print(f"Using sequence column: {seq_col}")
    print(f"Unique bc1s before collapsing: {combined_df['bc1'].nunique()}")
    if trafo_cols_present:
        print(f"trafo columns to preserve: {trafo_cols_present}")

    # Aggregate counts per bc1-donor-bc0 combination within each sample
    group_cols = ["sample_number", "outer_primers", "bc1", seq_col]
    if "bc0" in combined_df.columns and has_full_donor:
        group_cols.append("bc0")

    combined_df = combined_df.groupby(group_cols)["counts"].sum().reset_index()

    # Step 1: Collapse similar donors within each sample
    print("\nStep 1: Collapsing similar donors within samples...")
    collapse_group_cols = ["sample_number", "outer_primers", "bc1"]
    if "bc0" in combined_df.columns and has_full_donor:
        collapse_group_cols.append("bc0")

    collapsed_donors = collapse_sequences(
        combined_df,
        group_cols=collapse_group_cols,
        seq_col=seq_col,
        count_cols="counts",
        edit_distance_threshold=args.donor_edit_distance,
        n_jobs=args.n_jobs,
        preserve_cols=trafo_cols_present
    )

    # Step 2: Collapse similar bc1s within each sample
    print("\nStep 2: Collapsing similar bc1s within samples...")
    bc1_group_cols = ["sample_number", "outer_primers", seq_col]
    if "bc0" in combined_df.columns and has_full_donor:
        bc1_group_cols.append("bc0")

    collapsed_bc1 = collapse_sequences(
        collapsed_donors,
        group_cols=bc1_group_cols,
        seq_col="bc1",
        count_cols="counts",
        edit_distance_threshold=args.bc1_edit_distance,
        n_jobs=args.n_jobs,
        preserve_cols=trafo_cols_present
    )

    # Step 3: Collapse similar bc0s within each sample (if present)
    if "bc0" in collapsed_bc1.columns and has_full_donor:
        print("\nStep 3: Collapsing similar bc0s within samples...")
        collapsed_bc0 = collapse_sequences(
            collapsed_bc1,
            group_cols=["sample_number", "outer_primers", "bc1", seq_col],
            seq_col="bc0",
            count_cols="counts",
            edit_distance_threshold=args.bc0_edit_distance,
            n_jobs=args.n_jobs,
            preserve_cols=trafo_cols_present
        )
    else:
        collapsed_bc0 = collapsed_bc1

    # Create PCR replicate ID
    collapsed_bc0["PCR_replicate_id"] = (
        collapsed_bc0["sample_number"].astype(str) + "_" + collapsed_bc0["outer_primers"]
    )

    # Aggregate across all samples
    print("\nAggregating across samples...")
    agg_group_cols = ["bc1", seq_col]
    if "bc0" in collapsed_bc0.columns and has_full_donor:
        agg_group_cols.append("bc0")

    # Build aggregation dict
    agg_dict = {
        "counts": "sum",
        "PCR_replicate_id": "nunique"
    }

    # Include trafo columns in aggregation if present (take first non-null value per bc1)
    # trafo_ID should be consistent within each bc1, so first() is appropriate
    trafo_cols_in_collapsed = [col for col in trafo_cols_present if col in collapsed_bc0.columns]
    for col in trafo_cols_in_collapsed:
        agg_dict[col] = "first"

    all_samples_df = collapsed_bc0.groupby(agg_group_cols).agg(agg_dict).reset_index()
    all_samples_df = all_samples_df.rename(columns={"PCR_replicate_id": "num_PCR_replicates"})

    if trafo_cols_in_collapsed:
        print(f"trafo columns propagated: {trafo_cols_in_collapsed}")

    # Step 4: Collapse similar sequences across all samples
    print("\nStep 4: Collapsing across all samples...")

    # Determine which trafo columns to preserve
    trafo_cols_to_preserve = [col for col in trafo_cols_in_collapsed if col in all_samples_df.columns]

    collapsed_all_donors = collapse_sequences(
        all_samples_df,
        group_cols=["bc1"] + (["bc0"] if "bc0" in agg_group_cols else []),
        seq_col=seq_col,
        count_cols=["counts", "num_PCR_replicates"],
        edit_distance_threshold=args.donor_edit_distance,
        n_jobs=args.n_jobs,
        preserve_cols=trafo_cols_to_preserve
    )

    collapsed_all_bc1 = collapse_sequences(
        collapsed_all_donors,
        group_cols=[seq_col] + (["bc0"] if "bc0" in agg_group_cols else []),
        seq_col="bc1",
        count_cols=["counts", "num_PCR_replicates"],
        edit_distance_threshold=args.bc1_edit_distance,
        n_jobs=args.n_jobs,
        preserve_cols=trafo_cols_to_preserve
    )

    if "bc0" in agg_group_cols:
        collapsed_final = collapse_sequences(
            collapsed_all_bc1,
            group_cols=["bc1", seq_col],
            seq_col="bc0",
            count_cols=["counts", "num_PCR_replicates"],
            edit_distance_threshold=args.bc0_edit_distance,
            n_jobs=args.n_jobs,
            preserve_cols=trafo_cols_to_preserve
        )
    else:
        collapsed_final = collapsed_all_bc1

    # Sort by counts
    collapsed_final = collapsed_final.sort_values("counts", ascending=False)

    # Identify top donor per bc1
    top_per_bc1 = collapsed_final.drop_duplicates(subset=["bc1"], keep="first").copy()
    top_per_bc1 = top_per_bc1.rename(columns={
        seq_col: f"top_{seq_col}",
        "counts": "top_counts",
        "num_PCR_replicates": "top_num_PCR_replicates"
    })
    if "bc0" in top_per_bc1.columns:
        top_per_bc1 = top_per_bc1.rename(columns={"bc0": "top_bc0"})

    # Merge top info back
    merge_cols = ["bc1", f"top_{seq_col}", "top_counts", "top_num_PCR_replicates"]
    if "top_bc0" in top_per_bc1.columns:
        merge_cols.append("top_bc0")
    collapsed_final = collapsed_final.merge(top_per_bc1[merge_cols], on="bc1", how="left")

    # Filter for reproducible bc1s
    filtered_df = collapsed_final[
        collapsed_final["num_PCR_replicates"] >= args.min_pcr_replicates
    ].copy()

    print(f"\nFiltering for >= {args.min_pcr_replicates} PCR replicates:")
    print(f"  Before: {collapsed_final['bc1'].nunique()} unique bc1s")
    print(f"  After: {filtered_df['bc1'].nunique()} unique bc1s")

    # Add purity metrics
    filtered_df["total_bc1_counts"] = filtered_df.groupby("bc1")["counts"].transform("sum")
    filtered_df["bc1_purity"] = filtered_df["counts"] / filtered_df["total_bc1_counts"]
    filtered_df["num_distinct_donors"] = filtered_df.groupby("bc1")[seq_col].transform("nunique")

    # Calculate edit distance to top
    top_col = f"top_{seq_col}"
    filtered_df["edit_distance_to_top"] = filtered_df.apply(
        lambda row: distance(str(row[seq_col]), str(row[top_col]))
        if pd.notna(row[seq_col]) and pd.notna(row[top_col]) else None,
        axis=1
    )

    # Save collapsed data
    os.makedirs(os.path.dirname(args.output_collapsed), exist_ok=True)
    filtered_df.to_csv(args.output_collapsed, sep="\t", index=False)
    print(f"\nSaved collapsed data to: {args.output_collapsed}")

    # Generate purity summary
    purity_summary = filtered_df.groupby("bc1").agg({
        "counts": "sum",
        "num_PCR_replicates": "max",
        "bc1_purity": "max",
        "num_distinct_donors": "first"
    }).reset_index()
    purity_summary.columns = ["bc1", "total_counts", "max_PCR_replicates", "max_purity", "num_donors"]

    purity_summary.to_csv(args.output_purity, sep="\t", index=False)
    print(f"Saved purity summary to: {args.output_purity}")

    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total unique bc1s: {filtered_df['bc1'].nunique()}")
    print(f"Total rows: {len(filtered_df)}")
    print(f"Total counts: {filtered_df['counts'].sum()}")
    print(f"BC1s with single donor: {(purity_summary['num_donors'] == 1).sum()}")
    print(f"BC1s with multiple donors: {(purity_summary['num_donors'] > 1).sum()}")
    print(f"Mean purity: {purity_summary['max_purity'].mean():.3f}")
    print(f"Median purity: {purity_summary['max_purity'].median():.3f}")

    # Report trafo_ID distribution if present
    if "trafo_ID" in filtered_df.columns:
        print("\ntrafo_ID distribution:")
        trafo_counts = filtered_df.groupby("trafo_ID")["bc1"].nunique()
        for trafo_id, count in trafo_counts.items():
            print(f"  {trafo_id}: {count} unique bc1s")


if __name__ == "__main__":
    main()
