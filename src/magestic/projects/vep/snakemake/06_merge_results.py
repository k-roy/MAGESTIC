#!/usr/bin/env python3
"""
Merge results from all VEP models.

This script combines predictions from Evo2, Yorzoi, Shorkie, and ESM1-v
into a single comprehensive prediction table.

Usage:
    python 06_merge_results.py \
        --evo2 evo2_scores.tsv \
        --yorzoi yorzoi_predictions.tsv \
        --shorkie shorkie_predictions.tsv \
        --esm1v esm1v_scores.tsv \
        --output merged_predictions.tsv

Author: Kevin R. Roy
Date: 2026-02-20
"""

import argparse
import sys
from pathlib import Path
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge VEP predictions from all models"
    )
    parser.add_argument(
        "--evo2",
        type=Path,
        help="Evo2 scores TSV file"
    )
    parser.add_argument(
        "--yorzoi",
        type=Path,
        help="Yorzoi predictions TSV file"
    )
    parser.add_argument(
        "--shorkie",
        type=Path,
        help="Shorkie predictions TSV file"
    )
    parser.add_argument(
        "--esm1v",
        type=Path,
        help="ESM1-v scores TSV file"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Metadata TSV file to merge with results"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output merged TSV file"
    )
    parser.add_argument(
        "--join-on",
        default="test_id",
        help="Column to join on (default: test_id)"
    )
    return parser.parse_args()


def load_and_prefix(filepath: Path, prefix: str, join_col: str) -> pd.DataFrame:
    """Load TSV and prefix non-join columns with model name."""
    if not filepath or not filepath.exists():
        return None

    df = pd.read_csv(filepath, sep='\t')

    # Rename columns except the join column
    rename_map = {}
    for col in df.columns:
        if col != join_col:
            rename_map[col] = f"{prefix}_{col}"

    return df.rename(columns=rename_map)


def main():
    args = parse_args()

    print("=== Merging VEP Results ===")

    dfs_to_merge = []

    # Load each model's results
    models = [
        (args.evo2, "evo2"),
        (args.yorzoi, "yorzoi"),
        (args.shorkie, "shorkie"),
        (args.esm1v, "esm1v"),
    ]

    for filepath, prefix in models:
        if filepath and filepath.exists():
            df = load_and_prefix(filepath, prefix, args.join_on)
            if df is not None:
                print(f"  {prefix}: {len(df)} entries")
                dfs_to_merge.append(df)
        else:
            print(f"  {prefix}: not provided or not found")

    if not dfs_to_merge:
        print("ERROR: No model results found to merge")
        sys.exit(1)

    # Start with first dataframe
    merged = dfs_to_merge[0]

    # Merge remaining dataframes
    for df in dfs_to_merge[1:]:
        merged = merged.merge(df, on=args.join_on, how="outer")

    print(f"\nMerged: {len(merged)} entries")

    # Add metadata if provided
    if args.metadata and args.metadata.exists():
        metadata = pd.read_csv(args.metadata, sep='\t')
        if args.join_on in metadata.columns:
            merged = merged.merge(metadata, on=args.join_on, how="left")
            print(f"Added metadata columns: {list(metadata.columns)}")

    # Reorder columns: join_col first, then metadata, then model results
    cols = list(merged.columns)
    # Move join column to front
    cols.remove(args.join_on)
    cols = [args.join_on] + cols
    merged = merged[cols]

    # Save merged results
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep='\t', index=False)

    print(f"\nOutput: {args.output}")
    print(f"Columns: {list(merged.columns)}")

    # Summary statistics
    print("\n=== Summary ===")
    for filepath, prefix in models:
        score_col = f"{prefix}_score"
        if score_col in merged.columns:
            valid = merged[score_col].notna().sum()
            print(f"  {prefix}: {valid} scored ({100*valid/len(merged):.1f}%)")


if __name__ == "__main__":
    main()
