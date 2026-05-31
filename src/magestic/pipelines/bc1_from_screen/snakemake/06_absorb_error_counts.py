#!/usr/bin/env python3
"""
Step 03c: Absorb Error BC1 Counts into Parent BC1s

Takes the error BC1 mapping from 03b_error_detection.py and:
1. Adds error BC1 counts to their parent BC1
2. Removes error BC1 rows from the counts file
3. Outputs a cleaned long-format counts file

This consolidation should be done BEFORE tier classification so that
parent BC1s have the correct (higher) counts for tier assignment.

Input:
    - Combined long-format counts file (from 03_combine_counts.py)
    - Error BC1 mapping file (from 03b_error_detection.py)

Output:
    - Cleaned long-format counts file with error counts absorbed
    - Absorption log file

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple


def absorb_error_counts(
    counts_df: pd.DataFrame,
    error_mappings: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Absorb error BC1 counts into parent BC1s.

    Args:
        counts_df: Long-format counts (bc1, sample_name, counts)
        error_mappings: Error to parent mapping (error_bc1, parent_bc1, ...)

    Returns:
        (cleaned_counts_df, absorption_log_df)
    """
    if len(error_mappings) == 0:
        print("No error mappings provided, returning original counts")
        return counts_df.copy(), pd.DataFrame()

    # Create error -> parent mapping dictionary
    error_to_parent = dict(zip(error_mappings['error_bc1'], error_mappings['parent_bc1']))
    error_bc1_set = set(error_to_parent.keys())

    print(f"Error BC1s to absorb: {len(error_bc1_set):,}")

    # Split into error BC1 rows and non-error BC1 rows
    is_error = counts_df['bc1'].isin(error_bc1_set)
    error_rows = counts_df[is_error].copy()
    clean_rows = counts_df[~is_error].copy()

    print(f"Error BC1 count rows: {len(error_rows):,}")
    print(f"Non-error BC1 count rows: {len(clean_rows):,}")

    if len(error_rows) == 0:
        print("No error BC1 rows found in counts, returning original")
        return counts_df.copy(), pd.DataFrame()

    # Map error BC1s to their parents
    error_rows['parent_bc1'] = error_rows['bc1'].map(error_to_parent)

    # Group error rows by parent and sample
    error_contributions = error_rows.groupby(
        ['parent_bc1', 'sample_name']
    )['counts'].sum().reset_index()
    error_contributions.columns = ['bc1', 'sample_name', 'error_counts']

    print(f"Error contributions to aggregate: {len(error_contributions):,}")

    # Merge error contributions with clean rows
    merged = clean_rows.merge(
        error_contributions,
        on=['bc1', 'sample_name'],
        how='left'
    )
    merged['error_counts'] = merged['error_counts'].fillna(0).astype(int)
    merged['counts'] = merged['counts'] + merged['error_counts']

    # Also add new rows for parents that didn't exist in clean_rows for certain samples
    # (rare case where parent BC1 wasn't present in a sample but its error BC1 was)
    # Vectorized approach: use merge with indicator to find non-matching pairs
    check_merge = error_contributions.merge(
        clean_rows[['bc1', 'sample_name']].drop_duplicates(),
        on=['bc1', 'sample_name'],
        how='left',
        indicator=True
    )
    missing_parents = check_merge[check_merge['_merge'] == 'left_only'].copy()

    if len(missing_parents) > 0:
        print(f"New parent rows created (parent absent in sample): {len(missing_parents):,}")
        new_rows_df = missing_parents[['bc1', 'sample_name', 'error_counts']].copy()
        new_rows_df['counts'] = new_rows_df['error_counts']
        merged = pd.concat([merged, new_rows_df], ignore_index=True)

    # Create absorption log
    absorption_log = error_rows.groupby('parent_bc1').agg(
        n_error_bc1s=('bc1', 'nunique'),
        total_absorbed_counts=('counts', 'sum'),
        n_sample_contributions=('sample_name', 'count'),
        error_bc1s=('bc1', lambda x: ','.join(sorted(set(x))))
    ).reset_index()
    absorption_log.columns = ['parent_bc1', 'n_error_bc1s', 'total_absorbed_counts',
                              'n_sample_contributions', 'error_bc1s']

    # Select and order output columns
    result = merged[['bc1', 'sample_name', 'counts']].copy()

    # Summary statistics
    original_total = counts_df['counts'].sum()
    result_total = result['counts'].sum()
    absorbed_total = error_rows['counts'].sum()

    print(f"\nAbsorption complete:")
    print(f"  Original total counts: {original_total:,}")
    print(f"  Counts absorbed from errors: {absorbed_total:,}")
    print(f"  Result total counts: {result_total:,}")
    print(f"  Verification: {original_total} should equal {result_total} (diff: {abs(original_total - result_total)})")

    print(f"\nOriginal unique BC1s: {counts_df['bc1'].nunique():,}")
    print(f"Result unique BC1s: {result['bc1'].nunique():,}")
    print(f"BC1s removed (absorbed): {counts_df['bc1'].nunique() - result['bc1'].nunique():,}")

    return result, absorption_log


def main():
    parser = argparse.ArgumentParser(
        description="Absorb error BC1 counts into parent BC1s"
    )
    parser.add_argument(
        "--counts",
        type=Path,
        required=True,
        help="Combined long-format counts file from 03_combine_counts.py"
    )
    parser.add_argument(
        "--error-mappings",
        type=Path,
        required=True,
        help="Error BC1 mapping file from 03b_error_detection.py"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output path for cleaned counts file"
    )
    parser.add_argument(
        "--log-output",
        type=Path,
        help="Output path for absorption log (optional)"
    )

    args = parser.parse_args()

    # Load counts
    print(f"Loading counts from {args.counts}...")
    counts_df = pd.read_csv(args.counts, sep='\t')
    print(f"Loaded {len(counts_df):,} rows, {counts_df['bc1'].nunique():,} unique BC1s")

    # Load error mappings
    print(f"\nLoading error mappings from {args.error_mappings}...")
    error_mappings = pd.read_csv(args.error_mappings, sep='\t')
    print(f"Loaded {len(error_mappings):,} error BC1 mappings")

    # Absorb error counts
    print("\nAbsorbing error counts...")
    cleaned_counts, absorption_log = absorb_error_counts(counts_df, error_mappings)

    # Save cleaned counts
    print(f"\nSaving cleaned counts to {args.output}...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    cleaned_counts.to_csv(args.output, sep='\t', index=False)

    # Save absorption log if requested
    if args.log_output and len(absorption_log) > 0:
        print(f"Saving absorption log to {args.log_output}...")
        args.log_output.parent.mkdir(parents=True, exist_ok=True)
        absorption_log.to_csv(args.log_output, sep='\t', index=False)

    print("\nDone!")


if __name__ == "__main__":
    main()
