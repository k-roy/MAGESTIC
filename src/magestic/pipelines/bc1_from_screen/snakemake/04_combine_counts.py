#!/usr/bin/env python3
"""
Step 04: Combine Parsed BC1 Counts into Long Format

Combines all individual sample BC1 count files into a single long-format TSV file.

Input:
    - Directory containing parsed BC1 count files (*_bc1_counts.tsv)

Output:
    - Combined long-format counts file with columns: bc1, sample_name, counts

Author: Kevin R. Roy
"""

import argparse
import os
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import List


def read_count_file(filepath: Path) -> pd.DataFrame:
    """Read a single BC1 count file."""
    try:
        df = pd.read_csv(filepath, sep='\t')
        return df
    except Exception as e:
        print(f"Warning: Could not read {filepath}: {e}")
        return pd.DataFrame()


def combine_count_files(
    input_dir: Path,
    output_file: Path,
    n_workers: int = 4
) -> None:
    """
    Combine all BC1 count files in a directory into a single long-format file.

    Args:
        input_dir: Directory containing *_bc1_counts.tsv files
        output_file: Output path for combined file
        n_workers: Number of parallel workers for reading files
    """
    # Find all count files
    count_files = sorted(input_dir.glob("*_bc1_counts.tsv"))

    if not count_files:
        raise ValueError(f"No BC1 count files found in {input_dir}")

    print(f"Found {len(count_files)} BC1 count files to combine")

    # Read files in parallel
    print(f"Reading files with {n_workers} workers...")
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        dfs = list(executor.map(read_count_file, count_files))

    # Filter out empty dataframes
    dfs = [df for df in dfs if len(df) > 0]

    if not dfs:
        raise ValueError("All count files were empty or could not be read")

    print(f"Successfully read {len(dfs)} files")

    # Combine all dataframes
    print("Concatenating dataframes...")
    combined = pd.concat(dfs, ignore_index=True)

    # Ensure expected columns
    expected_cols = ['bc1', 'sample_name', 'counts']
    if not all(col in combined.columns for col in expected_cols):
        raise ValueError(f"Expected columns {expected_cols}, got {list(combined.columns)}")

    # Select and order columns
    combined = combined[expected_cols]

    # Print summary stats
    n_unique_bc1 = combined['bc1'].nunique()
    n_unique_samples = combined['sample_name'].nunique()
    total_counts = combined['counts'].sum()

    print(f"\nCombined statistics:")
    print(f"  Total rows: {len(combined):,}")
    print(f"  Unique BC1s: {n_unique_bc1:,}")
    print(f"  Unique samples: {n_unique_samples:,}")
    print(f"  Total counts: {total_counts:,}")

    # Save to output
    print(f"\nSaving to {output_file}...")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, sep='\t', index=False)

    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description="Combine parsed BC1 count files into long format"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing *_bc1_counts.tsv files"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output path for combined long-format file"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel workers for reading files (default: 4)"
    )

    args = parser.parse_args()

    combine_count_files(
        input_dir=args.input_dir,
        output_file=args.output,
        n_workers=args.workers
    )


if __name__ == "__main__":
    main()
