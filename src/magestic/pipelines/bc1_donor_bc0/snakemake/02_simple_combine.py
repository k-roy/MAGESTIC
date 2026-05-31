#!/usr/bin/env python3
"""
Simple combine script for bc1-donor-bc0 parsed files.

Combines all parsed TSV files from a directory into a single table,
aggregating counts by bc1+donor+sps+bc0.

Usage:
    python 02_simple_combine.py --parsed-dir <dir> --output <file>

Author: Kevin R. Roy
Date: 2026-03-22
"""

import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import os

def main():
    parser = argparse.ArgumentParser(description="Combine parsed bc1-donor-bc0 files")
    parser.add_argument("--parsed-dir", required=True, help="Directory with parsed TSV files")
    parser.add_argument("--output", required=True, help="Output combined TSV")
    parser.add_argument("--min-counts", type=int, default=2, help="Minimum counts to keep")
    
    args = parser.parse_args()
    
    parsed_dir = Path(args.parsed_dir)
    files = list(parsed_dir.glob("*_bc1_donor_bc0.tsv"))
    
    print(f"Found {len(files)} parsed files")
    
    all_data = []
    for f in tqdm(files, desc="Loading files"):
        try:
            df = pd.read_csv(f, sep="\t")
            if len(df) > 0:
                # Extract sample info from filename
                sample = f.stem.replace("_bc1_donor_bc0", "")
                df["sample"] = sample
                all_data.append(df)
        except Exception as e:
            print(f"Error loading {f}: {e}")
    
    if not all_data:
        print("No data loaded!")
        return
    
    # Combine all data
    combined = pd.concat(all_data, ignore_index=True)
    print(f"Combined: {len(combined)} rows from {len(all_data)} files")
    
    # Aggregate by bc1+donor+SPS+bc0
    # Use the column names from the parsed output
    group_cols = ["bc1", "donor", "SPS", "bc0"]
    available_cols = [c for c in group_cols if c in combined.columns]
    
    if "counts" in combined.columns:
        count_col = "counts"
    elif "count" in combined.columns:
        count_col = "count"
    else:
        # Create count column from row counts
        combined["counts"] = 1
        count_col = "counts"
    
    # Also include donor_bc0_fragment if present
    if "donor_bc0_fragment" in combined.columns:
        available_cols.append("donor_bc0_fragment")
    
    print(f"Aggregating by: {available_cols}")
    
    agg = combined.groupby(available_cols, as_index=False)[count_col].sum()
    agg = agg.rename(columns={count_col: "total_counts"})
    
    # Count number of samples per group
    sample_counts = combined.groupby(available_cols)["sample"].nunique().reset_index()
    sample_counts = sample_counts.rename(columns={"sample": "num_samples"})
    
    agg = agg.merge(sample_counts, on=available_cols)
    
    # Filter by minimum counts
    agg = agg[agg["total_counts"] >= args.min_counts]
    
    # Sort by total counts
    agg = agg.sort_values("total_counts", ascending=False)
    
    # Save
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    agg.to_csv(args.output, sep="\t", index=False)
    
    print(f"\nSaved {len(agg)} unique bc1-donor-bc0 combinations to {args.output}")
    print(f"Unique bc1s: {agg['bc1'].nunique()}")

if __name__ == "__main__":
    main()
