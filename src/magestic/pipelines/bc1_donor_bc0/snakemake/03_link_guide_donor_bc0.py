#!/usr/bin/env python3
"""
Link bc1-donor-bc0 data with guide-donor-bc0 data.

This script matches bc1 entries to their corresponding guides by:
1. Loading guide-donor-bc0 data (with donor_bc0_fragment)
2. Creating a 40bp linking key from both datasets:
   - bc1_donor_bc0: donor[-30:] + bc0[:10]
   - guide_donor_bc0: donor_bc0_fragment[:30] + donor_bc0_fragment[-10:]
3. Adding guide and SPS information to bc1 reference

Usage:
    python 03_link_guide_donor_bc0.py \\
        --input <bc1_data.tsv> \\
        --output <output.tsv> \\
        --guide-dir <guide_donor_bc0_dir>

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm


def create_linking_key_bc1(row):
    """
    Create 40bp linking key from bc1-donor-bc0 data.
    Key = last 30bp of donor + first 10bp of bc0
    """
    donor = row.get("donor") or row.get("top_donor")
    bc0 = row.get("bc0") or row.get("top_bc0")

    if pd.isna(donor) or pd.isna(bc0):
        return None

    donor = str(donor).upper()
    bc0 = str(bc0).upper()

    if len(donor) < 30 or len(bc0) < 10:
        return None

    return donor[-30:] + bc0[:10]


def create_linking_key_guide(row):
    """
    Create 40bp linking key from guide-donor-bc0 data.
    Key = first 30bp of donor_bc0_fragment + last 10bp of donor_bc0_fragment
    (This skips the 20bp SPS in the middle)
    """
    fragment = row.get("donor_bc0_fragment")

    if pd.isna(fragment) or fragment == "NA":
        return None

    fragment = str(fragment).upper()

    if len(fragment) < 60:
        return None

    # donor_bc0_fragment structure: donor[-30:] + SPS[20] + bc0[:10] = 60bp
    # Extract: first 30bp (donor part) + last 10bp (bc0 part)
    return fragment[:30] + fragment[-10:]


def main():
    parser = argparse.ArgumentParser(description="Link bc1 to guide data")
    parser.add_argument("--input", required=True, help="Input bc1 collapsed TSV")
    parser.add_argument("--output", required=True, help="Output linked TSV")
    parser.add_argument("--guide-dir", help="Directory with guide-donor-bc0 files")

    args = parser.parse_args()

    # Load bc1 data
    print("Loading bc1-donor-bc0 data...")
    bc1_df = pd.read_csv(args.input, sep="\t")
    print(f"Loaded {len(bc1_df)} rows with {bc1_df['bc1'].nunique()} unique bc1s")

    if args.guide_dir and Path(args.guide_dir).exists():
        guide_dir = Path(args.guide_dir)
        print(f"\nLoading guide-donor-bc0 data from: {guide_dir}")

        # Find guide count files (both purity filtered and regular)
        guide_files = list(guide_dir.glob("*_guide_donor_bc0_counts.tsv"))
        if not guide_files:
            guide_files = list(guide_dir.glob("*_guide_donor_bc0_purity_filtered.tsv"))

        if guide_files:
            guide_dfs = []
            for f in tqdm(guide_files, desc="Loading guide data"):
                df = pd.read_csv(f, sep="\t")
                # Filter out rows with NA values in key columns
                df = df[df["donor_bc0_fragment"].notna() & (df["donor_bc0_fragment"] != "NA")]
                guide_dfs.append(df)

            guide_df = pd.concat(guide_dfs, ignore_index=True)
            print(f"Loaded {len(guide_df)} valid guide-donor-bc0 entries")

            # Create linking keys
            print("\nCreating linking keys...")
            bc1_df["linking_key"] = bc1_df.apply(create_linking_key_bc1, axis=1)
            guide_df["linking_key"] = guide_df.apply(create_linking_key_guide, axis=1)

            # Check how many have valid keys
            bc1_valid = bc1_df["linking_key"].notna().sum()
            guide_valid = guide_df["linking_key"].notna().sum()
            print(f"  BC1 entries with valid linking keys: {bc1_valid}/{len(bc1_df)}")
            print(f"  Guide entries with valid linking keys: {guide_valid}/{len(guide_df)}")

            # Deduplicate guide data on linking key (keep highest count)
            count_col = "total_counts" if "total_counts" in guide_df.columns else "counts"
            guide_df = guide_df.sort_values(count_col, ascending=False)
            guide_dedup = guide_df.drop_duplicates(subset=["linking_key"], keep="first").copy()
            print(f"  Unique guide linking keys: {len(guide_dedup)}")

            # Select columns to transfer from guide data
            guide_cols = ["linking_key"]
            for col in ["guide", "SPS", "library_ID"]:
                if col in guide_dedup.columns:
                    guide_cols.append(col)

            # Merge on linking key
            print("\nMerging bc1 and guide data...")
            merged = bc1_df.merge(
                guide_dedup[guide_cols],
                on="linking_key",
                how="left",
                suffixes=("", "_from_guide")
            )

            matched = merged["guide"].notna().sum() if "guide" in merged.columns else 0
            print(f"Matched {matched}/{len(merged)} ({100*matched/len(merged):.1f}%) entries to guide data")

            # Report by library
            if "library_ID" in merged.columns:
                print("\nMatches by library_ID:")
                print(merged["library_ID"].value_counts(dropna=False).head(10))

            # Drop the linking key column before saving
            merged = merged.drop(columns=["linking_key"])

            merged.to_csv(args.output, sep="\t", index=False)
            print(f"\nSaved linked data to: {args.output}")
            return

    # If no guide data or no match columns, just copy input to output
    print("No guide data available or matching columns not found. Copying input to output.")
    bc1_df.to_csv(args.output, sep="\t", index=False)
    print(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()
