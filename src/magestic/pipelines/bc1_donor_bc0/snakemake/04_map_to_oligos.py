#!/usr/bin/env python3
"""
Map bc1-donor-bc0 data to designed oligo pools.

This script matches donor sequences to the designed oligo pool to:
1. Identify which designed oligo each bc1 is associated with
2. Extract guide and full donor information
3. Calculate match quality metrics

BC1 Donor Structure:
    The bc1 "donor" column from parsing includes:
    - ~25bp prefix before RE site
    - 11bp RE site (GTTTGAAGAGC)
    - 129bp actual designed donor

    This script extracts the 129bp designed donor portion for matching.

Usage:
    # Using harmonized oligo design (recommended - has pre-extracted donors)
    python 04_map_to_oligos.py \\
        --input <linked_tsv> [linked_tsv ...] \\
        --oligo-pools <harmonized_oligos.tsv> \\
        --output <output_tsv>

    # Using raw Twist oligo files (will extract donors)
    python 04_map_to_oligos.py \\
        --input <linked_tsv> [linked_tsv ...] \\
        --oligo-pools <2024_twist.tsv> <2021_twist.tsv> \\
        --output <output_tsv> \\
        --donor-start 51 --donor-end 180

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
from pathlib import Path
from tqdm import tqdm

# RE sites used to extract the actual donor from bc1 data
# HANDOFF_20260527 forward-port: LbCas12a/impLbCas12a donors use a distinct RE site
# that the SpCas9/SpG-only check missed → all V450/V454 donors returned None.
RE_SITE_SPCAS9_SPG = "GTTTGAAGAGC"   # SpCas9 / SpG
RE_SITE_LBCAS12A   = "TTTCGAAGAGC"   # LbCas12a / impLbCas12a
RE_SITES = [RE_SITE_SPCAS9_SPG, RE_SITE_LBCAS12A]
DESIGNED_DONOR_LENGTH = 129


def extract_actual_donor_from_bc1(donor_seq):
    """
    Extract the designed donor portion from bc1 donor sequence.

    The bc1 donor includes prefix + RE_site + designed_donor.
    Tries SpCas9/SpG RE site first, then LbCas12a RE site.
    """
    if pd.isna(donor_seq):
        return None
    upper = str(donor_seq).upper()
    for re_site in RE_SITES:
        if re_site in upper:
            pos = upper.find(re_site)
            actual_donor = upper[pos + len(re_site):]
            return actual_donor[:DESIGNED_DONOR_LENGTH] if len(actual_donor) >= DESIGNED_DONOR_LENGTH else actual_donor
    return None


def load_oligo_pools(pool_files, donor_start=51, donor_end=180):
    """
    Load and combine oligo pool design files.

    Handles multiple formats:
    - Harmonized design (has pre-extracted 'donor' column)
    - Raw Twist files (has 'oligo_seq' or 'oligo_sequence', extracts donor)
    """
    pools = []
    for f in pool_files:
        print(f"  Loading: {Path(f).name}")
        df = pd.read_csv(f, sep="\t", low_memory=False)
        df["pool_file"] = Path(f).name

        # Normalize oligo name column
        if "oligo_name" in df.columns and "Name" not in df.columns:
            df["Name"] = df["oligo_name"]

        # Check if donor is already extracted (harmonized format)
        if "donor" in df.columns:
            df["designed_donor"] = df["donor"].str.upper()
            print(f"    Using pre-extracted donor column ({df['designed_donor'].notna().sum()} donors)")
        else:
            # Extract donor from oligo sequence
            seq_col = None
            for col in ["Sequence", "sequence", "oligo_seq", "oligo_sequence"]:
                if col in df.columns:
                    seq_col = col
                    break

            if seq_col:
                df["designed_donor"] = df[seq_col].apply(
                    lambda x: str(x)[donor_start:donor_end].upper() if pd.notna(x) and len(str(x)) >= donor_end else None
                )
                print(f"    Extracted donor from {seq_col} positions [{donor_start}:{donor_end}]")
            else:
                print(f"    Warning: No sequence column found")
                df["designed_donor"] = None

        # Extract guide if available
        if "guide" in df.columns:
            df["designed_guide"] = df["guide"].str.upper()
        else:
            seq_col = None
            for col in ["Sequence", "sequence", "oligo_seq", "oligo_sequence"]:
                if col in df.columns:
                    seq_col = col
                    break
            if seq_col:
                df["designed_guide"] = df[seq_col].apply(
                    lambda x: str(x)[20:40].upper() if pd.notna(x) and len(str(x)) >= 40 else None
                )

        pools.append(df)

    combined = pd.concat(pools, ignore_index=True)
    print(f"Total: {len(combined)} oligos from {len(pool_files)} pool files")
    return combined


def main():
    parser = argparse.ArgumentParser(
        description="Map donors to designed oligos",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Using harmonized oligo design (recommended)
    python 04_map_to_oligos.py \\
        --input bc1_to_guide_donor_bc0.tsv \\
        --oligo-pools harmonized_oligos.tsv \\
        --output bc1_mapped_to_oligos.tsv
        """
    )
    parser.add_argument("--input", nargs="+", required=True, help="Input TSV file(s) from Step 03")
    parser.add_argument("--oligo-pools", nargs="+", required=True, help="Oligo pool design TSV file(s)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--guide-start", type=int, default=20, help="Guide start position in raw oligo")
    parser.add_argument("--guide-end", type=int, default=40, help="Guide end position in raw oligo")
    parser.add_argument("--donor-start", type=int, default=51, help="Donor start position in raw oligo")
    parser.add_argument("--donor-end", type=int, default=180, help="Donor end position in raw oligo")

    args = parser.parse_args()

    # Load input data
    print("Loading bc1-donor-bc0 data...")
    dfs = [pd.read_csv(f, sep="\t") for f in args.input]
    bc1_data = pd.concat(dfs, ignore_index=True)
    print(f"Loaded {len(bc1_data)} rows with {bc1_data['bc1'].nunique()} unique bc1s")

    # Load oligo pools
    print("\nLoading oligo pool designs...")
    oligo_pool = load_oligo_pools(args.oligo_pools, args.donor_start, args.donor_end)

    # Create donor lookup (uppercase for case-insensitive matching)
    donor_lookup = oligo_pool[["designed_donor", "designed_guide", "Name", "pool_file"]].dropna(subset=["designed_donor"])
    donor_lookup = donor_lookup.drop_duplicates(subset=["designed_donor"])
    print(f"\nCreated donor lookup with {len(donor_lookup)} unique designed donors")

    # Determine which donor column to use in bc1 data
    if "donor" in bc1_data.columns:
        donor_col = "donor"
    elif "top_donor" in bc1_data.columns:
        donor_col = "top_donor"
    else:
        donor_col = None
        print("Warning: No donor column found. Skipping oligo mapping.")

    if donor_col:
        print(f"\nExtracting actual donor from bc1 data (after RE site)...")
        bc1_data["actual_donor"] = bc1_data[donor_col].apply(extract_actual_donor_from_bc1)
        valid_donors = bc1_data["actual_donor"].notna().sum()
        print(f"Extracted actual donor for {valid_donors}/{len(bc1_data)} ({100*valid_donors/len(bc1_data):.1f}%) entries")

        print(f"\nMatching to designed oligos...")

        # Merge on actual_donor
        merged = bc1_data.merge(
            donor_lookup,
            left_on="actual_donor",
            right_on="designed_donor",
            how="left"
        )

        matched = merged["designed_donor"].notna().sum()
        total = len(merged)
        print(f"Matched {matched}/{total} ({100*matched/total:.1f}%) entries to designed oligos")

        # Rename columns for clarity
        merged = merged.rename(columns={
            "Name": "oligo_name",
            "designed_guide": "oligo_guide",
            "pool_file": "oligo_pool"
        })

        # Create matching quality columns
        merged["perfect_donor_match"] = merged["designed_donor"].notna()

        # Clean up temporary columns
        merged = merged.drop(columns=["actual_donor", "designed_donor"], errors="ignore")

        # Summary by oligo pool
        if "oligo_pool" in merged.columns:
            print("\nMatches by oligo pool:")
            print(merged["oligo_pool"].value_counts(dropna=False).to_string())

        # Summary by SPS/V_library if available
        if "SPS" in merged.columns:
            print("\nMatches by SPS (top 5):")
            for sps in merged["SPS"].value_counts().head(5).index:
                sps_df = merged[merged["SPS"] == sps]
                sps_matched = sps_df["perfect_donor_match"].sum()
                print(f"  {sps[:20]}: {sps_matched}/{len(sps_df)} ({100*sps_matched/len(sps_df):.1f}%)")

        # Save output
        merged.to_csv(args.output, sep="\t", index=False)
        print(f"\nSaved mapped data to: {args.output}")
        print(f"Output columns: {', '.join(merged.columns[:10])}...")
    else:
        # Just pass through the data
        bc1_data.to_csv(args.output, sep="\t", index=False)
        print(f"\nSaved data (without oligo mapping) to: {args.output}")


if __name__ == "__main__":
    main()
