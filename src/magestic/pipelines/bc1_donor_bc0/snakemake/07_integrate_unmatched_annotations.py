#!/usr/bin/env python3
"""
Integrate unmatched donor annotations into bc1 reference table.

This script merges the annotations from 06_annotate_unmatched_donors.py back
into the bc1 reference table, adding key columns for downstream analysis:
- closest_oligo_name: For Tier 2 aggregation when middle_ed = 0
- closest_middle_donor_edit_distance: Critical for attribution (0 = attributable)
- Extended donor information
- codon_variant_type: Harmonized variant classification

The updated reference table can then be used by bc1_from_screen pipeline
for proper tier assignment and aggregation.

Usage:
    python 07_integrate_unmatched_annotations.py \\
        --bc1-reference <bc1_reference_table.tsv> \\
        --unmatched-annotations <unmatched_donor_annotations.tsv> \\
        --output <output_bc1_reference_table.tsv>

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Integrate unmatched donor annotations into bc1 reference table"
    )

    parser.add_argument("--bc1-reference", "-r", required=True, type=Path,
                        help="Input bc1 reference table (TSV)")
    parser.add_argument("--unmatched-annotations", "-u", required=True, type=Path,
                        help="Unmatched donor annotations from 06_annotate (TSV)")
    parser.add_argument("--output", "-o", required=True, type=Path,
                        help="Output bc1 reference table with integrated annotations")

    args = parser.parse_args()

    print("=" * 70)
    print("INTEGRATING UNMATCHED DONOR ANNOTATIONS")
    print("=" * 70)

    # Load bc1 reference table
    print(f"\n1. Loading bc1 reference table...")
    print(f"   {args.bc1_reference}")
    bc1_df = pd.read_csv(args.bc1_reference, sep="\t", low_memory=False)
    print(f"   Loaded {len(bc1_df):,} rows, {len(bc1_df.columns)} columns")

    # Load unmatched annotations
    print(f"\n2. Loading unmatched annotations...")
    print(f"   {args.unmatched_annotations}")
    unmatched_df = pd.read_csv(args.unmatched_annotations, sep="\t", low_memory=False)
    print(f"   Loaded {len(unmatched_df):,} rows")

    # Columns to integrate from unmatched annotations
    cols_to_integrate = [
        "bc1",  # Key for merging
        # New closest donor info
        "closest_oligo_name",
        "closest_donor_edit_distance",
        "closest_middle_donor_edit_distance",
        "inferred_guide",
        # Extended donor info
        "is_extended",
        "extension_non_bc0_side_bp",
        "extension_bc0_side_bp",
        "observed_donor_length",
        # Chimera info
        "is_chimera",
        "chimera_left_parent",
        "chimera_right_parent",
        # Confidence tier
        "confidence_tier",
        "rescue_status",
        # Alignment and variant info
        "alignment_chrom",
        "alignment_start",
        "alignment_end",
        "observed_individual_variants",
        "designed_individual_variants",
        "delta_variants",
        "delta_variant_count",
        "codon_variant_type",
    ]

    # Filter to only columns that exist
    cols_to_integrate = [c for c in cols_to_integrate if c in unmatched_df.columns]
    unmatched_subset = unmatched_df[cols_to_integrate].copy()

    print(f"\n3. Integrating annotations...")
    print(f"   Columns to integrate: {len(cols_to_integrate) - 1}")  # -1 for bc1

    # Mark which BC1s are from unmatched annotation
    unmatched_bc1s = set(unmatched_subset["bc1"])

    # For each column, we want to:
    # 1. Add new columns that don't exist
    # 2. Update existing columns only for unmatched BC1s
    for col in cols_to_integrate:
        if col == "bc1":
            continue

        if col not in bc1_df.columns:
            bc1_df[col] = pd.NA
            print(f"   Added new column: {col}")

    # Merge using update logic
    bc1_df_indexed = bc1_df.set_index("bc1")
    unmatched_indexed = unmatched_subset.set_index("bc1")

    # Update only the unmatched BC1s
    bc1_df_indexed.update(unmatched_indexed, overwrite=True)

    # Reset index
    bc1_df = bc1_df_indexed.reset_index()

    # Verify integration
    n_with_closest = bc1_df["closest_oligo_name"].notna().sum()
    n_middle_ed_0 = 0
    if "closest_middle_donor_edit_distance" in bc1_df.columns:
        n_middle_ed_0 = (bc1_df["closest_middle_donor_edit_distance"] == 0).sum()
    n_extended = bc1_df["is_extended"].sum() if "is_extended" in bc1_df.columns else 0

    print(f"\n4. Integration summary:")
    print(f"   BC1s with closest_oligo_name: {n_with_closest:,}")
    print(f"   BC1s with middle ED = 0 (attributable): {n_middle_ed_0:,}")
    print(f"   Extended donors: {n_extended:,}")

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save updated reference table
    bc1_df.to_csv(args.output, sep="\t", index=False)
    print(f"\n5. Saved updated reference table:")
    print(f"   {args.output}")
    print(f"   {len(bc1_df):,} rows, {len(bc1_df.columns)} columns")

    # Create a summary of attributable BC1s
    print("\n" + "=" * 70)
    print("ATTRIBUTION SUMMARY")
    print("=" * 70)

    # Detect oligo_name column
    oligo_col = None
    for col in ["gdb_oligo_name", "bdb_oligo_name", "oligo_name"]:
        if col in bc1_df.columns:
            oligo_col = col
            break

    if oligo_col:
        # Previously matched
        previously_matched = bc1_df[oligo_col].notna()
        n_previously_matched = previously_matched.sum()

        # Newly attributable (middle_ed = 0 and not previously matched)
        newly_attributable = (
            (bc1_df["closest_middle_donor_edit_distance"] == 0) &
            (~previously_matched)
        )
        n_newly_attributable = newly_attributable.sum()

        # Total attributable for Tier 2
        total_tier2 = n_previously_matched + n_newly_attributable

        print(f"\nFor Tier 2 (oligo_name) aggregation:")
        print(f"  Previously matched ({oligo_col}): {n_previously_matched:,}")
        print(f"  Newly attributable (middle_ed=0): {n_newly_attributable:,}")
        print(f"  Total Tier 2 attributable: {total_tier2:,} ({100*total_tier2/len(bc1_df):.1f}%)")
    else:
        print("\nNo oligo_name column found for attribution analysis.")

    print("\n" + "=" * 70)
    print("INTEGRATION COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
