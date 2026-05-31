#!/usr/bin/env python3
"""
Step 05: Create DESeq2 Count Matrix with T0 Consolidation

Creates a wide-format count matrix combining abundance_high and abundance_medium BC1s
for a library, with configurable t0 sample consolidation strategies.

NOTE: Abundance tiers are distinct from confidence tiers in bc1_donor_bc0_pipeline.

T0 Consolidation Modes:
- screen_date: Pool all t0 samples from the same screen_date (yL406, yL437)
- pseudo_replicates: Create N pseudo-replicates from t0 samples (yL442)
- none: No consolidation

Usage:
    # Screen-date consolidation (yL437)
    python 05_create_count_matrix.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --tier-mapping tier_mapping.tsv \\
        --sample-key sample_key.tsv \\
        --library yL437 \\
        --t0-consolidation-mode screen_date \\
        --output count_matrix.tsv

    # Pseudo-replicate consolidation (yL442)
    python 05_create_count_matrix.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --tier-mapping tier_mapping.tsv \\
        --sample-key sample_key.tsv \\
        --library yL442 \\
        --t0-consolidation-mode pseudo_replicates \\
        --t0-pseudo-rep-count 5 \\
        --output count_matrix.tsv

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def consolidate_t0_by_screen_date(
    wide_matrix: pd.DataFrame,
    sample_key: pd.DataFrame,
    library: str,
    sample_name_col: str = "sample_name",
    generation_col: str = "generations",
    screen_date_col: str = "screen_date",
    library_col: str = "library_ID"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Consolidate t0 samples by screen_date.

    All t0 samples from the same screen_date are pooled together.
    This is appropriate when all conditions on the same date share
    the same starting population (e.g., yL406, yL437).

    Returns:
        (consolidated_matrix, metadata_df)
    """
    # Filter sample_key to library
    if library_col in sample_key.columns:
        lib_key = sample_key[sample_key[library_col] == library].copy()
    else:
        lib_key = sample_key.copy()

    # Identify t0 samples (generation == 0)
    t0_mask = lib_key[generation_col].astype(str).str.strip() == "0"
    t0_key = lib_key[t0_mask]

    if len(t0_key) == 0:
        logger.warning("No t0 samples found in sample_key")
        return wide_matrix, pd.DataFrame()

    # Get unique screen dates
    screen_dates = t0_key[screen_date_col].unique()
    logger.info(f"  Found {len(t0_key)} t0 samples across {len(screen_dates)} screen dates")

    # Create consolidated t0 columns (one per screen_date)
    consolidated_cols = {}
    metadata_rows = []

    for screen_date in screen_dates:
        date_mask = t0_key[screen_date_col] == screen_date
        date_samples = t0_key[date_mask][sample_name_col].tolist()

        # Filter to samples that exist in matrix
        date_samples = [s for s in date_samples if s in wide_matrix.columns]

        if len(date_samples) == 0:
            logger.warning(f"  No t0 samples found in matrix for screen_date {screen_date}")
            continue

        consolidated_name = f"t0_{screen_date}_consolidated"
        consolidated_cols[consolidated_name] = wide_matrix[date_samples].sum(axis=1)

        metadata_rows.append({
            "sample_name": consolidated_name,
            "screen_date": screen_date,
            "generation": 0,
            "n_samples_consolidated": len(date_samples),
            "source_samples": ",".join(date_samples)
        })

        logger.info(f"  Screen date {screen_date}: {len(date_samples)} t0 samples -> {consolidated_name}")

    # Remove original t0 samples, add consolidated
    all_t0_samples = t0_key[sample_name_col].tolist()
    non_t0_cols = [c for c in wide_matrix.columns if c not in all_t0_samples]

    result = pd.concat([
        wide_matrix[non_t0_cols],
        pd.DataFrame(consolidated_cols, index=wide_matrix.index)
    ], axis=1)

    return result, pd.DataFrame(metadata_rows)


def consolidate_t0_as_pseudo_replicates(
    wide_matrix: pd.DataFrame,
    sample_key: pd.DataFrame,
    library: str,
    n_pseudo_reps: int,
    seed: int = 42,
    sample_name_col: str = "sample_name",
    generation_col: str = "generations",
    library_col: str = "library_ID"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create pseudo-replicates from t0 samples to match treatment replicate counts.

    This is appropriate when all t0 samples are from the same screen_date
    but there are more t0 samples than condition replicates (e.g., yL442).

    Args:
        wide_matrix: BC1 x samples count matrix
        sample_key: Sample metadata
        library: Library name
        n_pseudo_reps: Number of pseudo-replicates to create
        seed: Random seed for reproducible sample assignment

    Returns:
        (consolidated_matrix, metadata_df)
    """
    # Filter sample_key to library
    if library_col in sample_key.columns:
        lib_key = sample_key[sample_key[library_col] == library].copy()
    else:
        lib_key = sample_key.copy()

    # Identify t0 samples (generation == 0)
    t0_mask = lib_key[generation_col].astype(str).str.strip() == "0"
    t0_key = lib_key[t0_mask]

    t0_samples = t0_key[sample_name_col].tolist()
    # Filter to samples that exist in matrix
    t0_samples = [s for s in t0_samples if s in wide_matrix.columns]

    if len(t0_samples) == 0:
        logger.warning("No t0 samples found")
        return wide_matrix, pd.DataFrame()

    logger.info(f"  Creating {n_pseudo_reps} pseudo-replicates from {len(t0_samples)} t0 samples")

    # Shuffle and partition
    np.random.seed(seed)
    shuffled = np.random.permutation(t0_samples)

    base_size = len(t0_samples) // n_pseudo_reps
    remainder = len(t0_samples) % n_pseudo_reps

    groups = []
    start = 0
    for i in range(n_pseudo_reps):
        # First 'remainder' groups get one extra sample
        size = base_size + (1 if i < remainder else 0)
        groups.append(shuffled[start:start+size].tolist())
        start += size

    # Create consolidated columns
    consolidated_cols = {}
    metadata_rows = []

    for i, group_samples in enumerate(groups, 1):
        pseudo_name = f"t0_pseudorep_{i}"
        consolidated_cols[pseudo_name] = wide_matrix[group_samples].sum(axis=1)

        metadata_rows.append({
            "sample_name": pseudo_name,
            "pseudo_rep_id": i,
            "generation": 0,
            "n_samples_consolidated": len(group_samples),
            "source_samples": ",".join(group_samples)
        })

        logger.info(f"  Pseudo-replicate {i}: {len(group_samples)} samples -> {pseudo_name}")

    # Remove original t0 samples, add pseudo-replicates
    non_t0_cols = [c for c in wide_matrix.columns if c not in t0_samples]

    result = pd.concat([
        wide_matrix[non_t0_cols],
        pd.DataFrame(consolidated_cols, index=wide_matrix.index)
    ], axis=1)

    return result, pd.DataFrame(metadata_rows)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create wide-format count matrix for DESeq2 with T0 consolidation"
    )
    parser.add_argument(
        "--bc1-counts", type=Path, required=True,
        help="Input BC1 counts long-format TSV"
    )
    parser.add_argument(
        "--tier-mapping", type=Path, required=True,
        help="BC1 tier mapping TSV"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition/generation info"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name to process"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output wide-format count matrix TSV"
    )
    parser.add_argument(
        "--metadata-output", type=Path,
        help="Output path for consolidation metadata (optional)"
    )
    parser.add_argument(
        "--t0-consolidation-mode", type=str,
        choices=["screen_date", "pseudo_replicates", "none"],
        default="screen_date",
        help="T0 consolidation strategy (default: screen_date)"
    )
    parser.add_argument(
        "--t0-pseudo-rep-count", type=int,
        help="Number of t0 pseudo-replicates to create (for pseudo_replicates mode)"
    )
    parser.add_argument(
        "--screen-dates", type=str, nargs="+",
        help="Screen dates for this library (optional filter)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for pseudo-replicate assignment (default: 42)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    logger.info(f"Step 05: Creating count matrix for library {args.library}")
    logger.info(f"  BC1 counts: {args.bc1_counts}")
    logger.info(f"  Tier mapping: {args.tier_mapping}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  T0 consolidation mode: {args.t0_consolidation_mode}")

    # Validate args
    if args.t0_consolidation_mode == "pseudo_replicates" and not args.t0_pseudo_rep_count:
        logger.error("--t0-pseudo-rep-count is required for pseudo_replicates mode")
        return 1

    # Load data
    logger.info("Loading BC1 counts...")
    bc1_counts = pd.read_csv(args.bc1_counts, sep="\t")
    logger.info(f"  Loaded {len(bc1_counts):,} rows")

    logger.info("Loading tier mapping...")
    tier_mapping = pd.read_csv(args.tier_mapping, sep="\t")
    logger.info(f"  Loaded {len(tier_mapping):,} BC1 tier assignments")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} sample entries")

    # Filter to relevant abundance tiers (high and medium)
    relevant_tiers = ["abundance_high", "abundance_medium"]
    tier_bc1s = tier_mapping[tier_mapping["tier"].isin(relevant_tiers)]["bc1"].unique()
    logger.info(f"  BC1s in abundance_high/abundance_medium: {len(tier_bc1s):,}")

    # Filter BC1 counts to relevant BC1s
    bc1_counts = bc1_counts[bc1_counts["bc1"].isin(tier_bc1s)]
    logger.info(f"  Filtered counts: {len(bc1_counts):,} rows")

    # Filter to library-specific samples if screen_dates provided
    if args.screen_dates:
        date_pattern = "|".join(args.screen_dates)
        if "sample_name" in bc1_counts.columns:
            bc1_counts = bc1_counts[bc1_counts["sample_name"].str.contains(date_pattern, na=False)]
        logger.info(f"  After screen date filter: {len(bc1_counts):,} rows")

    # Create wide-format matrix
    logger.info("Creating wide-format count matrix...")
    if "sample_name" in bc1_counts.columns and "counts" in bc1_counts.columns:
        wide_matrix = bc1_counts.pivot_table(
            index="bc1",
            columns="sample_name",
            values="counts",
            aggfunc="sum",
            fill_value=0
        )
    else:
        logger.error("Expected columns 'sample_name' and 'counts' in bc1_counts")
        return 1

    logger.info(f"  Matrix shape: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")

    # Apply T0 consolidation
    metadata_df = pd.DataFrame()

    if args.t0_consolidation_mode == "screen_date":
        logger.info("Consolidating t0 samples by screen_date...")
        wide_matrix, metadata_df = consolidate_t0_by_screen_date(
            wide_matrix, sample_key, args.library
        )
        logger.info(f"  Consolidated matrix: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")

    elif args.t0_consolidation_mode == "pseudo_replicates":
        logger.info(f"Creating {args.t0_pseudo_rep_count} t0 pseudo-replicates...")
        wide_matrix, metadata_df = consolidate_t0_as_pseudo_replicates(
            wide_matrix, sample_key, args.library,
            n_pseudo_reps=args.t0_pseudo_rep_count,
            seed=args.seed
        )
        logger.info(f"  Consolidated matrix: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")

    else:
        logger.info("No t0 consolidation applied")

    # Add tier information to output
    tier_info = tier_mapping[tier_mapping["bc1"].isin(wide_matrix.index)][["bc1", "tier"]].drop_duplicates()

    # Reset index for output
    wide_matrix = wide_matrix.reset_index()

    # Merge tier info
    wide_matrix = wide_matrix.merge(tier_info, on="bc1", how="left")

    # Reorder columns: bc1, tier, then samples
    cols = ["bc1", "tier"] + [c for c in wide_matrix.columns if c not in ["bc1", "tier"]]
    wide_matrix = wide_matrix[cols]

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    logger.info(f"Saving count matrix to {args.output}")
    wide_matrix.to_csv(args.output, sep="\t", index=False)

    # Save metadata if requested
    if args.metadata_output and len(metadata_df) > 0:
        logger.info(f"Saving consolidation metadata to {args.metadata_output}")
        args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
        metadata_df.to_csv(args.metadata_output, sep="\t", index=False)

    # Log summary by abundance tier
    for tier in ["abundance_high", "abundance_medium"]:
        n_tier = (wide_matrix["tier"] == tier).sum()
        logger.info(f"  {tier}: {n_tier:,} BC1s")

    logger.info("Step 05: Count matrix creation complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
