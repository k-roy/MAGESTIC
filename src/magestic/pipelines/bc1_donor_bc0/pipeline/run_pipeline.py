#!/usr/bin/env python3
"""
Main pipeline runner for BC1-Donor-BC0 processing.

This script orchestrates the full pipeline:
1. Load data from available sources (2x150, 2x300, gdb)
2. Merge data sources at the earliest stage
3. Apply unified sliding window matching for oligo assignment
4. Output reference table with bc1 -> oligo mappings

Usage:
    python run_pipeline.py --project-dir /path/to/project [options]

The pipeline automatically detects available data sources and adjusts
processing accordingly. It supports three scenarios:
- Only 2x150 data
- Only 2x300 data
- Both 2x150 and 2x300 data (merged for maximum bc1 recovery)
"""

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np
from tqdm import tqdm

from ..config import PipelineConfig
from ..readers.data_loader import load_all_data_sources
from ..core.merge import merge_data_sources
from ..core.matching import (
    load_and_build_lookup,
    match_dataframe,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='BC1-Donor-BC0 Processing Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Run with auto-detection of data sources
    python run_pipeline.py --project-dir /path/to/project

    # Run in debug mode with limited data
    python run_pipeline.py --project-dir /path/to/project --debug

    # Specify output directory
    python run_pipeline.py --project-dir /path/to/project --output-dir /path/to/output
        """
    )

    parser.add_argument(
        '--project-dir', '-p',
        type=Path,
        required=True,
        help='Root directory of the project'
    )

    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=None,
        help='Output directory (default: project_dir/processed_data/bc1_donor_bc0)'
    )

    parser.add_argument(
        '--common-dir', '-c',
        type=Path,
        default=Path('/path/to/common'),
        help='Directory containing common resources (oligo designs, etc.)'
    )

    parser.add_argument(
        '--min-reads',
        type=int,
        default=2,
        help='Minimum reads for a bc1-donor-bc0 to be included (default: 2)'
    )

    parser.add_argument(
        '--min-pcr-replicates',
        type=int,
        default=2,
        help='Minimum PCR replicates for confidence (default: 2)'
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help='Run in debug mode with limited data'
    )

    parser.add_argument(
        '--debug-n-samples',
        type=int,
        default=4,
        help='Number of samples to load in debug mode (default: 4)'
    )

    parser.add_argument(
        '--n-jobs',
        type=int,
        default=16,
        help='Number of parallel jobs (default: 16)'
    )

    return parser.parse_args()


def collapse_sequences_by_edit_distance(
    df: pd.DataFrame,
    group_cols: list,
    seq_col: str,
    count_cols: list,
    edit_distance_threshold: int = 3,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Collapse similar sequences within groups by edit distance.

    Args:
        df: Input DataFrame
        group_cols: Columns to group by
        seq_col: Column containing sequences to collapse
        count_cols: Columns to sum when collapsing
        edit_distance_threshold: Maximum edit distance for collapsing
        show_progress: Show progress bar

    Returns:
        DataFrame with collapsed sequences
    """
    from Levenshtein import distance

    def process_group(group):
        """Collapse sequences within a single group."""
        group = group.copy()
        result_rows = []

        while not group.empty:
            # Find row with highest count
            top_idx = group[count_cols].sum(axis=1).idxmax()
            top_seq = group.loc[top_idx, seq_col]

            # Find all sequences within edit distance
            group['_edit_dist'] = group[seq_col].apply(
                lambda x: distance(x, top_seq) if pd.notna(x) else 999
            )
            collapse_mask = group['_edit_dist'] <= edit_distance_threshold

            # Aggregate counts
            collapsed_counts = {
                col: group.loc[collapse_mask, col].sum()
                for col in count_cols
            }

            # Create result row
            result_rows.append({
                **{col: group.iloc[0][col] for col in group_cols},
                seq_col: top_seq,
                **collapsed_counts,
            })

            # Remove collapsed rows
            group = group.loc[~collapse_mask].drop(columns=['_edit_dist'])

        return result_rows

    # Process each group
    grouped = df.groupby(group_cols)
    all_results = []

    groups = list(grouped)
    iterator = tqdm(groups, desc=f"Collapsing {seq_col}") if show_progress else groups

    for _, group in iterator:
        all_results.extend(process_group(group))

    return pd.DataFrame(all_results)


def aggregate_and_filter(
    df: pd.DataFrame,
    min_reads: int = 2,
    min_pcr_replicates: int = 2,
    collapse_sequences: bool = True,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    Aggregate counts, collapse sequences, and apply quality filters.

    Args:
        df: Merged DataFrame with bc1, donor, bc0, counts columns
        min_reads: Minimum total reads
        min_pcr_replicates: Minimum PCR replicates
        collapse_sequences: Whether to collapse similar sequences
        show_progress: Show progress bars

    Returns:
        Filtered and aggregated DataFrame
    """
    logger.info("Aggregating and filtering data...")

    # Filter for entries with full donor
    if 'donor' in df.columns:
        df = df[df['donor'].notna()].copy()
        logger.info(f"  After filtering for full donor: {len(df)} rows")

    # HANDOFF_20260528 fix (Agent B, yL437/yL442): coalesce per-data-source counts
    # columns (counts_2x300 / counts_2x150) into a unified 'counts' column when the
    # merge produced data-source-suffixed names. Without this, 2x300-only screens
    # (yL406, yL437, yL442) crash with `KeyError: 'counts'` at the groupby below.
    if 'counts' not in df.columns:
        if 'counts_2x300' in df.columns and 'counts_2x150' in df.columns:
            df = df.copy()
            df['counts'] = df['counts_2x300'].fillna(0) + df['counts_2x150'].fillna(0)
        elif 'counts_2x300' in df.columns:
            df = df.copy()
            df['counts'] = df['counts_2x300']
        elif 'counts_2x150' in df.columns:
            df = df.copy()
            df['counts'] = df['counts_2x150']
        else:
            raise KeyError(
                "aggregate_and_filter: no 'counts', 'counts_2x300', or 'counts_2x150' "
                "column found in DataFrame. Upstream merge step is broken."
            )

    # Aggregate by bc1, donor, bc0
    agg_cols = {
        'counts': 'sum',
    }

    # Track PCR replicates and data sources
    if 'sample_number' in df.columns:
        agg_cols['sample_number'] = lambda x: len(set(x))

    if 'data_source' in df.columns:
        agg_cols['data_source'] = lambda x: ','.join(sorted(set(x)))

    if 'donor_bc0_fragment' in df.columns:
        agg_cols['donor_bc0_fragment'] = 'first'

    group_cols = ['bc1', 'donor', 'bc0']
    group_cols = [c for c in group_cols if c in df.columns]

    agg_df = df.groupby(group_cols).agg(agg_cols).reset_index()

    if 'sample_number' in agg_df.columns:
        agg_df.rename(columns={'sample_number': 'num_PCR_replicates'}, inplace=True)

    logger.info(f"  After aggregation: {len(agg_df)} rows, {agg_df['bc1'].nunique()} unique bc1s")

    # Optional sequence collapsing
    if collapse_sequences and len(agg_df) > 0:
        count_cols = ['counts']
        if 'num_PCR_replicates' in agg_df.columns:
            count_cols.append('num_PCR_replicates')

        # Collapse similar donors within bc1+bc0
        logger.info("  Collapsing similar donors...")
        agg_df = collapse_sequences_by_edit_distance(
            agg_df,
            group_cols=['bc1', 'bc0'],
            seq_col='donor',
            count_cols=count_cols,
            edit_distance_threshold=5,
            show_progress=show_progress
        )

        # Collapse similar bc1s within bc0+donor
        logger.info("  Collapsing similar bc1s...")
        agg_df = collapse_sequences_by_edit_distance(
            agg_df,
            group_cols=['bc0', 'donor'],
            seq_col='bc1',
            count_cols=count_cols,
            edit_distance_threshold=1,
            show_progress=show_progress
        )

    # Apply filters
    if 'counts' in agg_df.columns and min_reads > 0:
        agg_df = agg_df[agg_df['counts'] >= min_reads]
        logger.info(f"  After min reads filter: {len(agg_df)} rows")

    if 'num_PCR_replicates' in agg_df.columns and min_pcr_replicates > 0:
        agg_df = agg_df[agg_df['num_PCR_replicates'] >= min_pcr_replicates]
        logger.info(f"  After min PCR replicates filter: {len(agg_df)} rows")

    # Calculate bc1 purity metrics
    if 'bc1' in agg_df.columns and 'counts' in agg_df.columns:
        agg_df['total_bc1_counts'] = agg_df.groupby('bc1')['counts'].transform('sum')
        agg_df['bc1_purity'] = agg_df['counts'] / agg_df['total_bc1_counts']
        if 'donor' in agg_df.columns:
            agg_df['num_distinct_donors'] = agg_df.groupby('bc1')['donor'].transform('nunique')

    return agg_df


def run_pipeline(args):
    """Main pipeline execution."""
    logger.info("=" * 60)
    logger.info("BC1-DONOR-BC0 PROCESSING PIPELINE")
    logger.info("=" * 60)

    # Initialize configuration
    logger.info(f"\nProject directory: {args.project_dir}")
    config = PipelineConfig.from_project_dir(args.project_dir, args.common_dir)

    if args.output_dir:
        config.output_dir = args.output_dir

    config.min_reads = args.min_reads
    config.min_pcr_replicates = args.min_pcr_replicates

    # Validate configuration
    errors = config.validate()
    if errors:
        for err in errors:
            logger.error(err)
        return None

    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Log detected data sources
    logger.info("\nDetected data sources:")
    logger.info(f"  2x150 data: {config.has_2x150_data}")
    logger.info(f"  2x300 data: {config.has_2x300_data}")
    logger.info(f"  gdb data: {config.has_gdb_data}")

    if not any([config.has_2x150_data, config.has_2x300_data]):
        logger.error("No data sources found. Please check project directory structure.")
        return None

    # Step 1: Load all available data sources
    logger.info("\n" + "=" * 60)
    logger.info("STEP 1: Loading data sources")
    logger.info("=" * 60)

    data_sources = load_all_data_sources(
        config,
        debug_mode=args.debug,
        debug_n_samples=args.debug_n_samples
    )

    if not data_sources:
        logger.error("No data loaded from any source")
        return None

    for name, (df, info) in data_sources.items():
        logger.info(f"\n{name}: {info.num_rows} rows, {info.num_unique_bc1s} unique bc1s")

    # Step 2: Merge data sources
    logger.info("\n" + "=" * 60)
    logger.info("STEP 2: Merging data sources")
    logger.info("=" * 60)

    df_2x150 = data_sources.get('2x150', (None, None))[0]
    df_2x300 = data_sources.get('2x300', (None, None))[0]
    df_gdb = data_sources.get('gdb', (None, None))[0]

    merged_df = merge_data_sources(df_2x150, df_2x300, df_gdb)

    if merged_df.empty:
        logger.error("Merge resulted in empty DataFrame")
        return None

    logger.info(f"\nMerged dataset: {len(merged_df)} rows, {merged_df['bc1'].nunique()} unique bc1s")

    # Step 3: Aggregate and filter
    logger.info("\n" + "=" * 60)
    logger.info("STEP 3: Aggregating and filtering")
    logger.info("=" * 60)

    filtered_df = aggregate_and_filter(
        merged_df,
        min_reads=config.min_reads,
        min_pcr_replicates=config.min_pcr_replicates,
        collapse_sequences=not args.debug,  # Skip collapsing in debug mode for speed
        show_progress=True
    )

    logger.info(f"\nFiltered dataset: {len(filtered_df)} rows, {filtered_df['bc1'].nunique()} unique bc1s")

    # Step 4: Oligo matching
    logger.info("\n" + "=" * 60)
    logger.info("STEP 4: Matching to designed oligos")
    logger.info("=" * 60)

    # Build lookup tables from the harmonized oligo table (covers all four nuclease
    # classes -- SpG / SpCas9 / LbCas12a / impLbCas12a -- plus the cAMP/Ras/PKA
    # saturation oligos (V574-V597) added 2026-05-30). The previous two-file path
    # hard-coded `20210422_Twist_200mer_SpCas9_oligos_with_controls.tsv`, which has
    # never existed on /oak, so this step was a silent no-op (logged "Oligo design
    # files not found, skipping matching"). The harmonized table lives at the
    # canonical `common/annotation_files/` path resolved inside
    # `magestic.utils.oligo_annotations.HARMONIZED_OLIGO_TABLE`.
    logger.info("Building oligo lookup tables from harmonized table...")
    sublibrary_tables = load_and_build_lookup(use_harmonized_table=True)

    # Determine which columns to use for matching
    donor_col = 'donor' if 'donor' in filtered_df.columns else 'bdb_donor'
    guide_col = 'guide' if 'guide' in filtered_df.columns else 'gdb_guide'
    # Check for yeast_sublibrary column (or legacy "sublibrary" for backward compatibility)
    yeast_sublibrary_col = None
    if 'yeast_sublibrary' in filtered_df.columns:
        yeast_sublibrary_col = 'yeast_sublibrary'
    elif 'sublibrary' in filtered_df.columns:
        yeast_sublibrary_col = 'sublibrary'

    logger.info("Matching sequences to oligos...")
    matched_df = match_dataframe(
        filtered_df,
        sublibrary_tables,
        donor_col=donor_col,
        guide_col=guide_col if guide_col in filtered_df.columns else 'guide',
        prefix=''
    )

    # Step 5: Output results
    logger.info("\n" + "=" * 60)
    logger.info("STEP 5: Saving results")
    logger.info("=" * 60)

    timestamp = datetime.now().strftime("%Y%m%d")
    output_file = config.output_dir / f"{timestamp}_bc1_reference_table_min_{config.min_reads}_reads.tsv"

    matched_df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Saved reference table to: {output_file}")

    # Summary statistics
    logger.info("\n" + "=" * 60)
    logger.info("PIPELINE SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total unique bc1s: {matched_df['bc1'].nunique()}")
    logger.info(f"Total rows: {len(matched_df)}")

    if 'oligo_name' in matched_df.columns:
        matched = matched_df['oligo_name'].notna().sum()
        logger.info(f"Matched to oligos: {matched} ({100*matched/len(matched_df):.1f}%)")

    if 'data_source' in matched_df.columns:
        logger.info("\nData source breakdown:")
        for source in matched_df['data_source'].unique():
            count = len(matched_df[matched_df['data_source'].str.contains(source, na=False)])
            logger.info(f"  {source}: {count}")

    return matched_df


def main():
    """Entry point."""
    args = parse_args()
    result = run_pipeline(args)

    if result is not None:
        logger.info("\nPipeline completed successfully!")
        return 0
    else:
        logger.error("\nPipeline failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
