#!/usr/bin/env python3
"""
Step 06: Combine Technical Replicates

Combines counts from technical replicates (same biological sample amplified with
different primer pairs) before DESeq2 analysis.

Technical replicates are identified by the `sample_ID_with_inner_primers_combined`
column in the sample key - samples with the same ID differ only by primer pair
(e.g., KR1967-KR1966 vs KR1966-KR1967 are the same primers in reverse orientation).

This script:
1. Identifies technical replicate groups from sample key
2. Calculates R² correlation between replicates for QC
3. Flags outliers with low R² (potential amplification issues)
4. Sums counts across technical replicates
5. Outputs cleaned count matrix with one column per biological sample

Usage:
    python 06_combine_technical_replicates.py \\
        --count-matrix count_matrix.tsv \\
        --sample-key sample_key.tsv \\
        --output-matrix count_matrix_combined.tsv \\
        --qc-output-dir qc/technical_replicates/

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine technical replicates (fwd/rev primer pairs)"
    )
    parser.add_argument(
        "--count-matrix", type=Path, required=True,
        help="Input wide-format count matrix TSV"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with sample_ID_with_inner_primers_combined column"
    )
    parser.add_argument(
        "--output-matrix", type=Path, required=True,
        help="Output combined count matrix TSV"
    )
    parser.add_argument(
        "--qc-output-dir", type=Path, default=None,
        help="Directory for QC outputs (correlations, outliers, plots)"
    )
    parser.add_argument(
        "--r2-threshold", type=float, default=0.8,
        help="R² threshold below which to flag as outlier (default: 0.8)"
    )
    parser.add_argument(
        "--min-counts", type=int, default=10,
        help="Minimum counts in both replicates for correlation (default: 10)"
    )
    parser.add_argument(
        "--replicate-id-column", type=str,
        default="sample_ID_with_inner_primers_combined",
        help="Column name in sample key that groups technical replicates"
    )
    return parser.parse_args()


def calculate_log_correlation(counts1: pd.Series, counts2: pd.Series,
                               min_counts: int = 10) -> dict:
    """
    Calculate Pearson correlation on log10-transformed counts.

    Args:
        counts1: First count series (indexed by bc1)
        counts2: Second count series (indexed by bc1)
        min_counts: Minimum counts in both samples for inclusion

    Returns:
        Dictionary with r, r2, n_bc1s, and other stats
    """
    # Align by index (bc1)
    common_idx = counts1.index.intersection(counts2.index)
    c1 = counts1.loc[common_idx]
    c2 = counts2.loc[common_idx]

    # Filter to bc1s with sufficient counts in both
    mask = (c1 >= min_counts) & (c2 >= min_counts)
    c1_filt = c1[mask]
    c2_filt = c2[mask]

    if len(c1_filt) < 10:
        return {
            'r': np.nan,
            'r2': np.nan,
            'n_bc1s': len(c1_filt),
            'n_bc1s_total': len(common_idx),
            'pearson_p': np.nan,
            'status': 'insufficient_data'
        }

    # Log transform (add 1 for zeros)
    log_c1 = np.log10(c1_filt + 1)
    log_c2 = np.log10(c2_filt + 1)

    # Calculate correlation
    r, p = stats.pearsonr(log_c1, log_c2)

    return {
        'r': r,
        'r2': r ** 2,
        'n_bc1s': len(c1_filt),
        'n_bc1s_total': len(common_idx),
        'pearson_p': p,
        'status': 'ok'
    }


def identify_replicate_groups(sample_key: pd.DataFrame,
                               replicate_id_column: str) -> dict:
    """
    Group samples by their biological replicate ID.

    Args:
        sample_key: Sample key DataFrame
        replicate_id_column: Column name that groups technical replicates

    Returns:
        Dictionary: {replicate_id: [sample_name1, sample_name2, ...]}
    """
    if replicate_id_column not in sample_key.columns:
        logger.warning(f"Column '{replicate_id_column}' not in sample key. "
                      f"Available columns: {list(sample_key.columns)}")
        return {}

    groups = sample_key.groupby(replicate_id_column)['sample_name'].apply(list).to_dict()

    # Only keep groups with >1 sample (actual replicates)
    replicate_groups = {k: v for k, v in groups.items() if len(v) > 1}
    singleton_groups = {k: v for k, v in groups.items() if len(v) == 1}

    logger.info(f"  Found {len(replicate_groups)} replicate groups "
               f"({sum(len(v) for v in replicate_groups.values())} samples)")
    logger.info(f"  Found {len(singleton_groups)} singleton samples (no replicates)")

    return replicate_groups, singleton_groups


def combine_technical_replicates(count_matrix: pd.DataFrame,
                                  sample_key: pd.DataFrame,
                                  replicate_id_column: str,
                                  r2_threshold: float = 0.8,
                                  min_counts: int = 10) -> tuple:
    """
    Combine counts from technical replicates.

    Args:
        count_matrix: BC1 x samples count matrix
        sample_key: Sample key with replicate grouping column
        replicate_id_column: Column name for replicate grouping
        r2_threshold: R² threshold for outlier flagging
        min_counts: Minimum counts for correlation calculation

    Returns:
        Tuple of (combined_matrix, correlation_results, outliers)
    """
    replicate_groups, singleton_groups = identify_replicate_groups(
        sample_key, replicate_id_column
    )

    # Get sample columns from count matrix (exclude bc1 and metadata columns)
    meta_cols = ['bc1', 'tier', 'oligo_name']
    sample_cols = [c for c in count_matrix.columns if c not in meta_cols]

    # Build sample name to column mapping
    sample_to_col = {c: c for c in sample_cols}

    # Start with metadata columns
    combined = count_matrix[['bc1']].copy()
    if 'tier' in count_matrix.columns:
        combined['tier'] = count_matrix['tier']

    correlation_results = []
    outliers = []
    processed_samples = set()

    # Process replicate groups
    for group_id, samples in replicate_groups.items():
        # Get samples that exist in the count matrix
        available_samples = [s for s in samples if s in sample_to_col]

        if len(available_samples) == 0:
            logger.warning(f"No samples found for group {group_id}")
            continue
        elif len(available_samples) == 1:
            # Only one replicate available, use as-is
            combined[group_id] = count_matrix[available_samples[0]]
            processed_samples.add(available_samples[0])
            continue

        # Calculate pairwise correlations
        for i in range(len(available_samples)):
            for j in range(i + 1, len(available_samples)):
                s1, s2 = available_samples[i], available_samples[j]
                corr = calculate_log_correlation(
                    count_matrix[s1],
                    count_matrix[s2],
                    min_counts=min_counts
                )

                corr['group_id'] = group_id
                corr['sample_1'] = s1
                corr['sample_2'] = s2
                corr['is_outlier'] = corr['r2'] < r2_threshold if not np.isnan(corr['r2']) else True

                correlation_results.append(corr)

                if corr['is_outlier']:
                    outliers.append({
                        'group_id': group_id,
                        'sample_1': s1,
                        'sample_2': s2,
                        'r2': corr['r2'],
                        'n_bc1s': corr['n_bc1s'],
                        'reason': 'low_correlation' if not np.isnan(corr['r2']) else 'insufficient_data'
                    })

        # Sum counts across all replicates
        combined[group_id] = count_matrix[available_samples].sum(axis=1)

        for s in available_samples:
            processed_samples.add(s)

    # Add singleton samples (no replicates)
    for group_id, samples in singleton_groups.items():
        sample = samples[0]
        if sample in sample_to_col:
            combined[group_id] = count_matrix[sample]
            processed_samples.add(sample)

    # Check for any unprocessed samples
    unprocessed = set(sample_cols) - processed_samples
    if unprocessed:
        logger.warning(f"  {len(unprocessed)} samples not in sample key, kept as-is")
        for sample in unprocessed:
            combined[sample] = count_matrix[sample]

    return combined, pd.DataFrame(correlation_results), pd.DataFrame(outliers)


def save_qc_outputs(correlation_df: pd.DataFrame,
                    outliers_df: pd.DataFrame,
                    output_dir: Path):
    """Save QC outputs to files."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save correlations
    corr_file = output_dir / "technical_replicate_correlations.tsv"
    correlation_df.to_csv(corr_file, sep='\t', index=False)
    logger.info(f"  Saved correlations to {corr_file}")

    # Save outliers
    outlier_file = output_dir / "technical_replicate_outliers.tsv"
    outliers_df.to_csv(outlier_file, sep='\t', index=False)
    logger.info(f"  Saved outliers to {outlier_file}")

    # Print summary statistics
    if len(correlation_df) > 0:
        valid_r2 = correlation_df['r2'].dropna()
        if len(valid_r2) > 0:
            logger.info(f"  R² statistics: mean={valid_r2.mean():.3f}, "
                       f"median={valid_r2.median():.3f}, "
                       f"min={valid_r2.min():.3f}, max={valid_r2.max():.3f}")
            logger.info(f"  Outliers flagged: {len(outliers_df)}")

    # Create simple histogram plot if matplotlib available
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        if len(correlation_df) > 0:
            valid_r2 = correlation_df['r2'].dropna()
            if len(valid_r2) > 5:
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.hist(valid_r2, bins=30, edgecolor='black', alpha=0.7)
                ax.axvline(x=0.8, color='r', linestyle='--', label='Outlier threshold (0.8)')
                ax.set_xlabel('R² (log10 counts)')
                ax.set_ylabel('Number of replicate pairs')
                ax.set_title(f'Technical Replicate Correlation Distribution\n'
                            f'(n={len(valid_r2)}, median={valid_r2.median():.3f})')
                ax.legend()

                plot_file = output_dir / "correlation_distribution.png"
                plt.savefig(plot_file, dpi=150, bbox_inches='tight')
                plt.close()
                logger.info(f"  Saved plot to {plot_file}")
    except ImportError:
        logger.debug("matplotlib not available, skipping plot")


def main():
    args = parse_args()

    logger.info("Step 06: Combining Technical Replicates")
    logger.info(f"  Count matrix: {args.count_matrix}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  Output: {args.output_matrix}")
    logger.info(f"  R² threshold: {args.r2_threshold}")

    # Load data
    logger.info("Loading count matrix...")
    count_matrix = pd.read_csv(args.count_matrix, sep='\t')
    logger.info(f"  Loaded {len(count_matrix):,} BC1s x {len(count_matrix.columns)} columns")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep='\t')
    logger.info(f"  Loaded {len(sample_key):,} sample entries")

    # Check for required column
    if args.replicate_id_column not in sample_key.columns:
        logger.error(f"Column '{args.replicate_id_column}' not found in sample key")
        logger.error(f"Available columns: {list(sample_key.columns)}")
        return 1

    # Combine technical replicates
    logger.info("Combining technical replicates...")
    combined, correlations, outliers = combine_technical_replicates(
        count_matrix,
        sample_key,
        args.replicate_id_column,
        r2_threshold=args.r2_threshold,
        min_counts=args.min_counts
    )

    logger.info(f"  Output matrix: {len(combined):,} BC1s x {len(combined.columns)} columns")

    # Save combined matrix
    args.output_matrix.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output_matrix, sep='\t', index=False)
    logger.info(f"  Saved combined matrix to {args.output_matrix}")

    # Save QC outputs if directory specified
    if args.qc_output_dir:
        logger.info("Saving QC outputs...")
        save_qc_outputs(correlations, outliers, args.qc_output_dir)

    logger.info("Step 06: Complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
