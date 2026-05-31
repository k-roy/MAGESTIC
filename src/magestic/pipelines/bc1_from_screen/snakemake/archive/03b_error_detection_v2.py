#!/usr/bin/env python3
"""
Step 03b: Error BC1 Detection (v2 - Optimized)

Identifies BC1s that are likely PCR or sequencing errors based on:
1. Low count + sparse sample presence
2. Hamming distance 1-2 from a high-abundance BC1

OPTIMIZATION: Instead of generating variants for each of 2M potential errors,
we precompute all HD1/HD2 variants from the ~200k abundant BC1s once.
This reduces complexity from O(errors × variants) to O(abundant × variants).

For HD=2 with 20bp BC1s:
- Old approach: 2M errors × 1800 variants = 3.6 billion lookups
- New approach: 200k abundant × 1800 variants = 360M to build dict + 2M lookups

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor
import sys


def detect_count_column(df: pd.DataFrame) -> str:
    """Auto-detect the count column name."""
    if 'counts' in df.columns:
        return 'counts'
    elif 'sample_counts' in df.columns:
        return 'sample_counts'
    else:
        count_cols = [c for c in df.columns if 'count' in c.lower()]
        if count_cols:
            return count_cols[0]
        raise ValueError(f"Could not find count column. Available columns: {df.columns.tolist()}")


def calculate_bc1_statistics(counts_df: pd.DataFrame, count_col: str = 'counts') -> pd.DataFrame:
    """Calculate per-BC1 statistics from long-format counts."""
    stats = counts_df.groupby('bc1').agg(
        total_counts=(count_col, 'sum'),
        n_samples=(count_col, 'count'),
        max_count=(count_col, 'max'),
        mean_count=(count_col, 'mean')
    ).reset_index()
    return stats


def generate_variants_for_parent(parent_bc1: str, parent_counts: int, max_hamming: int = 2) -> list:
    """
    Generate all HD1 (and optionally HD2) variants for a parent BC1.

    Returns list of tuples: (variant, parent_bc1, hamming_distance, parent_counts)
    """
    bases = 'ACGT'
    results = []
    seq_len = len(parent_bc1)

    # HD=1 variants
    for i in range(seq_len):
        for base in bases:
            if base != parent_bc1[i]:
                variant = parent_bc1[:i] + base + parent_bc1[i+1:]
                results.append((variant, parent_bc1, 1, parent_counts))

    # HD=2 variants (if requested)
    if max_hamming >= 2:
        for i in range(seq_len):
            for j in range(i+1, seq_len):
                for base_i in bases:
                    for base_j in bases:
                        if base_i != parent_bc1[i] or base_j != parent_bc1[j]:
                            variant = parent_bc1[:i] + base_i + parent_bc1[i+1:j] + base_j + parent_bc1[j+1:]
                            if variant != parent_bc1:
                                results.append((variant, parent_bc1, 2, parent_counts))

    return results


def build_variant_lookup(high_abundance_bc1s: list, high_abundance_counts: dict,
                         max_hamming: int = 2, n_workers: int = 8) -> Dict[str, Tuple[str, int, int]]:
    """
    Build a dictionary mapping variant sequences to their parent BC1.

    For collisions (multiple parents could produce the same variant),
    we keep the parent with the highest count.

    Returns:
        dict: variant_seq -> (parent_bc1, hamming_distance, parent_counts)
    """
    print(f"Building variant lookup table from {len(high_abundance_bc1s):,} abundant BC1s...")
    print(f"  Max Hamming distance: {max_hamming}")

    # Estimate size
    bc1_len = len(high_abundance_bc1s[0]) if high_abundance_bc1s else 20
    hd1_per_bc1 = bc1_len * 3  # 3 alternative bases per position
    hd2_per_bc1 = bc1_len * (bc1_len - 1) // 2 * 9  # C(n,2) * 9 combinations

    if max_hamming == 1:
        estimated_variants = len(high_abundance_bc1s) * hd1_per_bc1
    else:
        estimated_variants = len(high_abundance_bc1s) * (hd1_per_bc1 + hd2_per_bc1)

    print(f"  Estimated variants to generate: {estimated_variants:,}")

    variant_lookup = {}

    # Process in chunks with progress reporting
    chunk_size = 10000
    total = len(high_abundance_bc1s)

    for chunk_start in range(0, total, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total)
        chunk = high_abundance_bc1s[chunk_start:chunk_end]

        for parent_bc1 in chunk:
            parent_counts = high_abundance_counts[parent_bc1]
            variants = generate_variants_for_parent(parent_bc1, parent_counts, max_hamming)

            for variant, parent, hd, counts in variants:
                # If variant already exists, keep the higher-count parent
                # Also prefer HD=1 over HD=2 if counts are equal
                if variant in variant_lookup:
                    existing_parent, existing_hd, existing_counts = variant_lookup[variant]
                    if counts > existing_counts or (counts == existing_counts and hd < existing_hd):
                        variant_lookup[variant] = (parent, hd, counts)
                else:
                    variant_lookup[variant] = (parent, hd, counts)

        if (chunk_end % 50000 == 0) or chunk_end == total:
            print(f"  Processed {chunk_end:,}/{total:,} abundant BC1s, lookup table size: {len(variant_lookup):,}")

    print(f"  Final lookup table size: {len(variant_lookup):,} variants")

    # Memory estimate
    mem_bytes = sys.getsizeof(variant_lookup)
    for k, v in list(variant_lookup.items())[:1000]:
        mem_bytes += sys.getsizeof(k) + sys.getsizeof(v)
    mem_estimate_gb = (mem_bytes / 1000) * len(variant_lookup) / 1e9
    print(f"  Estimated memory usage: ~{mem_estimate_gb:.1f} GB")

    return variant_lookup


def identify_error_bc1s_optimized(
    bc1_stats: pd.DataFrame,
    min_count_ratio: float = 100.0,
    max_error_count: int = 50,
    max_error_samples: int = 3,
    max_hamming: int = 2,
    n_workers: int = 8
) -> pd.DataFrame:
    """
    Identify BC1s that are likely errors using precomputed variant lookup.
    """
    print(f"Total BC1s to analyze: {len(bc1_stats):,}")

    # Split into potential errors and high-abundance BC1s
    potential_errors = bc1_stats[
        (bc1_stats['total_counts'] <= max_error_count) &
        (bc1_stats['n_samples'] <= max_error_samples)
    ].copy()

    high_abundance = bc1_stats[bc1_stats['total_counts'] > max_error_count].copy()

    print(f"Potential error BC1s (count<={max_error_count}, samples<={max_error_samples}): {len(potential_errors):,}")
    print(f"High-abundance BC1s: {len(high_abundance):,}")

    if len(potential_errors) == 0 or len(high_abundance) == 0:
        print("No potential errors or no high-abundance BC1s found")
        return pd.DataFrame(columns=['error_bc1', 'parent_bc1', 'hamming_distance',
                                     'error_counts', 'parent_counts', 'count_ratio'])

    # Build precomputed lookup table (the key optimization)
    high_abundance_bc1s = high_abundance['bc1'].tolist()
    high_abundance_counts = dict(zip(high_abundance['bc1'], high_abundance['total_counts']))

    variant_lookup = build_variant_lookup(
        high_abundance_bc1s, high_abundance_counts, max_hamming, n_workers
    )

    # Now just do O(1) lookups for each potential error
    print(f"\nSearching for parent BC1s using precomputed lookup...")

    error_mappings = []
    potential_error_list = potential_errors['bc1'].tolist()
    error_counts_dict = dict(zip(potential_errors['bc1'], potential_errors['total_counts']))

    chunk_size = 100000
    total = len(potential_error_list)

    for chunk_start in range(0, total, chunk_size):
        chunk_end = min(chunk_start + chunk_size, total)
        chunk = potential_error_list[chunk_start:chunk_end]

        for error_bc1 in chunk:
            if error_bc1 in variant_lookup:
                parent_bc1, hamming_dist, parent_counts = variant_lookup[error_bc1]
                error_count = error_counts_dict[error_bc1]
                count_ratio = parent_counts / error_count if error_count > 0 else float('inf')

                if count_ratio >= min_count_ratio:
                    error_mappings.append({
                        'error_bc1': error_bc1,
                        'parent_bc1': parent_bc1,
                        'hamming_distance': hamming_dist,
                        'error_counts': error_count,
                        'parent_counts': parent_counts,
                        'count_ratio': count_ratio
                    })

        if (chunk_end % 500000 == 0) or chunk_end == total:
            print(f"  Processed {chunk_end:,}/{total:,} potential errors, found {len(error_mappings):,} mappings")

    result_df = pd.DataFrame(error_mappings)

    print(f"\nIdentified {len(result_df):,} error BC1s")
    if len(result_df) > 0:
        print(f"  Hamming distance 1: {(result_df['hamming_distance'] == 1).sum():,}")
        print(f"  Hamming distance 2: {(result_df['hamming_distance'] == 2).sum():,}")

    return result_df


def main():
    parser = argparse.ArgumentParser(
        description="Identify error BC1s by Hamming distance analysis (optimized v2)"
    )
    parser.add_argument("--counts", type=Path, required=True,
                        help="Combined long-format counts file")
    parser.add_argument("--output", type=Path, required=True,
                        help="Output path for error BC1 mapping file")
    parser.add_argument("--stats-output", type=Path,
                        help="Output path for BC1 statistics (optional)")
    parser.add_argument("--max-error-count", type=int, default=50,
                        help="Maximum total counts for error BC1 (default: 50)")
    parser.add_argument("--max-error-samples", type=int, default=3,
                        help="Maximum samples for error BC1 (default: 3)")
    parser.add_argument("--min-count-ratio", type=float, default=100.0,
                        help="Minimum ratio of parent to error counts (default: 100)")
    parser.add_argument("--max-hamming", type=int, default=2, choices=[1, 2],
                        help="Maximum Hamming distance (default: 2)")
    parser.add_argument("--workers", type=int, default=8,
                        help="Number of parallel workers (default: 8)")

    args = parser.parse_args()

    # Load counts
    print(f"Loading counts from {args.counts}...")
    counts_df = pd.read_csv(args.counts, sep='\t')
    print(f"Loaded {len(counts_df):,} rows")

    # Auto-detect count column
    count_col = detect_count_column(counts_df)
    print(f"Using count column: '{count_col}'")

    # Calculate BC1 statistics
    print("\nCalculating BC1 statistics...")
    bc1_stats = calculate_bc1_statistics(counts_df, count_col=count_col)
    print(f"Unique BC1s: {len(bc1_stats):,}")
    print(f"Total counts: {bc1_stats['total_counts'].sum():,}")

    # Save stats if requested
    if args.stats_output:
        print(f"\nSaving BC1 statistics to {args.stats_output}...")
        args.stats_output.parent.mkdir(parents=True, exist_ok=True)
        bc1_stats.to_csv(args.stats_output, sep='\t', index=False)

    # Identify error BC1s using optimized algorithm
    print("\nIdentifying error BC1s (optimized algorithm)...")
    error_mappings = identify_error_bc1s_optimized(
        bc1_stats,
        min_count_ratio=args.min_count_ratio,
        max_error_count=args.max_error_count,
        max_error_samples=args.max_error_samples,
        max_hamming=args.max_hamming,
        n_workers=args.workers
    )

    # Save error mappings
    print(f"\nSaving error mappings to {args.output}...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    error_mappings.to_csv(args.output, sep='\t', index=False)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total BC1s analyzed: {len(bc1_stats):,}")
    print(f"Error BC1s identified: {len(error_mappings):,}")
    if len(error_mappings) > 0:
        total_error_counts = error_mappings['error_counts'].sum()
        print(f"Total counts in error BC1s: {total_error_counts:,}")
        print(f"These counts will be absorbed by parent BC1s in the next step")

    print("\nDone!")


if __name__ == "__main__":
    main()
