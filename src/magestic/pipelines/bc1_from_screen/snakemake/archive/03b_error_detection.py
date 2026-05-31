#!/usr/bin/env python3
"""
Step 03b: Error BC1 Detection

Identifies BC1s that are likely PCR or sequencing errors based on:
1. Low count + sparse sample presence
2. Hamming distance 1-2 from a high-abundance BC1

Error BC1s can have their counts absorbed by their parent BC1 in the next step,
which improves tier classification by consolidating fragmented counts.

Input:
    - Combined long-format counts file (from 03_combine_counts.py)
    - BC1 reference table (for additional context)

Output:
    - Error BC1 mapping file: error_bc1 -> parent_bc1 with metadata
    - Summary statistics file

Author: Kevin R. Roy
"""

import argparse
import os
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Set, Optional
from concurrent.futures import ProcessPoolExecutor
import itertools


def calculate_bc1_statistics(counts_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate per-BC1 statistics from long-format counts.

    Returns DataFrame with columns:
        bc1, total_counts, n_samples, max_count, mean_count
    """
    stats = counts_df.groupby('bc1').agg(
        total_counts=('counts', 'sum'),
        n_samples=('counts', 'count'),
        max_count=('counts', 'max'),
        mean_count=('counts', 'mean')
    ).reset_index()

    return stats


def generate_hamming_1_variants(seq: str) -> Set[str]:
    """Generate all sequences that are Hamming distance 1 from input."""
    variants = set()
    bases = 'ACGT'

    for i, char in enumerate(seq):
        for base in bases:
            if base != char:
                variant = seq[:i] + base + seq[i+1:]
                variants.add(variant)

    return variants


def generate_hamming_2_variants(seq: str) -> Set[str]:
    """Generate all sequences that are Hamming distance 2 from input."""
    variants = set()
    bases = 'ACGT'
    seq_len = len(seq)

    for i in range(seq_len):
        for j in range(i+1, seq_len):
            for base_i in bases:
                for base_j in bases:
                    if base_i != seq[i] or base_j != seq[j]:
                        variant = seq[:i] + base_i + seq[i+1:j] + base_j + seq[j+1:]
                        if variant != seq:  # Exclude original
                            variants.add(variant)

    return variants


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two equal-length strings."""
    if len(s1) != len(s2):
        return -1
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def find_parent_for_error_bc1(
    error_bc1: str,
    high_abundance_set: Set[str],
    high_abundance_counts: Dict[str, int],
    max_hamming: int = 2
) -> Optional[Tuple[str, int, int]]:
    """
    Find the most likely parent BC1 for an error BC1.

    Returns:
        (parent_bc1, hamming_distance, parent_total_counts) or None if no parent found
    """
    best_parent = None
    best_distance = max_hamming + 1
    best_counts = 0

    # Check Hamming distance 1 first
    for variant in generate_hamming_1_variants(error_bc1):
        if variant in high_abundance_set:
            parent_counts = high_abundance_counts[variant]
            # Prefer higher-abundance parent if multiple matches
            if best_distance > 1 or (best_distance == 1 and parent_counts > best_counts):
                best_parent = variant
                best_distance = 1
                best_counts = parent_counts

    # If no distance-1 match and we allow distance-2
    if best_parent is None and max_hamming >= 2:
        for variant in generate_hamming_2_variants(error_bc1):
            if variant in high_abundance_set:
                parent_counts = high_abundance_counts[variant]
                if parent_counts > best_counts:
                    best_parent = variant
                    best_distance = 2
                    best_counts = parent_counts

    if best_parent:
        return (best_parent, best_distance, best_counts)
    return None


def identify_error_bc1s(
    bc1_stats: pd.DataFrame,
    min_count_ratio: float = 100.0,
    max_error_count: int = 50,
    max_error_samples: int = 3,
    max_hamming: int = 2,
    n_workers: int = 8
) -> pd.DataFrame:
    """
    Identify BC1s that are likely errors.

    Args:
        bc1_stats: DataFrame with bc1 statistics
        min_count_ratio: Minimum ratio of parent to error counts
        max_error_count: Maximum total counts for a BC1 to be considered an error
        max_error_samples: Maximum samples for a BC1 to be considered an error
        max_hamming: Maximum Hamming distance to consider (1 or 2)
        n_workers: Number of parallel workers

    Returns:
        DataFrame with error->parent mappings
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

    # Create lookup structures
    high_abundance_set = set(high_abundance['bc1'].values)
    high_abundance_counts = dict(zip(high_abundance['bc1'], high_abundance['total_counts']))

    # Find parents for potential errors
    error_mappings = []
    potential_error_list = potential_errors['bc1'].tolist()
    error_counts = dict(zip(potential_errors['bc1'], potential_errors['total_counts']))

    print(f"Searching for parent BC1s (max Hamming distance: {max_hamming})...")

    # Process in chunks for progress reporting
    chunk_size = 10000
    n_chunks = (len(potential_error_list) + chunk_size - 1) // chunk_size

    for i, start in enumerate(range(0, len(potential_error_list), chunk_size)):
        chunk = potential_error_list[start:start + chunk_size]

        for error_bc1 in chunk:
            result = find_parent_for_error_bc1(
                error_bc1, high_abundance_set, high_abundance_counts, max_hamming
            )
            if result:
                parent_bc1, distance, parent_counts = result
                error_count = error_counts[error_bc1]
                count_ratio = parent_counts / error_count if error_count > 0 else float('inf')

                if count_ratio >= min_count_ratio:
                    error_mappings.append({
                        'error_bc1': error_bc1,
                        'parent_bc1': parent_bc1,
                        'hamming_distance': distance,
                        'error_counts': error_count,
                        'parent_counts': parent_counts,
                        'count_ratio': count_ratio
                    })

        if (i + 1) % 10 == 0 or i == n_chunks - 1:
            print(f"  Processed {min(start + chunk_size, len(potential_error_list)):,}/{len(potential_error_list):,} potential errors, found {len(error_mappings):,} mappings")

    result_df = pd.DataFrame(error_mappings)

    print(f"\nIdentified {len(result_df):,} error BC1s")
    if len(result_df) > 0:
        print(f"  Hamming distance 1: {(result_df['hamming_distance'] == 1).sum():,}")
        print(f"  Hamming distance 2: {(result_df['hamming_distance'] == 2).sum():,}")

    return result_df


def main():
    parser = argparse.ArgumentParser(
        description="Identify error BC1s by Hamming distance analysis"
    )
    parser.add_argument(
        "--counts",
        type=Path,
        required=True,
        help="Combined long-format counts file from 03_combine_counts.py"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output path for error BC1 mapping file"
    )
    parser.add_argument(
        "--stats-output",
        type=Path,
        help="Output path for BC1 statistics (optional)"
    )
    parser.add_argument(
        "--max-error-count",
        type=int,
        default=50,
        help="Maximum total counts for a BC1 to be considered an error (default: 50)"
    )
    parser.add_argument(
        "--max-error-samples",
        type=int,
        default=3,
        help="Maximum number of samples for a BC1 to be considered an error (default: 3)"
    )
    parser.add_argument(
        "--min-count-ratio",
        type=float,
        default=100.0,
        help="Minimum ratio of parent to error counts (default: 100)"
    )
    parser.add_argument(
        "--max-hamming",
        type=int,
        default=2,
        choices=[1, 2],
        help="Maximum Hamming distance to consider (default: 2)"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help="Number of parallel workers (default: 8)"
    )

    args = parser.parse_args()

    # Load counts
    print(f"Loading counts from {args.counts}...")
    counts_df = pd.read_csv(args.counts, sep='\t')
    print(f"Loaded {len(counts_df):,} rows")

    # Calculate BC1 statistics
    print("\nCalculating BC1 statistics...")
    bc1_stats = calculate_bc1_statistics(counts_df)
    print(f"Unique BC1s: {len(bc1_stats):,}")
    print(f"Total counts: {bc1_stats['total_counts'].sum():,}")

    # Save stats if requested
    if args.stats_output:
        print(f"\nSaving BC1 statistics to {args.stats_output}...")
        args.stats_output.parent.mkdir(parents=True, exist_ok=True)
        bc1_stats.to_csv(args.stats_output, sep='\t', index=False)

    # Identify error BC1s
    print("\nIdentifying error BC1s...")
    error_mappings = identify_error_bc1s(
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
