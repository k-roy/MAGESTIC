#!/usr/bin/env python3
"""
Step 03b: Error BC1 Detection (v3 - Reference-based)

Uses the BC1 reference table as ground truth to identify errors:

1. Build HD=1 lookup from BC1 reference table (~1M designed BC1s)
2. For each observed BC1:
   - Exact match to ref → real BC1
   - HD=1 match to ref → error of known BC1
3. Remaining non-matching BC1s: check HD=1 among themselves
   (catches errors of novel/unexpected BC1s)

This is more biologically meaningful than pure abundance-based detection.

Author: Kevin R. Roy
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Set, Tuple, Optional
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


def analyze_bc1_length_distribution(bc1_stats: pd.DataFrame, expected_length: int = 20) -> dict:
    """
    Analyze BC1 length distribution and return statistics.

    BC1s are expected to be 20bp, but synthesis/sequencing errors can produce
    19bp or 21bp variants. These length variants cannot be matched by HD=1
    to 20bp reference BC1s.

    Returns dict with length distribution stats.
    """
    bc1_stats = bc1_stats.copy()
    bc1_stats['bc1_length'] = bc1_stats['bc1'].str.len()

    length_counts = bc1_stats.groupby('bc1_length').agg(
        n_unique_bc1s=('bc1', 'count'),
        total_reads=('total_counts', 'sum')
    ).reset_index()

    total_unique = len(bc1_stats)
    total_reads = bc1_stats['total_counts'].sum()

    print(f"\n=== BC1 Length Distribution ===")
    print(f"{'Length':>8} {'Unique BC1s':>15} {'% BC1s':>10} {'Total Reads':>15} {'% Reads':>10}")
    print("-" * 60)

    length_stats = {}
    for _, row in length_counts.iterrows():
        length = int(row['bc1_length'])
        n_bc1s = int(row['n_unique_bc1s'])
        n_reads = int(row['total_reads'])
        pct_bc1s = 100 * n_bc1s / total_unique if total_unique > 0 else 0
        pct_reads = 100 * n_reads / total_reads if total_reads > 0 else 0

        marker = " <-- expected" if length == expected_length else ""
        print(f"{length:>6}bp {n_bc1s:>15,} {pct_bc1s:>9.1f}% {n_reads:>15,} {pct_reads:>9.1f}%{marker}")

        length_stats[length] = {
            'unique_bc1s': n_bc1s,
            'total_reads': n_reads,
            'pct_bc1s': pct_bc1s,
            'pct_reads': pct_reads
        }

    # Count length variants
    length_variant_bc1s = bc1_stats[bc1_stats['bc1_length'] != expected_length]
    n_length_variants = len(length_variant_bc1s)
    length_variant_reads = length_variant_bc1s['total_counts'].sum()

    print(f"\nSummary:")
    print(f"  BC1s with expected length ({expected_length}bp): {total_unique - n_length_variants:,} ({100 * (total_unique - n_length_variants) / total_unique:.1f}%)")
    print(f"  Length variant BC1s: {n_length_variants:,} ({100 * n_length_variants / total_unique:.1f}%)")
    print(f"  Length variant reads: {length_variant_reads:,} ({100 * length_variant_reads / total_reads:.1f}%)")
    print(f"\nNote: Length variant BC1s cannot be matched by HD=1 to reference.")
    print(f"      They will appear in 'remaining_non_matches' below.")

    return {
        'expected_length': expected_length,
        'by_length': length_stats,
        'total_unique_bc1s': total_unique,
        'total_reads': total_reads,
        'length_variant_bc1s': n_length_variants,
        'length_variant_reads': int(length_variant_reads)
    }


def generate_hd1_variants(seq: str) -> list:
    """Generate all HD=1 variants of a sequence."""
    variants = []
    bases = 'ACGT'
    for i in range(len(seq)):
        for base in bases:
            if base != seq[i]:
                variants.append(seq[:i] + base + seq[i+1:])
    return variants


def build_ref_hd1_lookup(ref_bc1s: Set[str]) -> Dict[str, str]:
    """
    Build HD=1 lookup from reference BC1s.

    Returns dict: hd1_variant -> reference_bc1

    For collisions (rare: two ref BC1s have overlapping HD=1 neighborhoods),
    we keep the first one encountered.
    """
    print(f"Building HD=1 lookup from {len(ref_bc1s):,} reference BC1s...")

    lookup = {}
    ref_list = list(ref_bc1s)

    chunk_size = 50000
    for chunk_start in range(0, len(ref_list), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(ref_list))

        for ref_bc1 in ref_list[chunk_start:chunk_end]:
            for variant in generate_hd1_variants(ref_bc1):
                # Don't overwrite - keep first (arbitrary but consistent)
                if variant not in lookup:
                    lookup[variant] = ref_bc1

        if chunk_end % 100000 == 0 or chunk_end == len(ref_list):
            print(f"  Processed {chunk_end:,}/{len(ref_list):,} reference BC1s, lookup size: {len(lookup):,}")

    print(f"  Final HD=1 lookup: {len(lookup):,} variants")
    return lookup


def identify_errors_reference_based(
    bc1_stats: pd.DataFrame,
    ref_bc1s: Set[str],
    min_count_ratio: float = 100.0,
    max_error_count: int = 50,
    max_error_samples: int = 3
) -> Tuple[pd.DataFrame, dict]:
    """
    Reference-based error detection.

    Returns:
        - DataFrame with error->parent mappings
        - dict with statistics
    """
    observed_bc1s = set(bc1_stats['bc1'].values)
    bc1_counts = dict(zip(bc1_stats['bc1'], bc1_stats['total_counts']))
    bc1_samples = dict(zip(bc1_stats['bc1'], bc1_stats['n_samples']))

    print(f"\n=== Reference-Based Error Detection ===")
    print(f"Observed BC1s: {len(observed_bc1s):,}")
    print(f"Reference BC1s: {len(ref_bc1s):,}")

    # Step 1: Classify observed BC1s
    exact_matches = observed_bc1s & ref_bc1s
    non_matches = observed_bc1s - ref_bc1s

    print(f"\nStep 1: Exact match classification")
    print(f"  Exact matches to reference: {len(exact_matches):,}")
    print(f"  Non-matching BC1s: {len(non_matches):,}")

    # Step 2: Build HD=1 lookup from reference and check non-matches
    ref_hd1_lookup = build_ref_hd1_lookup(ref_bc1s)

    print(f"\nStep 2: Checking non-matches against reference HD=1...")

    error_mappings = []
    hd1_to_ref_matches = set()
    remaining_non_matches = set()

    for bc1 in non_matches:
        if bc1 in ref_hd1_lookup:
            parent_bc1 = ref_hd1_lookup[bc1]
            error_count = bc1_counts[bc1]
            parent_count = bc1_counts.get(parent_bc1, 0)
            n_samples = bc1_samples[bc1]

            # Apply error criteria
            if error_count <= max_error_count and n_samples <= max_error_samples:
                if parent_count > 0:
                    count_ratio = parent_count / error_count if error_count > 0 else float('inf')
                    if count_ratio >= min_count_ratio:
                        error_mappings.append({
                            'error_bc1': bc1,
                            'parent_bc1': parent_bc1,
                            'hamming_distance': 1,
                            'error_counts': error_count,
                            'parent_counts': parent_count,
                            'count_ratio': count_ratio,
                            'error_type': 'ref_hd1'
                        })
                        hd1_to_ref_matches.add(bc1)
                        continue

            hd1_to_ref_matches.add(bc1)
        else:
            remaining_non_matches.add(bc1)

    print(f"  HD=1 matches to reference: {len(hd1_to_ref_matches):,}")
    print(f"  Errors identified (ref HD=1): {len(error_mappings):,}")
    print(f"  Remaining non-matches: {len(remaining_non_matches):,}")

    # Step 3: Check HD=1 among remaining non-matches
    # These might be novel BC1s with their own errors
    print(f"\nStep 3: Checking HD=1 among remaining non-matches...")

    if len(remaining_non_matches) > 0:
        # Split into potential errors (low count) and potential parents (higher count)
        remaining_stats = bc1_stats[bc1_stats['bc1'].isin(remaining_non_matches)].copy()

        potential_errors = remaining_stats[
            (remaining_stats['total_counts'] <= max_error_count) &
            (remaining_stats['n_samples'] <= max_error_samples)
        ]
        potential_parents = remaining_stats[remaining_stats['total_counts'] > max_error_count]

        print(f"  Potential errors among remaining: {len(potential_errors):,}")
        print(f"  Potential parents among remaining: {len(potential_parents):,}")

        if len(potential_parents) > 0 and len(potential_errors) > 0:
            # Build HD=1 lookup from potential parents
            parent_bc1s = set(potential_parents['bc1'].values)
            parent_counts = dict(zip(potential_parents['bc1'], potential_parents['total_counts']))
            parent_hd1_lookup = {}

            for parent_bc1 in parent_bc1s:
                for variant in generate_hd1_variants(parent_bc1):
                    if variant not in parent_hd1_lookup:
                        parent_hd1_lookup[variant] = parent_bc1

            # Check potential errors against this lookup
            novel_errors = 0
            for _, row in potential_errors.iterrows():
                bc1 = row['bc1']
                if bc1 in parent_hd1_lookup:
                    parent_bc1 = parent_hd1_lookup[bc1]
                    error_count = row['total_counts']
                    parent_count = parent_counts[parent_bc1]
                    count_ratio = parent_count / error_count if error_count > 0 else float('inf')

                    if count_ratio >= min_count_ratio:
                        error_mappings.append({
                            'error_bc1': bc1,
                            'parent_bc1': parent_bc1,
                            'hamming_distance': 1,
                            'error_counts': error_count,
                            'parent_counts': parent_count,
                            'count_ratio': count_ratio,
                            'error_type': 'novel_hd1'
                        })
                        novel_errors += 1

            print(f"  Novel errors identified (HD=1 among non-ref): {novel_errors:,}")

    result_df = pd.DataFrame(error_mappings)

    stats = {
        'observed_bc1s': len(observed_bc1s),
        'reference_bc1s': len(ref_bc1s),
        'exact_matches': len(exact_matches),
        'hd1_to_ref': len(hd1_to_ref_matches),
        'remaining_non_matches': len(remaining_non_matches),
        'total_errors': len(result_df),
        'ref_hd1_errors': len(result_df[result_df['error_type'] == 'ref_hd1']) if len(result_df) > 0 else 0,
        'novel_hd1_errors': len(result_df[result_df['error_type'] == 'novel_hd1']) if len(result_df) > 0 else 0
    }

    return result_df, stats


def main():
    parser = argparse.ArgumentParser(
        description="Reference-based error BC1 detection (v3)"
    )
    parser.add_argument("--counts", type=Path, required=True,
                        help="Combined long-format counts file")
    parser.add_argument("--bc1-reference", type=Path, required=True,
                        help="BC1 reference table with designed BC1s")
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

    args = parser.parse_args()

    # Load counts
    print(f"Loading counts from {args.counts}...")
    counts_df = pd.read_csv(args.counts, sep='\t')
    print(f"Loaded {len(counts_df):,} rows")

    # Auto-detect count column
    count_col = detect_count_column(counts_df)
    print(f"Using count column: '{count_col}'")

    # Load reference BC1s
    print(f"\nLoading BC1 reference from {args.bc1_reference}...")
    ref_df = pd.read_csv(args.bc1_reference, sep='\t')
    ref_bc1s = set(ref_df['bc1'].values)
    print(f"Loaded {len(ref_bc1s):,} reference BC1s")

    # Calculate BC1 statistics
    print("\nCalculating BC1 statistics...")
    bc1_stats = calculate_bc1_statistics(counts_df, count_col=count_col)
    print(f"Unique BC1s: {len(bc1_stats):,}")
    print(f"Total counts: {bc1_stats['total_counts'].sum():,}")

    # Analyze BC1 length distribution (important for understanding non-matches)
    length_stats = analyze_bc1_length_distribution(bc1_stats, expected_length=20)

    # Save stats if requested
    if args.stats_output:
        print(f"\nSaving BC1 statistics to {args.stats_output}...")
        args.stats_output.parent.mkdir(parents=True, exist_ok=True)
        bc1_stats.to_csv(args.stats_output, sep='\t', index=False)

    # Reference-based error detection
    error_mappings, stats = identify_errors_reference_based(
        bc1_stats,
        ref_bc1s,
        min_count_ratio=args.min_count_ratio,
        max_error_count=args.max_error_count,
        max_error_samples=args.max_error_samples
    )

    # Save error mappings
    print(f"\nSaving error mappings to {args.output}...")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    error_mappings.to_csv(args.output, sep='\t', index=False)

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Observed BC1s: {stats['observed_bc1s']:,}")
    print(f"Reference BC1s: {stats['reference_bc1s']:,}")
    print(f"Exact matches to reference: {stats['exact_matches']:,}")
    print(f"HD=1 matches to reference: {stats['hd1_to_ref']:,}")
    print(f"Remaining non-matches: {stats['remaining_non_matches']:,}")

    # Add context about length variants
    if length_stats['length_variant_bc1s'] > 0:
        print(f"\n  Note: {length_stats['length_variant_bc1s']:,} length variant BC1s ({length_stats['length_variant_reads']:,} reads)")
        print(f"        contribute to remaining non-matches (HD=1 requires equal length)")

    print(f"\nTotal errors identified: {stats['total_errors']:,}")
    print(f"  - Errors of reference BC1s: {stats['ref_hd1_errors']:,}")
    print(f"  - Errors of novel BC1s: {stats['novel_hd1_errors']:,}")

    if stats['total_errors'] > 0:
        total_error_counts = error_mappings['error_counts'].sum()
        print(f"\nTotal counts in error BC1s: {total_error_counts:,}")

    print("\nDone!")


if __name__ == "__main__":
    main()
