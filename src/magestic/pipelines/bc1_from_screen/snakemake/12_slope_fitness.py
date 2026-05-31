#!/usr/bin/env python3
"""
Step 12: Multi-Timepoint Slope Fitness Calculation
===================================================

This script calculates fitness as the slope of log2(relative abundance) vs generation,
using all available timepoints instead of just t0 vs endpoint.

Key Features:
    1. Weighted linear regression across timepoints
    2. Linearity assessment (R², curvature detection)
    3. Robust aggregation from BC1-level to variant-level
    4. Comparison with endpoint-based LFC

Input:
    - Count matrix with samples at multiple timepoints
    - Sample key mapping samples to generations
    - BC1 reference table (for BC1->oligo->variant mapping)
    - Variance model from step E1 (optional, for weighting)

Output:
    - BC1-level slope fitness results
    - Variant-level aggregated slope fitness
    - Summary statistics and QC metrics

Usage:
    python 12_slope_fitness.py \
        --count-matrix counts.tsv \
        --sample-key sample_key.tsv \
        --bc1-reference bc1_reference.tsv \
        --output-bc1 bc1_slope_fitness.tsv \
        --output-variant variant_slope_fitness.tsv \
        --output-summary slope_fitness_summary.tsv \
        --library yL406 \
        --condition "Caffeine" \
        [--variance-model a,b] \
        [--min-timepoints 3] \
        [--min-bc1s 2]

Author: Kevin R. Roy
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd
import numpy as np

from magestic.pipelines.bc1_from_screen.core.slope_fitness import (
    calculate_slope_fitness_dataframe,
    compare_slope_vs_endpoint_fitness
)


def load_count_matrix(filepath: Path) -> pd.DataFrame:
    """Load count matrix (BC1s as rows, samples as columns)."""
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    return df


def parse_sample_key(
    sample_key_path: Path,
    condition: str,
    library: str
) -> Dict[str, int]:
    """
    Parse sample key to extract sample->generation mapping for a specific condition.

    Returns dict mapping sample names to generation numbers.
    """
    sample_key = pd.read_csv(sample_key_path, sep='\t')

    # Filter to specified condition and library
    if 'condition' in sample_key.columns:
        mask = sample_key['condition'] == condition
    elif 'media' in sample_key.columns:
        mask = sample_key['media'] == condition
    else:
        mask = pd.Series([True] * len(sample_key))

    if 'library' in sample_key.columns:
        mask &= sample_key['library'] == library
    elif 'yL' in sample_key.columns:
        mask &= sample_key['yL'] == library

    filtered = sample_key[mask]

    # Extract generation mapping
    generation_mapping = {}

    # Common column names for generation/timepoint
    gen_cols = ['generations', 'generation', 'timepoint', 'gen', 't']
    gen_col = None
    for col in gen_cols:
        if col in filtered.columns:
            gen_col = col
            break

    if gen_col is None:
        raise ValueError(f"Could not find generation column. Available: {filtered.columns.tolist()}")

    # Common column names for sample ID
    sample_cols = ['sample_name', 'sample_id', 'sample', 'sample_number']
    sample_col = None
    for col in sample_cols:
        if col in filtered.columns:
            sample_col = col
            break

    if sample_col is None:
        raise ValueError(f"Could not find sample column. Available: {filtered.columns.tolist()}")

    for _, row in filtered.iterrows():
        sample_name = str(row[sample_col])
        gen = row[gen_col]

        # Handle string generation values like "0", "4g", "t0", etc.
        if isinstance(gen, str):
            gen = gen.lower().replace('g', '').replace('t', '')
            if gen == 'none' or gen == '':
                continue
            gen = int(gen)
        else:
            gen = int(gen)

        generation_mapping[sample_name] = gen

    return generation_mapping


def load_bc1_reference(filepath: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Load BC1 reference table and extract mappings.

    Returns:
        bc1_to_oligo: Dict mapping BC1 sequence to oligo ID
        oligo_to_variant: Dict mapping oligo ID to variant ID
    """
    ref = pd.read_csv(filepath, sep='\t')

    bc1_to_oligo = {}
    oligo_to_variant = {}

    # Identify columns
    bc1_col = None
    for col in ['bc1', 'BC1', 'barcode', 'bc1_sequence']:
        if col in ref.columns:
            bc1_col = col
            break

    oligo_col = None
    for col in ['oligo_id', 'oligo', 'oligo_name', 'sublibrary_oligo_combo']:
        if col in ref.columns:
            oligo_col = col
            break

    variant_col = None
    for col in ['variant_id', 'variant', 'sublibrary_variant_combo', 'variant_name']:
        if col in ref.columns:
            variant_col = col
            break

    if bc1_col and oligo_col:
        bc1_to_oligo = dict(zip(ref[bc1_col], ref[oligo_col]))

    if oligo_col and variant_col:
        # Get unique oligo->variant mappings
        oligo_to_variant = dict(zip(ref[oligo_col], ref[variant_col]))

    return bc1_to_oligo, oligo_to_variant


def parse_variance_model(model_str: str) -> Tuple[float, float]:
    """Parse variance model from string 'a,b'."""
    parts = model_str.split(',')
    if len(parts) != 2:
        raise ValueError(f"Variance model must be 'a,b', got: {model_str}")
    return float(parts[0]), float(parts[1])


def calculate_summary_statistics(
    bc1_fitness: pd.DataFrame,
    variant_fitness: pd.DataFrame
) -> pd.DataFrame:
    """Calculate summary statistics for QC."""
    summary = []

    # BC1-level stats
    if len(bc1_fitness) > 0:
        summary.append({'metric': 'n_bc1s_analyzed', 'value': len(bc1_fitness)})
        summary.append({'metric': 'bc1_mean_slope', 'value': bc1_fitness['slope'].mean()})
        summary.append({'metric': 'bc1_median_slope', 'value': bc1_fitness['slope'].median()})
        summary.append({'metric': 'bc1_std_slope', 'value': bc1_fitness['slope'].std()})
        summary.append({'metric': 'bc1_mean_r_squared', 'value': bc1_fitness['r_squared'].mean()})
        summary.append({'metric': 'bc1_mean_linearity', 'value': bc1_fitness['linearity_score'].mean()})
        summary.append({'metric': 'bc1_pct_monotonic', 'value': bc1_fitness['is_monotonic'].mean() * 100})
        summary.append({'metric': 'bc1_mean_timepoints', 'value': bc1_fitness['n_timepoints'].mean()})

        # Correlation between slope and endpoint LFC
        if 'endpoint_lfc' in bc1_fitness.columns:
            corr = bc1_fitness['slope'].corr(bc1_fitness['endpoint_lfc'])
            summary.append({'metric': 'bc1_slope_endpoint_correlation', 'value': corr})

    # Variant-level stats
    if len(variant_fitness) > 0:
        summary.append({'metric': 'n_variants_analyzed', 'value': len(variant_fitness)})
        summary.append({'metric': 'variant_mean_slope', 'value': variant_fitness['slope_fitness'].mean()})
        summary.append({'metric': 'variant_median_slope', 'value': variant_fitness['slope_fitness'].median()})
        summary.append({'metric': 'variant_std_slope', 'value': variant_fitness['slope_fitness'].std()})
        summary.append({'metric': 'variant_mean_concordance', 'value': variant_fitness['concordance'].mean()})
        summary.append({'metric': 'variant_mean_n_bc1s', 'value': variant_fitness['n_bc1s'].mean()})

        # Significant variants (|z| > 3)
        n_sig = (np.abs(variant_fitness['z_score']) > 3).sum()
        summary.append({'metric': 'variant_n_significant_z3', 'value': n_sig})
        summary.append({'metric': 'variant_pct_significant_z3', 'value': n_sig / len(variant_fitness) * 100})

    return pd.DataFrame(summary)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate multi-timepoint slope fitness'
    )
    parser.add_argument(
        '--count-matrix', type=Path, required=True,
        help='Count matrix (BC1s x samples)'
    )
    parser.add_argument(
        '--sample-key', type=Path, required=True,
        help='Sample key with generation information'
    )
    parser.add_argument(
        '--bc1-reference', type=Path, required=True,
        help='BC1 reference table with oligo/variant mappings'
    )
    parser.add_argument(
        '--output-bc1', type=Path, required=True,
        help='Output BC1-level slope fitness results'
    )
    parser.add_argument(
        '--output-variant', type=Path, required=True,
        help='Output variant-level slope fitness results'
    )
    parser.add_argument(
        '--output-summary', type=Path, required=True,
        help='Output summary statistics'
    )
    parser.add_argument(
        '--library', type=str, required=True,
        help='Library name (e.g., yL406)'
    )
    parser.add_argument(
        '--condition', type=str, required=True,
        help='Condition to analyze (e.g., Caffeine)'
    )
    parser.add_argument(
        '--variance-model', type=str, default=None,
        help='Variance model as "a,b" for Var = a + b/sqrt(count)'
    )
    parser.add_argument(
        '--min-timepoints', type=int, default=3,
        help='Minimum timepoints for slope calculation (default: 3)'
    )
    parser.add_argument(
        '--min-bc1s', type=int, default=2,
        help='Minimum BC1s for variant aggregation (default: 2)'
    )
    parser.add_argument(
        '--endpoint-results', type=Path, default=None,
        help='Optional: endpoint-based DESeq2 results for comparison'
    )

    args = parser.parse_args()

    print(f"=== Step E2: Multi-Timepoint Slope Fitness ===")
    print(f"Library: {args.library}")
    print(f"Condition: {args.condition}")
    print(f"Min timepoints: {args.min_timepoints}")
    print(f"Min BC1s: {args.min_bc1s}")

    # Load data
    print("\nLoading count matrix...")
    count_matrix = load_count_matrix(args.count_matrix)
    print(f"  Loaded {count_matrix.shape[0]} BC1s x {count_matrix.shape[1]} samples")

    print("\nParsing sample key for generation mapping...")
    generation_mapping = parse_sample_key(
        args.sample_key,
        condition=args.condition,
        library=args.library
    )
    print(f"  Found {len(generation_mapping)} samples for {args.condition}")
    print(f"  Generations: {sorted(set(generation_mapping.values()))}")

    print("\nLoading BC1 reference table...")
    bc1_to_oligo, oligo_to_variant = load_bc1_reference(args.bc1_reference)
    print(f"  BC1->oligo mappings: {len(bc1_to_oligo)}")
    print(f"  Oligo->variant mappings: {len(oligo_to_variant)}")

    # Parse variance model if provided
    variance_model = None
    if args.variance_model:
        variance_model = parse_variance_model(args.variance_model)
        print(f"\nUsing variance model: Var = {variance_model[0]:.4f} + {variance_model[1]:.4f}/sqrt(count)")

    # Calculate slope fitness
    print("\nCalculating slope fitness...")
    bc1_fitness, variant_fitness = calculate_slope_fitness_dataframe(
        count_matrix=count_matrix,
        generation_mapping=generation_mapping,
        variance_model=variance_model,
        bc1_to_oligo=bc1_to_oligo,
        oligo_to_variant=oligo_to_variant,
        min_timepoints=args.min_timepoints,
        min_bc1s=args.min_bc1s
    )

    print(f"  BC1-level results: {len(bc1_fitness)}")
    print(f"  Variant-level results: {len(variant_fitness)}")

    # Calculate summary statistics
    print("\nCalculating summary statistics...")
    summary = calculate_summary_statistics(bc1_fitness, variant_fitness)

    # Compare with endpoint results if provided
    if args.endpoint_results and args.endpoint_results.exists():
        print("\nComparing with endpoint-based results...")
        endpoint_results = pd.read_csv(args.endpoint_results, sep='\t')
        comparison = compare_slope_vs_endpoint_fitness(
            variant_fitness, endpoint_results,
            merge_on='variant_id'
        )
        if len(comparison) > 0 and 'sign_concordant' in comparison.columns:
            concordance_rate = comparison['sign_concordant'].mean() * 100
            summary = pd.concat([summary, pd.DataFrame([
                {'metric': 'slope_endpoint_sign_concordance_pct', 'value': concordance_rate}
            ])], ignore_index=True)
            print(f"  Sign concordance with endpoint: {concordance_rate:.1f}%")

    # Create output directories
    args.output_bc1.parent.mkdir(parents=True, exist_ok=True)
    args.output_variant.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)

    # Save results
    print("\nSaving results...")
    bc1_fitness.to_csv(args.output_bc1, sep='\t', index=False)
    print(f"  BC1 results: {args.output_bc1}")

    variant_fitness.to_csv(args.output_variant, sep='\t', index=False)
    print(f"  Variant results: {args.output_variant}")

    summary.to_csv(args.output_summary, sep='\t', index=False)
    print(f"  Summary: {args.output_summary}")

    # Print key metrics
    print("\n=== Summary ===")
    for _, row in summary.iterrows():
        print(f"  {row['metric']}: {row['value']:.4f}" if isinstance(row['value'], float) else f"  {row['metric']}: {row['value']}")

    print("\nDone!")


if __name__ == '__main__':
    main()
