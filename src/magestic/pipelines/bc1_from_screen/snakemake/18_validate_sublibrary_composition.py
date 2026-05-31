#!/usr/bin/env python3
"""
Step 14: Validate Yeast Sublibrary Composition

Validates BC1 yeast_sublibrary assignments and verifies that consolidated libraries
(e.g., yL406) match expected mixing proportions.

Note: "yeast_sublibrary" refers to transformation batch (trafo_ID level), which defines
the nuclease-PAM combination (e.g., SpG_NGNG, SpCas9_NGG). This is distinct from:
- "guide_donor_oligo_subpool": oligo design pools with different SPS sequences
- "step_1_plasmid_library": V guide-donor-bc0 libraries (e.g., V629, V631, V634)

This script:
1. Loads BC1 counts from individual yeast_sublibrary samples
2. Loads BC1 counts from the consolidated library (e.g., yL406)
3. Cross-references with BC1 reference table assignments
4. Validates that BC1s appear in their assigned yeast_sublibraries
5. Compares observed vs expected mixing proportions

Usage:
    python 14_validate_sublibrary_composition.py \\
        --bc1-counts-long bc1_counts_long_format.tsv \\
        --sample-key sample_key.tsv \\
        --bc1-reference bc1_reference_table.tsv \\
        --consolidated-library yL406 \\
        --yeast-sublibraries yL386,yL387,yL391,yL398,yL401 \\
        --expected-proportions 3,3,3,1,1 \\
        --output-dir validation_output/

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate yeast_sublibrary composition and mixing proportions"
    )
    parser.add_argument(
        "--bc1-counts-long", type=Path, required=True,
        help="Long-format BC1 counts TSV (bc1, sample_name, counts columns)"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with yL column for yeast_sublibrary identification"
    )
    parser.add_argument(
        "--bc1-reference", type=Path, required=True,
        help="BC1 reference table with yeast_sublibrary assignments"
    )
    parser.add_argument(
        "--consolidated-library", type=str, required=True,
        help="Name of consolidated library (e.g., yL406)"
    )
    parser.add_argument(
        "--yeast-sublibraries", type=str, required=True,
        help="Comma-separated list of component yeast_sublibraries (e.g., yL386,yL387,yL391,yL398,yL401)"
    )
    parser.add_argument(
        "--expected-proportions", type=str, required=True,
        help="Comma-separated expected proportions (e.g., 3,3,3,1,1 for volumes)"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for validation results"
    )
    parser.add_argument(
        "--filter-gDNA-plate", type=str, default=None,
        help="Optional: filter to specific gDNA plate (e.g., 20240706_MAGESTIC...)"
    )
    parser.add_argument(
        "--filter-condition", type=str, default=None,
        help="Optional: filter to specific condition (e.g., 'YPD+G418+1 mM 5FC outgrowth')"
    )
    parser.add_argument(
        "--min-counts", type=int, default=5,
        help="Minimum counts to consider a BC1 present (default: 5)"
    )
    return parser.parse_args()


def load_yeast_sublibrary_bc1s(bc1_counts: pd.DataFrame,
                                sample_key: pd.DataFrame,
                                yeast_sublibraries: List[str],
                                min_counts: int = 5,
                                filter_gDNA_plate: str = None,
                                filter_condition: str = None) -> Dict[str, set]:
    """
    Load BC1s present in each yeast_sublibrary sample.

    Returns:
        {yeast_sublibrary: set of bc1s present in that yeast_sublibrary}
    """
    # Get sample names for each yeast_sublibrary
    sample_key_filt = sample_key.copy()

    if filter_gDNA_plate:
        if 'gDNA_plate_name' in sample_key_filt.columns:
            sample_key_filt = sample_key_filt[
                sample_key_filt['gDNA_plate_name'].str.contains(filter_gDNA_plate, na=False)
            ]

    if filter_condition:
        if 'condition' in sample_key_filt.columns:
            sample_key_filt = sample_key_filt[
                sample_key_filt['condition'].str.contains(filter_condition, na=False)
            ]

    yeast_sublibrary_bc1s = {}

    for yeast_sublib in yeast_sublibraries:
        # Get samples for this yeast_sublibrary
        if 'yL' in sample_key_filt.columns:
            yeast_sublib_samples = sample_key_filt[sample_key_filt['yL'] == yeast_sublib]['sample_name'].tolist()
        else:
            logger.warning(f"No 'yL' column in sample key, cannot filter by yeast_sublibrary")
            yeast_sublib_samples = []

        if not yeast_sublib_samples:
            logger.warning(f"  No samples found for yeast_sublibrary {yeast_sublib}")
            yeast_sublibrary_bc1s[yeast_sublib] = set()
            continue

        # Get BC1s with sufficient counts in these samples
        counts_filt = bc1_counts[bc1_counts['sample_name'].isin(yeast_sublib_samples)]

        # Aggregate counts per BC1
        bc1_totals = counts_filt.groupby('bc1')['counts'].sum()
        bc1s_present = set(bc1_totals[bc1_totals >= min_counts].index)

        yeast_sublibrary_bc1s[yeast_sublib] = bc1s_present
        logger.info(f"  {yeast_sublib}: {len(yeast_sublib_samples)} samples, {len(bc1s_present):,} BC1s present")

    return yeast_sublibrary_bc1s


def load_consolidated_bc1s(bc1_counts: pd.DataFrame,
                            sample_key: pd.DataFrame,
                            library: str,
                            min_counts: int = 5,
                            filter_gDNA_plate: str = None,
                            filter_condition: str = None) -> Tuple[set, pd.Series]:
    """
    Load BC1s present in the consolidated library.

    Returns:
        (set of bc1s, Series of bc1 -> total counts)
    """
    sample_key_filt = sample_key.copy()

    if filter_gDNA_plate:
        if 'gDNA_plate_name' in sample_key_filt.columns:
            sample_key_filt = sample_key_filt[
                sample_key_filt['gDNA_plate_name'].str.contains(filter_gDNA_plate, na=False)
            ]

    if filter_condition:
        if 'condition' in sample_key_filt.columns:
            sample_key_filt = sample_key_filt[
                sample_key_filt['condition'].str.contains(filter_condition, na=False)
            ]

    # Get samples for consolidated library
    if 'yL' in sample_key_filt.columns:
        lib_samples = sample_key_filt[sample_key_filt['yL'] == library]['sample_name'].tolist()
    else:
        lib_samples = []

    logger.info(f"  {library}: {len(lib_samples)} samples found")

    if not lib_samples:
        return set(), pd.Series(dtype=float)

    # Get BC1 counts
    counts_filt = bc1_counts[bc1_counts['sample_name'].isin(lib_samples)]
    bc1_totals = counts_filt.groupby('bc1')['counts'].sum()

    bc1s_present = set(bc1_totals[bc1_totals >= min_counts].index)
    logger.info(f"  {library}: {len(bc1s_present):,} BC1s with counts >= {min_counts}")

    return bc1s_present, bc1_totals


def validate_assignments(yeast_sublibrary_bc1s: Dict[str, set],
                          consolidated_bc1s: set,
                          bc1_reference: pd.DataFrame,
                          yeast_sublibraries: List[str]) -> pd.DataFrame:
    """
    Validate BC1 yeast_sublibrary assignments.

    For each BC1 in the consolidated library:
    - Check which yeast_sublibrary(ies) it appears in
    - Compare with assigned yeast_sublibrary from reference table
    - Flag mismatches and multi-yeast_sublibrary BC1s
    """
    results = []

    # Get yeast_sublibrary column from reference table
    yeast_sublib_col = None
    for col in ['yeast_sublibrary', 'sublibrary', 'yL', 'library']:
        if col in bc1_reference.columns:
            yeast_sublib_col = col
            break

    if yeast_sublib_col is None:
        logger.warning("No yeast_sublibrary column found in reference table")
        return pd.DataFrame()

    # Create lookup
    bc1_to_assigned = bc1_reference.set_index('bc1')[yeast_sublib_col].to_dict() if 'bc1' in bc1_reference.columns else {}

    for bc1 in consolidated_bc1s:
        # Check which yeast_sublibraries this BC1 appears in
        observed_yeast_sublibs = [sl for sl in yeast_sublibraries if bc1 in yeast_sublibrary_bc1s.get(sl, set())]

        # Get assigned yeast_sublibrary
        assigned = bc1_to_assigned.get(bc1, None)

        # Determine status
        if len(observed_yeast_sublibs) == 0:
            status = "NOT_IN_YEAST_SUBLIBRARIES"
        elif len(observed_yeast_sublibs) == 1:
            if assigned == observed_yeast_sublibs[0]:
                status = "MATCH"
            elif assigned is None:
                status = "UNASSIGNED"
            else:
                status = "MISMATCH"
        else:
            status = "MULTI_YEAST_SUBLIBRARY"

        results.append({
            'bc1': bc1,
            'assigned_yeast_sublibrary': assigned,
            'observed_yeast_sublibraries': ','.join(observed_yeast_sublibs) if observed_yeast_sublibs else 'none',
            'n_observed': len(observed_yeast_sublibs),
            'status': status
        })

    return pd.DataFrame(results)


def calculate_proportions(yeast_sublibrary_bc1s: Dict[str, set],
                           consolidated_counts: pd.Series,
                           yeast_sublibraries: List[str],
                           expected_proportions: List[float]) -> pd.DataFrame:
    """
    Calculate observed vs expected yeast_sublibrary proportions.

    Uses counts from the consolidated library, grouped by which yeast_sublibrary
    each BC1 belongs to.
    """
    # Calculate total counts per yeast_sublibrary
    observed_counts = {}
    for yeast_sublib in yeast_sublibraries:
        bc1s = yeast_sublibrary_bc1s.get(yeast_sublib, set())
        # Sum counts for BC1s that belong to this yeast_sublibrary
        yeast_sublib_counts = consolidated_counts[consolidated_counts.index.isin(bc1s)].sum()
        observed_counts[yeast_sublib] = yeast_sublib_counts

    total_observed = sum(observed_counts.values())

    # Calculate proportions
    total_expected = sum(expected_proportions)
    expected_fractions = [p / total_expected for p in expected_proportions]

    results = []
    for i, yeast_sublib in enumerate(yeast_sublibraries):
        observed = observed_counts.get(yeast_sublib, 0)
        observed_frac = observed / total_observed if total_observed > 0 else 0
        expected_frac = expected_fractions[i]

        # Calculate deviation
        deviation = observed_frac - expected_frac
        fold_change = observed_frac / expected_frac if expected_frac > 0 else np.inf

        results.append({
            'yeast_sublibrary': yeast_sublib,
            'expected_proportion': expected_proportions[i],
            'expected_fraction': expected_frac,
            'observed_counts': observed,
            'observed_fraction': observed_frac,
            'deviation': deviation,
            'fold_change': fold_change,
            'n_bc1s': len(yeast_sublibrary_bc1s.get(yeast_sublib, set()))
        })

    return pd.DataFrame(results)


def create_proportion_plot(proportions_df: pd.DataFrame, output_path: Path):
    """Create bar chart comparing observed vs expected proportions."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(10, 6))

        x = range(len(proportions_df))
        width = 0.35

        bars1 = ax.bar([i - width/2 for i in x],
                       proportions_df['expected_fraction'] * 100,
                       width, label='Expected', color='steelblue', alpha=0.7)
        bars2 = ax.bar([i + width/2 for i in x],
                       proportions_df['observed_fraction'] * 100,
                       width, label='Observed', color='coral', alpha=0.7)

        ax.set_ylabel('Proportion (%)')
        ax.set_xlabel('Yeast Sublibrary')
        ax.set_title('Yeast Sublibrary Composition: Expected vs Observed')
        ax.set_xticks(x)
        ax.set_xticklabels(proportions_df['yeast_sublibrary'])
        ax.legend()

        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            ax.annotate(f'{height:.1f}%',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=8)

        for bar in bars2:
            height = bar.get_height()
            ax.annotate(f'{height:.1f}%',
                       xy=(bar.get_x() + bar.get_width() / 2, height),
                       xytext=(0, 3),
                       textcoords="offset points",
                       ha='center', va='bottom', fontsize=8)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved plot to {output_path}")

    except ImportError:
        logger.warning("matplotlib not available, skipping plot")


def main():
    args = parse_args()

    # Parse yeast_sublibrary and proportion lists
    yeast_sublibraries = [s.strip() for s in args.yeast_sublibraries.split(',')]
    expected_proportions = [float(p.strip()) for p in args.expected_proportions.split(',')]

    if len(yeast_sublibraries) != len(expected_proportions):
        logger.error(f"Number of yeast_sublibraries ({len(yeast_sublibraries)}) must match "
                    f"number of proportions ({len(expected_proportions)})")
        return 1

    logger.info("Step 14: Validating Yeast Sublibrary Composition")
    logger.info(f"  Consolidated library: {args.consolidated_library}")
    logger.info(f"  Component yeast_sublibraries: {yeast_sublibraries}")
    logger.info(f"  Expected proportions: {expected_proportions}")

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    logger.info("Loading BC1 counts (long format)...")
    bc1_counts = pd.read_csv(args.bc1_counts_long, sep='\t')

    # Normalize column names
    col_map = {}
    for col in bc1_counts.columns:
        if col.lower() in ['bc1', 'barcode']:
            col_map[col] = 'bc1'
        elif col.lower() in ['sample_name', 'sample']:
            col_map[col] = 'sample_name'
        elif col.lower() in ['counts', 'count', 'sample_counts']:
            col_map[col] = 'counts'
    bc1_counts = bc1_counts.rename(columns=col_map)

    logger.info(f"  Loaded {len(bc1_counts):,} count records")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep='\t')
    logger.info(f"  Loaded {len(sample_key):,} sample entries")

    logger.info("Loading BC1 reference table...")
    bc1_reference = pd.read_csv(args.bc1_reference, sep='\t')
    logger.info(f"  Loaded {len(bc1_reference):,} BC1 reference entries")

    # Load BC1s for each yeast_sublibrary
    logger.info("Loading BC1s for each yeast_sublibrary...")
    yeast_sublibrary_bc1s = load_yeast_sublibrary_bc1s(
        bc1_counts, sample_key, yeast_sublibraries,
        min_counts=args.min_counts,
        filter_gDNA_plate=args.filter_gDNA_plate,
        filter_condition=args.filter_condition
    )

    # Load BC1s for consolidated library
    logger.info(f"Loading BC1s for consolidated library {args.consolidated_library}...")
    consolidated_bc1s, consolidated_counts = load_consolidated_bc1s(
        bc1_counts, sample_key, args.consolidated_library,
        min_counts=args.min_counts,
        filter_gDNA_plate=args.filter_gDNA_plate,
        filter_condition=args.filter_condition
    )

    # Validate assignments
    logger.info("Validating BC1 yeast_sublibrary assignments...")
    validation_df = validate_assignments(
        yeast_sublibrary_bc1s, consolidated_bc1s, bc1_reference, yeast_sublibraries
    )

    if len(validation_df) > 0:
        # Summary statistics
        status_counts = validation_df['status'].value_counts()
        logger.info("Assignment validation summary:")
        for status, count in status_counts.items():
            pct = 100 * count / len(validation_df)
            logger.info(f"  {status}: {count:,} ({pct:.1f}%)")

        # Save validation results
        validation_file = args.output_dir / "yeast_sublibrary_assignment_validation.tsv"
        validation_df.to_csv(validation_file, sep='\t', index=False)
        logger.info(f"  Saved validation to {validation_file}")

    # Calculate proportions
    logger.info("Calculating yeast_sublibrary proportions...")
    proportions_df = calculate_proportions(
        yeast_sublibrary_bc1s, consolidated_counts, yeast_sublibraries, expected_proportions
    )

    logger.info("Proportion comparison:")
    for _, row in proportions_df.iterrows():
        exp = row['expected_fraction'] * 100
        obs = row['observed_fraction'] * 100
        fc = row['fold_change']
        logger.info(f"  {row['yeast_sublibrary']}: expected={exp:.1f}%, observed={obs:.1f}%, fold_change={fc:.2f}")

    # Save proportions
    proportions_file = args.output_dir / "yeast_sublibrary_proportions.tsv"
    proportions_df.to_csv(proportions_file, sep='\t', index=False)
    logger.info(f"  Saved proportions to {proportions_file}")

    # Create plot
    plot_file = args.output_dir / "yeast_sublibrary_proportions.png"
    create_proportion_plot(proportions_df, plot_file)

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("VALIDATION SUMMARY")
    logger.info("=" * 60)

    # Check if proportions match expectations (within 20% relative deviation)
    max_deviation = proportions_df['fold_change'].apply(lambda x: abs(x - 1)).max()
    if max_deviation < 0.2:
        logger.info("Proportions match expectations (all within 20% of expected)")
    else:
        worst_yeast_sublib = proportions_df.loc[
            proportions_df['fold_change'].apply(lambda x: abs(x - 1)).idxmax(), 'yeast_sublibrary'
        ]
        logger.warning(f"Proportions deviate from expectations (worst: {worst_yeast_sublib}, {max_deviation:.1%} deviation)")

    if len(validation_df) > 0:
        match_rate = (validation_df['status'] == 'MATCH').sum() / len(validation_df) * 100
        logger.info(f"Assignment match rate: {match_rate:.1f}%")

    logger.info("Step 14: Complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
