"""
Guide-based disambiguation for multi-PAM oligo assignments.

PROBLEM:
When multiple oligos target the same variant with different PAM sites, they have
IDENTICAL donor sequences. The donor-only matching cannot distinguish between them
and arbitrarily picks the first match. This results in "imperfect guide" assignments
where the guide sequence doesn't match the assigned oligo.

SOLUTION:
Use the guide sequence to correctly identify which PAM variant was actually used.
This script:
1. Identifies BC1s with misassigned oligos (imperfect guide but perfect donor)
2. Finds alternative oligos with the same donor
3. Reassigns to the oligo whose guide best matches the observed guide

IMPACT:
In the 20250912_SpG_QTL_screen, this fixes ~10,400 BC1s (47% of "imperfect guide" cases),
reducing guide mismatches from 10-16 down to 1-2.

Author: Kevin Roy / Claude
Date: 2026-02-04
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
from pathlib import Path
import pandas as pd
import numpy as np
import logging
from collections import defaultdict

logger = logging.getLogger(__name__)


@dataclass
class ReassignmentResult:
    """Result of attempting to reassign an oligo based on guide."""
    bc1: str
    original_oligo: str
    original_guide_distance: int

    reassigned: bool = False
    new_oligo: Optional[str] = None
    new_guide_distance: int = -1

    # Additional info
    num_alternatives: int = 0
    same_variant: bool = True  # All alternatives target same variant


@dataclass
class MultiPamOligoLookup:
    """
    Lookup tables for multi-PAM oligo disambiguation.

    Groups oligos by their middle donor sequence (which is identical for
    oligos targeting the same variant with different PAMs).
    """
    # donor (middle 60bp) -> list of (oligo_name, guide)
    donor_to_oligo_guides: Dict[str, List[Tuple[str, str]]] = field(default_factory=dict)

    # oligo_name -> guide (for reference)
    oligo_to_guide: Dict[str, str] = field(default_factory=dict)

    # Statistics
    num_oligos: int = 0
    num_multi_pam_donors: int = 0
    num_oligos_in_multi_pam_groups: int = 0


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings of equal length."""
    if len(s1) != len(s2):
        return max(len(s1), len(s2))  # Penalize length mismatch
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def load_oligo_designs_for_disambiguation(
    oligo_files: List[Path],
    middle_donor_start: int = 30,
    middle_donor_end: int = 90,
) -> MultiPamOligoLookup:
    """
    Load oligo design files and build lookup for multi-PAM disambiguation.

    Args:
        oligo_files: List of paths to oligo design TSV files
        middle_donor_start: Start position of middle donor in 129bp donor (default 30)
        middle_donor_end: End position of middle donor (default 90, giving 60bp)

    Returns:
        MultiPamOligoLookup with donor -> oligo+guide mapping
    """
    lookup = MultiPamOligoLookup()

    GUIDE_START = 20
    GUIDE_END = 40
    DONOR_START = 51
    DONOR_END = 180

    for filepath in oligo_files:
        if not filepath.exists():
            logger.warning(f"Oligo file not found: {filepath}")
            continue

        logger.info(f"Loading oligo designs from {filepath.name}")
        df = pd.read_csv(filepath, sep='\t')

        for _, row in df.iterrows():
            oligo_name = row['oligo_name']
            oligo_seq = row['oligo_seq']

            # Extract guide and donor
            guide = oligo_seq[GUIDE_START:GUIDE_END].upper()
            donor = oligo_seq[DONOR_START:DONOR_END].upper()

            # Get middle donor for grouping
            middle_donor = donor[middle_donor_start:middle_donor_end]

            # Add to lookups
            lookup.oligo_to_guide[oligo_name] = guide

            if middle_donor not in lookup.donor_to_oligo_guides:
                lookup.donor_to_oligo_guides[middle_donor] = []
            lookup.donor_to_oligo_guides[middle_donor].append((oligo_name, guide))

            lookup.num_oligos += 1

    # Calculate statistics
    for donor, oligos in lookup.donor_to_oligo_guides.items():
        if len(oligos) > 1:
            lookup.num_multi_pam_donors += 1
            lookup.num_oligos_in_multi_pam_groups += len(oligos)

    logger.info(f"  Loaded {lookup.num_oligos:,} oligos")
    logger.info(f"  Multi-PAM donors: {lookup.num_multi_pam_donors:,}")
    logger.info(f"  Oligos in multi-PAM groups: {lookup.num_oligos_in_multi_pam_groups:,}")

    return lookup


def find_best_oligo_by_guide(
    observed_guide: str,
    middle_donor: str,
    lookup: MultiPamOligoLookup,
    max_guide_distance: int = 2,
) -> Tuple[Optional[str], int, int]:
    """
    Find the best matching oligo based on guide sequence.

    Args:
        observed_guide: The observed guide sequence (20bp)
        middle_donor: The middle donor sequence (60bp) to identify candidate oligos
        lookup: MultiPamOligoLookup with donor -> oligo+guide mapping
        max_guide_distance: Maximum allowed Hamming distance to consider a match

    Returns:
        Tuple of (best_oligo_name, best_distance, num_alternatives)
        Returns (None, -1, 0) if no match found
    """
    if middle_donor not in lookup.donor_to_oligo_guides:
        return None, -1, 0

    alternatives = lookup.donor_to_oligo_guides[middle_donor]

    if len(alternatives) == 0:
        return None, -1, 0

    # Find best matching guide
    # Note: Use deterministic tie-breaking (alphabetical oligo name) for reproducibility
    best_oligo = None
    best_distance = 100

    for oligo_name, designed_guide in alternatives:
        dist = hamming_distance(observed_guide, designed_guide)
        # Deterministic tie-breaking: prefer lower distance, then alphabetically first oligo name
        if dist < best_distance or (dist == best_distance and (best_oligo is None or oligo_name < best_oligo)):
            best_distance = dist
            best_oligo = oligo_name

    # Only return if within threshold
    if best_distance <= max_guide_distance:
        return best_oligo, best_distance, len(alternatives)

    return None, best_distance, len(alternatives)


def reassign_oligos_by_guide(
    df: pd.DataFrame,
    lookup: MultiPamOligoLookup,
    oligo_col: str = 'bdb_oligo_name',
    guide_col: str = 'guide',
    donor_col: str = 'bdb_donor',
    perfect_guide_col: str = 'perfect_guide',
    max_guide_distance: int = 2,
    middle_donor_start: int = 30,
    middle_donor_end: int = 90,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reassign oligos for BC1s with imperfect guides using guide-based disambiguation.

    This function:
    1. Identifies BC1s with assigned oligo but imperfect guide
    2. For each, checks if an alternative oligo (same donor, different PAM) has better guide match
    3. Reassigns if the new assignment reduces guide distance to ≤ max_guide_distance

    Args:
        df: DataFrame with BC1 reference data
        lookup: MultiPamOligoLookup for guide disambiguation
        oligo_col: Column name for current oligo assignment
        guide_col: Column name for observed guide sequence
        donor_col: Column name for observed donor sequence
        perfect_guide_col: Column name for perfect guide flag
        max_guide_distance: Maximum guide distance for reassignment
        middle_donor_start: Start of middle donor region
        middle_donor_end: End of middle donor region

    Returns:
        Tuple of (updated_df, reassignment_log_df)
    """
    logger.info("Starting guide-based oligo reassignment...")

    # Make a copy to avoid modifying original
    df = df.copy()

    # Add columns for tracking
    df['guide_reassigned'] = False
    df['guide_reassigned_from'] = None
    df['guide_distance_before'] = -1
    df['guide_distance_after'] = -1

    # Identify candidates for reassignment
    # - Has oligo assignment
    # - Has guide sequence
    # - Guide is not perfect
    has_oligo = df[oligo_col].notna() & (df[oligo_col] != '')
    has_guide = df[guide_col].notna() & (df[guide_col] != '')
    imperfect_guide = df[perfect_guide_col] == False

    candidates = df[has_oligo & has_guide & imperfect_guide].copy()
    logger.info(f"  Found {len(candidates):,} candidates for reassignment")

    reassignment_log = []
    reassigned_count = 0
    checked_count = 0

    for idx, row in candidates.iterrows():
        observed_guide = row[guide_col]
        current_oligo = row[oligo_col]
        donor = row.get(donor_col)

        # Skip if no donor
        if pd.isna(donor) or len(donor) < middle_donor_end:
            continue

        checked_count += 1
        middle_donor = donor[middle_donor_start:middle_donor_end].upper()

        # Calculate current guide distance
        current_designed_guide = lookup.oligo_to_guide.get(current_oligo)
        if current_designed_guide:
            current_distance = hamming_distance(observed_guide, current_designed_guide)
        else:
            current_distance = -1

        # Find best alternative
        best_oligo, best_distance, num_alternatives = find_best_oligo_by_guide(
            observed_guide.upper(),
            middle_donor,
            lookup,
            max_guide_distance=max_guide_distance
        )

        # Reassign if better match found
        if best_oligo and best_oligo != current_oligo and best_distance < current_distance:
            df.at[idx, oligo_col] = best_oligo
            df.at[idx, 'guide_reassigned'] = True
            df.at[idx, 'guide_reassigned_from'] = current_oligo
            df.at[idx, 'guide_distance_before'] = current_distance
            df.at[idx, 'guide_distance_after'] = best_distance

            # Update perfect_guide if now perfect
            if best_distance == 0:
                df.at[idx, perfect_guide_col] = True

            reassigned_count += 1

            reassignment_log.append({
                'bc1': row.get('bc1', idx),
                'original_oligo': current_oligo,
                'new_oligo': best_oligo,
                'guide_distance_before': current_distance,
                'guide_distance_after': best_distance,
                'num_alternatives': num_alternatives,
            })

    logger.info(f"  Checked {checked_count:,} BC1s with multi-PAM donors")
    logger.info(f"  Reassigned {reassigned_count:,} BC1s ({100*reassigned_count/max(checked_count,1):.1f}%)")

    # Create log DataFrame
    log_df = pd.DataFrame(reassignment_log) if reassignment_log else pd.DataFrame()

    return df, log_df


def run_guide_disambiguation(
    reference_table_path: Path,
    oligo_design_files: List[Path],
    output_path: Path,
    oligo_col: str = 'bdb_oligo_name',
    guide_col: str = 'guide',
    donor_col: str = 'bdb_donor',
    perfect_guide_col: str = 'perfect_guide',
    max_guide_distance: int = 2,
) -> Dict:
    """
    Main function to run guide-based disambiguation on a reference table.

    Args:
        reference_table_path: Path to input reference table TSV
        oligo_design_files: List of oligo design file paths
        output_path: Path for output reference table
        oligo_col: Column name for oligo assignment
        guide_col: Column name for observed guide
        donor_col: Column name for observed donor
        perfect_guide_col: Column name for perfect guide flag
        max_guide_distance: Max guide distance for reassignment

    Returns:
        Dict with statistics about the reassignment
    """
    logger.info(f"Loading reference table from {reference_table_path}")
    df = pd.read_csv(reference_table_path, sep='\t', low_memory=False)
    logger.info(f"  Loaded {len(df):,} rows")

    # Build lookup
    lookup = load_oligo_designs_for_disambiguation(oligo_design_files)

    # Run reassignment
    df_updated, log_df = reassign_oligos_by_guide(
        df,
        lookup,
        oligo_col=oligo_col,
        guide_col=guide_col,
        donor_col=donor_col,
        perfect_guide_col=perfect_guide_col,
        max_guide_distance=max_guide_distance,
    )

    # Save results
    logger.info(f"Saving updated reference table to {output_path}")
    df_updated.to_csv(output_path, sep='\t', index=False)

    # Save log
    log_path = output_path.parent / (output_path.stem + '_reassignment_log.tsv')
    if len(log_df) > 0:
        logger.info(f"Saving reassignment log to {log_path}")
        log_df.to_csv(log_path, sep='\t', index=False)

    # Calculate statistics
    stats = {
        'total_rows': len(df),
        'reassigned': df_updated['guide_reassigned'].sum(),
        'multi_pam_donors': lookup.num_multi_pam_donors,
        'oligos_in_multi_pam_groups': lookup.num_oligos_in_multi_pam_groups,
    }

    # Guide distance improvement
    if len(log_df) > 0:
        stats['mean_distance_before'] = log_df['guide_distance_before'].mean()
        stats['mean_distance_after'] = log_df['guide_distance_after'].mean()

    return stats


if __name__ == '__main__':
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

    parser = argparse.ArgumentParser(
        description='Reassign oligos based on guide sequence for multi-PAM disambiguation'
    )
    parser.add_argument('--input', '-i', required=True, help='Input reference table TSV')
    parser.add_argument('--output', '-o', required=True, help='Output reference table TSV')
    parser.add_argument('--oligo-files', '-f', nargs='+', required=True,
                        help='Oligo design files')
    parser.add_argument('--max-distance', '-d', type=int, default=2,
                        help='Maximum guide distance for reassignment (default: 2)')

    args = parser.parse_args()

    stats = run_guide_disambiguation(
        reference_table_path=Path(args.input),
        oligo_design_files=[Path(f) for f in args.oligo_files],
        output_path=Path(args.output),
        max_guide_distance=args.max_distance,
    )

    print("\n=== Reassignment Statistics ===")
    for key, value in stats.items():
        print(f"  {key}: {value}")
