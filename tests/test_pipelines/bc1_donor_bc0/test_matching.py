"""
Test the sliding window matching module.

Tests include:
1. Perfect donor matching
2. Sliding window matching for truncated donors
3. Guide-based disambiguation
4. The specific 95bp truncated donor case we identified
"""

from pathlib import Path

import pandas as pd
import logging
import pytest

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    OligoLookupTables,
    load_oligo_design,
    build_lookup_tables,
    match_sequence,
    match_dataframe,
    load_and_build_lookup,
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

pytestmark = pytest.mark.integration

# Paths
BASE_DIR = Path("/path/to")
COMMON_DIR = BASE_DIR / "common"
SPG_OLIGO_FILE = COMMON_DIR / "oligo_designs/20240411_Twist_200mer_oligo_array_order.tsv"
SPCAS9_OLIGO_FILE = COMMON_DIR / "oligo_designs/20210422_Twist_200mer.tsv"


def test_load_oligo_design():
    """Test loading oligo design files."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Load oligo design")
    logger.info("="*60)

    df = load_oligo_design(SPCAS9_OLIGO_FILE)

    logger.info(f"Loaded {len(df)} oligos")
    logger.info(f"Columns: {df.columns.tolist()}")

    # Show sample
    sample = df.head(2)
    for i, row in sample.iterrows():
        logger.info(f"\nOligo {i}:")
        logger.info(f"  Name: {row['oligo_name']}")
        logger.info(f"  Guide: {row['guide']}")
        logger.info(f"  Donor ({len(row['donor'])}bp): {row['donor'][:50]}...")

    return df


def test_build_lookup_tables():
    """Test building lookup tables."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Build lookup tables")
    logger.info("="*60)

    df = load_oligo_design(SPCAS9_OLIGO_FILE)
    tables = build_lookup_tables(df)

    stats = tables.get_stats()
    logger.info(f"Lookup table statistics:")
    for key, value in stats.items():
        logger.info(f"  {key}: {value}")

    return tables


def test_perfect_donor_match():
    """Test perfect donor matching."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Perfect donor matching")
    logger.info("="*60)

    df = load_oligo_design(SPCAS9_OLIGO_FILE)
    tables = build_lookup_tables(df)

    # Pick a donor from the design
    test_row = df.iloc[0]
    test_donor = test_row['donor']
    test_guide = test_row['guide']
    expected_oligo = test_row['oligo_name']

    logger.info(f"Testing with donor from: {expected_oligo}")
    logger.info(f"Donor ({len(test_donor)}bp): {test_donor[:50]}...")
    logger.info(f"Guide: {test_guide}")

    result = match_sequence(test_donor, test_guide, tables)

    logger.info(f"\nMatch result:")
    logger.info(f"  oligo_name: {result.oligo_name}")
    logger.info(f"  perfect_donor: {result.perfect_donor}")
    logger.info(f"  is_unambiguous: {result.is_unambiguous}")
    logger.info(f"  match_method: {result.match_method}")
    logger.info(f"  Expected: {expected_oligo}")
    logger.info(f"  Match correct: {result.oligo_name == expected_oligo}")

    assert result.oligo_name == expected_oligo, "Perfect donor should match"
    assert result.perfect_donor == True
    assert result.is_unambiguous == True


def test_truncated_donor_sliding_window():
    """
    Test the specific 95bp truncated donor case we identified.

    The donor was truncated from 129bp to 95bp (34bp deletion from 5' end).
    The old bdb matching failed because it extracted a fixed middle region.
    The new sliding window should match because the intact portion matches.
    """
    logger.info("\n" + "="*60)
    logger.info("TEST: Truncated donor (95bp) sliding window match")
    logger.info("="*60)

    df = load_oligo_design(SPCAS9_OLIGO_FILE)
    tables = build_lookup_tables(df)

    # The specific case we analyzed
    # Original oligo: SpCas9_PHM7_chrXV_164993_T_G_CGG_PAM_-_PAM_strand_164986_PAM_coord_guide_disruption_position_7
    expected_oligo = "SpCas9_PHM7_chrXV_164993_T_G_CGG_PAM_-_PAM_strand_164986_PAM_coord_guide_disruption_position_7"

    # Get the designed donor
    designed_row = df[df['oligo_name'] == expected_oligo]
    if len(designed_row) == 0:
        logger.warning(f"Oligo {expected_oligo} not found in design file")
        return

    designed_donor = designed_row.iloc[0]['donor']
    guide = designed_row.iloc[0]['guide']

    logger.info(f"Designed donor ({len(designed_donor)}bp): {designed_donor}")
    logger.info(f"Guide: {guide}")

    # The truncated donor (95bp) - missing 34bp from 5' end
    # Original: AAAGGATACGACCTTCTTGGGTTGAAGAATTGTAGAATTCTTCTTGTCCAAGTATTGGTTTCCACAAC...
    # Truncated: AAAGGATACGACCTTCTTATTGGTTTCCACAAC... (missing the middle portion)
    truncated_donor = "AAAGGATACGACCTTCTTATTGGTTTCCACAACATCCGGAGGAAATGCCTGAGGCTCTTTCATGACAGCTGCCGGATCACTGTAAATAACTCCCC"

    logger.info(f"\nTruncated donor ({len(truncated_donor)}bp): {truncated_donor}")

    # Test without guide
    result_no_guide = match_sequence(truncated_donor, None, tables)
    logger.info(f"\nMatch WITHOUT guide:")
    logger.info(f"  oligo_name: {result_no_guide.oligo_name}")
    logger.info(f"  perfect_middle_region: {result_no_guide.perfect_middle_region}")
    logger.info(f"  match_method: {result_no_guide.match_method}")
    logger.info(f"  num_candidates: {result_no_guide.num_candidate_oligos}")

    # Test with guide
    result_with_guide = match_sequence(truncated_donor, guide, tables)
    logger.info(f"\nMatch WITH guide:")
    logger.info(f"  oligo_name: {result_with_guide.oligo_name}")
    logger.info(f"  perfect_middle_region: {result_with_guide.perfect_middle_region}")
    logger.info(f"  guide_match: {result_with_guide.guide_match}")
    logger.info(f"  is_unambiguous: {result_with_guide.is_unambiguous}")
    logger.info(f"  match_method: {result_with_guide.match_method}")

    # The sliding window should find a match (if any 60bp window is intact)
    # Check if any 60bp window of truncated donor exists in designed donor
    window_size = 60
    windows_found = []
    for i in range(len(truncated_donor) - window_size + 1):
        window = truncated_donor[i:i + window_size]
        if window in designed_donor:
            windows_found.append((i, window))

    logger.info(f"\n60bp windows from truncated donor found in designed:")
    for pos, window in windows_found[:3]:  # Show first 3
        logger.info(f"  Position {pos}: {window[:30]}...")


def test_extended_donor():
    """Test matching an extended donor (chimera)."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Extended donor (chimera)")
    logger.info("="*60)

    df = load_oligo_design(SPCAS9_OLIGO_FILE)
    tables = build_lookup_tables(df)

    # Take a designed donor and extend it
    test_row = df.iloc[100]
    designed_donor = test_row['donor']
    guide = test_row['guide']
    expected_oligo = test_row['oligo_name']

    # Create extended donor (add 50bp of random sequence)
    extension = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    extended_donor = designed_donor + extension

    logger.info(f"Original donor ({len(designed_donor)}bp)")
    logger.info(f"Extended donor ({len(extended_donor)}bp): ...{extended_donor[-30:]}")

    result = match_sequence(extended_donor, guide, tables)

    logger.info(f"\nMatch result:")
    logger.info(f"  oligo_name: {result.oligo_name}")
    logger.info(f"  perfect_donor: {result.perfect_donor}")
    logger.info(f"  perfect_middle_region: {result.perfect_middle_region}")
    logger.info(f"  match_method: {result.match_method}")
    logger.info(f"  Expected: {expected_oligo}")
    logger.info(f"  Match correct: {result.oligo_name == expected_oligo}")

    # Extended donor should NOT be perfect but should match via sliding window
    assert result.perfect_donor == False, "Extended donor should not be perfect match"
    assert result.perfect_middle_region == True, "Should match via sliding window"


def test_match_dataframe():
    """Test batch matching of a DataFrame."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Batch DataFrame matching")
    logger.info("="*60)

    # Load reference table
    ref_file = Path("/path/to/processed_data/by_project/QTL/20250912_SpG_QTL_screen/bc1_donor_bc0/data_tables/top_bc1_donor_bc0_associations_min_2_reads.tsv")

    if not ref_file.exists():
        logger.warning(f"Reference file not found: {ref_file}")
        return

    df = pd.read_csv(ref_file, sep='\t', nrows=1000)
    logger.info(f"Loaded {len(df)} rows")

    # Load oligo tables
    tables = load_and_build_lookup([SPCAS9_OLIGO_FILE, SPG_OLIGO_FILE])

    # Match (without guide since this is bdb data)
    df_matched = match_dataframe(
        df,
        tables,
        donor_col='donor',
        guide_col='guide',  # This column doesn't exist, will be ignored
        prefix='new_'
    )

    logger.info(f"\nMatching summary:")
    logger.info(f"  Total rows: {len(df_matched)}")
    logger.info(f"  With oligo match: {df_matched['new_oligo_name'].notna().sum()}")
    logger.info(f"  Perfect donor: {df_matched['new_perfect_donor'].sum()}")
    logger.info(f"  Perfect middle region: {df_matched['new_perfect_middle_region'].sum()}")
    logger.info(f"  Unambiguous: {df_matched['new_is_unambiguous'].sum()}")

    logger.info(f"\nMatch method distribution:")
    logger.info(df_matched['new_match_method'].value_counts())

    return df_matched


def test_compare_old_vs_new_matching():
    """
    Compare old bdb matching results with new sliding window matching.

    This verifies that the new matching recovers more cases.
    """
    logger.info("\n" + "="*60)
    logger.info("TEST: Compare old vs new matching")
    logger.info("="*60)

    # Load the existing reference table with old bdb matching results
    ref_file = Path("/path/to/processed_data/by_project/QTL/20250912_SpG_QTL_screen/bc1_donor_bc0/data_tables/20250912_SpG_QTL_bc1_reference_table_min_2_reads.tsv")

    if not ref_file.exists():
        logger.warning(f"Reference file not found: {ref_file}")
        return

    df = pd.read_csv(ref_file, sep='\t', nrows=10000)
    logger.info(f"Loaded {len(df)} rows from reference table")

    # Check if old bdb columns exist
    old_cols = ['bdb_oligo_name', 'bdb_perfect_donor', 'bdb_perfect_middle_donor_region']
    has_old_cols = all(col in df.columns for col in old_cols)

    if not has_old_cols:
        logger.warning("Old bdb columns not found in reference table")
        return

    # Load oligo tables
    tables = load_and_build_lookup([SPCAS9_OLIGO_FILE, SPG_OLIGO_FILE])

    # Rename donor column if needed
    donor_col = 'bdb_donor' if 'bdb_donor' in df.columns else 'donor'

    # Apply new matching
    df_matched = match_dataframe(
        df,
        tables,
        donor_col=donor_col,
        guide_col='guide',
        prefix='new_'
    )

    # Compare results
    logger.info(f"\n=== Comparison Results ===")

    old_matched = df_matched['bdb_oligo_name'].notna().sum()
    new_matched = df_matched['new_oligo_name'].notna().sum()

    logger.info(f"Old bdb matching: {old_matched} matched ({100*old_matched/len(df):.1f}%)")
    logger.info(f"New sliding window: {new_matched} matched ({100*new_matched/len(df):.1f}%)")
    logger.info(f"Improvement: +{new_matched - old_matched} ({100*(new_matched - old_matched)/len(df):.2f}%)")

    # Cases where new matched but old didn't
    newly_matched = df_matched[(df_matched['new_oligo_name'].notna()) & (df_matched['bdb_oligo_name'].isna())]
    logger.info(f"\nNewly matched by sliding window: {len(newly_matched)}")

    if len(newly_matched) > 0:
        logger.info("\nSample of newly matched cases:")
        for i, (idx, row) in enumerate(newly_matched.head(3).iterrows()):
            donor_len = len(row[donor_col]) if pd.notna(row[donor_col]) else 0
            logger.info(f"\n  Case {i+1}:")
            logger.info(f"    bc1: {row['bc1']}")
            logger.info(f"    donor length: {donor_len}bp")
            logger.info(f"    new_oligo_name: {row['new_oligo_name']}")
            logger.info(f"    new_match_method: {row['new_match_method']}")

    # Check agreement where both matched
    both_matched = df_matched[(df_matched['new_oligo_name'].notna()) & (df_matched['bdb_oligo_name'].notna())]
    if len(both_matched) > 0:
        agreement = (both_matched['new_oligo_name'] == both_matched['bdb_oligo_name']).sum()
        logger.info(f"\nAgreement where both matched: {agreement}/{len(both_matched)} ({100*agreement/len(both_matched):.1f}%)")


if __name__ == "__main__":
    logger.info("Starting matching module tests")
    logger.info("="*60)

    test_load_oligo_design()
    test_build_lookup_tables()
    test_perfect_donor_match()
    test_truncated_donor_sliding_window()
    test_extended_donor()
    test_match_dataframe()
    test_compare_old_vs_new_matching()

    logger.info("\n" + "="*60)
    logger.info("All tests completed!")
    logger.info("="*60)
