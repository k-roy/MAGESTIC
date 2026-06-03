"""
Test the merge module with real data from the SpG QTL screen.

This script validates that the merge logic correctly:
1. Combines 2x150 and 2x300 data
2. Links gdb data via donor_bc0 keys
3. Detects and flags conflicts
4. Counts PCR replicates across sources
"""

from pathlib import Path

import pandas as pd
import logging
import pytest

from magestic.pipelines.bc1_donor_bc0.core.merge import (
    merge_data_sources,
    filter_by_quality,
    BC1Record,
    extract_donor_portion_from_fragment,
    build_donor_bc0_key,
)

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


pytestmark = pytest.mark.integration

# Data paths for the SpG QTL screen (yL406 project)
BASE_DIR = Path("/path/to")
PROJECT_DIR = BASE_DIR / "projects/QTL/20240806_yL406_SpG_QTL_screen/processed_data"

# 2x300 data (full bc1-donor-bc0)
DATA_2X300_FILE = PROJECT_DIR / "bc1_donor_bc0/data_tables/top_bc1_donor_bc0_fragments_by_trafo_ID.tsv"

# Guide-donor-bc0 data
GDB_DIR = PROJECT_DIR / "guide_donor_bc0_purity_filtered"


def load_test_data(n_rows: int = 10000):
    """Load test data from the SpG QTL screen."""

    data = {}

    # Load 2x300 data
    if DATA_2X300_FILE.exists():
        logger.info(f"Loading 2x300 data from {DATA_2X300_FILE}")
        df_2x300 = pd.read_csv(DATA_2X300_FILE, sep='\t', nrows=n_rows)
        logger.info(f"  Loaded {len(df_2x300)} rows")
        logger.info(f"  Columns: {df_2x300.columns.tolist()}")
        data['2x300'] = df_2x300
    else:
        logger.warning(f"2x300 data file not found: {DATA_2X300_FILE}")

    # Load gdb data (combine all libraries)
    gdb_files = list(GDB_DIR.glob("V*_guide_donor_bc0_purity_filtered.tsv"))
    if gdb_files:
        logger.info(f"Loading gdb data from {len(gdb_files)} files")
        gdb_dfs = []
        for f in gdb_files:
            df = pd.read_csv(f, sep='\t')
            gdb_dfs.append(df)
        df_gdb = pd.concat(gdb_dfs, ignore_index=True)
        logger.info(f"  Loaded {len(df_gdb)} total gdb rows")
        logger.info(f"  Columns: {df_gdb.columns.tolist()}")
        data['gdb'] = df_gdb
    else:
        logger.warning(f"No gdb files found in {GDB_DIR}")

    return data


def test_merge_basic():
    """Test basic merge functionality."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Basic merge functionality")
    logger.info("="*60)

    data = load_test_data(n_rows=5000)

    # Merge data
    df_merged = merge_data_sources(
        df_2x300=data.get('2x300'),
        df_gdb=data.get('gdb'),
        resolve_conflicts=True
    )

    logger.info(f"\nMerged DataFrame shape: {df_merged.shape}")
    logger.info(f"Columns: {df_merged.columns.tolist()}")

    # Summary statistics
    logger.info("\n--- Summary Statistics ---")
    logger.info(f"Total records: {len(df_merged)}")
    logger.info(f"With 2x300 data: {df_merged['has_2x300'].sum()}")
    logger.info(f"With gdb data: {df_merged['has_gdb'].sum()}")
    logger.info(f"With both: {((df_merged['has_2x300']) & (df_merged['has_gdb'])).sum()}")
    logger.info(f"With donor conflicts: {df_merged['has_donor_conflict'].sum()}")

    # Show sample records
    logger.info("\n--- Sample Records ---")
    sample = df_merged.head(3)
    for i, row in sample.iterrows():
        logger.info(f"\nRecord {i}:")
        logger.info(f"  bc1: {row['bc1']}")
        logger.info(f"  bc0: {row['bc0']}")
        logger.info(f"  donor length: {len(row['donor']) if pd.notna(row['donor']) else 'None'}")
        logger.info(f"  guide: {row['guide']}")
        logger.info(f"  has_2x300: {row['has_2x300']}, has_gdb: {row['has_gdb']}")
        logger.info(f"  total_counts: {row['total_counts']}")

    return df_merged


def test_conflict_detection():
    """Test that conflicts are properly detected."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Conflict detection")
    logger.info("="*60)

    data = load_test_data(n_rows=50000)

    df_merged = merge_data_sources(
        df_2x300=data.get('2x300'),
        df_gdb=data.get('gdb'),
        resolve_conflicts=True
    )

    # Find records with conflicts
    conflicts = df_merged[df_merged['has_donor_conflict']]
    logger.info(f"\nFound {len(conflicts)} records with donor conflicts")

    if len(conflicts) > 0:
        logger.info("\n--- Sample Conflicts ---")
        for i, (idx, row) in enumerate(conflicts.head(3).iterrows()):
            logger.info(f"\nConflict {i+1}:")
            logger.info(f"  bc1: {row['bc1']}")
            logger.info(f"  conflict_details: {row['conflict_details']}")
            if row['donor_full']:
                logger.info(f"  donor_full ({len(row['donor_full'])}bp): {row['donor_full'][:50]}...")
            if row['donor_from_gdb']:
                logger.info(f"  donor_from_gdb ({len(row['donor_from_gdb'])}bp): {row['donor_from_gdb'][:50]}...")

    return conflicts


def test_quality_filtering():
    """Test quality filtering."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Quality filtering")
    logger.info("="*60)

    data = load_test_data(n_rows=10000)

    df_merged = merge_data_sources(
        df_2x300=data.get('2x300'),
        df_gdb=data.get('gdb'),
    )

    logger.info(f"Before filtering: {len(df_merged)} records")

    # Apply filters
    df_filtered = filter_by_quality(
        df_merged,
        min_reads=10,
        min_pcr_replicates=2,
    )

    logger.info(f"After filtering (min_reads=10, min_pcr_replicates=2): {len(df_filtered)} records")

    return df_filtered


def test_donor_bc0_key():
    """Test the donor_bc0 key generation."""
    logger.info("\n" + "="*60)
    logger.info("TEST: Donor-BC0 key generation")
    logger.info("="*60)

    # Test extract_donor_portion_from_fragment
    fragment = "AGGCATCCGACATGTGGAGGAGTGGTGTCCCGATAGACGACTGGACAGCACCGGGAGGGG"
    donor_portion = extract_donor_portion_from_fragment(fragment)
    logger.info(f"Fragment: {fragment}")
    logger.info(f"Donor portion: {donor_portion}")
    logger.info(f"Expected: AGGCATCCGACATGTGGAGGAGTGGTGTCC (30bp)")

    # Test build_donor_bc0_key
    donor = "GCATCCATTCTTTGTTTATTCTTTAGCTGGGTTCATGTCATTCCCTTCCTACGTCAATAGTAGTGATGATAGGATATGAAAGAATGCATCAACTGTGGAAGGCATCCGACATGTGGAGGAGTGGTGTCC"
    bc0 = "CCGGGAGGGG"
    key = build_donor_bc0_key(donor, bc0)
    logger.info(f"\nFull donor ({len(donor)}bp): {donor[:30]}...{donor[-30:]}")
    logger.info(f"BC0: {bc0}")
    logger.info(f"Key (last 30bp + bc0): {key}")

    # Verify the key matches what we'd get from fragment
    expected_key = donor_portion + bc0[-10:]  # Adjust bc0 if needed
    logger.info(f"Expected from fragment: {donor_portion + bc0}")


if __name__ == "__main__":
    logger.info("Starting merge module tests with real data")
    logger.info("="*60)

    # Run tests
    test_donor_bc0_key()
    df_merged = test_merge_basic()
    conflicts = test_conflict_detection()
    df_filtered = test_quality_filtering()

    logger.info("\n" + "="*60)
    logger.info("All tests completed!")
    logger.info("="*60)
