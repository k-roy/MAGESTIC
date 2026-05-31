"""
Coordinate conversion utilities for MAGESTIC and S288c reference genomes.

The MAGESTIC background strain has coordinate offsets from the S288c (R64-2-1)
reference genome due to strain-specific insertions and deletions.

For plotting and display purposes, S288c coordinates are preferred since
they match published reference databases and literature.

Author: Kevin R. Roy
"""

from pathlib import Path
from typing import Optional

import pandas as pd


# Default coordinate offset file path
DEFAULT_COORD_OFFSET_FILE = Path(
    "/path/to/scripts_and_keyfiles/guide_donor_design"
    "/20240404_SpG_natural_variants/MAGESTIC_background_strain_coord_offset_from_R64-2-1.tsv"
)

# Module-level cache for coordinate offset data
_COORD_OFFSET_DF: Optional[pd.DataFrame] = None


def load_coordinate_offsets(offset_file: Path = DEFAULT_COORD_OFFSET_FILE) -> pd.DataFrame:
    """
    Load coordinate offset data for MAGESTIC to S288c (R64-2-1) conversion.

    The offset file contains cumulative coordinate offsets per chromosome
    that account for strain-specific insertions/deletions.

    Args:
        offset_file: Path to TSV file with offset data.
                    Columns: chrom, coord, cumulative_offset

    Returns:
        DataFrame with offset information, or empty DataFrame if file not found.
    """
    global _COORD_OFFSET_DF

    if _COORD_OFFSET_DF is not None:
        return _COORD_OFFSET_DF

    if not offset_file.exists():
        print(f"Warning: Coordinate offset file not found: {offset_file}")
        return pd.DataFrame()

    _COORD_OFFSET_DF = pd.read_csv(offset_file, sep='\t')
    return _COORD_OFFSET_DF


def magestic_to_s288c_coord(
    chrom: str,
    magestic_coord: int,
    offset_df: Optional[pd.DataFrame] = None
) -> int:
    """
    Convert MAGESTIC coordinate to S288c (R64-2-1) coordinate.

    MAGESTIC data uses coordinates relative to our lab strain background,
    which has offsets from the S288c (R64-2-1) reference used by Bloom et al.
    and most published databases.

    **PLOTTING STANDARD**: All displayed coordinates should use S288c.
    Use this function to convert MAGESTIC coords before displaying.

    Args:
        chrom: Chromosome name (e.g., 'chrX')
        magestic_coord: Coordinate in MAGESTIC background strain
        offset_df: DataFrame with offset information (loaded if None)

    Returns:
        Coordinate in S288c (R64-2-1) reference for display
    """
    if offset_df is None:
        offset_df = load_coordinate_offsets()

    if offset_df.empty:
        return magestic_coord

    chrom_offsets = offset_df[offset_df['chrom'] == chrom]
    if chrom_offsets.empty:
        return magestic_coord

    # Find the applicable offset (last offset position before our coord)
    applicable = chrom_offsets[chrom_offsets['coord'] <= magestic_coord]
    if applicable.empty:
        return magestic_coord

    offset = applicable.iloc[-1]['cumulative_offset']
    return magestic_coord - int(offset)


# Alias for backward compatibility
magestic_to_bloom_coord = magestic_to_s288c_coord


def s288c_to_magestic_coord(
    chrom: str,
    s288c_coord: int,
    offset_df: Optional[pd.DataFrame] = None
) -> int:
    """
    Convert S288c (R64-2-1) coordinate to MAGESTIC coordinate.

    Inverse of magestic_to_s288c_coord. Use for internal data lookups
    when querying MAGESTIC-aligned data using S288c coordinates.

    Args:
        chrom: Chromosome name
        s288c_coord: Coordinate in S288c (R64-2-1) reference
        offset_df: DataFrame with offset information

    Returns:
        Coordinate in MAGESTIC background strain
    """
    if offset_df is None:
        offset_df = load_coordinate_offsets()

    if offset_df.empty:
        return s288c_coord

    chrom_offsets = offset_df[offset_df['chrom'] == chrom]
    if chrom_offsets.empty:
        return s288c_coord

    # This is approximate - find offset at closest position
    applicable = chrom_offsets[chrom_offsets['coord'] <= s288c_coord]
    if applicable.empty:
        return s288c_coord

    offset = applicable.iloc[-1]['cumulative_offset']
    return s288c_coord + int(offset)


# Alias for backward compatibility
bloom_to_magestic_coord = s288c_to_magestic_coord


def clear_offset_cache() -> None:
    """Clear the cached offset DataFrame. Useful for testing."""
    global _COORD_OFFSET_DF
    _COORD_OFFSET_DF = None
