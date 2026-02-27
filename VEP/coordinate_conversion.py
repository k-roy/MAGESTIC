"""
Coordinate Conversion Utilities for VEP Pipeline

Consolidated functions for converting between S288C reference coordinates
and strain-specific coordinates using offset tables.

Offset tables are JSON files that map S288C positions to cumulative offsets
caused by insertions and deletions in strain genomes relative to S288C.

Functions:
- load_offset_table: Load offset table from JSON file
- save_offset_table: Save offset table to JSON file
- get_offset_at_position: Get cumulative offset at S288C position
- s288c_to_strain_pos: Convert S288C position to strain position
- strain_to_s288c_pos: Convert strain position to S288C position (approximate)
- build_offset_table_from_variants: Build offset table from VCF variants

Author: Kevin R. Roy
Date: 2026-02-26
"""

import bisect
import json
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Tuple, NamedTuple

# =============================================================================
# Constants
# =============================================================================

BASE_DIR = Path("/path/to")

# Default offset tables directory
DEFAULT_OFFSET_DIR = (
    BASE_DIR / "projects/QTL/20260210_variant_effect_prediction"
    / "processed_data/20260205_in_silico_QTN_dissection/strain_genomes/offset_tables"
)


# =============================================================================
# Type Definitions
# =============================================================================

class OffsetTable(NamedTuple):
    """Offset table with sorted keys for efficient lookup."""
    positions: List[int]   # Sorted S288C positions where offset changes
    offsets: List[int]     # Cumulative offset at each position


# =============================================================================
# Core Functions
# =============================================================================

def load_offset_table(
    strain: str,
    chrom: str,
    offset_dir: Optional[Path] = None,
    use_cache: bool = True
) -> Dict[int, int]:
    """Load offset table for a strain/chromosome pair.

    Offset tables map S288C positions to cumulative offsets. These offsets
    represent the net insertion/deletion between S288C and the strain at
    each variant position.

    Args:
        strain: Strain name (e.g., 'RMx', 'BYa')
        chrom: Chromosome name (e.g., 'chrVII')
        offset_dir: Directory containing offset JSON files.
                   Default: DEFAULT_OFFSET_DIR
        use_cache: If True, cache loaded tables (recommended)

    Returns:
        Dictionary mapping S288C position -> cumulative offset

    Example:
        >>> offsets = load_offset_table('RMx', 'chrVII')
        >>> offsets[100000]  # Offset at position 100000
        -5  # RMx has 5 fewer bases than S288C up to this point
    """
    if offset_dir is None:
        offset_dir = DEFAULT_OFFSET_DIR

    offset_dir = Path(offset_dir)
    offset_file = offset_dir / f"{strain}_{chrom}_offsets.json"

    if not offset_file.exists():
        return {}

    if use_cache:
        return _load_offset_table_cached(offset_file)
    else:
        return _load_offset_table_impl(offset_file)


@lru_cache(maxsize=256)
def _load_offset_table_cached(offset_file: Path) -> Dict[int, int]:
    """Cached implementation of load_offset_table."""
    return _load_offset_table_impl(offset_file)


def _load_offset_table_impl(offset_file: Path) -> Dict[int, int]:
    """Load offset table from JSON file."""
    with open(offset_file, 'r', encoding='utf-8') as f:
        raw_offsets = json.load(f)
    # Convert string keys to int (JSON serializes dict keys as strings)
    return {int(k): v for k, v in raw_offsets.items()}


def save_offset_table(
    offset_table: Dict[int, int],
    strain: str,
    chrom: str,
    offset_dir: Optional[Path] = None
) -> Path:
    """Save offset table to JSON file.

    Args:
        offset_table: Dictionary mapping S288C position -> cumulative offset
        strain: Strain name
        chrom: Chromosome name
        offset_dir: Output directory. Default: DEFAULT_OFFSET_DIR

    Returns:
        Path to saved file.
    """
    if offset_dir is None:
        offset_dir = DEFAULT_OFFSET_DIR

    offset_dir = Path(offset_dir)
    offset_dir.mkdir(parents=True, exist_ok=True)

    output_file = offset_dir / f"{strain}_{chrom}_offsets.json"

    # Sort by position for consistency
    sorted_offsets = {str(k): v for k, v in sorted(offset_table.items())}

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(sorted_offsets, f, indent=2)

    return output_file


def get_offset_at_position(
    offset_table: Dict[int, int],
    s288c_pos: int
) -> int:
    """Get cumulative offset at an S288C position.

    The offset represents the net insertion/deletion between S288C and
    the strain at positions <= s288c_pos.

    For efficient repeated lookups, consider using make_offset_lookup().

    Args:
        offset_table: Dictionary mapping S288C position -> offset
        s288c_pos: Position in S288C coordinates (1-based)

    Returns:
        Cumulative offset at this position.
        Positive = strain has extra bases (insertions)
        Negative = strain has fewer bases (deletions)
    """
    if not offset_table:
        return 0

    # Find the largest key <= s288c_pos
    offset = 0
    for pos in sorted(offset_table.keys()):
        if pos <= s288c_pos:
            offset = offset_table[pos]
        else:
            break

    return offset


def make_offset_lookup(offset_table: Dict[int, int]) -> OffsetTable:
    """Create efficient lookup structure from offset table.

    Use this when you need to make many lookups on the same table.
    Uses bisect for O(log n) lookups instead of O(n).

    Args:
        offset_table: Dictionary mapping S288C position -> offset

    Returns:
        OffsetTable namedtuple for efficient lookups
    """
    if not offset_table:
        return OffsetTable([], [])

    sorted_items = sorted(offset_table.items())
    positions = [item[0] for item in sorted_items]
    offsets = [item[1] for item in sorted_items]

    return OffsetTable(positions, offsets)


def get_offset_fast(lookup: OffsetTable, s288c_pos: int) -> int:
    """Fast offset lookup using binary search.

    Args:
        lookup: OffsetTable from make_offset_lookup()
        s288c_pos: Position in S288C coordinates

    Returns:
        Cumulative offset at this position
    """
    if not lookup.positions:
        return 0

    idx = bisect.bisect_right(lookup.positions, s288c_pos)
    if idx == 0:
        return 0
    return lookup.offsets[idx - 1]


def s288c_to_strain_pos(
    s288c_pos: int,
    offset_table: Dict[int, int],
    locus_start: int = 1
) -> int:
    """Convert S288C position to strain sequence position.

    Args:
        s288c_pos: Position in S288C coordinates (1-based)
        offset_table: Offset table for the strain/chromosome
        locus_start: S288C position of locus start (default 1 for chromosome)

    Returns:
        Position in strain coordinates (0-based relative to locus start)
    """
    offset = get_offset_at_position(offset_table, s288c_pos)
    return (s288c_pos - locus_start) + offset


def strain_to_s288c_pos(
    strain_pos: int,
    offset_table: Dict[int, int],
    locus_start: int = 1
) -> int:
    """Convert strain position back to S288C coordinates (approximate).

    Note: This is an approximation because multiple S288C positions may
    map to the same strain position (within deletions).

    Args:
        strain_pos: Position in strain sequence (0-based relative to locus)
        offset_table: Offset table for the strain/chromosome
        locus_start: S288C position of locus start

    Returns:
        Approximate S288C position (1-based)
    """
    # Simple approximation: walk through offset table
    if not offset_table:
        return locus_start + strain_pos

    # Target: find s288c_pos such that s288c_to_strain_pos(s288c_pos) == strain_pos
    # Use binary search
    sorted_positions = sorted(offset_table.keys())

    low = locus_start
    high = max(sorted_positions) + 1000  # Extend beyond last position

    while low < high:
        mid = (low + high) // 2
        mapped = s288c_to_strain_pos(mid, offset_table, locus_start)
        if mapped < strain_pos:
            low = mid + 1
        else:
            high = mid

    return low


def build_offset_table_from_variants(
    variants: List[Tuple[int, str, str]],
    chrom_length: int
) -> Dict[int, int]:
    """Build offset table from a list of variants.

    Variants must be sorted by position. Each variant is a tuple of
    (position, ref_allele, alt_allele).

    Args:
        variants: List of (pos, ref, alt) tuples, sorted by position
        chrom_length: Length of the chromosome

    Returns:
        Offset table mapping S288C position -> cumulative offset
    """
    offset_table = {}
    cumulative_offset = 0

    for pos, ref, alt in sorted(variants, key=lambda x: x[0]):
        # Offset change from this variant
        delta = len(alt) - len(ref)
        cumulative_offset += delta
        offset_table[pos] = cumulative_offset

    return offset_table


def build_reverse_mapping(
    offset_table: Dict[int, int],
    chrom_length: int
) -> Dict[int, int]:
    """Build reverse mapping from strain position to S288C position.

    For each strain position, find the corresponding S288C position.

    Args:
        offset_table: Forward offset table (S288C -> offset)
        chrom_length: Length of the chromosome

    Returns:
        Dictionary mapping strain position -> S288C position
    """
    reverse_map = {}

    # For each S288C position, calculate strain position
    for s288c_pos in range(1, chrom_length + 1):
        strain_pos = s288c_to_strain_pos(s288c_pos, offset_table, locus_start=1)
        if strain_pos >= 0:
            reverse_map[strain_pos] = s288c_pos

    return reverse_map


# =============================================================================
# Test
# =============================================================================

if __name__ == "__main__":
    # Test basic functionality with mock data
    print("Testing coordinate_conversion module...")

    # Create test offset table
    test_offsets = {
        100: -2,   # 2 bp deletion at position 100
        200: -2,   # No change (still -2)
        300: 1,    # 3 bp insertion at position 300 (net: -2 + 3 = +1)
        400: 1,    # No change
    }

    # Test get_offset_at_position
    print("\nTesting get_offset_at_position...")
    assert get_offset_at_position(test_offsets, 50) == 0, "Before any variant"
    assert get_offset_at_position(test_offsets, 100) == -2, "At first variant"
    assert get_offset_at_position(test_offsets, 150) == -2, "Between variants"
    assert get_offset_at_position(test_offsets, 300) == 1, "At insertion"
    assert get_offset_at_position(test_offsets, 500) == 1, "After all variants"
    print("  get_offset_at_position: PASS")

    # Test s288c_to_strain_pos
    print("\nTesting s288c_to_strain_pos...")
    assert s288c_to_strain_pos(50, test_offsets, locus_start=1) == 49
    assert s288c_to_strain_pos(150, test_offsets, locus_start=1) == 147  # 150 - 1 + (-2)
    print("  s288c_to_strain_pos: PASS")

    # Test make_offset_lookup and get_offset_fast
    print("\nTesting fast offset lookup...")
    lookup = make_offset_lookup(test_offsets)
    assert get_offset_fast(lookup, 50) == 0
    assert get_offset_fast(lookup, 100) == -2
    assert get_offset_fast(lookup, 150) == -2
    assert get_offset_fast(lookup, 300) == 1
    print("  fast offset lookup: PASS")

    # Test build_offset_table_from_variants
    print("\nTesting build_offset_table_from_variants...")
    test_variants = [
        (100, 'ATG', 'A'),    # 2 bp deletion (-2)
        (300, 'G', 'GCTA'),   # 3 bp insertion (+3)
    ]
    built_offsets = build_offset_table_from_variants(test_variants, 1000)
    assert built_offsets[100] == -2
    assert built_offsets[300] == 1  # cumulative: -2 + 3 = +1
    print("  build_offset_table_from_variants: PASS")

    # Test with real data if available
    print("\nTesting with real offset table...")
    test_file = DEFAULT_OFFSET_DIR / "RMx_chrVII_offsets.json"
    if test_file.exists():
        offsets = load_offset_table('RMx', 'chrVII')
        print(f"  Loaded {len(offsets)} offset entries")

        # Find a position with non-zero offset
        for pos in sorted(offsets.keys())[:5]:
            offset = offsets[pos]
            strain_pos = s288c_to_strain_pos(pos, offsets, locus_start=1)
            print(f"  S288C pos {pos}: offset={offset}, strain_pos={strain_pos}")

        print("  load_offset_table: PASS")
    else:
        print(f"  Skipping - offset file not found: {test_file}")

    print("\nAll coordinate_conversion tests completed!")
