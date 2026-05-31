"""
Strain utilities for the Bloom et al. 16-strain yeast panel.

This module provides utilities for working with the 16-strain panel used in
Bloom et al. 2019 for QTL mapping. The strains are arranged in a round-robin
cross design where each strain is crossed with its left and right neighbors.

Key concepts:
- Shared parent: The strain that two adjacent crosses share
- Left/right neighbors: The two strains adjacent to a shared parent
- Strain origin classification: Categorizes variants by which strains carry them

Reference:
    Bloom JS, et al. (2019). Rare variants contribute disproportionately
    to quantitative trait variation in yeast. eLife 8:e49212.

Author: Kevin R. Roy
"""

from typing import Optional, Set, Tuple

import pandas as pd


# Round-robin cross strain order (Bloom et al.)
STRAIN_ORDER = (
    'RMx', 'BYa', 'M22', 'CLIB219x', 'CBS2888a', 'YJM981x',
    '273614xa', 'PW5a', 'Y10x', 'I14a', 'YPS1009x',
    'YJM454a', 'YJM978x', 'CLIB413a', 'YJM145x', 'YPS163a'
)

# Extended for wrap-around lookups (adds first two at end)
STRAIN_ORDER_EXTENDED = STRAIN_ORDER + ('RMx', 'BYa')

# Strain colors for visualization
STRAIN_COLORS = {
    'RMx': '#e41a1c', 'BYa': '#377eb8', 'M22': '#4daf4a',
    'CLIB219x': '#984ea3', 'CBS2888a': '#ff7f00', 'YJM981x': '#a65628',
    '273614xa': '#f781bf', 'PW5a': '#999999', 'Y10x': '#66c2a5',
    'I14a': '#fc8d62', 'YPS1009x': '#8da0cb', 'YJM454a': '#e78ac3',
    'YJM978x': '#a6d854', 'CLIB413a': '#ffd92f', 'YJM145x': '#e5c494',
    'YPS163a': '#b3b3b3',
}


def get_strain_neighbors(shared_parent: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Get the left and right neighbor strains for a shared parent in round-robin.

    In the Bloom et al. design, each strain participates in two crosses:
    - One with its left neighbor (previous in STRAIN_ORDER)
    - One with its right neighbor (next in STRAIN_ORDER)

    Args:
        shared_parent: The shared parent strain name

    Returns:
        Tuple of (left_neighbor, right_neighbor), or (None, None) if not found
    """
    if shared_parent not in STRAIN_ORDER:
        return None, None
    idx = STRAIN_ORDER.index(shared_parent)
    left = STRAIN_ORDER_EXTENDED[idx - 1] if idx > 0 else STRAIN_ORDER[-1]
    right = STRAIN_ORDER_EXTENDED[idx + 1]
    return left, right


def classify_strain_origin(samples: Set[str], shared_parent: str) -> str:
    """
    Classify variant based on which strains carry it.

    This classification is useful for understanding variant segregation
    patterns in the round-robin cross design.

    Args:
        samples: Set of strain names that carry the variant
        shared_parent: The shared parent strain for the cross of interest

    Returns:
        Classification string:
        - 'all_three': Present in shared parent and both neighbors
        - 'shared_only': Only in shared parent
        - 'shared_plus_left': In shared parent and left neighbor
        - 'shared_plus_right': In shared parent and right neighbor
        - 'both_neighbors': In both neighbors but not shared parent
        - 'left_only': Only in left neighbor
        - 'right_only': Only in right neighbor
        - 'other_strains': Other combinations
    """
    left, right = get_strain_neighbors(shared_parent)
    if left is None:
        return 'other_strains'

    has_shared = shared_parent in samples
    has_left = left in samples
    has_right = right in samples

    if has_shared and has_left and has_right:
        return 'all_three'
    elif has_shared and has_left:
        return 'shared_plus_left'
    elif has_shared and has_right:
        return 'shared_plus_right'
    elif has_shared:
        return 'shared_only'
    elif has_left and has_right:
        return 'both_neighbors'
    elif has_left:
        return 'left_only'
    elif has_right:
        return 'right_only'
    else:
        return 'other_strains'


def is_variant_filtered_by_shared_logic(samples: Set[str], shared_parent: str) -> bool:
    """
    Determine if a variant would be FILTERED (not designed) by shared parent logic.

    In the MAGESTIC design, variants are selected based on segregation patterns
    in the round-robin crosses.

    Variants are KEPT (not filtered) if:
    1. Variant is in shared parent AND NOT in both left and right neighbors
    2. Variant is in BOTH left and right neighbors AND NOT in shared parent

    Variants are FILTERED if:
    - In all three strains (shared + both neighbors)
    - Other combinations that don't meet the segregation criteria

    Args:
        samples: Set of strain names that carry the variant
        shared_parent: The shared parent strain

    Returns:
        True if variant would be filtered, False if kept
    """
    left, right = get_strain_neighbors(shared_parent)
    if left is None:
        return True  # Can't determine, mark as filtered

    has_shared = shared_parent in samples
    has_left = left in samples
    has_right = right in samples

    # Kept: in shared parent but NOT in both neighbors
    if has_shared and not (has_left and has_right):
        return False  # Not filtered (kept)

    # Kept: in both neighbors but NOT in shared parent
    if (has_left and has_right) and not has_shared:
        return False  # Not filtered (kept)

    # Everything else is filtered
    return True


def parse_samples(samples_str: str) -> Set[str]:
    """
    Parse comma-separated sample string into set of strain names.

    Args:
        samples_str: Comma-separated string of strain names (e.g., "RMx,BYa,M22")

    Returns:
        Set of strain names, empty set if input is NA or empty
    """
    if pd.isna(samples_str) or not samples_str:
        return set()
    return set(s.strip() for s in str(samples_str).split(','))


def get_cross_id(shared_parent: str) -> str:
    """
    Get the cross ID for a shared parent.

    Cross IDs follow the pattern "cross_NN" based on the position
    of the shared parent in STRAIN_ORDER.

    Args:
        shared_parent: The shared parent strain name

    Returns:
        Cross ID string (e.g., "cross_01"), or empty string if not found
    """
    if shared_parent not in STRAIN_ORDER:
        return ""
    idx = STRAIN_ORDER.index(shared_parent)
    return f"cross_{idx + 1:02d}"


def get_shared_parent_from_cross(cross_id: str) -> Optional[str]:
    """
    Get the shared parent strain from a cross ID.

    Args:
        cross_id: Cross ID string (e.g., "cross_01")

    Returns:
        Shared parent strain name, or None if invalid cross ID
    """
    if not cross_id.startswith("cross_"):
        return None
    try:
        idx = int(cross_id.split("_")[1]) - 1
        if 0 <= idx < len(STRAIN_ORDER):
            return STRAIN_ORDER[idx]
    except (ValueError, IndexError):
        pass
    return None
