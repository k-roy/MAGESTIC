"""
Anchor-based sequence extraction utilities for MAGESTIC amplicon parsing.

This module provides robust extraction of sequence elements using constant
sequence anchors (msd, RE_site) rather than fixed position offsets.

Key concept: The 3' end of amplicons is more consistent than the 5' end:
- msd_fragment (AGGAAAACAGAC) is constant
- bc0 (10bp) + SPS (20bp) are fixed-length before msd
- donor length can vary (extended/retracted donors)

The end-slice approach extracts a fixed-length fragment from the 3' end,
providing a consistent key for linking guide_donor_bc0 and bc1_donor_bc0 data.

Amplicon structures:
    guide_donor_bc0:  [leader][guide:20bp][RE_site][donor][SPS:20bp][bc0:10bp][msd:12bp]
    bc1_donor_bc0:    [bc1:20bp][linker][RE_site][donor][SPS:20bp][bc0:10bp][msd:12bp]

Author: Kevin R. Roy
"""

from dataclasses import dataclass
from typing import Optional, Tuple, Dict
from rapidfuzz import fuzz

# =============================================================================
# Constants
# =============================================================================

@dataclass(frozen=True)
class AnchorSequences:
    """Constant sequence anchors for different Cas types."""
    RE_site: str
    msd: str
    RE_site_length: int = None
    msd_length: int = None

    def __post_init__(self):
        # Set lengths as derived values
        object.__setattr__(self, 'RE_site_length', len(self.RE_site))
        object.__setattr__(self, 'msd_length', len(self.msd))


# Anchor sequences by Cas type
ANCHORS: Dict[str, AnchorSequences] = {
    "SpCas9": AnchorSequences(
        RE_site="GTTTGAAGAGC",
        msd="AGGAAAACAGAC",
    ),
    "SpG": AnchorSequences(
        RE_site="GTTTGAAGAGC",  # Same as SpCas9
        msd="AGGAAAACAGAC",
    ),
    "LbCas12a": AnchorSequences(
        RE_site="TTTCGAAGAGC",  # One less G at start
        msd="AGGAAAACAGAC",
    ),
}

# Default lengths
DEFAULT_GUIDE_LENGTH = 20
DEFAULT_SPS_LENGTH = 20
DEFAULT_BC0_LENGTH = 10
DEFAULT_BC1_LENGTH = 20
DEFAULT_DONOR_BC0_FRAGMENT_LENGTH = 60  # For end-slice matching

# BC1 length tolerance - allow 19-21bp for synthesis/sequencing errors
# This is used when extracting BC1 from anchor-based parsing
BC1_LENGTH_RANGE = (19, 21)  # min, max inclusive


# ============================================================================
# BC1 Length Handling
# ============================================================================

def validate_bc1_length(bc1: str, expected: int = DEFAULT_BC1_LENGTH,
                        tolerance: int = 1) -> tuple:
    """
    Validate BC1 length and return (is_valid, length, message).

    Args:
        bc1: BC1 sequence
        expected: Expected length (default 20)
        tolerance: Allowed deviation from expected (default 1 = 19-21bp)

    Returns:
        Tuple of (is_valid, actual_length, message)
    """
    actual = len(bc1)
    min_len = expected - tolerance
    max_len = expected + tolerance

    if actual < min_len:
        return False, actual, f"BC1 too short: {actual}bp (min {min_len}bp)"
    elif actual > max_len:
        return False, actual, f"BC1 too long: {actual}bp (max {max_len}bp)"
    elif actual != expected:
        return True, actual, f"BC1 length {actual}bp (expected {expected}bp, within tolerance)"
    else:
        return True, actual, "OK"


def normalize_bc1_length(
    bc1: str,
    target_length: int = DEFAULT_BC1_LENGTH,
    truncate_side: str = "3prime",
    pad_char: str = "N",
) -> str:
    """
    Normalize BC1 to target length by truncating or padding.

    Use this when downstream analysis requires fixed-length BC1s.
    Prefer NOT normalizing when possible (anchor-based extraction + variable-length handling).

    Args:
        bc1: BC1 sequence
        target_length: Target length (default 20)
        truncate_side: Which side to truncate if too long ("3prime" or "5prime")
        pad_char: Character to use for padding (default "N")

    Returns:
        Normalized BC1 of target_length
    """
    actual = len(bc1)

    if actual == target_length:
        return bc1
    elif actual > target_length:
        # Truncate
        excess = actual - target_length
        if truncate_side == "3prime":
            return bc1[:target_length]
        else:  # 5prime
            return bc1[excess:]
    else:
        # Pad
        deficit = target_length - actual
        if truncate_side == "3prime":
            return bc1 + (pad_char * deficit)
        else:  # 5prime
            return (pad_char * deficit) + bc1


def extract_bc1_between_anchors(
    seq: str,
    prefix: str,
    suffix: str,
    expected_length: int = DEFAULT_BC1_LENGTH,
    length_tolerance: int = 1,
    use_fuzzy: bool = False,
    min_fuzzy_score: float = 85.0,
) -> Optional[Tuple[str, int, int]]:
    """
    Extract BC1 by finding prefix and suffix anchors.

    This is the recommended approach as it handles length variations naturally.

    Args:
        seq: Sequence to search
        prefix: BC1_PREFIX sequence (immediately before BC1)
        suffix: BC1_SUFFIX sequence (immediately after BC1)
        expected_length: Expected BC1 length
        length_tolerance: Allowed deviation from expected length
        use_fuzzy: Use fuzzy matching for anchors
        min_fuzzy_score: Minimum score for fuzzy matching

    Returns:
        Tuple of (bc1, start_pos, end_pos) or None if not found
    """
    seq = seq.upper()
    prefix = prefix.upper()
    suffix = suffix.upper()

    # Find prefix
    if use_fuzzy:
        prefix_pos = find_anchor_fuzzy(seq, prefix, min_score=min_fuzzy_score)
    else:
        prefix_pos = find_anchor_exact(seq, prefix)

    if prefix_pos is None:
        return None

    bc1_start = prefix_pos + len(prefix)

    # Find suffix (must be after BC1)
    min_suffix_pos = bc1_start + expected_length - length_tolerance
    max_suffix_pos = bc1_start + expected_length + length_tolerance + len(suffix)

    search_region = seq[min_suffix_pos:max_suffix_pos + len(suffix)]

    if use_fuzzy:
        suffix_offset = find_anchor_fuzzy(search_region, suffix, min_score=min_fuzzy_score)
    else:
        suffix_offset = search_region.find(suffix)
        suffix_offset = suffix_offset if suffix_offset >= 0 else None

    if suffix_offset is None:
        # Suffix not found - try extracting fixed length
        bc1_end = bc1_start + expected_length
        if bc1_end <= len(seq):
            bc1 = seq[bc1_start:bc1_end]
            return bc1, bc1_start, bc1_end
        return None

    bc1_end = min_suffix_pos + suffix_offset
    bc1 = seq[bc1_start:bc1_end]

    # Validate length
    bc1_len = len(bc1)
    min_len = expected_length - length_tolerance
    max_len = expected_length + length_tolerance

    if min_len <= bc1_len <= max_len:
        return bc1, bc1_start, bc1_end
    else:
        return None


# =============================================================================
# Anchor Finding Functions
# =============================================================================

def find_anchor_exact(
    seq: str,
    pattern: str,
    start: int = 0,
    end: Optional[int] = None,
) -> Optional[int]:
    """
    Find exact match of pattern in sequence.

    Args:
        seq: Sequence to search
        pattern: Pattern to find
        start: Start index for search window
        end: End index for search window (None = end of sequence)

    Returns:
        Start index of match, or None if not found
    """
    if end is None:
        end = len(seq)

    search_region = seq[start:end]
    idx = search_region.find(pattern)

    if idx != -1:
        return start + idx
    return None


def find_anchor_fuzzy(
    seq: str,
    pattern: str,
    start: int = 0,
    end: Optional[int] = None,
    min_score: float = 85.0,
) -> Optional[int]:
    """
    Find pattern in sequence using rapidfuzz fuzzy matching.

    Args:
        seq: Sequence to search
        pattern: Pattern to find
        start: Start index for search window
        end: End index for search window (None = end of sequence)
        min_score: Minimum fuzz ratio to accept match (0-100)

    Returns:
        Start index of best match, or None if no match >= min_score
    """
    if end is None:
        end = len(seq)

    search_region = seq[start:end]
    pattern_len = len(pattern)

    if len(search_region) < pattern_len:
        return None

    best_score = 0
    best_pos = None

    for i in range(len(search_region) - pattern_len + 1):
        candidate = search_region[i:i + pattern_len]
        score = fuzz.ratio(candidate, pattern)
        if score > best_score:
            best_score = score
            best_pos = i

    if best_score >= min_score:
        return start + best_pos
    return None


def find_anchor(
    seq: str,
    pattern: str,
    start: int = 0,
    end: Optional[int] = None,
    use_fuzzy: bool = True,
    min_score: float = 85.0,
) -> Optional[int]:
    """
    Find pattern in sequence, trying exact match first then fuzzy.

    Args:
        seq: Sequence to search
        pattern: Pattern to find
        start: Start index for search window
        end: End index for search window
        use_fuzzy: Whether to try fuzzy matching if exact fails
        min_score: Minimum fuzz ratio for fuzzy matching

    Returns:
        Start index of match, or None if not found
    """
    # Try exact match first
    pos = find_anchor_exact(seq, pattern, start, end)
    if pos is not None:
        return pos

    # Try fuzzy match if enabled
    if use_fuzzy:
        return find_anchor_fuzzy(seq, pattern, start, end, min_score)

    return None


def find_anchor_from_end(
    seq: str,
    pattern: str,
    max_distance_from_end: int = 50,
    use_fuzzy: bool = True,
    min_score: float = 85.0,
) -> Optional[int]:
    """
    Find pattern searching from end of sequence.

    Useful for finding msd_fragment which is always at the 3' end.

    Args:
        seq: Sequence to search
        pattern: Pattern to find
        max_distance_from_end: How far from end to search
        use_fuzzy: Whether to try fuzzy matching
        min_score: Minimum fuzz ratio for fuzzy matching

    Returns:
        Start index of match, or None if not found
    """
    start = max(0, len(seq) - max_distance_from_end - len(pattern))
    return find_anchor(seq, pattern, start=start, use_fuzzy=use_fuzzy, min_score=min_score)


# =============================================================================
# End-Slice Extraction
# =============================================================================

def extract_donor_bc0_fragment(
    seq: str,
    msd_pos: Optional[int] = None,
    msd_pattern: str = "AGGAAAACAGAC",
    fragment_length: int = DEFAULT_DONOR_BC0_FRAGMENT_LENGTH,
    bc0_length: int = DEFAULT_BC0_LENGTH,
    sps_length: int = DEFAULT_SPS_LENGTH,
) -> Optional[str]:
    """
    Extract a fixed-length donor_bc0 fragment from the 3' end.

    This provides a consistent key for linking guide_donor_bc0 and bc1_donor_bc0
    data, regardless of 5' trimming variations or full donor length.

    The fragment includes: [partial_donor][SPS:20bp][bc0:10bp]

    Args:
        seq: Full sequence
        msd_pos: Position of msd (if already found), or None to find it
        msd_pattern: msd sequence to search for
        fragment_length: Total length of fragment to extract
        bc0_length: Length of bc0 barcode
        sps_length: Length of SPS

    Returns:
        Fixed-length donor_bc0 fragment, or None if msd not found
    """
    # Find msd position if not provided
    if msd_pos is None:
        msd_pos = find_anchor_from_end(seq, msd_pattern)
        if msd_pos is None:
            return None

    # Extract fragment ending at msd
    fragment_start = msd_pos - fragment_length
    if fragment_start < 0:
        return None

    return seq[fragment_start:msd_pos]


def extract_sps_bc0(
    seq: str,
    msd_pos: Optional[int] = None,
    msd_pattern: str = "AGGAAAACAGAC",
    bc0_length: int = DEFAULT_BC0_LENGTH,
    sps_length: int = DEFAULT_SPS_LENGTH,
) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract SPS and bc0 from end of sequence.

    Args:
        seq: Full sequence
        msd_pos: Position of msd (if already found)
        msd_pattern: msd sequence to search for
        bc0_length: Length of bc0 barcode
        sps_length: Length of SPS

    Returns:
        Tuple of (SPS, bc0), or (None, None) if extraction fails
    """
    if msd_pos is None:
        msd_pos = find_anchor_from_end(seq, msd_pattern)
        if msd_pos is None:
            return None, None

    bc0_start = msd_pos - bc0_length
    sps_start = bc0_start - sps_length

    if sps_start < 0:
        return None, None

    sps = seq[sps_start:bc0_start]
    bc0 = seq[bc0_start:msd_pos]

    return sps, bc0


# =============================================================================
# Full Feature Extraction
# =============================================================================

@dataclass
class ExtractedFeatures:
    """Container for extracted sequence features."""
    guide: Optional[str] = None
    donor: Optional[str] = None
    sps: Optional[str] = None
    bc0: Optional[str] = None
    msd: Optional[str] = None
    donor_bc0_fragment: Optional[str] = None
    re_site_pos: Optional[int] = None
    msd_pos: Optional[int] = None
    orientation: str = "forward"


def extract_features_anchored(
    seq: str,
    cas_type: str = "SpCas9",
    guide_length: int = DEFAULT_GUIDE_LENGTH,
    sps_length: int = DEFAULT_SPS_LENGTH,
    bc0_length: int = DEFAULT_BC0_LENGTH,
    donor_bc0_fragment_length: int = DEFAULT_DONOR_BC0_FRAGMENT_LENGTH,
    use_fuzzy: bool = True,
) -> Optional[ExtractedFeatures]:
    """
    Extract all sequence features using anchor-based approach.

    Args:
        seq: Full sequence string
        cas_type: Cas type for anchor sequences ("SpCas9", "SpG", "LbCas12a")
        guide_length: Expected guide length
        sps_length: Expected SPS length
        bc0_length: Expected bc0 length
        donor_bc0_fragment_length: Length of end-slice fragment
        use_fuzzy: Whether to use fuzzy matching for anchors

    Returns:
        ExtractedFeatures object, or None if critical anchors not found
    """
    if cas_type not in ANCHORS:
        raise ValueError(f"Unknown cas_type: {cas_type}. Valid: {list(ANCHORS.keys())}")

    anchors = ANCHORS[cas_type]

    # Find msd from end (most reliable anchor)
    msd_pos = find_anchor_from_end(seq, anchors.msd, use_fuzzy=use_fuzzy)
    if msd_pos is None:
        return None

    # Find RE_site from start
    re_site_pos = find_anchor(seq, anchors.RE_site, end=min(200, len(seq)), use_fuzzy=use_fuzzy)
    if re_site_pos is None:
        return None

    # Calculate positions
    re_site_end = re_site_pos + anchors.RE_site_length
    bc0_end = msd_pos
    bc0_start = bc0_end - bc0_length
    sps_end = bc0_start
    sps_start = sps_end - sps_length
    donor_end = sps_start
    donor_start = re_site_end
    guide_start = re_site_pos - guide_length

    # Validate positions
    if guide_start < 0 or donor_end <= donor_start:
        return None

    # Extract features
    features = ExtractedFeatures(
        guide=seq[guide_start:re_site_pos],
        donor=seq[donor_start:donor_end],
        sps=seq[sps_start:sps_end],
        bc0=seq[bc0_start:bc0_end],
        msd=seq[msd_pos:msd_pos + anchors.msd_length],
        re_site_pos=re_site_pos,
        msd_pos=msd_pos,
    )

    # Extract end-slice donor_bc0_fragment
    features.donor_bc0_fragment = extract_donor_bc0_fragment(
        seq,
        msd_pos=msd_pos,
        msd_pattern=anchors.msd,
        fragment_length=donor_bc0_fragment_length,
        bc0_length=bc0_length,
        sps_length=sps_length,
    )

    return features


def extract_features_with_orientation(
    seq: str,
    cas_type: str = "SpCas9",
    **kwargs,
) -> Optional[ExtractedFeatures]:
    """
    Extract features, trying both orientations if needed.

    Args:
        seq: Full sequence string
        cas_type: Cas type for anchor sequences
        **kwargs: Additional arguments passed to extract_features_anchored

    Returns:
        ExtractedFeatures object with orientation set, or None if extraction fails
    """
    # Try forward orientation first
    features = extract_features_anchored(seq, cas_type, **kwargs)
    if features is not None:
        features.orientation = "forward"
        return features

    # Try reverse complement
    from .sequence_utils import reverse_complement
    rc_seq = reverse_complement(seq)

    features = extract_features_anchored(rc_seq, cas_type, **kwargs)
    if features is not None:
        features.orientation = "reverse"
        return features

    return None
