"""
Auto-inference of amplicon structure from oligo design file.

Author: Kevin R. Roy

This module provides automatic detection of:
- Nuclease type (SpCas9/SpG vs LbCas12a) based on oligo structure
- Read orientation (forward vs reverse complement)
- Appropriate anchor sequences for parsing

User supplies:
  - collapsed FASTQ files
  - oligo design TSV (with oligo_sequence column)

Script automatically detects:
  - Nuclease type and associated structure parameters
  - Read orientation from sequencing
  - Appropriate anchor sequences for parsing
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import pysam

from magestic.utils.sequence_utils import rev_comp

logger = logging.getLogger(__name__)


# ============================================================================
# Nuclease-specific anchor definitions
# ============================================================================

@dataclass
class NucleaseAnchors:
    """Anchor sequences and structure parameters for a nuclease type."""
    name: str
    pad_5p: str              # 5' pad/cloning adapter sequence
    re_site: str             # Internal restriction enzyme cloning site
    guide_start: int         # Guide start position in oligo
    guide_end: int           # Guide end position in oligo
    donor_start: int         # Donor start position in oligo
    donor_end: int = 180     # Donor end position (usually same for all)
    sps_start: int = 180     # SPS start position
    sps_end: int = 200       # SPS end position

    @property
    def guide_length(self) -> int:
        return self.guide_end - self.guide_start

    @property
    def donor_length(self) -> int:
        return self.donor_end - self.donor_start

    @property
    def pad_5p_length(self) -> int:
        return len(self.pad_5p)


# Define anchors for each nuclease type
NUCLEASE_ANCHORS: Dict[str, NucleaseAnchors] = {
    'SpCas9_SpG': NucleaseAnchors(
        name='SpCas9_SpG',
        pad_5p='ctgcgattggcaggcgcgcc',      # 20bp
        re_site='gtttgaagagc',               # 11bp
        guide_start=20,
        guide_end=40,
        donor_start=51,
        donor_end=180,
    ),
    'LbCas12a': NucleaseAnchors(
        name='LbCas12a',
        pad_5p='caaagattaaataattcccactaagtgtgg',  # 30bp (longer!)
        re_site='tttcgaagagc',               # 11bp (different - starts with ttt)
        guide_start=30,
        guide_end=50,
        donor_start=61,
        donor_end=180,
    ),
}

# SPS → nuclease/library mapping
# Built from known library designs
SPS_TO_LIBRARY: Dict[str, Dict] = {
    'cgatagacgactggacagca': {'nuclease': 'SpCas9_NGG', 'library': 'V416', 'anchors': 'SpCas9_SpG'},
    'cttgactcgacagtgagagc': {'nuclease': 'SpG_NGNG', 'library': 'V520', 'anchors': 'SpCas9_SpG'},
    'cgtatggagatgacgtcagg': {'nuclease': 'SpG_NGNH', 'library': 'V521', 'anchors': 'SpCas9_SpG'},
    'gaactacgtgtctcgtcagg': {'nuclease': 'LbCas12a_TTTV', 'library': 'V419', 'anchors': 'LbCas12a'},
    'gttatgctggtcctaggtcg': {'nuclease': 'impLbCas12a_nonTTTV', 'library': 'V419', 'anchors': 'LbCas12a'},
}


# ============================================================================
# Structure inference result
# ============================================================================

@dataclass
class StructureInferenceResult:
    """Result of amplicon structure inference."""
    nuclease_type: str                   # e.g., 'SpCas9_SpG' or 'LbCas12a'
    anchors: NucleaseAnchors             # Anchor sequences and positions
    orientation: str                      # 'forward' or 'reverse'
    confidence: float                     # 0.0 to 1.0
    sps_to_library: Dict[str, Dict]      # SPS → library/nuclease mapping

    # Detection statistics
    total_sampled: int = 0
    detected_forward: int = 0
    detected_reverse: int = 0
    undetected: int = 0

    def __str__(self) -> str:
        return (
            f"StructureInferenceResult(\n"
            f"  nuclease_type={self.nuclease_type},\n"
            f"  orientation={self.orientation},\n"
            f"  confidence={self.confidence:.2%},\n"
            f"  guide_positions=[{self.anchors.guide_start}:{self.anchors.guide_end}],\n"
            f"  donor_positions=[{self.anchors.donor_start}:{self.anchors.donor_end}],\n"
            f"  sampled={self.total_sampled}, "
            f"fwd={self.detected_forward}, rev={self.detected_reverse}, "
            f"undetected={self.undetected}\n"
            f")"
        )


# ============================================================================
# Utility functions
# ============================================================================


def fuzzy_match(seq: str, pattern: str, threshold: float = 0.85) -> bool:
    """
    Check if seq matches pattern with at least threshold similarity.

    Uses simple character overlap ratio (not edit distance).
    """
    if len(seq) != len(pattern):
        # Compare up to shorter length
        min_len = min(len(seq), len(pattern))
        seq = seq[:min_len]
        pattern = pattern[:min_len]

    if len(seq) == 0:
        return False

    matches = sum(1 for a, b in zip(seq.lower(), pattern.lower()) if a == b)
    ratio = matches / len(pattern)
    return ratio >= threshold


def sample_reads_from_fastq(
    fastq_path: Path,
    n_sample: int = 5000
) -> List[Tuple[str, str]]:
    """
    Sample first N reads from a collapsed FASTQ file.

    Returns:
        List of (sequence, quality) tuples
    """
    reads = []
    try:
        with pysam.FastxFile(str(fastq_path)) as fh:
            for i, entry in enumerate(fh):
                if i >= n_sample:
                    break
                reads.append((entry.sequence, entry.quality))
    except Exception as e:
        logger.warning(f"Error reading {fastq_path}: {e}")

    return reads


# ============================================================================
# Main inference functions
# ============================================================================

def build_sps_mapping_from_oligo_design(oligo_design_path: Path) -> Dict[str, Dict]:
    """
    Build SPS → library/nuclease mapping from oligo design file.

    Args:
        oligo_design_path: Path to oligo design TSV with oligo_sequence column

    Returns:
        Dict mapping SPS sequence to {nuclease, library, anchors}
    """
    import pandas as pd

    # Start with known mappings
    sps_mapping = dict(SPS_TO_LIBRARY)

    try:
        df = pd.read_csv(oligo_design_path, sep='\t')

        # Check for required columns
        if 'oligo_sequence' not in df.columns:
            logger.warning(f"oligo_sequence column not found in {oligo_design_path}")
            return sps_mapping

        # Extract SPS from each oligo (last 20bp)
        df['SPS'] = df['oligo_sequence'].str[-20:].str.lower()

        # Try to determine nuclease type from oligo structure
        for sps in df['SPS'].unique():
            if sps in sps_mapping:
                continue  # Already known

            # Get a sample oligo with this SPS
            sample_oligo = df[df['SPS'] == sps]['oligo_sequence'].iloc[0]

            # Detect nuclease type from structure
            nuclease_type = detect_nuclease_from_oligo(sample_oligo)
            if nuclease_type:
                sps_mapping[sps] = {
                    'nuclease': 'unknown',
                    'library': 'unknown',
                    'anchors': nuclease_type
                }
                logger.info(f"Detected {nuclease_type} structure for SPS {sps[:8]}...")

    except Exception as e:
        logger.warning(f"Error building SPS mapping from {oligo_design_path}: {e}")

    return sps_mapping


def detect_nuclease_from_oligo(oligo_seq: str) -> Optional[str]:
    """
    Detect nuclease type from a designed oligo sequence.

    Args:
        oligo_seq: Full 200bp oligo sequence

    Returns:
        Nuclease anchor type ('SpCas9_SpG' or 'LbCas12a') or None
    """
    oligo_lower = oligo_seq.lower()

    # Check for LbCas12a markers (check first - more distinctive)
    if NUCLEASE_ANCHORS['LbCas12a'].re_site in oligo_lower:
        # Verify with 5' pad (should be longer)
        pad_5p = NUCLEASE_ANCHORS['LbCas12a'].pad_5p
        if oligo_lower.startswith(pad_5p) or fuzzy_match(oligo_lower[:35], pad_5p):
            return 'LbCas12a'

    # Check for SpCas9/SpG markers
    if NUCLEASE_ANCHORS['SpCas9_SpG'].re_site in oligo_lower:
        pad_5p = NUCLEASE_ANCHORS['SpCas9_SpG'].pad_5p
        if oligo_lower.startswith(pad_5p) or fuzzy_match(oligo_lower[:25], pad_5p):
            return 'SpCas9_SpG'

    return None


def detect_structure_from_read(
    seq: str,
    anchors_dict: Dict[str, NucleaseAnchors],
    min_length: int = 100
) -> Tuple[Optional[str], str]:
    """
    Detect nuclease type and orientation from a single read.

    Args:
        seq: Read sequence
        anchors_dict: Dict of nuclease type → NucleaseAnchors
        min_length: Minimum read length to consider (skip short reads)

    Returns:
        Tuple of (nuclease_type or None, orientation)
    """
    # Skip very short reads (likely failed merges or contamination)
    if len(seq) < min_length:
        return None, 'too_short'

    seq_lower = seq.lower()
    seq_rc = rev_comp(seq).lower()

    # Additional anchor patterns to check
    # PERFECT_LEADER is the last 11bp of SpCas9/SpG 5' pad
    perfect_leader_spg = 'gcaggcgcgcc'

    for nuclease_type, anchors in anchors_dict.items():
        pad_5p = anchors.pad_5p.lower()
        re_site = anchors.re_site.lower()

        # Check forward orientation
        # Perfect match to 5' pad
        if seq_lower.startswith(pad_5p):
            return nuclease_type, 'forward'

        # Perfect match to last portion of 5' pad (PERFECT_LEADER for SpCas9/SpG)
        if nuclease_type == 'SpCas9_SpG' and seq_lower.startswith(perfect_leader_spg):
            return nuclease_type, 'forward'

        # Perfect match to RE site (if 5' pad was trimmed)
        if re_site in seq_lower[:80]:
            return nuclease_type, 'forward'

        # Check reverse complement orientation
        if seq_rc.startswith(pad_5p):
            return nuclease_type, 'reverse'

        # Check for PERFECT_LEADER at start of RC (common pattern)
        if nuclease_type == 'SpCas9_SpG' and seq_rc.startswith(perfect_leader_spg):
            return nuclease_type, 'reverse'

        if re_site in seq_rc[:80]:
            return nuclease_type, 'reverse'

    # Try fuzzy matching if perfect matching failed
    for nuclease_type, anchors in anchors_dict.items():
        pad_5p = anchors.pad_5p.lower()
        pad_len = len(pad_5p)

        # Fuzzy match forward
        if fuzzy_match(seq_lower[:pad_len + 5], pad_5p, threshold=0.80):
            return nuclease_type, 'forward'

        # Fuzzy match reverse
        if fuzzy_match(seq_rc[:pad_len + 5], pad_5p, threshold=0.80):
            return nuclease_type, 'reverse'

        # Check for partial 5' pad match (last 11bp - PERFECT_LEADER)
        if nuclease_type == 'SpCas9_SpG':
            if fuzzy_match(seq_lower[:15], perfect_leader_spg, threshold=0.80):
                return nuclease_type, 'forward'
            if fuzzy_match(seq_rc[:15], perfect_leader_spg, threshold=0.80):
                return nuclease_type, 'reverse'

    return None, 'unknown'


def infer_amplicon_structure(
    collapsed_fastq_path: Path,
    oligo_design_path: Optional[Path] = None,
    n_sample: int = 5000,
    min_confidence: float = 0.5
) -> StructureInferenceResult:
    """
    Infer amplicon structure from first N reads by matching to known anchors.

    This function auto-detects:
    - Nuclease type (SpCas9/SpG vs LbCas12a)
    - Read orientation (forward vs reverse complement)
    - Appropriate anchor sequences for parsing

    Args:
        collapsed_fastq_path: Path to collapsed FASTQ file
        oligo_design_path: Optional path to oligo design file (for SPS mapping)
        n_sample: Number of reads to sample for inference
        min_confidence: Minimum confidence threshold (0-1)

    Returns:
        StructureInferenceResult with detected structure

    Raises:
        ValueError: If structure cannot be determined with sufficient confidence
    """
    logger.info(f"Inferring amplicon structure from {collapsed_fastq_path.name}")
    logger.info(f"Sampling first {n_sample} reads...")

    # Build SPS mapping if oligo design provided
    if oligo_design_path and oligo_design_path.exists():
        sps_mapping = build_sps_mapping_from_oligo_design(oligo_design_path)
    else:
        sps_mapping = dict(SPS_TO_LIBRARY)

    # Sample reads
    reads = sample_reads_from_fastq(collapsed_fastq_path, n_sample)
    if not reads:
        raise ValueError(f"No reads found in {collapsed_fastq_path}")

    # Count detections
    nuclease_counts: Counter = Counter()
    orientation_counts: Dict[str, Counter] = {
        'SpCas9_SpG': Counter(),
        'LbCas12a': Counter(),
    }
    undetected = 0

    for seq, _ in reads:
        nuclease_type, orientation = detect_structure_from_read(seq, NUCLEASE_ANCHORS)

        if nuclease_type:
            nuclease_counts[nuclease_type] += 1
            orientation_counts[nuclease_type][orientation] += 1
        else:
            undetected += 1

    # Determine dominant nuclease type
    if not nuclease_counts:
        raise ValueError(
            f"Could not detect nuclease type in any of {len(reads)} sampled reads. "
            "Check that reads contain expected anchor sequences."
        )

    dominant_nuclease, dominant_count = nuclease_counts.most_common(1)[0]
    confidence = dominant_count / len(reads)

    # Determine dominant orientation for the dominant nuclease
    orient_counts = orientation_counts[dominant_nuclease]
    if orient_counts['reverse'] > orient_counts['forward']:
        dominant_orientation = 'reverse'
    else:
        dominant_orientation = 'forward'

    # Log results
    logger.info(f"Detected nuclease type: {dominant_nuclease} "
                f"({confidence:.1%} confidence)")
    logger.info(f"Detected orientation: {dominant_orientation}")
    logger.info(f"Forward: {orient_counts['forward']}, "
                f"Reverse: {orient_counts['reverse']}, "
                f"Undetected: {undetected}")

    # Check confidence threshold
    if confidence < min_confidence:
        raise ValueError(
            f"Low detection confidence ({confidence:.1%}). "
            f"Expected >= {min_confidence:.0%}. "
            f"Check that reads match expected oligo structure."
        )

    result = StructureInferenceResult(
        nuclease_type=dominant_nuclease,
        anchors=NUCLEASE_ANCHORS[dominant_nuclease],
        orientation=dominant_orientation,
        confidence=confidence,
        sps_to_library=sps_mapping,
        total_sampled=len(reads),
        detected_forward=orient_counts['forward'],
        detected_reverse=orient_counts['reverse'],
        undetected=undetected,
    )

    return result


def infer_structure_from_multiple_fastqs(
    fastq_paths: List[Path],
    oligo_design_path: Optional[Path] = None,
    n_sample_per_file: int = 1000,
) -> StructureInferenceResult:
    """
    Infer amplicon structure from multiple FASTQ files.

    Useful when you have multiple samples and want to pool detections.

    Args:
        fastq_paths: List of collapsed FASTQ paths
        oligo_design_path: Optional path to oligo design file
        n_sample_per_file: Number of reads to sample per file

    Returns:
        StructureInferenceResult with detected structure
    """
    all_nuclease_counts: Counter = Counter()
    all_orientation_counts: Dict[str, Counter] = {
        'SpCas9_SpG': Counter(),
        'LbCas12a': Counter(),
    }
    total_undetected = 0
    total_sampled = 0

    for fastq_path in fastq_paths:
        if not fastq_path.exists():
            logger.warning(f"Skipping non-existent file: {fastq_path}")
            continue

        reads = sample_reads_from_fastq(fastq_path, n_sample_per_file)
        total_sampled += len(reads)

        for seq, _ in reads:
            nuclease_type, orientation = detect_structure_from_read(seq, NUCLEASE_ANCHORS)

            if nuclease_type:
                all_nuclease_counts[nuclease_type] += 1
                all_orientation_counts[nuclease_type][orientation] += 1
            else:
                total_undetected += 1

    if not all_nuclease_counts:
        raise ValueError("Could not detect nuclease type in any sampled reads.")

    dominant_nuclease, dominant_count = all_nuclease_counts.most_common(1)[0]
    confidence = dominant_count / total_sampled if total_sampled > 0 else 0

    orient_counts = all_orientation_counts[dominant_nuclease]
    dominant_orientation = 'reverse' if orient_counts['reverse'] > orient_counts['forward'] else 'forward'

    # Build SPS mapping
    if oligo_design_path and oligo_design_path.exists():
        sps_mapping = build_sps_mapping_from_oligo_design(oligo_design_path)
    else:
        sps_mapping = dict(SPS_TO_LIBRARY)

    return StructureInferenceResult(
        nuclease_type=dominant_nuclease,
        anchors=NUCLEASE_ANCHORS[dominant_nuclease],
        orientation=dominant_orientation,
        confidence=confidence,
        sps_to_library=sps_mapping,
        total_sampled=total_sampled,
        detected_forward=orient_counts['forward'],
        detected_reverse=orient_counts['reverse'],
        undetected=total_undetected,
    )


# ============================================================================
# Convenience functions for pipeline integration
# ============================================================================

def get_structure_for_nuclease(nuclease_type: str) -> NucleaseAnchors:
    """
    Get structure parameters for a known nuclease type.

    Args:
        nuclease_type: Either 'SpCas9_SpG' or 'LbCas12a'

    Returns:
        NucleaseAnchors dataclass

    Raises:
        ValueError: If nuclease type is unknown
    """
    if nuclease_type not in NUCLEASE_ANCHORS:
        raise ValueError(
            f"Unknown nuclease type: {nuclease_type}. "
            f"Known types: {list(NUCLEASE_ANCHORS.keys())}"
        )
    return NUCLEASE_ANCHORS[nuclease_type]


def get_nuclease_from_sps(sps: str) -> Optional[str]:
    """
    Get nuclease anchor type from SPS sequence.

    Args:
        sps: 20bp SPS sequence

    Returns:
        Nuclease anchor type or None if unknown
    """
    sps_lower = sps.lower()
    if sps_lower in SPS_TO_LIBRARY:
        return SPS_TO_LIBRARY[sps_lower]['anchors']
    return None


# ============================================================================
# Module-level exports
# ============================================================================

__all__ = [
    'NucleaseAnchors',
    'NUCLEASE_ANCHORS',
    'SPS_TO_LIBRARY',
    'StructureInferenceResult',
    'infer_amplicon_structure',
    'infer_structure_from_multiple_fastqs',
    'get_structure_for_nuclease',
    'get_nuclease_from_sps',
    'rev_comp',
    'fuzzy_match',
    'detect_nuclease_from_oligo',
]
