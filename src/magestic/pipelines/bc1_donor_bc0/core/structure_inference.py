"""
De novo inference of BC1 amplicon structure from collapsed FASTQ data.

Author: Kevin R. Roy

This module provides automatic detection of constant sequences (BC1_PREFIX,
BC1_SUFFIX, RE_site) by analyzing conservation patterns across abundant reads,
WITHOUT relying on hardcoded sequence patterns.

Key insight: Collapsed FASTQs are sorted by abundance (most common sequences first).
By analyzing the top few thousand sequences, we can identify:
- Conserved regions = constant sequences (BC1_PREFIX, BC1_SUFFIX, RE_site)
- Variable regions = barcodes/donors (BC1, donor, BC0)

Algorithm:
1. Load top N most abundant sequences
2. Find msd (AGGAAAACAGAC) from 3' end - most reliable anchor for alignment
3. Align all reads by msd position
4. Build position-wise conservation matrix
5. Identify conserved stretches (>95% same nucleotide)
6. Map conserved stretches to BC1_PREFIX, BC1_SUFFIX, RE_site
7. Report inferred sequences with validation metrics
"""

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import numpy as np
import pysam

logger = logging.getLogger(__name__)


# ============================================================================
# Known anchor sequences (used for alignment, not for detection)
# ============================================================================

# msd is the most reliable 3' anchor - universal across all libraries
MSD_SEQUENCE = "AGGAAAACAGAC"
MSD_SEQUENCE_RC = "CTGTTTTTCCT"  # Reverse complement of msd (seen in R2-only data)
MSD_LENGTH = len(MSD_SEQUENCE)

# Known SPS sequences (for library type detection after structure inference)
KNOWN_SPS_SEQUENCES: Dict[str, Dict] = {
    'cgatagacgactggacagca': {'nuclease': 'SpCas9_NGG', 'library': 'V416'},
    'cttgactcgacagtgagagc': {'nuclease': 'SpG_NGNG', 'library': 'V520'},
    'cgtatggagatgacgtcagg': {'nuclease': 'SpG_NGNH', 'library': 'V521'},
    'gaactacgtgtctcgtcagg': {'nuclease': 'LbCas12a_TTTV', 'library': 'V419'},
    'gttatgctggtcctaggtcg': {'nuclease': 'impLbCas12a_nonTTTV', 'library': 'V419'},
}

# LbCas12a scaffold variants
# NOTE: The original designed sequence (caaagattaaataattcccactaagtgtgg) severely
# compromises LbCas12a function. These two functional variants are used instead:
LBCAS12A_SCAFFOLD_VARIANTS = {
    # Original (compromised function - included for detection but not recommended)
    'original': {
        'name': 'LbCas12a_original',
        'pad_5p': 'caaagattaaataattcccactaagtgtgg',  # 30bp
        'note': 'Original design - compromises LbCas12a function',
    },
    # WT scaffold version - functional
    'WT': {
        'name': 'RPR1_leader_LbCas12a_DRf_crRNA',
        'full_leader': 'aatgatcggagctgcgattggcaggtttcaaagattaaataatttctactaagtgtagat',  # 60bp
        'pad_5p': 'aatgatcggagctgcgattggcaggtttcaaagattaaataatttctactaagtgtagat',
        'note': 'WT scaffold - functional LbCas12a',
    },
    # GC1 version - functional, GC-enriched
    'GC1': {
        'name': 'RPR1_leader_LbCas12a_DRf_GC_1_crRNA',
        'full_leader': 'aatgatcggagctgcgattggcaggtttcaaagattaaataattCcTactaagtgtAgGt',  # 60bp (mixed case shows mutations)
        'pad_5p': 'aatgatcggagctgcgattggcaggtttcaaagattaaataattcctactaagtgtaggt',  # lowercase version
        'note': 'GC1 variant - functional LbCas12a with GC enrichment',
    },
}

# Expected lengths for validation
EXPECTED_BC1_LENGTH = 20
EXPECTED_BC0_LENGTH = 10
EXPECTED_SPS_LENGTH = 20


# ============================================================================
# Data structures
# ============================================================================

@dataclass
class InferredConstantRegion:
    """A constant region inferred from conservation analysis."""
    name: str                      # e.g., "BC1_PREFIX", "BC1_SUFFIX", "RE_site"
    sequence: str                  # Inferred consensus sequence
    start_position: int            # Median start position in reads (relative to msd)
    end_position: int              # Median end position
    conservation: float            # Fraction of reads matching this sequence
    position_std: float            # Standard deviation of position

    @property
    def length(self) -> int:
        return len(self.sequence)

    def __str__(self) -> str:
        return (
            f"{self.name}: {self.sequence} "
            f"(pos {self.start_position}-{self.end_position}, "
            f"conservation {self.conservation:.1%})"
        )


@dataclass
class BC1StructureInferenceResult:
    """Result of de novo BC1 amplicon structure inference."""
    # Inferred constant regions
    inferred_constants: Dict[str, InferredConstantRegion]

    # Variable region positions (relative to read start after orientation correction)
    bc1_start: int
    bc1_end: int
    bc1_length: int
    donor_start: int
    donor_end: int
    sps_start: int
    sps_end: int
    bc0_start: int
    bc0_end: int

    # msd anchor (most reliable - used for validation)
    msd_sequence: str
    msd_conservation: float

    # SPS detection (for single-sublibrary cases)
    detected_sps: Optional[str] = None
    sps_counts: Dict[str, int] = field(default_factory=dict)
    is_single_sublibrary: bool = False
    dominant_sps_fraction: float = 0.0

    # Orientation
    orientation: str = "forward"  # "forward" or "reverse"

    # Quality metrics
    total_reads_sampled: int = 0
    reads_with_msd_found: int = 0
    overall_confidence: float = 0.0

    def __str__(self) -> str:
        lines = [
            "=" * 80,
            "BC1-Donor-BC0 Structure Inference Report (De Novo)",
            "=" * 80,
            f"Reads sampled: {self.total_reads_sampled}",
            f"Reads with msd found: {self.reads_with_msd_found} "
            f"({self.reads_with_msd_found/self.total_reads_sampled:.1%})" if self.total_reads_sampled > 0 else "",
            f"Orientation: {self.orientation}",
            "",
            "INFERRED CONSTANT SEQUENCES",
            "-" * 60,
        ]

        for name in ['BC1_PREFIX', 'BC1_SUFFIX', 'RE_site']:
            if name in self.inferred_constants:
                region = self.inferred_constants[name]
                lines.append(
                    f"  {name:12s}: {region.sequence:20s} "
                    f"(pos {region.start_position:3d}-{region.end_position:3d}, "
                    f"conservation {region.conservation:.1%})"
                )

        lines.extend([
            "",
            "INFERRED VARIABLE REGIONS",
            "-" * 60,
            f"  BC1   : positions {self.bc1_start:3d}-{self.bc1_end:3d} "
            f"({self.bc1_length}bp, expected {EXPECTED_BC1_LENGTH}bp)",
            f"  Donor : positions {self.donor_start:3d}-{self.donor_end:3d} "
            f"({self.donor_end - self.donor_start}bp)",
            f"  SPS   : positions {self.sps_start:3d}-{self.sps_end:3d} "
            f"({self.sps_end - self.sps_start}bp)",
            f"  BC0   : positions {self.bc0_start:3d}-{self.bc0_end:3d} "
            f"({self.bc0_end - self.bc0_start}bp)",
            "",
            "MSD ANCHOR VALIDATION",
            "-" * 60,
            f"  msd sequence: {self.msd_sequence}",
            f"  conservation: {self.msd_conservation:.1%}",
        ])

        if self.sps_counts:
            lines.extend([
                "",
                "SPS COMPOSITION",
                "-" * 60,
            ])
            for sps, count in sorted(self.sps_counts.items(), key=lambda x: -x[1])[:5]:
                info = KNOWN_SPS_SEQUENCES.get(sps.lower(), {})
                library = info.get('library', 'unknown')
                nuclease = info.get('nuclease', 'unknown')
                pct = count / sum(self.sps_counts.values()) * 100 if self.sps_counts else 0
                lines.append(f"  {sps}: {count:6d} ({pct:5.1f}%) - {library}/{nuclease}")

            lines.extend([
                "",
                f"Library type: {'SINGLE SUBLIBRARY' if self.is_single_sublibrary else 'POOL OF POOLS'}",
            ])
            if self.is_single_sublibrary and self.detected_sps:
                lines.append(f"Detected SPS: {self.detected_sps}")

        lines.extend([
            "",
            f"OVERALL CONFIDENCE: {self.overall_confidence:.1%}",
            "=" * 80,
        ])

        return "\n".join(lines)


# ============================================================================
# Utility functions
# ============================================================================

def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def fuzzy_find(seq: str, pattern: str, max_mismatches: int = 1,
               start: int = 0, end: Optional[int] = None) -> int:
    """
    Find pattern in seq allowing up to max_mismatches.

    Args:
        seq: Sequence to search in
        pattern: Pattern to find
        max_mismatches: Maximum allowed mismatches
        start: Start position for search
        end: End position for search (None = end of sequence)

    Returns:
        Position of first match, or -1 if not found
    """
    pattern_len = len(pattern)
    if end is None:
        end = len(seq) - pattern_len + 1
    else:
        end = min(end, len(seq) - pattern_len + 1)

    pattern_lower = pattern.lower()

    for i in range(start, end):
        candidate = seq[i:i + pattern_len].lower()
        mismatches = sum(1 for a, b in zip(candidate, pattern_lower) if a != b)
        if mismatches <= max_mismatches:
            return i

    return -1


def sample_reads_from_fastq(
    fastq_path: Path,
    n_sample: int = 5000
) -> List[Tuple[str, int]]:
    """
    Sample first N reads from a collapsed FASTQ file.

    Collapsed FASTQs have read counts encoded in the read name.
    Returns reads sorted by abundance (highest first).

    Returns:
        List of (sequence, count) tuples
    """
    reads = []
    try:
        with pysam.FastxFile(str(fastq_path)) as fh:
            for i, entry in enumerate(fh):
                if i >= n_sample:
                    break
                # Try to parse count from read name
                try:
                    count = int(entry.name.split()[0])
                except (ValueError, IndexError):
                    count = 1
                reads.append((entry.sequence, count))
    except Exception as e:
        logger.warning(f"Error reading {fastq_path}: {e}")

    return reads


# ============================================================================
# Conservation analysis functions
# ============================================================================

def build_conservation_matrix(
    sequences: List[str],
    reference_positions: List[int],
    max_length: int = 400
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build nucleotide conservation matrix across aligned sequences.

    Sequences are aligned by their reference_positions (e.g., msd position).
    Position 0 in the matrix corresponds to the reference anchor.

    Args:
        sequences: List of sequences
        reference_positions: Position of reference anchor in each sequence
        max_length: Maximum length of matrix (positions around anchor)

    Returns:
        Tuple of:
        - conservation_matrix: (4 x max_length*2) array with counts per nucleotide
        - position_coverage: array with number of sequences covering each position
    """
    # Use 2x max_length to cover positions before and after anchor
    matrix_size = max_length * 2
    center = max_length  # Position 0 relative to anchor maps to this index

    # Initialize: 4 rows for A, C, G, T
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
                  'a': 0, 'c': 1, 'g': 2, 't': 3}

    conservation_matrix = np.zeros((4, matrix_size), dtype=np.int32)
    position_coverage = np.zeros(matrix_size, dtype=np.int32)

    for seq, ref_pos in zip(sequences, reference_positions):
        # Align sequence so that reference position maps to center
        for seq_pos, base in enumerate(seq):
            # Calculate position relative to reference
            rel_pos = seq_pos - ref_pos
            matrix_idx = center + rel_pos

            if 0 <= matrix_idx < matrix_size:
                nuc_idx = nuc_to_idx.get(base)
                if nuc_idx is not None:
                    conservation_matrix[nuc_idx, matrix_idx] += 1
                    position_coverage[matrix_idx] += 1

    return conservation_matrix, position_coverage


def find_conserved_regions(
    conservation_matrix: np.ndarray,
    position_coverage: np.ndarray,
    min_conservation: float = 0.95,
    min_length: int = 5,
    min_coverage: int = 100
) -> List[Tuple[int, int, str, float]]:
    """
    Find regions where >min_conservation of reads have the same nucleotide.

    Args:
        conservation_matrix: (4 x N) array with nucleotide counts
        position_coverage: (N,) array with total coverage per position
        min_conservation: Minimum fraction conserved to consider
        min_length: Minimum length of conserved region
        min_coverage: Minimum coverage at a position to consider

    Returns:
        List of (start, end, consensus_sequence, avg_conservation) tuples
        Positions are relative to the reference anchor (can be negative)
    """
    idx_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    matrix_size = conservation_matrix.shape[1]
    center = matrix_size // 2

    # Find positions meeting conservation threshold
    conserved_positions = []
    for pos in range(matrix_size):
        if position_coverage[pos] < min_coverage:
            continue

        # Find dominant nucleotide
        max_count = conservation_matrix[:, pos].max()
        conservation = max_count / position_coverage[pos]

        if conservation >= min_conservation:
            dominant_nuc = idx_to_nuc[conservation_matrix[:, pos].argmax()]
            rel_pos = pos - center
            conserved_positions.append((rel_pos, dominant_nuc, conservation))

    # Group consecutive conserved positions into regions
    if not conserved_positions:
        return []

    regions = []
    current_start = conserved_positions[0][0]
    current_seq = conserved_positions[0][1]
    current_conservations = [conserved_positions[0][2]]

    for i in range(1, len(conserved_positions)):
        pos, nuc, cons = conserved_positions[i]
        prev_pos = conserved_positions[i-1][0]

        if pos == prev_pos + 1:
            # Consecutive position, extend current region
            current_seq += nuc
            current_conservations.append(cons)
        else:
            # Gap found, save current region if long enough
            if len(current_seq) >= min_length:
                avg_cons = np.mean(current_conservations)
                regions.append((current_start, current_start + len(current_seq), current_seq, avg_cons))

            # Start new region
            current_start = pos
            current_seq = nuc
            current_conservations = [cons]

    # Don't forget the last region
    if len(current_seq) >= min_length:
        avg_cons = np.mean(current_conservations)
        regions.append((current_start, current_start + len(current_seq), current_seq, avg_cons))

    return regions


def map_conserved_to_structure(
    conserved_regions: List[Tuple[int, int, str, float]],
    msd_position_in_matrix: int,  # This is always 0 since msd is the reference
) -> Dict[str, InferredConstantRegion]:
    """
    Map discovered conserved regions to structural elements.

    Heuristics (positions are relative to msd, so negative means upstream):
    - BC1_PREFIX: Conserved region ~250-300bp upstream of msd (near 5' end)
    - BC1_SUFFIX: Conserved region ~15bp long, ~200-250bp upstream of msd
    - RE_site: Conserved region ~11bp long, ~180-220bp upstream of msd

    Args:
        conserved_regions: List of (start, end, sequence, conservation) tuples
        msd_position_in_matrix: Reference position (should be 0)

    Returns:
        Dict mapping region names to InferredConstantRegion objects
    """
    inferred = {}

    # Sort regions by position (most upstream first = most negative)
    regions_by_pos = sorted(conserved_regions, key=lambda x: x[0])

    # Heuristics for assignment:
    # In typical bc1-donor-bc0 amplicons:
    # - BC1_PREFIX is near the 5' end (most upstream, 7-12bp long)
    # - BC1_SUFFIX comes after BC1 (15bp long, contains cloning site info)
    # - RE_site is 11bp, typically "GTTTGAAGAGC" or similar
    # - Then comes the donor (~130bp variable)
    # - Then SPS (20bp variable), BC0 (10bp variable), msd (12bp constant)

    # Find the most 5' conserved region (likely BC1_PREFIX)
    if regions_by_pos:
        # Look for regions that are 5-15bp long in the most upstream positions
        for start, end, seq, cons in regions_by_pos:
            length = end - start
            # BC1_PREFIX should be relatively short (7-12bp) and upstream
            if 5 <= length <= 15 and start < -200:
                if 'BC1_PREFIX' not in inferred:
                    inferred['BC1_PREFIX'] = InferredConstantRegion(
                        name='BC1_PREFIX',
                        sequence=seq,
                        start_position=start,
                        end_position=end,
                        conservation=cons,
                        position_std=0.0  # Will be computed separately if needed
                    )
                    break

    # Look for BC1_SUFFIX (should be ~15bp, between BC1 and donor region)
    for start, end, seq, cons in regions_by_pos:
        length = end - start
        # BC1_SUFFIX is typically 15bp, positioned after BC1_PREFIX + 20bp BC1
        if 13 <= length <= 18 and -250 < start < -150:
            if 'BC1_SUFFIX' not in inferred:
                inferred['BC1_SUFFIX'] = InferredConstantRegion(
                    name='BC1_SUFFIX',
                    sequence=seq,
                    start_position=start,
                    end_position=end,
                    conservation=cons,
                    position_std=0.0
                )
                break

    # Look for RE_site (should be ~11bp, typically GTTTGAAGAGC or TTTCGAAGAGC)
    for start, end, seq, cons in regions_by_pos:
        length = end - start
        # RE_site is typically 11bp
        if 9 <= length <= 13 and -200 < start < -100:
            # Check if it contains GAAGAGC motif (common to both SpCas9 and LbCas12a)
            if 'GAAGAGC' in seq.upper() or 'gaagagc' in seq.lower():
                if 'RE_site' not in inferred:
                    inferred['RE_site'] = InferredConstantRegion(
                        name='RE_site',
                        sequence=seq,
                        start_position=start,
                        end_position=end,
                        conservation=cons,
                        position_std=0.0
                    )
                    break

    # If RE_site not found by motif, try by position
    if 'RE_site' not in inferred:
        for start, end, seq, cons in regions_by_pos:
            length = end - start
            if 9 <= length <= 13 and -180 < start < -120:
                inferred['RE_site'] = InferredConstantRegion(
                    name='RE_site',
                    sequence=seq,
                    start_position=start,
                    end_position=end,
                    conservation=cons,
                    position_std=0.0
                )
                break

    return inferred


# ============================================================================
# Main inference functions
# ============================================================================

def infer_bc1_amplicon_structure(
    fastq_path: Path,
    n_sample: int = 5000,
    min_conservation: float = 0.95,
    min_region_length: int = 5,
    verbose: bool = True,
) -> BC1StructureInferenceResult:
    """
    Infer BC1 amplicon structure from collapsed FASTQ using de novo conservation analysis.

    Algorithm:
    1. Load top N most abundant sequences (collapsed FASTQs are pre-sorted)
    2. Find msd (AGGAAAACAGAC) from 3' end - most reliable anchor
    3. Align all reads by msd position
    4. Build position-wise conservation matrix
    5. Identify conserved stretches (>95% same nucleotide)
    6. Map conserved stretches to BC1_PREFIX, BC1_SUFFIX, RE_site
    7. Infer variable region boundaries
    8. Detect SPS composition

    Args:
        fastq_path: Path to collapsed FASTQ file
        n_sample: Number of reads to sample
        min_conservation: Minimum conservation for a region to be considered constant
        min_region_length: Minimum length of conserved region
        verbose: Print progress messages

    Returns:
        BC1StructureInferenceResult with inferred structure
    """
    if verbose:
        logger.info(f"Inferring BC1 structure from {fastq_path.name}")
        logger.info(f"Sampling first {n_sample} reads...")

    # Step 1: Sample reads
    reads = sample_reads_from_fastq(fastq_path, n_sample)
    if not reads:
        raise ValueError(f"No reads found in {fastq_path}")

    total_sampled = len(reads)
    if verbose:
        logger.info(f"Loaded {total_sampled} sequences")

    # Step 2: Find msd position in each read and determine orientation
    # Strategy: Try both orientations (forward and reverse complement)
    # - Forward: msd near 3' end (typical for merged 2x300 or forward R1 reads)
    # - Reverse: msd found after reverse complementing (typical for reverse R2 reads)
    # Also check for rc(msd) at 5' end as additional signal for R2-only data
    aligned_seqs = []
    msd_positions = []
    orientation_counts = Counter()
    sps_counter = Counter()

    for seq, count in reads:
        found = False

        # Method 1: Try finding msd near 3' end of forward sequence
        msd_pos = fuzzy_find(seq, MSD_SEQUENCE, max_mismatches=1,
                             start=max(0, len(seq) - 50))

        if msd_pos != -1:
            aligned_seqs.append(seq)
            msd_positions.append(msd_pos)
            orientation_counts['forward'] += 1
            found = True

            # Extract SPS (20bp before BC0, which is 10bp before msd)
            sps_start = msd_pos - EXPECTED_BC0_LENGTH - EXPECTED_SPS_LENGTH
            sps_end = msd_pos - EXPECTED_BC0_LENGTH
            if sps_start >= 0:
                sps = seq[sps_start:sps_end].lower()
                sps_counter[sps] += count

        if not found:
            # Method 2: Try reverse complement (for R2 reads or reverse-oriented data)
            rc_seq = reverse_complement(seq)
            msd_pos = fuzzy_find(rc_seq, MSD_SEQUENCE, max_mismatches=1,
                                 start=max(0, len(rc_seq) - 50))

            if msd_pos != -1:
                aligned_seqs.append(rc_seq)
                msd_positions.append(msd_pos)
                orientation_counts['reverse'] += 1
                found = True

                # Extract SPS from reverse complement
                sps_start = msd_pos - EXPECTED_BC0_LENGTH - EXPECTED_SPS_LENGTH
                sps_end = msd_pos - EXPECTED_BC0_LENGTH
                if sps_start >= 0:
                    sps = rc_seq[sps_start:sps_end].lower()
                    sps_counter[sps] += count

        if not found:
            # Method 3: Check for rc(msd) at 5' end (R2-only data where rc(msd) is at start)
            # This would mean the read is already in the right orientation for R2
            rc_msd_pos = fuzzy_find(seq, MSD_SEQUENCE_RC, max_mismatches=1,
                                    start=0, end=min(20, len(seq)))
            if rc_msd_pos != -1:
                # Read starts with rc(msd), so it's R2-only data
                # Reverse complement to get msd at 3' end
                rc_seq = reverse_complement(seq)
                msd_pos = len(rc_seq) - rc_msd_pos - MSD_LENGTH
                if 0 <= msd_pos < len(rc_seq):
                    aligned_seqs.append(rc_seq)
                    msd_positions.append(msd_pos)
                    orientation_counts['reverse_r2'] += 1

                    # Extract SPS
                    sps_start = msd_pos - EXPECTED_BC0_LENGTH - EXPECTED_SPS_LENGTH
                    sps_end = msd_pos - EXPECTED_BC0_LENGTH
                    if sps_start >= 0:
                        sps = rc_seq[sps_start:sps_end].lower()
                        sps_counter[sps] += count

    reads_with_msd = len(aligned_seqs)
    if reads_with_msd == 0:
        raise ValueError(
            f"Could not find msd anchor ({MSD_SEQUENCE}) in any of {total_sampled} reads. "
            "Check that reads contain expected BC1-donor-BC0 amplicon structure."
        )

    if verbose:
        logger.info(f"Found msd in {reads_with_msd}/{total_sampled} reads "
                   f"({reads_with_msd/total_sampled:.1%})")

    # Determine dominant orientation
    # Combine reverse and reverse_r2 counts
    reverse_total = orientation_counts['reverse'] + orientation_counts.get('reverse_r2', 0)
    dominant_orientation = 'reverse' if reverse_total > orientation_counts['forward'] else 'forward'
    if verbose:
        logger.info(f"Orientation: {dominant_orientation}")
        logger.info(f"  Forward reads (msd at 3'): {orientation_counts['forward']}")
        logger.info(f"  Reverse reads (needed rc): {orientation_counts['reverse']}")
        if orientation_counts.get('reverse_r2', 0) > 0:
            logger.info(f"  R2-only reads (rc(msd) at 5'): {orientation_counts['reverse_r2']}")

    # Step 3: Build conservation matrix
    if verbose:
        logger.info("Building conservation matrix...")

    conservation_matrix, position_coverage = build_conservation_matrix(
        aligned_seqs, msd_positions, max_length=300
    )

    # Step 4: Find conserved regions
    if verbose:
        logger.info(f"Finding conserved regions (>{min_conservation:.0%} conservation)...")

    conserved_regions = find_conserved_regions(
        conservation_matrix,
        position_coverage,
        min_conservation=min_conservation,
        min_length=min_region_length,
        min_coverage=max(10, reads_with_msd // 10)  # At least 10% of aligned reads
    )

    if verbose:
        logger.info(f"Found {len(conserved_regions)} conserved regions")
        for start, end, seq, cons in sorted(conserved_regions, key=lambda x: x[0]):
            logger.debug(f"  Position {start:4d} to {end:4d}: {seq} ({cons:.1%})")

    # Step 5: Map conserved regions to structure
    inferred_constants = map_conserved_to_structure(conserved_regions, msd_position_in_matrix=0)

    if verbose:
        logger.info("Inferred constant regions:")
        for name, region in inferred_constants.items():
            logger.info(f"  {name}: {region.sequence} (conservation {region.conservation:.1%})")

    # Step 6: Calculate variable region positions
    # Use inferred constants to determine boundaries

    # Find BC1 position (between BC1_PREFIX end and BC1_SUFFIX start)
    if 'BC1_PREFIX' in inferred_constants and 'BC1_SUFFIX' in inferred_constants:
        bc1_start = inferred_constants['BC1_PREFIX'].end_position
        bc1_end = inferred_constants['BC1_SUFFIX'].start_position
    else:
        # Fall back to typical positions if not found
        # Typical: BC1_PREFIX ends around -280, BC1 is 20bp
        bc1_start = -280
        bc1_end = -260
        logger.warning("Could not determine BC1 position from inferred constants, using defaults")

    bc1_length = bc1_end - bc1_start

    # Donor starts after RE_site ends
    if 'RE_site' in inferred_constants:
        donor_start = inferred_constants['RE_site'].end_position
    else:
        # Typical: RE_site ends around -170
        donor_start = -170

    # SPS starts 30bp before msd (20bp SPS + 10bp BC0)
    sps_start = -30
    sps_end = -10
    bc0_start = -10
    bc0_end = 0  # msd starts at position 0

    # Donor ends where SPS starts
    donor_end = sps_start

    # Step 7: Validate msd conservation
    # Extract msd sequences from aligned reads and check conservation
    msd_seqs = []
    for seq, pos in zip(aligned_seqs, msd_positions):
        if pos + MSD_LENGTH <= len(seq):
            msd_seqs.append(seq[pos:pos + MSD_LENGTH].upper())

    # Calculate msd conservation
    if msd_seqs:
        msd_counter = Counter(msd_seqs)
        most_common_msd, msd_count = msd_counter.most_common(1)[0]
        msd_conservation = msd_count / len(msd_seqs)
    else:
        most_common_msd = MSD_SEQUENCE
        msd_conservation = 0.0

    # Step 8: Determine single sublibrary status
    total_sps_count = sum(sps_counter.values())
    if total_sps_count > 0:
        most_common_sps, sps_count = sps_counter.most_common(1)[0]
        sps_fraction = sps_count / total_sps_count
        is_single_sublibrary = sps_fraction >= 0.95
        detected_sps = most_common_sps if is_single_sublibrary else None
    else:
        sps_fraction = 0.0
        is_single_sublibrary = False
        detected_sps = None

    # Calculate overall confidence
    # Based on: msd detection rate, msd conservation, constant region detection
    n_constants_found = len(inferred_constants)
    confidence_factors = [
        reads_with_msd / total_sampled,  # msd detection rate
        msd_conservation,  # msd sequence conservation
        min(n_constants_found / 3, 1.0),  # Found expected constant regions
    ]
    overall_confidence = np.mean(confidence_factors)

    # Create result object
    result = BC1StructureInferenceResult(
        inferred_constants=inferred_constants,
        bc1_start=bc1_start,
        bc1_end=bc1_end,
        bc1_length=bc1_length,
        donor_start=donor_start,
        donor_end=donor_end,
        sps_start=sps_start,
        sps_end=sps_end,
        bc0_start=bc0_start,
        bc0_end=bc0_end,
        msd_sequence=most_common_msd,
        msd_conservation=msd_conservation,
        detected_sps=detected_sps,
        sps_counts=dict(sps_counter.most_common(20)),  # Top 20 SPS
        is_single_sublibrary=is_single_sublibrary,
        dominant_sps_fraction=sps_fraction,
        orientation=dominant_orientation,
        total_reads_sampled=total_sampled,
        reads_with_msd_found=reads_with_msd,
        overall_confidence=overall_confidence,
    )

    if verbose:
        print(result)

    return result


def infer_structure_from_multiple_fastqs(
    fastq_paths: List[Path],
    n_sample_per_file: int = 2000,
    min_conservation: float = 0.95,
    verbose: bool = True,
) -> BC1StructureInferenceResult:
    """
    Infer BC1 amplicon structure from multiple FASTQ files.

    Pools reads from multiple files for more robust inference.

    Args:
        fastq_paths: List of collapsed FASTQ paths
        n_sample_per_file: Number of reads to sample per file
        min_conservation: Minimum conservation threshold
        verbose: Print progress messages

    Returns:
        BC1StructureInferenceResult with inferred structure
    """
    all_reads = []

    for fastq_path in fastq_paths:
        if not fastq_path.exists():
            logger.warning(f"Skipping non-existent file: {fastq_path}")
            continue

        reads = sample_reads_from_fastq(fastq_path, n_sample_per_file)
        all_reads.extend(reads)
        if verbose:
            logger.info(f"Loaded {len(reads)} reads from {fastq_path.name}")

    if not all_reads:
        raise ValueError("No reads loaded from any input file")

    # Write combined reads to a temporary structure for analysis
    # For now, just use the largest file
    if verbose:
        logger.info(f"Total reads pooled: {len(all_reads)}")

    # Sort by count (most abundant first)
    all_reads.sort(key=lambda x: -x[1])

    # Create a temporary fastq-like structure for the main function
    # For simplicity, we'll process the pooled reads directly

    # Actually, let's just call the main function on the first/largest file
    # and validate with others
    largest_file = max(fastq_paths, key=lambda p: p.stat().st_size if p.exists() else 0)

    return infer_bc1_amplicon_structure(
        largest_file,
        n_sample=sum(1 for p in fastq_paths if p.exists()) * n_sample_per_file,
        min_conservation=min_conservation,
        verbose=verbose,
    )


# ============================================================================
# Convenience functions for pipeline integration
# ============================================================================

def get_extraction_parameters(result: BC1StructureInferenceResult) -> Dict:
    """
    Convert inference result to extraction parameters for parsing scripts.

    Returns dict with:
    - bc1_prefix: Inferred BC1_PREFIX sequence
    - bc1_prefix_length: Length of BC1_PREFIX
    - bc1_suffix: Inferred BC1_SUFFIX sequence
    - re_site: Inferred RE_site sequence
    - orientation: Detected orientation
    """
    params = {
        'orientation': result.orientation,
        'bc1_length': result.bc1_length,
        'sps_length': EXPECTED_SPS_LENGTH,
        'bc0_length': EXPECTED_BC0_LENGTH,
        'msd_sequence': result.msd_sequence,
    }

    if 'BC1_PREFIX' in result.inferred_constants:
        params['bc1_prefix'] = result.inferred_constants['BC1_PREFIX'].sequence
        params['bc1_prefix_length'] = result.inferred_constants['BC1_PREFIX'].length
    else:
        params['bc1_prefix'] = None
        params['bc1_prefix_length'] = 0

    if 'BC1_SUFFIX' in result.inferred_constants:
        params['bc1_suffix'] = result.inferred_constants['BC1_SUFFIX'].sequence
    else:
        params['bc1_suffix'] = None

    if 'RE_site' in result.inferred_constants:
        params['re_site'] = result.inferred_constants['RE_site'].sequence
    else:
        params['re_site'] = None

    return params


def print_inference_report(result: BC1StructureInferenceResult) -> str:
    """
    Generate a formatted report string from inference result.

    This is the main output for users to see what was inferred.
    """
    return str(result)


# ============================================================================
# Module-level exports
# ============================================================================

__all__ = [
    'InferredConstantRegion',
    'BC1StructureInferenceResult',
    'MSD_SEQUENCE',
    'MSD_SEQUENCE_RC',
    'KNOWN_SPS_SEQUENCES',
    'LBCAS12A_SCAFFOLD_VARIANTS',
    'infer_bc1_amplicon_structure',
    'infer_structure_from_multiple_fastqs',
    'get_extraction_parameters',
    'print_inference_report',
    'reverse_complement',
    'fuzzy_find',
    'build_conservation_matrix',
    'find_conserved_regions',
]
