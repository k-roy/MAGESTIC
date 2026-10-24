#!/usr/bin/env python3
"""
Step 02: Parse BC1 Barcodes from Collapsed FASTQ Files

Extracts 20bp BC1 barcodes from collapsed amplicon sequencing data by
searching for flanking sequences in both forward and reverse-complement.

Features:
- Auto-inference of BC1_PREFIX and BC1_SUFFIX from read conservation analysis
- Optional validation against expected sequences from bc1_donor_bc0 pipeline
- Handles both forward and reverse-complement orientations

Usage:
    # With auto-inference (recommended)
    python 02_parse_bc1_reads.py \\
        --input collapsed.fastq.gz \\
        --output bc1_counts.tsv \\
        --auto-infer

    # With explicit sequences (legacy)
    python 02_parse_bc1_reads.py \\
        --input collapsed.fastq.gz \\
        --output bc1_counts.tsv \\
        --bc1-prefix CATGCTC \\
        --bc1-suffix CCCTAGG

    # With validation against expected library
    python 02_parse_bc1_reads.py \\
        --input collapsed.fastq.gz \\
        --output bc1_counts.tsv \\
        --auto-infer \\
        --expected-prefix CATGCTC \\
        --expected-suffix CCCTAGG

Author: Kevin R. Roy
"""

import argparse
import gzip
import logging
import sys
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
from magestic.utils.bc1_extraction import (
    extract_bc1 as _canonical_extract_bc1, ANCHOR_TIERS, EMITTING_TIERS)
from magestic.utils.collapsed_fastq import (
    iter_collapsed_fastq, COLLAPSED_HEADER_FORMATS)

logger = logging.getLogger(__name__)


# ============================================================================
# Auto-Inference Functions
# ============================================================================

def sample_reads_from_collapsed_fastq(
    fastq_path: Path,
    n_sample: int = 5000
) -> List[Tuple[str, int]]:
    """
    Sample first N reads from a collapsed FASTQ file.

    Collapsed FASTQs have read counts encoded in the read name. Both producer
    formats are handled by magestic.utils.collapsed_fastq -- see there, and
    planning/210_bc1_recall_code_gate.md, for why this is not inlined.

    Returns:
        List of (sequence, count) tuples
    """
    return [
        (sequence, count)
        for sequence, count, _qual, _fmt in iter_collapsed_fastq(
            fastq_path, limit=n_sample)
    ]


def build_position_conservation(
    sequences: List[str],
    weights: Optional[List[int]] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build nucleotide conservation matrix for aligned sequences.

    Assumes all sequences are similar length and roughly aligned.

    Returns:
        Tuple of (conservation_matrix, coverage)
        - conservation_matrix: (4 x max_len) array with weighted counts
        - coverage: (max_len,) array with total weight per position
    """
    max_len = max(len(s) for s in sequences)
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
                  'a': 0, 'c': 1, 'g': 2, 't': 3}

    matrix = np.zeros((4, max_len), dtype=np.float64)
    coverage = np.zeros(max_len, dtype=np.float64)

    if weights is None:
        weights = [1] * len(sequences)

    for seq, weight in zip(sequences, weights):
        for pos, base in enumerate(seq):
            idx = nuc_to_idx.get(base)
            if idx is not None:
                matrix[idx, pos] += weight
                coverage[pos] += weight

    return matrix, coverage


def find_conserved_flanking_regions(
    conservation_matrix: np.ndarray,
    coverage: np.ndarray,
    min_conservation: float = 0.95,
    min_length: int = 5,
    min_coverage_frac: float = 0.5
) -> List[Tuple[int, int, str, float]]:
    """
    Find conserved regions in the sequence.

    Returns:
        List of (start, end, consensus_seq, avg_conservation)
    """
    idx_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    max_coverage = coverage.max()
    min_coverage = max_coverage * min_coverage_frac

    # Find positions meeting conservation threshold
    conserved = []
    for pos in range(len(coverage)):
        if coverage[pos] < min_coverage:
            continue

        max_count = conservation_matrix[:, pos].max()
        cons_frac = max_count / coverage[pos]

        if cons_frac >= min_conservation:
            nuc = idx_to_nuc[conservation_matrix[:, pos].argmax()]
            conserved.append((pos, nuc, cons_frac))

    if not conserved:
        return []

    # Group into contiguous regions
    regions = []
    start = conserved[0][0]
    seq = conserved[0][1]
    cons_vals = [conserved[0][2]]

    for i in range(1, len(conserved)):
        pos, nuc, cons = conserved[i]
        prev_pos = conserved[i - 1][0]

        if pos == prev_pos + 1:
            seq += nuc
            cons_vals.append(cons)
        else:
            if len(seq) >= min_length:
                regions.append((start, start + len(seq), seq, np.mean(cons_vals)))
            start = pos
            seq = nuc
            cons_vals = [cons]

    # Last region
    if len(seq) >= min_length:
        regions.append((start, start + len(seq), seq, np.mean(cons_vals)))

    return regions


# Known BC1_PREFIX patterns from various libraries
# These are constant sequences that immediately precede the 20bp BC1
KNOWN_BC1_PREFIXES = [
    "CATGCTC",           # Standard SpCas9/SpG (7bp)
    "CAGGTCATGCTC",      # Extended prefix from 2x150 reads (12bp)
    "GGTCATGCTC",        # Alternative (10bp)
]

# Known BC1_SUFFIX patterns (immediately after BC1)
KNOWN_BC1_SUFFIXES = [
    "CCCTAGG",           # Standard SpCas9/SpG (7bp)
    "CCCTAGGATACG",      # Extended suffix (12bp)
]

# Reverse complement patterns for detecting RC orientation
# In RC reads: RC(SUFFIX) appears BEFORE BC1, RC(PREFIX) appears AFTER BC1
# RC(CCCTAGG) = CCTAGGG, RC(CATGCTC) = GAGCATG
KNOWN_RC_PREFIXES = [
    "CCTAGGG",           # RC of CCCTAGG - appears before BC1 in RC reads
    "CGTATCCTAGGG",      # Extended RC prefix
]

KNOWN_RC_SUFFIXES = [
    "GAGCATG",           # RC of CATGCTC - appears after BC1 in RC reads
    "GAGCATGACC",        # Extended RC suffix
]


def find_known_bc1_prefix(
    sequences: List[str],
    weights: List[int],
    known_prefixes: List[str] = None,
) -> Optional[Tuple[str, List[int], float]]:
    """
    Search for known BC1_PREFIX patterns in reads.

    Returns:
        Tuple of (prefix_found, positions_per_read, match_rate) or None
    """
    if known_prefixes is None:
        known_prefixes = KNOWN_BC1_PREFIXES

    total_weight = sum(weights)

    for prefix in known_prefixes:
        positions = []
        matched_weight = 0

        for seq, weight in zip(sequences, weights):
            pos = seq.upper().find(prefix.upper())
            positions.append(pos)
            if pos >= 0:
                matched_weight += weight

        match_rate = matched_weight / total_weight
        if match_rate >= 0.5:  # At least 50% of reads match
            return prefix, positions, match_rate

    return None


def find_known_rc_prefix(
    sequences: List[str],
    weights: List[int],
    known_rc_prefixes: List[str] = None,
) -> Optional[Tuple[str, List[int], float]]:
    """
    Search for known RC BC1_PREFIX patterns (these appear BEFORE BC1 in RC reads).

    In RC orientation, the suffix (CCCTAGG) becomes RC(CCCTAGG)=CCTAGGG and
    appears before the BC1.

    Returns:
        Tuple of (rc_prefix_found, positions_per_read, match_rate) or None
    """
    if known_rc_prefixes is None:
        known_rc_prefixes = KNOWN_RC_PREFIXES

    total_weight = sum(weights)

    for rc_prefix in known_rc_prefixes:
        positions = []
        matched_weight = 0

        for seq, weight in zip(sequences, weights):
            pos = seq.upper().find(rc_prefix.upper())
            positions.append(pos)
            if pos >= 0:
                matched_weight += weight

        match_rate = matched_weight / total_weight
        if match_rate >= 0.5:  # At least 50% of reads match
            return rc_prefix, positions, match_rate

    return None


def find_common_anchor(
    sequences: List[str],
    weights: List[int],
    min_length: int = 5,
    max_length: int = 15,
    min_frequency: float = 0.5,
    prefer_5prime: bool = True,
) -> Optional[Tuple[str, List[int]]]:
    """
    Find a common anchor sequence that appears in most reads.

    Uses k-mer counting to find the most frequent subsequence,
    with preference for anchors near the 5' end.

    Returns:
        Tuple of (anchor_sequence, list_of_positions_per_read) or None
    """
    # Count k-mers of various lengths
    for k in range(max_length, min_length - 1, -1):
        kmer_counts = Counter()
        kmer_positions = defaultdict(list)
        kmer_median_pos = {}  # Track median position for 5' preference

        for i, (seq, weight) in enumerate(zip(sequences, weights)):
            seen_in_read = set()
            for pos in range(len(seq) - k + 1):
                kmer = seq[pos:pos + k]
                if kmer not in seen_in_read:
                    kmer_counts[kmer] += weight
                    seen_in_read.add(kmer)
                kmer_positions[kmer].append((i, pos))

        # Calculate median position for each k-mer
        total_weight = sum(weights)
        candidates = []

        for kmer, count in kmer_counts.most_common(50):
            freq = count / total_weight
            if freq >= min_frequency:
                # Calculate median position
                all_positions = [pos for _, pos in kmer_positions[kmer]]
                median_pos = sorted(all_positions)[len(all_positions) // 2]

                # Score: prefer high frequency AND 5' position
                if prefer_5prime:
                    # Penalize k-mers found far from 5' end
                    # Position 0-30 is best, penalize positions > 50
                    position_score = max(0, 1 - (median_pos / 100))
                    combined_score = freq * (0.5 + 0.5 * position_score)
                else:
                    combined_score = freq

                candidates.append((kmer, count, median_pos, combined_score))

        # Sort by combined score
        candidates.sort(key=lambda x: -x[3])

        if candidates:
            best_kmer = candidates[0][0]
            positions = [-1] * len(sequences)
            for read_idx, pos in kmer_positions[best_kmer]:
                if positions[read_idx] == -1:
                    positions[read_idx] = pos
            return best_kmer, positions

    return None


def infer_bc1_flanking_sequences(
    fastq_path: Path,
    bc1_length: int = 20,
    n_sample: int = 5000,
    min_conservation: float = 0.90,
    verbose: bool = True
) -> Dict[str, str]:
    """
    Auto-infer BC1_PREFIX and BC1_SUFFIX from collapsed FASTQ data.

    Algorithm:
    1. Sample top N reads (collapsed FASTQs are sorted by abundance)
    2. Find common anchor sequence (likely BC1_PREFIX or known constant)
    3. Align reads by anchor position
    4. Build position-wise conservation matrix
    5. Find conserved regions flanking the variable BC1 region

    Returns:
        Dict with 'bc1_prefix', 'bc1_suffix', 'bc1_start', 'bc1_end',
        'orientation', 'confidence'
    """
    if verbose:
        logger.info(f"Auto-inferring BC1 flanking sequences from {fastq_path.name}")
        logger.info(f"Sampling first {n_sample} reads...")

    # Step 1: Sample reads
    reads = sample_reads_from_collapsed_fastq(fastq_path, n_sample)
    if not reads:
        raise ValueError(f"No reads found in {fastq_path}")

    sequences = [r[0] for r in reads]
    weights = [r[1] for r in reads]

    if verbose:
        logger.info(f"Loaded {len(sequences)} sequences")

    # Step 2: First try known BC1_PREFIX patterns (most reliable)
    known_result = find_known_bc1_prefix(sequences, weights)

    # Track orientation for correct BC1 extraction
    is_rc_orientation = False

    if known_result:
        anchor_seq, anchor_positions, match_rate = known_result
        if verbose:
            n_found = sum(1 for p in anchor_positions if p >= 0)
            logger.info(f"Found known BC1_PREFIX: {anchor_seq} (in {n_found}/{len(sequences)} reads, {match_rate:.1%} match rate)")
    else:
        if verbose:
            logger.info("No known BC1_PREFIX found, checking for RC orientation...")

        # Step 2b: Try known RC_PREFIX patterns (reverse complement orientation)
        # In RC reads, RC(CCCTAGG)=CCTAGGG appears BEFORE the BC1
        rc_result = find_known_rc_prefix(sequences, weights)

        if rc_result:
            anchor_seq, anchor_positions, match_rate = rc_result
            is_rc_orientation = True
            if verbose:
                n_found = sum(1 for p in anchor_positions if p >= 0)
                logger.info(f"Found RC orientation with prefix: {anchor_seq} (in {n_found}/{len(sequences)} reads, {match_rate:.1%} match rate)")
        else:
            if verbose:
                logger.info("No known RC_PREFIX found, searching for common anchor...")

            # Fall back to finding a common anchor
            anchor_result = find_common_anchor(sequences, weights, min_length=6, max_length=12, prefer_5prime=True)

            if anchor_result is None:
                if verbose:
                    logger.warning("No common anchor found, using position-based analysis")
                anchor_seq = None
                anchor_positions = None
            else:
                anchor_seq, anchor_positions = anchor_result
                if verbose:
                    n_found = sum(1 for p in anchor_positions if p >= 0)
                    logger.info(f"Found anchor: {anchor_seq} (in {n_found}/{len(sequences)} reads)")

    # Step 3: Align reads by anchor and build conservation matrix
    if anchor_positions:
        # Filter to reads with anchor and align
        aligned_sequences = []
        aligned_weights = []

        for seq, weight, pos in zip(sequences, weights, anchor_positions):
            if pos >= 0:
                # Align so anchor is at position 0
                # We want to see what's BEFORE the anchor and AFTER
                aligned_sequences.append(seq)
                aligned_weights.append(weight)

        if not aligned_sequences:
            raise ValueError("No reads matched the anchor sequence")

        # Build aligned conservation matrix (relative to anchor)
        # The anchor itself is the BC1_PREFIX candidate
        # Look for what comes after anchor + BC1 (20bp)

        # Extract sequences after anchor to find BC1_SUFFIX
        suffix_candidates = Counter()
        for seq, weight, pos in zip(sequences, weights, anchor_positions):
            if pos >= 0:
                bc1_end = pos + len(anchor_seq) + bc1_length
                if bc1_end + 10 <= len(seq):
                    # Extract up to 10bp after BC1
                    suffix = seq[bc1_end:bc1_end + 10]
                    suffix_candidates[suffix] += weight

        # Find most common suffix
        if suffix_candidates:
            bc1_suffix_full, _ = suffix_candidates.most_common(1)[0]
            # Trim to a reasonable length (use first 7-10bp)
            bc1_suffix = bc1_suffix_full[:7] if len(bc1_suffix_full) >= 7 else bc1_suffix_full
        else:
            bc1_suffix = None

        bc1_prefix = anchor_seq

        # Set orientation based on whether we detected RC orientation
        if is_rc_orientation:
            orientation = "reverse"
            # For RC orientation, the "prefix" we found is actually RC(original_suffix)
            # and the "suffix" we found is actually RC(original_prefix)
            # We need to swap and RC them for the output
            # But for extraction, we use them as-is (extract after anchor, then RC the BC1)
            if verbose:
                logger.info(f"RC orientation detected: will reverse-complement extracted BC1s")
        else:
            orientation = "forward"

        # Calculate positions (relative to typical read structure)
        # BC1 starts right after the anchor
        bc1_start_typical = len(anchor_seq)
        bc1_end_typical = bc1_start_typical + bc1_length

    else:
        # Fallback: position-based analysis
        matrix, coverage = build_position_conservation(sequences, weights)
        regions = find_conserved_flanking_regions(
            matrix, coverage,
            min_conservation=min_conservation,
            min_length=5
        )

        if verbose:
            logger.info(f"Found {len(regions)} conserved regions")
            for start, end, seq, cons in regions:
                logger.info(f"  Position {start:3d}-{end:3d}: {seq} ({cons:.1%})")

        # Look for gap of ~20bp
        bc1_prefix = None
        bc1_suffix = None
        orientation = "forward"

        regions_sorted = sorted(regions, key=lambda x: x[0])

        for i in range(len(regions_sorted) - 1):
            r1_start, r1_end, r1_seq, r1_cons = regions_sorted[i]
            r2_start, r2_end, r2_seq, r2_cons = regions_sorted[i + 1]

            gap = r2_start - r1_end

            if abs(gap - bc1_length) <= 3:
                bc1_prefix = r1_seq
                bc1_suffix = r2_seq
                bc1_start_typical = r1_end
                bc1_end_typical = r2_start
                break

        # Try reverse complement if not found
        if bc1_prefix is None:
            if verbose:
                logger.info("Trying reverse complement orientation...")

            rc_sequences = [reverse_complement(s) for s in sequences]
            matrix, coverage = build_position_conservation(rc_sequences, weights)
            regions = find_conserved_flanking_regions(
                matrix, coverage,
                min_conservation=min_conservation,
                min_length=5
            )

            regions_sorted = sorted(regions, key=lambda x: x[0])

            for i in range(len(regions_sorted) - 1):
                r1_start, r1_end, r1_seq, r1_cons = regions_sorted[i]
                r2_start, r2_end, r2_seq, r2_cons = regions_sorted[i + 1]

                gap = r2_start - r1_end

                if abs(gap - bc1_length) <= 3:
                    bc1_prefix = r1_seq
                    bc1_suffix = r2_seq
                    bc1_start_typical = r1_end
                    bc1_end_typical = r2_start
                    orientation = "reverse"
                    break

    if bc1_prefix is None:
        raise ValueError(
            f"Could not identify BC1 flanking sequences. "
            f"Try providing explicit --bc1-prefix and --bc1-suffix."
        )

    # Calculate confidence based on how consistently we found these sequences
    confidence = 0.95 if anchor_positions else 0.8

    # When RC orientation is detected, we found the RC patterns (e.g., CCTAGGG instead of CATGCTC)
    # Convert back to ORIGINAL patterns so the extraction function's RC logic works correctly
    if is_rc_orientation:
        # bc1_prefix is currently CCTAGGG (RC of CCCTAGG) - this is RC of original SUFFIX
        # bc1_suffix is currently GAGCATG (RC of CATGCTC) - this is RC of original PREFIX
        # Swap and RC to get original patterns
        original_prefix = reverse_complement(bc1_suffix) if bc1_suffix else "CATGCTC"  # GAGCATG -> CATGCTC
        original_suffix = reverse_complement(bc1_prefix)  # CCTAGGG -> CCCTAGG
        bc1_prefix = original_prefix
        bc1_suffix = original_suffix
        if verbose:
            logger.info(f"Converted RC patterns back to original: PREFIX={bc1_prefix}, SUFFIX={bc1_suffix}")

    result = {
        'bc1_prefix': bc1_prefix,
        'bc1_suffix': bc1_suffix if bc1_suffix else "UNKNOWN",
        'bc1_start': bc1_start_typical if 'bc1_start_typical' in dir() else 0,
        'bc1_end': bc1_end_typical if 'bc1_end_typical' in dir() else bc1_length,
        'orientation': orientation,
        'confidence': confidence,
    }

    if verbose:
        logger.info("=" * 60)
        logger.info("AUTO-INFERENCE RESULT")
        logger.info("=" * 60)
        logger.info(f"  BC1_PREFIX: {result['bc1_prefix']}")
        logger.info(f"  BC1_SUFFIX: {result['bc1_suffix']}")
        logger.info(f"  Orientation: {result['orientation']}")
        logger.info(f"  Confidence: {result['confidence']:.1%}")
        logger.info("=" * 60)

    return result


def validate_inferred_sequences(
    inferred: Dict[str, str],
    expected_prefix: Optional[str] = None,
    expected_suffix: Optional[str] = None,
    max_mismatches: int = 1
) -> Tuple[bool, List[str]]:
    """
    Validate inferred sequences against expected values.

    Returns:
        Tuple of (valid, list_of_warnings)
    """
    warnings = []
    valid = True

    def count_mismatches(seq1: str, seq2: str) -> int:
        if len(seq1) != len(seq2):
            return max(len(seq1), len(seq2))
        return sum(1 for a, b in zip(seq1.upper(), seq2.upper()) if a != b)

    if expected_prefix:
        mm = count_mismatches(inferred['bc1_prefix'], expected_prefix)
        if mm > max_mismatches:
            warnings.append(
                f"BC1_PREFIX mismatch: inferred {inferred['bc1_prefix']} vs "
                f"expected {expected_prefix} ({mm} mismatches)"
            )
            valid = False
        elif mm > 0:
            warnings.append(
                f"BC1_PREFIX has {mm} mismatch(es): inferred {inferred['bc1_prefix']} "
                f"vs expected {expected_prefix}"
            )

    if expected_suffix:
        mm = count_mismatches(inferred['bc1_suffix'], expected_suffix)
        if mm > max_mismatches:
            warnings.append(
                f"BC1_SUFFIX mismatch: inferred {inferred['bc1_suffix']} vs "
                f"expected {expected_suffix} ({mm} mismatches)"
            )
            valid = False
        elif mm > 0:
            warnings.append(
                f"BC1_SUFFIX has {mm} mismatch(es): inferred {inferred['bc1_suffix']} "
                f"vs expected {expected_suffix}"
            )

    return valid, warnings


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq.upper()))


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse BC1 barcodes from collapsed FASTQ",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-infer flanking sequences (recommended)
  python 02_parse_bc1_reads.py --input collapsed.fastq.gz --output bc1_counts.tsv --auto-infer

  # Auto-infer with validation against expected library sequences
  python 02_parse_bc1_reads.py --input collapsed.fastq.gz --output bc1_counts.tsv \\
      --auto-infer --expected-prefix CATGCTC --expected-suffix CCCTAGG

  # Use explicit sequences (legacy mode)
  python 02_parse_bc1_reads.py --input collapsed.fastq.gz --output bc1_counts.tsv \\
      --bc1-prefix CATGCTC --bc1-suffix CCCTAGG
        """
    )
    parser.add_argument(
        "--input", type=Path, required=True,
        help="Input collapsed FASTQ file (gzipped)"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output BC1 counts TSV"
    )
    parser.add_argument(
        "--bc1-length", type=int, default=20,
        help="BC1 barcode length (default: 20)"
    )
    parser.add_argument(
        "--length-tolerance", type=int, default=0,
        help="Accept a bc1 this many bp either side of --bc1-length. DEFAULT 0 "
             "= the exact-N20 filter. Reads rejected on length are reported in "
             "the log and the provenance sidecar, never silently dropped. Was "
             "an un-overridable +/-1 before 2026-09-19, which admitted 19 and "
             "21 bp calls into every downstream statistic "
             "(memory filter-bc1-to-exact-n20-before-any-stat)."
    )

    # Auto-inference arguments
    parser.add_argument(
        "--auto-infer", action="store_true",
        help="Auto-infer BC1_PREFIX and BC1_SUFFIX from read conservation"
    )
    parser.add_argument(
        "--inference-sample-size", type=int, default=5000,
        help="Number of reads to sample for auto-inference (default: 5000)"
    )
    parser.add_argument(
        "--min-conservation", type=float, default=0.90,
        help="Minimum conservation for flanking sequence detection (default: 0.90)"
    )

    # Explicit sequence arguments (used if --auto-infer not specified)
    parser.add_argument(
        "--bc1-prefix", type=str, default=None,
        help="Sequence immediately before BC1 (required if not using --auto-infer)"
    )
    parser.add_argument(
        "--bc1-suffix", type=str, default=None,
        help="Sequence immediately after BC1 (required if not using --auto-infer)"
    )

    # Validation arguments
    parser.add_argument(
        "--expected-prefix", type=str, default=None,
        help="Expected BC1_PREFIX for validation (optional)"
    )
    parser.add_argument(
        "--expected-suffix", type=str, default=None,
        help="Expected BC1_SUFFIX for validation (optional)"
    )
    parser.add_argument(
        "--strict-validation", action="store_true",
        help="Fail if inferred sequences don't match expected (default: warn only)"
    )

    parser.add_argument(
        "--min-counts", type=int, default=1,
        help="Minimum read count to output (default: 1)"
    )
    parser.add_argument(
        "--sample-name", type=str, default=None,
        help="Sample name (derived from filename if not provided)"
    )
    return parser.parse_args()


# BC1 length tolerance - allow 19-21bp for synthesis/sequencing errors
BC1_LENGTH_TOLERANCE = 1  # +/- 1bp from expected length


def extract_bc1_from_sequence(
    seq: str,
    bc1_prefix: str = None,
    bc1_suffix: str = None,
    bc1_length: int = 20,
    length_tolerance: int = None,
    k: int = 1,
    established_offset: int = None,
) -> tuple:
    """Extract BC1 -- DELEGATES to the canonical extractor.

    2026-09-09: the previous body was an exact ``str.find`` on the prefix (zero
    error tolerance) followed by a suffix check that degraded to no check at
    all, returning an UNLABELLED guess. Measured on real reads, with the
    shipped defaults it produced ZERO confirmed extractions on three of the four
    known layouts, and on the 2021 layout emitted 46.53% of reads as
    unvalidated barcodes. See planning/pipeline_provenance_rebuild/
    22_bc1_extraction_review.md.

    There is now ONE definition of a bc1, in magestic.utils.bc1_extraction,
    shared by the screen path, the bc1_donor_bc0 path and WGS.

    ``bc1_prefix`` / ``bc1_suffix`` / ``length_tolerance`` are accepted and
    IGNORED -- kept so existing callers do not break. The canonical anchors are
    layout-independent, so a per-era suffix is no longer a parameter: the
    universal downstream anchor's OBSERVED OFFSET is the layout call.

    Returns:
        (bc1, orientation, anchor_tier, actual_length). The third element was
        previously named ``primer_pair`` and carried "inner"/"outer"/"no_suffix"
        -- extraction-branch values in a field whose name collides with the real
        inner/outer PCR primers. It now carries the anchor tier:
        confirmed | offset_fallback | prefix_only | length_mismatch | none
    """
    call = _canonical_extract_bc1(
        seq, k=k, bc1_length=bc1_length, established_offset=established_offset
    )
    return call.bc1, call.orientation, call.anchor_tier, call.bc1_length

def _modal_offset(path, bc1_length=20, cap=20000):
    """This sample's modal universal-anchor offset, from STRICTLY confirmed reads.

    The offset IS the layout call (18=A, 1=B/D, 10=C). Establishing it enables the
    offset-constrained tiers -- without it a read that is truncated, carries one
    substitution in the anchor, or carries an indel cannot be recovered, and lands
    in prefix_only. Measured 2026-09-09: this is the difference between 97.0% and
    99.3% of reads yielding a bc1 on a real 2021 file.

    Sweep-wide (23,217 samples) the modal offset accounts for a median 97.4% of a
    sample's confirmed reads, so a per-sample modal call is safe.
    """
    from collections import Counter
    import gzip
    offs = Counter()
    op = gzip.open if str(path).endswith('.gz') else open
    n = 0
    with op(path, 'rt', errors='replace') as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            c = _canonical_extract_bc1(line.strip().upper(), k=1)
            if c.anchor_tier == 'confirmed':
                offs[c.anchor_offset] += 1
            n += 1
            if n >= cap:
                break
    if not offs:
        return None, 0.0
    off, hits = offs.most_common(1)[0]
    return off, hits / sum(offs.values())


def parse_collapsed_fastq(
    fastq_path: Path,
    bc1_prefix: str,
    bc1_suffix: str,
    bc1_length: int,
    length_tolerance: int = BC1_LENGTH_TOLERANCE,
    established_offset: int = None
) -> tuple:
    """
    Parse collapsed FASTQ and extract BC1 counts.

    Returns:
        A dict. ``n_records`` is the number of collapsed RECORDS read;
        ``total_reads`` is the number of READS they represent. The two differ by
        the collapse factor and conflating them was the defect fixed on
        2026-09-19 (planning/210_bc1_recall_code_gate.md).

        bc1_counts        {bc1: reads}
        bc1_tiers         {bc1: Counter(anchor_tier -> reads)}
        bc1_orientations  {bc1: Counter(orientation -> reads)}
        total_reads, matched_reads, n_records
        length_stats      Counter(accepted bc1 length -> reads)
        reject_lengths    Counter(REJECTED bc1 length -> reads) -- never silently
                          dropped, per handoff §4.1
        n_rejected_n_base reads whose bc1 contained an N
        tier_stats        Counter(anchor_tier -> reads), all reads
        header_formats    Counter(read-name format -> records)
    """
    bc1_counts = defaultdict(int)
    bc1_tiers = defaultdict(Counter)
    bc1_orientations = defaultdict(Counter)
    total_reads = 0
    matched_reads = 0
    length_stats = Counter()
    reject_lengths = Counter()
    n_rejected_n_base = 0
    tier_stats = Counter()
    header_formats = Counter()
    n_records = 0

    min_bc1_len = bc1_length - length_tolerance
    max_bc1_len = bc1_length + length_tolerance

    # Read counts come from magestic.utils.collapsed_fastq, which resolves BOTH
    # producer formats per record and RAISES on an unrecognized read name. The
    # previous inline block recognized only ``count=`` and silently weighted
    # every legacy ``@12`` record as one read -- see
    # planning/210_bc1_recall_code_gate.md for the measurement.
    for sequence, count, _quality, _fmt in iter_collapsed_fastq(
            fastq_path, format_counter=header_formats):
        n_records += 1
        total_reads += count

        # Extract BC1 (now returns 4 values)
        bc1, orientation, anchor_tier, actual_length = extract_bc1_from_sequence(
            sequence, bc1_prefix, bc1_suffix, bc1_length, length_tolerance,
            established_offset=established_offset
        )
        tier_stats[anchor_tier] += count

        if not bc1:
            continue
        if "N" in bc1:
            n_rejected_n_base += count
            continue
        if min_bc1_len <= actual_length <= max_bc1_len:
            bc1_counts[bc1] += count
            bc1_tiers[bc1][anchor_tier] += count
            bc1_orientations[bc1][orientation] += count
            matched_reads += count
            length_stats[actual_length] += count
        else:
            reject_lengths[actual_length] += count

    return dict(
        bc1_counts=bc1_counts, bc1_tiers=bc1_tiers,
        bc1_orientations=bc1_orientations, total_reads=total_reads,
        matched_reads=matched_reads, n_records=n_records,
        length_stats=length_stats, reject_lengths=reject_lengths,
        n_rejected_n_base=n_rejected_n_base, tier_stats=tier_stats,
        header_formats=header_formats,
    )


def main():
    args = parse_args()

    # Derive sample name from filename if not provided
    sample_name = args.sample_name
    if sample_name is None:
        sample_name = args.input.stem.replace("_collapsed", "").replace(".fastq", "")

    logger.info(f"Step 02: Parsing BC1 from {args.input}")
    logger.info(f"  Sample: {sample_name}")
    logger.info(f"  BC1 length: {args.bc1_length}")

    # Check input exists
    if not args.input.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1

    # Determine BC1 flanking sequences
    if args.auto_infer:
        # Auto-inference mode
        logger.info("=" * 60)
        logger.info("AUTO-INFERENCE MODE")
        logger.info("=" * 60)

        try:
            inferred = infer_bc1_flanking_sequences(
                args.input,
                bc1_length=args.bc1_length,
                n_sample=args.inference_sample_size,
                min_conservation=args.min_conservation,
                verbose=True
            )

            bc1_prefix = inferred['bc1_prefix']
            bc1_suffix = inferred['bc1_suffix']

            # Validate if expected sequences provided
            if args.expected_prefix or args.expected_suffix:
                valid, warnings = validate_inferred_sequences(
                    inferred,
                    expected_prefix=args.expected_prefix,
                    expected_suffix=args.expected_suffix
                )

                for w in warnings:
                    if valid or not args.strict_validation:
                        logger.warning(w)
                    else:
                        logger.error(w)

                if not valid and args.strict_validation:
                    logger.error("Validation failed with --strict-validation. Exiting.")
                    return 1

        except Exception as e:
            # A FAILED INFERENCE IS NOT FATAL, and since 2026-09-09 it is barely
            # even relevant: the inferred prefix/suffix are accepted and IGNORED,
            # because the canonical extractor calls the layout per read from the
            # universal downstream anchor's observed offset. Returning 1 here
            # aborted the run BEFORE the real extractor ever executed.
            #
            # Measured in Phase A of the re-call (2026-09-19): 44 of 50 unit
            # failures came through this path -- 35 on zero-byte collapsed files
            # and 9 on near-empty wells holding a handful of off-target reads.
            # Every one of those is an empty-well RESULT that belongs in the
            # output as an empty table with a sidecar, not a task failure.
            logger.warning(
                "Auto-inference could not determine flanking sequences (%s). "
                "Continuing anyway: the inferred values are IGNORED by the "
                "canonical extractor, so this does not affect the call. The "
                "anchor-tier counts below are the evidence to read.", e)
            bc1_prefix = args.bc1_prefix or "CATGCTC"
            bc1_suffix = args.bc1_suffix or "CCCTAGG"
    else:
        # Explicit sequences mode
        if not args.bc1_prefix or not args.bc1_suffix:
            # Use defaults if not provided
            bc1_prefix = args.bc1_prefix or "CATGCTC"
            bc1_suffix = args.bc1_suffix or "CCCTAGG"
            logger.warning(f"Using default BC1_PREFIX={bc1_prefix} BC1_SUFFIX={bc1_suffix}")
            
        else:
            bc1_prefix = args.bc1_prefix
            bc1_suffix = args.bc1_suffix

    logger.info(
        "  Anchors: CANONICAL (magestic.utils.bc1_extraction). "
        "--bc1-prefix/--bc1-suffix are accepted but IGNORED since 2026-09-09: "
        "the layout is identified from the universal downstream anchor's "
        "observed offset, not configured per era."
    )

    # Parse FASTQ
    logger.info("Parsing collapsed FASTQ...")
    # PASS 1 -- establish this sample's layout from its own reads.
    est_off, off_purity = _modal_offset(args.input, args.bc1_length)
    layout = {18: "A_1.0/2.5", 1: "B_or_D_no_spacer", 10: "C_3.0"}.get(est_off)
    if est_off is None:
        logger.warning("  No universal anchor found in this sample; offset-constrained "
                       "recovery is DISABLED and reads will fall to prefix_only.")
    else:
        logger.info(f"  Layout call: offset +{est_off}"
                    f"{' = ' + layout if layout else ' (UNRECOGNIZED offset)'}, "
                    f"purity {off_purity:.4f}")
        if off_purity < 0.90:
            logger.warning(f"  Offset purity {off_purity:.4f} < 0.90 -- this sample is "
                           f"MIXED. The layout call is least trustworthy here.")

    # PASS 2 -- extract with that offset established.
    stats = parse_collapsed_fastq(
        args.input,
        bc1_prefix,
        bc1_suffix,
        args.bc1_length,
        length_tolerance=args.length_tolerance,
        established_offset=est_off
    )
    bc1_counts = stats["bc1_counts"]
    total_reads = stats["total_reads"]
    matched_reads = stats["matched_reads"]
    length_stats = stats["length_stats"]
    tier_stats = stats["tier_stats"]
    n_records = stats["n_records"]
    header_formats = stats["header_formats"]

    match_rate = 100 * matched_reads / total_reads if total_reads > 0 else 0
    # READS and RECORDS are different quantities. Log both, always: a run whose
    # collapse factor is 1.00 on a legacy tree is the signature of the count
    # defect fixed 2026-09-19 (planning/210_bc1_recall_code_gate.md).
    collapse_factor = total_reads / n_records if n_records else 0
    logger.info(f"  Collapsed records: {n_records:,}")
    logger.info(f"  Total reads: {total_reads:,} "
                f"(collapse factor {collapse_factor:.2f} reads/record)")
    logger.info(f"  Read-name formats: "
                f"{dict(header_formats.most_common())}")
    logger.info(f"  Matched reads: {matched_reads:,} ({match_rate:.1f}%)")
    logger.info(f"  Unique BC1s: {len(bc1_counts):,}")

    # Report BC1 length statistics
    if length_stats:
        logger.info("  BC1 length distribution:")
        for length in sorted(length_stats.keys()):
            count = length_stats[length]
            pct = 100 * count / matched_reads if matched_reads > 0 else 0
            marker = " <-- expected" if length == args.bc1_length else ""
            logger.info(f"    {length}bp: {count:,} ({pct:.1f}%){marker}")

    # Anchor-tier breakdown. The old code discarded this entirely: every read
    # was counted the same whether both anchors were verified or neither was.
    if tier_stats:
        logger.info("  Anchor tier distribution (read-weighted):")
        for tier in ANCHOR_TIERS:
            n = tier_stats.get(tier, 0)
            if n:
                logger.info(f"    {tier}: {n:,} ({100*n/total_reads:.2f}%)")
        accounted = sum(tier_stats.get(x, 0) for x in ANCHOR_TIERS)
        if accounted != total_reads:
            logger.error(
                f"  TIER ACCOUNTING BROKEN: {accounted:,} of {total_reads:,} reads "
                f"land in a known tier. Unknown tiers: "
                f"{sorted(set(tier_stats) - set(ANCHOR_TIERS))}. "
                f"ANCHOR_TIERS is the single source of truth -- add the tier there."
            )
        unconfirmed = tier_stats.get("prefix_only", 0)
        if unconfirmed and total_reads and unconfirmed / total_reads > 0.05:
            logger.warning(
                f"  {100*unconfirmed/total_reads:.1f}% of reads yielded a bc1 with NO "
                f"downstream anchor confirmation (tier=prefix_only). This is the "
                f"failure mode that made the pre-2026-09-09 extractor emit unvalidated "
                f"barcodes on 2021-era layouts -- check the layout call."
            )

    # Reads whose bc1 was the WRONG LENGTH are reported, never silently dropped
    # (handoff §4.1 / memory filter-bc1-to-exact-n20-before-any-stat). With the
    # default --length-tolerance 0 this is the exact-N20 filter.
    reject_lengths = stats["reject_lengths"]
    n_rejected_length = sum(reject_lengths.values())
    logger.info(f"  Length filter: tolerance +/-{args.length_tolerance} around "
                f"{args.bc1_length} bp")
    if n_rejected_length:
        logger.info(f"  Rejected on length: {n_rejected_length:,} reads "
                    f"({100*n_rejected_length/total_reads:.3f}%) -- distribution:")
        for length in sorted(reject_lengths):
            n = reject_lengths[length]
            logger.info(f"    {length}bp: {n:,} ({100*n/total_reads:.3f}%)")
    else:
        logger.info("  Rejected on length: 0 reads")
    if stats["n_rejected_n_base"]:
        logger.info(f"  Rejected for an N base: {stats['n_rejected_n_base']:,} reads")

    # Singleton tail and reads-per-lineage, per sample. Handoff §4.2: richness
    # moves the WRONG WAY across an outgrowth without an error-aware view of
    # this tail, so it is reported on every sample rather than derived later.
    n_lineages = len(bc1_counts)
    n_singletons = sum(1 for v in bc1_counts.values() if v == 1)
    singleton_fraction = n_singletons / n_lineages if n_lineages else 0
    reads_per_lineage = matched_reads / n_lineages if n_lineages else 0
    logger.info(f"  Singletons: {n_singletons:,} of {n_lineages:,} lineages "
                f"({100*singleton_fraction:.1f}%)")
    logger.info(f"  Reads per lineage: {reads_per_lineage:.2f}")

    # Filter by minimum counts
    n_before_min_counts = len(bc1_counts)
    if args.min_counts > 1:
        bc1_counts = {k: v for k, v in bc1_counts.items() if v >= args.min_counts}
        logger.info(f"  After min_counts filter: {len(bc1_counts):,}")

    # Create output DataFrame. The anchor tier, the orientation and the layout
    # call were computed and then DISCARDED before 2026-09-19, which is why the
    # legacy count tables cannot be compared to each other -- a confirmed read
    # and an unconfirmed prefix_only guess were the same row. All seven emitting
    # tiers are written always, even when zero, so files stay joinable.
    import pandas as pd

    bc1_tiers = stats["bc1_tiers"]
    bc1_orientations = stats["bc1_orientations"]
    # The schema is declared, not inferred from the rows. A sample that yields
    # ZERO barcodes is a real result -- an empty or failed well -- and it must
    # still produce a valid table with the right header, so downstream joins and
    # the counts==sum(n_<tier>) check see a schema rather than a KeyError.
    # pd.DataFrame([]) has no columns at all, which crashed sort_values("counts")
    # on 4 genuinely empty 2021 samples during the first Phase A run.
    OUT_COLUMNS = (["bc1", "sample_name", "counts"]
                   + [f"n_{t}" for t in EMITTING_TIERS] + ["n_fwd", "n_rev"])
    rows = []
    for bc1, count in bc1_counts.items():
        tiers = bc1_tiers[bc1]
        orients = bc1_orientations[bc1]
        row = {"bc1": bc1, "sample_name": sample_name, "counts": count}
        for tier in EMITTING_TIERS:
            row[f"n_{tier}"] = tiers.get(tier, 0)
        # The orientation vocabulary is the module's: fwd | rev (Bc1Call
        # .orientation, bc1_extraction.py:107). Spelling it "forward"/"reverse"
        # here would make both columns silently zero.
        row["n_fwd"] = orients.get("fwd", 0)
        row["n_rev"] = orients.get("rev", 0)
        rows.append(row)
    output_df = pd.DataFrame(rows, columns=OUT_COLUMNS)
    if not len(output_df):
        logger.warning(
            "  ZERO barcodes passed extraction for this sample. Writing an EMPTY "
            "table with the full header plus its provenance sidecar -- an empty "
            "well is a result, not an error, and the sidecar records why "
            "(total_reads=%d, matched_reads=%d).", total_reads, matched_reads)

    # Sort by counts descending
    output_df = output_df.sort_values("counts", ascending=False)

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    logger.info(f"Saving {len(output_df):,} BC1s to {args.output}")
    output_df.to_csv(args.output, sep="\t", index=False)

    # counts == sum of the per-tier columns AND == n_fwd + n_rev, on every row.
    # If either fails the table is not internally consistent and nothing
    # downstream should trust it. The orientation check is what would have caught
    # the "forward"/"fwd" vocabulary slip that made both columns read zero.
    if len(output_df):
        tier_sum = output_df[[f"n_{t}" for t in EMITTING_TIERS]].sum(axis=1)
        n_bad = int((tier_sum != output_df["counts"]).sum())
        if n_bad:
            logger.error(f"  TIER COLUMN ACCOUNTING BROKEN on {n_bad:,} rows: "
                         f"counts != sum(n_<tier>).")
        orient_sum = output_df["n_fwd"] + output_df["n_rev"]
        n_bad_or = int((orient_sum != output_df["counts"]).sum())
        if n_bad_or:
            logger.error(f"  ORIENTATION COLUMN ACCOUNTING BROKEN on "
                         f"{n_bad_or:,} rows: counts != n_fwd + n_rev.")
        logger.info(f"  Orientation: fwd {int(output_df['n_fwd'].sum()):,} / "
                    f"rev {int(output_df['n_rev'].sum()):,} reads")

    # Provenance sidecar (planning/policies/PIPELINE_CORRECTNESS_POLICY.md §2).
    # Without this the outputs of two runs cannot be told apart, which is exactly
    # why the legacy counts are not comparable (handoff §4.5).
    try:
        from magestic.utils.provenance import StepProvenance
        prov = StepProvenance(
            pipeline="bc1_from_screen", step="02_parse_bc1_reads",
            output=args.output)
        prov.add_input(args.input, role="collapsed_fastq")
        prov.set_parameters(
            sample_name=sample_name,
            bc1_length=args.bc1_length,
            length_tolerance=args.length_tolerance,
            min_counts=args.min_counts,
            auto_infer=args.auto_infer,
            established_offset=est_off,
            layout_call=layout if est_off is not None else None,
            offset_purity=round(off_purity, 4) if est_off is not None else None,
            anchors="canonical (magestic.utils.bc1_extraction)",
            reported_inference_prefix=bc1_prefix,
            reported_inference_suffix=bc1_suffix,
        )
        prov.note(
            "--bc1-prefix/--bc1-suffix and any --auto-infer result are ACCEPTED "
            "AND IGNORED since 2026-09-09: the layout is called per read from "
            "the universal downstream anchor's observed offset. The inference "
            "values recorded here are what the log printed, NOT what ran."
        )
        prov.set_counts(
            collapsed_records=n_records,
            total_reads=total_reads,
            collapse_factor=round(collapse_factor, 4),
            matched_reads=matched_reads,
            rejected_on_length=n_rejected_length,
            rejected_length_distribution={str(k): v for k, v in
                                          sorted(reject_lengths.items())},
            rejected_n_base=stats["n_rejected_n_base"],
            read_name_formats=dict(header_formats),
            anchor_tiers={t: tier_stats.get(t, 0) for t in ANCHOR_TIERS},
            lineages=n_lineages,
            lineages_before_min_counts=n_before_min_counts,
            lineages_written=len(output_df),
            singletons=n_singletons,
            singleton_fraction=round(singleton_fraction, 6),
            reads_per_lineage=round(reads_per_lineage, 4),
        )
        prov.add_output(args.output, rows=len(output_df))
        sidecar = prov.write()
        logger.info(f"  Provenance sidecar: {sidecar}")
    except Exception as e:                       # retain data, refuse completion
        logger.error(f"  FAILED to write the provenance sidecar: {e!r}. The "
                     f"counts were written but are NOT provenance-backed.")
        return 1

    logger.info("Step 02: Parsing complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
