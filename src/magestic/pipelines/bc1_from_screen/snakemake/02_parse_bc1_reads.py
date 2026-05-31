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

    Collapsed FASTQs have read counts encoded in the read name.

    Returns:
        List of (sequence, count) tuples
    """
    reads = []
    open_fn = gzip.open if str(fastq_path).endswith(".gz") else open

    with open_fn(fastq_path, "rt") as f:
        while len(reads) < n_sample:
            header = f.readline().strip()
            if not header:
                break

            sequence = f.readline().strip()
            f.readline()  # plus line
            f.readline()  # quality line

            # Parse count from header
            count = 1
            if "count=" in header:
                try:
                    count = int(header.split("count=")[1].split()[0])
                except (ValueError, IndexError):
                    pass
            elif "-" in header:
                try:
                    count = int(header.split("-")[-1])
                except ValueError:
                    pass

            reads.append((sequence, count))

    return reads


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
    bc1_prefix: str,
    bc1_suffix: str,
    bc1_length: int,
    length_tolerance: int = BC1_LENGTH_TOLERANCE,
) -> tuple:
    """
    Extract BC1 barcode from a sequence.

    Supports BC1 lengths of +/- length_tolerance from expected.
    For example, with bc1_length=20 and tolerance=1, accepts 19-21bp BC1s.

    Returns:
        (bc1, orientation, primer_pair, actual_length) or (None, None, None, None) if not found
    """
    seq = seq.upper()
    bc1_prefix = bc1_prefix.upper()
    bc1_suffix = bc1_suffix.upper()

    min_bc1_len = bc1_length - length_tolerance
    max_bc1_len = bc1_length + length_tolerance

    # Try forward orientation
    prefix_idx = seq.find(bc1_prefix)
    if prefix_idx >= 0:
        bc1_start = prefix_idx + len(bc1_prefix)

        # PASS 1: Look for EXACT suffix match at any valid BC1 length
        # This prevents fuzzy match at 19bp from triggering when 20bp is the real match
        for try_len in range(min_bc1_len, max_bc1_len + 1):
            bc1_end = bc1_start + try_len
            if bc1_end > len(seq):
                continue

            suffix_start = bc1_end
            suffix_region = seq[suffix_start:suffix_start + len(bc1_suffix)]

            if suffix_region == bc1_suffix:
                bc1 = seq[bc1_start:bc1_end]
                return bc1, "forward", "inner", len(bc1)

        # PASS 2: If no exact match, try fuzzy match (partial suffix) at expected length only
        # This avoids false positives from finding suffix at wrong offset
        bc1_end = bc1_start + bc1_length
        if bc1_end <= len(seq):
            suffix_start = bc1_end
            # Check if first 3 chars of suffix appear at exactly the expected position
            suffix_region = seq[suffix_start:suffix_start + len(bc1_suffix)]
            if suffix_region.startswith(bc1_suffix[:3]):
                bc1 = seq[bc1_start:bc1_end]
                return bc1, "forward", "outer", len(bc1)

        # PASS 3: If still no match, return expected length without suffix validation
        if bc1_end <= len(seq):
            bc1 = seq[bc1_start:bc1_end]
            return bc1, "forward", "no_suffix", len(bc1)

    # Try reverse complement orientation
    rc_prefix = reverse_complement(bc1_suffix)
    rc_suffix = reverse_complement(bc1_prefix)

    prefix_idx = seq.find(rc_prefix)
    if prefix_idx >= 0:
        bc1_start = prefix_idx + len(rc_prefix)

        # PASS 1: Look for EXACT suffix match at any valid BC1 length
        for try_len in range(min_bc1_len, max_bc1_len + 1):
            bc1_end = bc1_start + try_len
            if bc1_end > len(seq):
                continue

            suffix_start = bc1_end
            suffix_region = seq[suffix_start:suffix_start + len(rc_suffix)]

            if suffix_region == rc_suffix:
                bc1_rc = seq[bc1_start:bc1_end]
                bc1 = reverse_complement(bc1_rc)
                return bc1, "reverse", "inner", len(bc1)

        # PASS 2: Try fuzzy match at expected length only
        bc1_end = bc1_start + bc1_length
        if bc1_end <= len(seq):
            suffix_start = bc1_end
            suffix_region = seq[suffix_start:suffix_start + len(rc_suffix)]
            if suffix_region.startswith(rc_suffix[:3]):
                bc1_rc = seq[bc1_start:bc1_end]
                bc1 = reverse_complement(bc1_rc)
                return bc1, "reverse", "outer", len(bc1)

        # PASS 3: Return expected length without suffix validation
        if bc1_end <= len(seq):
            bc1_rc = seq[bc1_start:bc1_end]
            bc1 = reverse_complement(bc1_rc)
            return bc1, "reverse", "no_suffix", len(bc1)

    return None, None, None, None


def parse_collapsed_fastq(
    fastq_path: Path,
    bc1_prefix: str,
    bc1_suffix: str,
    bc1_length: int,
    length_tolerance: int = BC1_LENGTH_TOLERANCE,
) -> tuple:
    """
    Parse collapsed FASTQ and extract BC1 counts.

    Returns:
        Tuple of (bc1_counts, total_reads, matched_reads, length_stats)
        where length_stats is a Counter of BC1 lengths
    """
    bc1_counts = defaultdict(int)
    total_reads = 0
    matched_reads = 0
    length_stats = Counter()

    min_bc1_len = bc1_length - length_tolerance
    max_bc1_len = bc1_length + length_tolerance

    open_fn = gzip.open if str(fastq_path).endswith(".gz") else open

    with open_fn(fastq_path, "rt") as f:
        while True:
            # Read 4 lines (FASTQ format)
            header = f.readline().strip()
            if not header:
                break

            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            # Extract count from collapsed FASTQ header (format: @seq_123_count=456)
            count = 1
            if "count=" in header:
                try:
                    count = int(header.split("count=")[1].split()[0])
                except (ValueError, IndexError):
                    pass
            elif "-" in header:
                # Alternative format: @seq-456 where 456 is count
                try:
                    count = int(header.split("-")[-1])
                except ValueError:
                    pass

            total_reads += count

            # Extract BC1 (now returns 4 values)
            bc1, orientation, primer_pair, actual_length = extract_bc1_from_sequence(
                sequence, bc1_prefix, bc1_suffix, bc1_length, length_tolerance
            )

            if bc1 and "N" not in bc1:
                # Validate length is within tolerance
                if min_bc1_len <= actual_length <= max_bc1_len:
                    bc1_counts[bc1] += count
                    matched_reads += count
                    length_stats[actual_length] += count

    return bc1_counts, total_reads, matched_reads, length_stats


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
            logger.error(f"Auto-inference failed: {e}")
            if args.bc1_prefix and args.bc1_suffix:
                logger.warning("Falling back to provided --bc1-prefix and --bc1-suffix")
                bc1_prefix = args.bc1_prefix
                bc1_suffix = args.bc1_suffix
            else:
                return 1
    else:
        # Explicit sequences mode
        if not args.bc1_prefix or not args.bc1_suffix:
            # Use defaults if not provided
            bc1_prefix = args.bc1_prefix or "CATGCTC"
            bc1_suffix = args.bc1_suffix or "CCCTAGG"
            logger.warning(f"Using default BC1_PREFIX={bc1_prefix} BC1_SUFFIX={bc1_suffix}")
            logger.warning("Consider using --auto-infer for library-specific detection")
        else:
            bc1_prefix = args.bc1_prefix
            bc1_suffix = args.bc1_suffix

    logger.info(f"  BC1 prefix: {bc1_prefix}")
    logger.info(f"  BC1 suffix: {bc1_suffix}")

    # Parse FASTQ
    logger.info("Parsing collapsed FASTQ...")
    bc1_counts, total_reads, matched_reads, length_stats = parse_collapsed_fastq(
        args.input,
        bc1_prefix,
        bc1_suffix,
        args.bc1_length
    )

    match_rate = 100 * matched_reads / total_reads if total_reads > 0 else 0
    logger.info(f"  Total reads: {total_reads:,}")
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

    # Filter by minimum counts
    if args.min_counts > 1:
        bc1_counts = {k: v for k, v in bc1_counts.items() if v >= args.min_counts}
        logger.info(f"  After min_counts filter: {len(bc1_counts):,}")

    # Create output DataFrame
    import pandas as pd

    output_df = pd.DataFrame([
        {"bc1": bc1, "sample_name": sample_name, "counts": count}
        for bc1, count in bc1_counts.items()
    ])

    # Sort by counts descending
    output_df = output_df.sort_values("counts", ascending=False)

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    logger.info(f"Saving {len(output_df):,} BC1s to {args.output}")
    output_df.to_csv(args.output, sep="\t", index=False)

    logger.info("Step 02: Parsing complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
