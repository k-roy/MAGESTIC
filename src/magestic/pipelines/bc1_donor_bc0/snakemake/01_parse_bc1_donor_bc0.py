#!/usr/bin/env python3
"""
Parse bc1-donor-bc0 sequences from paired FASTQ files with arbitrary read lengths.

This script handles any R1/R2 length configuration by:
1. Auto-detecting read lengths from the first N reads
2. Determining optimal extraction strategy based on amplicon structure
3. Using anchor-based extraction for robustness

Supported configurations:
- 2x300 bp: Full donor captured in merged reads
- 2x150 bp: Partial donor in R2, full donor via separate donor-bc0 amplicon
- R1=228, R2=90: Full donor in R1, bc0/SPS/msd from R2
- Any other configuration: Dynamically adjusts extraction

Amplicon structure (bc1-donor-bc0):
    5' [PS1][bc1_prefix][bc1:20bp][bc1_suffix][SaCas9_SceI][donor:~129bp][SPS:20bp][bc0:10bp][msd] 3'

R2 is reverse complement, reading from 3' end:
    R2: [msd_rc][bc0_rc][SPS_rc][partial_donor_rc]...

Author: Kevin R. Roy
Date: 2026-03-21
"""

import argparse
import sys
import os
import gzip
import timeit
from collections import Counter
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Try to import pysam and tqdm
try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False
    logger.warning("pysam not available, using fallback FASTQ parser")

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable


# =============================================================================
# Constants
# =============================================================================

MAX_READS = int(10e10)

# Known BC1 prefixes
KNOWN_BC1_PREFIXES = [
    "CAGGTCATGCTC",  # 12bp - common in 2x150 data
    "GGTCATGCTC",    # 10bp
    "CATGCTC",       # 7bp - common in 2x300 collapsed data
]

# BC1 constants
DEFAULT_BC1_PREFIX = "CATGCTC"
DEFAULT_BC1_SUFFIX = "CCCTAGGATACGTGC"
DEFAULT_BC1_LENGTH = 20
BC1_LENGTH_TOLERANCE = 1

# SaCas9/SceI fragment
SACAS9_SCEI_FRAGMENT = "CCCTAGGATACGTGC"

# MSD sequence (for anchor-based extraction)
MSD_SEQUENCE = "AGGAAAACAGAC"

# Expected component lengths
SPS_LENGTH = 20
BC0_LENGTH = 10
MSD_LENGTH = 12
DONOR_BC0_FRAGMENT_LENGTH = 60  # 30bp partial donor + 20bp SPS + 10bp bc0


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class ReadConfig:
    """Configuration derived from read length analysis."""
    r1_length: int
    r2_length: int
    bc1_prefix: str = DEFAULT_BC1_PREFIX
    bc1_suffix: str = DEFAULT_BC1_SUFFIX
    bc1_length: int = DEFAULT_BC1_LENGTH

    # Primer orientation: "standard" (R1=BC1, R2=msd) or "reversed" (R2=BC1, R1=msd)
    orientation: str = "standard"

    # Extraction strategy
    can_extract_full_donor_from_bc1_read: bool = False
    can_extract_bc0_sps_from_msd_read: bool = False

    def __post_init__(self):
        """Determine extraction strategy based on read lengths and orientation."""
        # Determine which read has BC1 vs msd based on orientation
        if self.orientation == "reversed":
            bc1_read_length = self.r2_length
            msd_read_length = self.r1_length
        else:  # standard
            bc1_read_length = self.r1_length
            msd_read_length = self.r2_length

        # BC1 read structure: bc1_prefix + bc1 + bc1_suffix + ... (optionally more)
        # Minimum to extract BC1: ~12 + 20 + 15 = 47bp
        # For full donor from BC1 read: ~12 + 20 + 15 + 15 + 129 + 20 = 211bp
        min_bc1_read_for_full_donor = 200

        # MSD read structure (RC): msd + bc0 + SPS + partial_donor
        # Minimum to extract bc0 and SPS: 12 + 10 + 20 = 42bp
        min_msd_read_for_bc0_sps = 42
        # For full donor from msd read (RC): ~12 + 10 + 20 + 129 = 171bp
        min_msd_read_for_full_donor = 170

        self.can_extract_bc0_sps_from_msd_read = msd_read_length >= min_msd_read_for_bc0_sps
        # Can extract donor from either read depending on length
        self.can_extract_full_donor_from_bc1_read = bc1_read_length >= min_bc1_read_for_full_donor
        self.can_extract_full_donor_from_msd_read = msd_read_length >= min_msd_read_for_full_donor

    @property
    def bc1_read(self) -> str:
        """Which read (R1 or R2) contains BC1."""
        return "R2" if self.orientation == "reversed" else "R1"

    @property
    def msd_read(self) -> str:
        """Which read (R1 or R2) contains msd/bc0/SPS."""
        return "R1" if self.orientation == "reversed" else "R2"

    @property
    def strategy(self) -> str:
        """Return extraction strategy description."""
        parts = []
        parts.append(f"orientation:{self.orientation}")
        parts.append(f"bc1_from_{self.bc1_read}")
        if self.can_extract_full_donor_from_msd_read:
            parts.append(f"donor_from_{self.msd_read}")
        elif self.can_extract_full_donor_from_bc1_read:
            parts.append(f"donor_from_{self.bc1_read}")
        else:
            parts.append("partial_donor")
        return "_".join(parts)


# =============================================================================
# Utility Functions
# =============================================================================

def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


def read_fastq_entries(filepath: str, n_entries: int = None):
    """Read FASTQ entries from a file (gzipped or not)."""
    open_fn = gzip.open if filepath.endswith('.gz') else open
    count = 0

    with open_fn(filepath, 'rt') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()

            yield (header, sequence, plus, quality)
            count += 1

            if n_entries and count >= n_entries:
                break


def detect_read_lengths(r1_path: str, r2_path: str, n_sample: int = 1000) -> Tuple[int, int]:
    """Detect R1 and R2 read lengths from FASTQ files.

    Returns (0, 0) if files are empty.
    """
    r1_lengths = []
    r2_lengths = []

    for i, (_, seq, _, _) in enumerate(read_fastq_entries(r1_path, n_sample)):
        r1_lengths.append(len(seq))
        if i >= n_sample - 1:
            break

    for i, (_, seq, _, _) in enumerate(read_fastq_entries(r2_path, n_sample)):
        r2_lengths.append(len(seq))
        if i >= n_sample - 1:
            break

    # Handle empty files
    if not r1_lengths or not r2_lengths:
        return 0, 0

    r1_median = sorted(r1_lengths)[len(r1_lengths)//2]
    r2_median = sorted(r2_lengths)[len(r2_lengths)//2]

    return r1_median, r2_median


def find_bc1_prefix(sequences: List[str], known_prefixes: List[str] = None) -> Optional[str]:
    """Find the BC1_PREFIX in sequences."""
    if known_prefixes is None:
        known_prefixes = KNOWN_BC1_PREFIXES

    best_prefix = None
    best_count = 0

    for prefix in known_prefixes:
        count = sum(1 for seq in sequences if prefix in seq[:50])
        if count > best_count:
            best_count = count
            best_prefix = prefix

    return best_prefix if best_count > len(sequences) * 0.3 else None


# MSD and its reverse complement for orientation detection
MSD_RC = rev_comp(MSD_SEQUENCE)  # GTCTGTTTTCCT


def detect_orientation(r1_seqs: List[str], r2_seqs: List[str], n_sample: int = 1000) -> str:
    """
    Detect primer orientation by checking which read contains BC1_PREFIX vs MSD_RC.

    Returns:
        "standard" if R1 has BC1_PREFIX (R1=BC1, R2=msd)
        "reversed" if R2 has BC1_PREFIX (R2=BC1, R1=msd)
    """
    # Check for BC1_PREFIX in first 50bp of each read
    r1_bc1_count = sum(1 for seq in r1_seqs[:n_sample]
                       if any(prefix in seq[:50] for prefix in KNOWN_BC1_PREFIXES))
    r2_bc1_count = sum(1 for seq in r2_seqs[:n_sample]
                       if any(prefix in seq[:50] for prefix in KNOWN_BC1_PREFIXES))

    # Check for MSD_RC in first 50bp of each read
    r1_msd_count = sum(1 for seq in r1_seqs[:n_sample] if MSD_RC in seq[:50])
    r2_msd_count = sum(1 for seq in r2_seqs[:n_sample] if MSD_RC in seq[:50])

    logger.info(f"Orientation detection:")
    logger.info(f"  R1 BC1_PREFIX count: {r1_bc1_count}/{n_sample}")
    logger.info(f"  R2 BC1_PREFIX count: {r2_bc1_count}/{n_sample}")
    logger.info(f"  R1 MSD_RC count: {r1_msd_count}/{n_sample}")
    logger.info(f"  R2 MSD_RC count: {r2_msd_count}/{n_sample}")

    # Reversed orientation: R2 has BC1, R1 has MSD_RC
    if r2_bc1_count > r1_bc1_count and r1_msd_count > r2_msd_count:
        return "reversed"
    # Standard orientation: R1 has BC1, R2 has MSD_RC
    elif r1_bc1_count > r2_bc1_count:
        return "standard"
    # If MSD_RC is in R1 but no clear BC1 pattern, assume reversed
    elif r1_msd_count > r2_msd_count:
        return "reversed"
    else:
        return "standard"  # Default to standard


# =============================================================================
# R1 Processing (BC1 and Donor extraction)
# =============================================================================

def extract_bc1_from_r1(seq: str, config: ReadConfig) -> Tuple[str, str, int, bool]:
    """
    Extract BC1 from R1 sequence.

    Returns:
        Tuple of (bc1, PS1, bc1_length, prefix_found)
    """
    bc1_prefix = config.bc1_prefix
    bc1_suffix = config.bc1_suffix
    bc1_length = config.bc1_length

    prefix_idx = seq.find(bc1_prefix)

    if prefix_idx >= 0:
        bc1_start = prefix_idx + len(bc1_prefix)
        PS1 = seq[:prefix_idx + len(bc1_prefix)]

        # Try to find exact suffix match at each valid BC1 length
        for try_len in range(bc1_length - BC1_LENGTH_TOLERANCE, bc1_length + BC1_LENGTH_TOLERANCE + 1):
            bc1_end = bc1_start + try_len
            if bc1_end > len(seq):
                continue

            suffix_region = seq[bc1_end:bc1_end + len(bc1_suffix)]
            if suffix_region == bc1_suffix:
                bc1 = seq[bc1_start:bc1_end]
                return bc1, PS1, len(bc1), True

        # Fall back to expected length
        bc1 = seq[bc1_start:bc1_start + bc1_length]
        return bc1, PS1, bc1_length, True
    else:
        # Default position if prefix not found
        bc1 = seq[16:16 + bc1_length]
        PS1 = seq[:16]
        return bc1, PS1, bc1_length, False


def extract_donor_from_r1(seq: str, config: ReadConfig) -> Optional[Tuple[str, str]]:
    """
    Extract full donor and SPS from R1 if read is long enough.

    Uses anchor-based extraction: find MSD from end, work backwards.

    Returns:
        Tuple of (donor, SPS) or None if extraction fails
    """
    # Find the CTGCTCTTCAAAC pattern (part of MAGESTIC 3.0 sequence near 3' end)
    magestic_3_marker = "CTGCTCTTCAAAC"
    marker_pos = seq.rfind(magestic_3_marker)

    if marker_pos < 0:
        return None

    # The SPS ends just before this marker
    # Structure: ...donor][SPS][bc0][msd/MAGESTIC_marker]
    # But in R1, we may not reach msd, so use the marker position

    # Estimate SPS position (20bp before marker)
    sps_end = marker_pos
    sps_start = sps_end - SPS_LENGTH

    if sps_start < 50:  # Too close to start, likely wrong
        return None

    sps = seq[sps_start:sps_end]

    # Find donor start (after SaCas9_SceI/RE site)
    re_site = "GTTTGAAGAGC"
    re_pos = seq.find(re_site)

    if re_pos < 0 or re_pos > sps_start - 50:  # Donor should be at least 50bp
        return None

    donor_start = re_pos + len(re_site)
    donor_end = sps_start
    donor = seq[donor_start:donor_end]

    return donor, sps


# =============================================================================
# MSD Read Processing (BC0, SPS, MSD, donor extraction from msd-side read)
# =============================================================================

def process_msd_read(seq: str, config: ReadConfig, is_reverse_complement: bool = True) -> Dict[str, str]:
    """
    Process the read containing msd/bc0/SPS/donor.

    For reversed orientation (228x90): R1 contains msd region in reverse complement.
    The R1 read structure when sequenced from 3' end is:
        [3' flanking][msd_rc][bc0_rc][SPS_rc][donor_rc][bc1_suffix_rc]...

    For standard orientation (2x150): R2 contains msd region similarly.

    Args:
        seq: The sequence to process
        config: ReadConfig with extraction parameters
        is_reverse_complement: If True, sequences need RC to get forward orientation

    Returns:
        Dict with extracted components (in forward orientation)
    """
    read_length = len(seq)
    result = {}

    # First, find MSD anchor (or its RC)
    msd_target = MSD_RC if is_reverse_complement else MSD_SEQUENCE
    msd_pos = seq.find(msd_target)

    if msd_pos >= 0:
        # MSD found - use anchor-based extraction
        if is_reverse_complement:
            # Structure: [3' flanking][msd_rc][bc0_rc][SPS_rc][donor_rc]...
            # msd_rc starts at msd_pos, bc0_rc comes AFTER it
            msd_end = msd_pos + MSD_LENGTH
            bc0_start = msd_end
            bc0_end = bc0_start + BC0_LENGTH
            sps_start = bc0_end
            sps_end = sps_start + SPS_LENGTH

            # Extract components (need to RC to get forward orientation)
            if bc0_end <= read_length:
                result['bc0'] = rev_comp(seq[bc0_start:bc0_end])

            if sps_end <= read_length:
                result['SPS'] = rev_comp(seq[sps_start:sps_end])

            # Donor is after SPS
            if sps_end < read_length:
                donor_rc = seq[sps_end:]
                result['donor'] = rev_comp(donor_rc)
                result['donor_source'] = 'msd_read'

            # donor_bc0_fragment: last 30bp of donor + SPS + bc0 = 60bp
            # In the RC read, this is: [30bp donor_rc][SPS_rc][bc0_rc]
            # We need to extract from the end of donor region backwards
            if sps_end < read_length:
                donor_end = read_length
                frag_start = bc0_start
                frag_end = min(sps_end + 30, donor_end)  # bc0_rc + SPS_rc + 30bp donor_rc
                if frag_end - frag_start >= 40:  # At least partial
                    result['donor_bc0_fragment'] = rev_comp(seq[frag_start:frag_end])

            result['msd'] = MSD_SEQUENCE  # We found it, so it's correct
        else:
            # Forward orientation - msd is at start
            msd_start = msd_pos
            msd_end = msd_start + MSD_LENGTH
            bc0_start = msd_end
            bc0_end = bc0_start + BC0_LENGTH
            sps_start = bc0_end
            sps_end = sps_start + SPS_LENGTH

            if read_length >= msd_end:
                result['msd'] = seq[msd_start:msd_end]
            if read_length >= bc0_end:
                result['bc0'] = seq[bc0_start:bc0_end]
            if read_length >= sps_end:
                result['SPS'] = seq[sps_start:sps_end]
            if read_length > sps_end:
                result['donor'] = seq[sps_end:]
                result['donor_source'] = 'msd_read'
    else:
        # MSD not found, try position-based extraction assuming msd_rc near start
        if is_reverse_complement:
            # Check if MSD_RC might be at a slightly different position (allow fuzzy)
            # For now, assume it's near the start (positions 0-40)
            for offset in range(0, min(40, read_length - MSD_LENGTH)):
                if seq[offset:offset + MSD_LENGTH] == MSD_RC:
                    msd_pos = offset
                    break
            else:
                # Still not found, skip extraction
                return result

            # Use the found position
            msd_end = msd_pos + MSD_LENGTH
            bc0_start = msd_end
            bc0_end = bc0_start + BC0_LENGTH
            sps_start = bc0_end
            sps_end = sps_start + SPS_LENGTH

            if bc0_end <= read_length:
                result['bc0'] = rev_comp(seq[bc0_start:bc0_end])
            if sps_end <= read_length:
                result['SPS'] = rev_comp(seq[sps_start:sps_end])
            if sps_end < read_length:
                result['donor'] = rev_comp(seq[sps_end:])
                result['donor_source'] = 'msd_read_positional'

    return result


def process_r2_seq_dynamic(seq: str, config: ReadConfig) -> Dict[str, str]:
    """
    Process R2 sequence with dynamic positions based on read length.
    Legacy wrapper for backward compatibility - delegates to process_msd_read.
    """
    return process_msd_read(seq, config, is_reverse_complement=True)


# =============================================================================
# Combined Processing
# =============================================================================

def process_read_pair(
    r1_seq: str,
    r2_seq: str,
    config: ReadConfig,
) -> Optional[Dict[str, Any]]:
    """
    Process a paired R1/R2 read to extract all bc1-donor-bc0 components.

    Handles both standard and reversed primer orientations:
    - Standard: R1 has BC1, R2 has msd/bc0/SPS/donor (in RC)
    - Reversed: R2 has BC1, R1 has msd/bc0/SPS/donor (in RC)

    Returns:
        Dict with all extracted components or None if extraction fails
    """
    result = {}

    # Determine which read to use for BC1 vs msd based on orientation
    if config.orientation == "reversed":
        bc1_seq = r2_seq
        msd_seq = r1_seq
    else:  # standard
        bc1_seq = r1_seq
        msd_seq = r2_seq

    # Extract BC1 from the BC1-containing read
    bc1, ps1, bc1_length, bc1_prefix_found = extract_bc1_from_r1(bc1_seq, config)

    if not bc1:
        return None

    result['bc1'] = bc1
    result['PS1'] = ps1
    result['bc1_length'] = bc1_length
    result['bc1_prefix_found'] = bc1_prefix_found

    # Extract bc0, SPS, msd, and donor from the msd-containing read
    msd_components = process_msd_read(msd_seq, config, is_reverse_complement=True)
    result.update(msd_components)

    # Build donor_bc0_fragment if we have the components but it wasn't already built
    if 'donor_bc0_fragment' not in result:
        if 'donor' in result and 'SPS' in result and 'bc0' in result:
            # Use last 30bp of donor + SPS + bc0
            donor = result['donor']
            if len(donor) >= 30:
                result['donor_bc0_fragment'] = donor[-30:] + result['SPS'] + result['bc0']

    return result


# =============================================================================
# File Processing
# =============================================================================

def process_files(
    r1_path: str,
    r2_path: str,
    output_path: str,
    config: ReadConfig = None,
    overwrite: bool = True,
):
    """Process paired FASTQ files and write results."""
    if not overwrite and os.path.exists(output_path):
        logger.info(f"Output file {output_path} exists. Skipping.")
        return

    # Auto-detect read lengths and orientation if config not provided
    if config is None:
        r1_len, r2_len = detect_read_lengths(r1_path, r2_path)
        logger.info(f"Detected read lengths: R1={r1_len}bp, R2={r2_len}bp")

        # Handle empty files - write empty output and exit
        if r1_len == 0 or r2_len == 0:
            logger.warning(f"Empty FASTQ files detected. Writing empty output.")
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, "w") as f:
                f.write("bc1\tdonor\tsps\tbc0\tcount\n")
            return

        # Sample sequences for BC1_PREFIX and orientation detection
        r1_seqs = [seq for _, seq, _, _ in read_fastq_entries(r1_path, 5000)]
        r2_seqs = [seq for _, seq, _, _ in read_fastq_entries(r2_path, 5000)]

        # Detect orientation
        orientation = detect_orientation(r1_seqs, r2_seqs, n_sample=1000)
        logger.info(f"Detected orientation: {orientation}")

        # Detect BC1_PREFIX from the BC1-containing read
        if orientation == "reversed":
            bc1_prefix = find_bc1_prefix(r2_seqs)
        else:
            bc1_prefix = find_bc1_prefix(r1_seqs)

        if bc1_prefix:
            logger.info(f"Detected BC1_PREFIX: {bc1_prefix}")
        else:
            bc1_prefix = DEFAULT_BC1_PREFIX
            logger.info(f"Using default BC1_PREFIX: {bc1_prefix}")

        config = ReadConfig(
            r1_length=r1_len,
            r2_length=r2_len,
            bc1_prefix=bc1_prefix,
            orientation=orientation,
        )

    logger.info(f"Extraction strategy: {config.strategy}")
    logger.info(f"  Orientation: {config.orientation}")
    logger.info(f"  BC1 read: {config.bc1_read}")
    logger.info(f"  MSD read: {config.msd_read}")
    logger.info(f"  Can extract full donor from msd read: {config.can_extract_full_donor_from_msd_read}")

    # Process reads
    results = Counter()
    length_stats = Counter()
    reads_processed = 0
    reads_successful = 0

    start_time = timeit.default_timer()

    if HAS_PYSAM:
        with pysam.FastxFile(r1_path) as r1_file, pysam.FastxFile(r2_path) as r2_file:
            for r1_read, r2_read in tqdm(zip(r1_file, r2_file), desc="Processing reads"):
                reads_processed += 1
                if reads_processed > MAX_READS:
                    break

                parsed = process_read_pair(r1_read.sequence, r2_read.sequence, config)

                if parsed and 'bc1' in parsed:
                    # Create key tuple for counting
                    key = (
                        parsed.get('bc1', ''),
                        parsed.get('donor', ''),
                        parsed.get('bc0', ''),
                        parsed.get('SPS', ''),
                        parsed.get('donor_bc0_fragment', ''),
                    )
                    results[key] += 1
                    reads_successful += 1

                    if 'bc1_length' in parsed:
                        length_stats[parsed['bc1_length']] += 1
    else:
        r1_iter = read_fastq_entries(r1_path)
        r2_iter = read_fastq_entries(r2_path)

        for (_, r1_seq, _, _), (_, r2_seq, _, _) in tqdm(zip(r1_iter, r2_iter), desc="Processing reads"):
            reads_processed += 1
            if reads_processed > MAX_READS:
                break

            parsed = process_read_pair(r1_seq, r2_seq, config)

            if parsed and 'bc1' in parsed:
                key = (
                    parsed.get('bc1', ''),
                    parsed.get('donor', ''),
                    parsed.get('bc0', ''),
                    parsed.get('SPS', ''),
                    parsed.get('donor_bc0_fragment', ''),
                )
                results[key] += 1
                reads_successful += 1

                if 'bc1_length' in parsed:
                    length_stats[parsed['bc1_length']] += 1

    elapsed = timeit.default_timer() - start_time

    # Log statistics
    logger.info(f"Processing complete: {elapsed:.2f}s")
    logger.info(f"  Reads processed: {reads_processed:,}")
    logger.info(f"  Successful extractions: {reads_successful:,} ({100*reads_successful/reads_processed:.1f}%)")
    logger.info(f"  Unique bc1-donor-bc0 combinations: {len(results):,}")

    if length_stats:
        logger.info("  BC1 length distribution:")
        total = sum(length_stats.values())
        for length in sorted(length_stats.keys()):
            count = length_stats[length]
            pct = 100 * count / total
            logger.info(f"    {length}bp: {count:,} ({pct:.1f}%)")

    # Write output
    with open(output_path, 'w') as f:
        header = "\t".join(['counts', 'bc1', 'donor', 'bc0', 'SPS', 'donor_bc0_fragment'])
        f.write(header + "\n")

        for key, count in sorted(results.items(), key=lambda x: -x[1]):
            bc1, donor, bc0, sps, donor_bc0_fragment = key
            f.write(f"{count}\t{bc1}\t{donor}\t{bc0}\t{sps}\t{donor_bc0_fragment}\n")

    logger.info(f"Wrote {len(results)} unique sequences to {output_path}")


# =============================================================================
# Main
# =============================================================================

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse bc1-donor-bc0 from paired FASTQ with arbitrary read lengths.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect read lengths
  python 01c_parse_arbitrary_length_fastq.py R1.fastq.gz R2.fastq.gz output.tsv

  # Specify read lengths explicitly
  python 01c_parse_arbitrary_length_fastq.py R1.fastq.gz R2.fastq.gz output.tsv --r1-length 228 --r2-length 90
        """
    )

    parser.add_argument("r1_fastq", help="R1 FASTQ file (gzipped or not)")
    parser.add_argument("r2_fastq", help="R2 FASTQ file (gzipped or not)")
    parser.add_argument("output_tsv", help="Output TSV file")

    parser.add_argument("--r1-length", type=int, help="R1 read length (auto-detected if not specified)")
    parser.add_argument("--r2-length", type=int, help="R2 read length (auto-detected if not specified)")
    parser.add_argument("--bc1-prefix", type=str, help=f"BC1 prefix (default: auto-detect or {DEFAULT_BC1_PREFIX})")
    parser.add_argument("--overwrite", action="store_true", default=True, help="Overwrite existing output")
    parser.add_argument("--no-overwrite", action="store_true", help="Skip if output exists")

    return parser.parse_args()


def main():
    """Entry point."""
    args = parse_args()

    config = None
    if args.r1_length and args.r2_length:
        config = ReadConfig(
            r1_length=args.r1_length,
            r2_length=args.r2_length,
            bc1_prefix=args.bc1_prefix or DEFAULT_BC1_PREFIX,
        )

    overwrite = not args.no_overwrite

    process_files(
        args.r1_fastq,
        args.r2_fastq,
        args.output_tsv,
        config=config,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    main()
