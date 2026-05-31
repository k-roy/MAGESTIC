#!/usr/bin/env python3
"""
Step 02: Combine Collapsed FASTQs and Parse Guide-Donor-BC0 Sequences

This script processes collapsed FASTQ files from different sequencing runs and
primer orientations to extract guide, donor, and BC0 sequences.

Processing:
- Reads collapsed FASTQs (where read name = count of identical reads)
- Reverse complements reads from reverse primer pairs to match forward orientation
- Parses sequence components using anchor-based extraction (msd from 3' end)
- Combines counts across sequencing replicates
- Outputs per-library TSV with counts by run and direction

Key feature: Uses end-slice approach for donor_bc0_fragment extraction, which is
consistent regardless of 5' trimming variations or full donor length.

Input: Collapsed FASTQ files from Step 01
       Format: {date}_{plate}_{sample}_collapsed.fastq
       Read names contain counts of identical sequences

Output: Per-library TSV files with columns:
        library_ID, total_counts, total_fwd_counts, total_rev_counts,
        fwd_MMDD, rev_MMDD (per sequencing run), sequence components

Usage:
    # Use with environment variables (set by Snakefile)
    python 02_combine_step_1_collapsed_fastq.py

    # Or specify paths directly
    python 02_combine_step_1_collapsed_fastq.py \
        --project-dir /path/to/project \
        --keyfile step_1_sample_key_unified.tsv

Author: Kevin R. Roy
"""

import argparse
import os
import sys
import timeit
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set

import pysam
from tqdm import tqdm

# Import anchor-based extraction utilities
from magestic.utils.sequence_anchors import (
    find_anchor,
    find_anchor_from_end,
    extract_donor_bc0_fragment,
    DEFAULT_DONOR_BC0_FRAGMENT_LENGTH,
)


# ============================================================================
# Configuration
# ============================================================================

def get_config_from_env_or_defaults() -> dict:
    """
    Get configuration from environment variables or defaults.

    Environment variables (typically set by Snakefile):
        PROJECT_DIR: Project directory path
        DATATYPE: Subproject name (default: guide_donor_bc0)
        SAMPLE_KEYFILE: Path to sample keyfile (relative to keyfiles dir or absolute)
    """
    BASE_DIR = Path(os.environ.get(
        "BASE_DIR",
        "/path/to"
    ))

    # Default project (can be overridden)
    DEFAULT_PROJECT_DIR = BASE_DIR / "projects/QTL/20240806_yL406_SpG_QTL_screen"

    project_dir = Path(os.environ.get("PROJECT_DIR", DEFAULT_PROJECT_DIR))
    datatype = os.environ.get("DATATYPE", "guide_donor_bc0")

    return {
        "base_dir": BASE_DIR,
        "project_dir": project_dir,
        "datatype": datatype,
        "processed_data_dir": project_dir / "processed_data" / datatype,
        "keyfile_dir": project_dir / "keyfiles" / datatype,
        "sample_keyfile": os.environ.get("SAMPLE_KEYFILE", "step_1_sample_key_unified.tsv"),
    }


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine collapsed FASTQs and parse guide-donor-BC0 sequences"
    )
    parser.add_argument(
        "--project-dir", "-p",
        type=Path,
        help="Project directory path (overrides PROJECT_DIR env var)"
    )
    parser.add_argument(
        "--keyfile", "-k",
        type=str,
        help="Sample keyfile name or path (overrides SAMPLE_KEYFILE env var)"
    )
    parser.add_argument(
        "--datatype", "-d",
        type=str,
        default="guide_donor_bc0",
        help="Subproject/datatype name (default: guide_donor_bc0)"
    )
    parser.add_argument(
        "--min-quality", "-q",
        type=int,
        default=0,
        help="Minimum average Phred quality score to keep a read (default: 0, no filtering)"
    )
    return parser.parse_args()


def calculate_avg_phred(qual_string: str) -> float:
    """
    Calculate average Phred quality score from quality string.

    Args:
        qual_string: ASCII-encoded quality string

    Returns:
        Average Phred quality score
    """
    if not qual_string:
        return 0.0
    return sum(ord(c) - 33 for c in qual_string) / len(qual_string)


# ============================================================================
# Sequence Constants
# ============================================================================

@dataclass(frozen=True)
class SequenceConstants:
    """Fixed sequence lengths and patterns for parsing reads."""
    guide_length: int = 20
    RE_site: str = "GTTTGAAGAGC"  # SpCas9/SpG
    RE_site_lbcas12a: str = "TTTCGAAGAGC"  # LbCas12a (one less G)
    RE_site_length: int = 11
    bc0_length: int = 10
    SPS_length: int = 20
    msd: str = "AGGAAAACAGAC"  # Constant 3' anchor
    msd_length: int = 12
    donor_bc0_fragment_length: int = 60  # For end-slice matching


SEQ_CONST = SequenceConstants()


# ============================================================================
# Utility Functions
# ============================================================================

def rev_comp(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


def setup_directories(dirs: List[Path]) -> None:
    """Create necessary output directories."""
    for directory in dirs:
        directory.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Sample Keyfile Loading
# ============================================================================

@dataclass
class SampleInfo:
    """Information about a single sample from the keyfile."""
    sample_number: str
    library_ID: str
    inner_primers: str
    outer_primers: str
    sequencing_dates: List[str]
    polymerase: str


def load_sample_keyfile(keyfile_path: Path) -> Tuple[Dict[str, List[SampleInfo]], Set[str]]:
    """
    Load sample keyfile and organize by library_ID.

    Keyfile format (6 columns, tab-separated):
    sample_number  library_ID  inner_primers  outer_primers  sequencing_dates  polymerase

    Returns:
        Tuple of:
        - Dict mapping library_ID to list of SampleInfo objects
        - Set of all sequencing dates found in keyfile
    """
    library_to_samples: Dict[str, List[SampleInfo]] = {}
    all_dates: Set[str] = set()

    with open(keyfile_path) as infile:
        header = infile.readline().strip().split("\t")
        expected_cols = ["sample_number", "library_ID", "inner_primers",
                         "outer_primers", "sequencing_dates", "polymerase"]

        if len(header) < 5:
            raise ValueError(
                f"Expected at least 5 columns in keyfile, got {len(header)}. "
                f"Expected: {expected_cols}"
            )

        for line_num, line in enumerate(infile, start=2):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 5:
                print(f"  WARNING: Line {line_num} has {len(fields)} fields, expected at least 5")
                continue

            sample_number = fields[0]
            library_ID = fields[1]
            inner_primers = fields[2]
            outer_primers = fields[3]
            dates_str = fields[4]
            polymerase = fields[5] if len(fields) > 5 else "Q5"

            # Parse comma-separated dates
            sequencing_dates = [d.strip() for d in dates_str.split(",") if d.strip()]
            all_dates.update(sequencing_dates)

            sample = SampleInfo(
                sample_number=sample_number,
                library_ID=library_ID,
                inner_primers=inner_primers,
                outer_primers=outer_primers,
                sequencing_dates=sequencing_dates,
                polymerase=polymerase,
            )

            library_to_samples.setdefault(library_ID, []).append(sample)

    return library_to_samples, all_dates


# ============================================================================
# Sequence Processing
# ============================================================================

def parse_sequence(seq: str, cas_type: str = "SpCas9") -> Tuple[str, str, str, str, str, str, str]:
    """
    Parse a sequence to extract guide, donor, BC0, and other components.

    Uses anchor-based extraction:
    - Find msd from 3' end (most reliable anchor)
    - Find RE site from 5' end
    - Extract features relative to anchors

    The end-slice donor_bc0_fragment is consistent regardless of:
    - 5' trimming variations
    - Full donor length

    Args:
        seq: Full sequence string
        cas_type: "SpCas9", "SpG", or "LbCas12a" (determines RE site pattern)

    Returns:
        Tuple of (leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment)
        Returns "NA" for all fields if anchors not found
    """
    # Find msd from end (most reliable anchor)
    msd_pos = find_anchor_from_end(seq, SEQ_CONST.msd, max_distance_from_end=30)
    if msd_pos is None:
        return ("NA",) * 7

    # Find RE site - try SpCas9/SpG first, then LbCas12a
    RE_site_pattern = SEQ_CONST.RE_site if cas_type != "LbCas12a" else SEQ_CONST.RE_site_lbcas12a
    RE_site_pos = find_anchor(seq, RE_site_pattern, end=min(150, len(seq)))

    # If SpCas9 pattern not found, try LbCas12a
    if RE_site_pos is None and cas_type != "LbCas12a":
        RE_site_pos = find_anchor(seq, SEQ_CONST.RE_site_lbcas12a, end=min(150, len(seq)))

    if RE_site_pos is None:
        return ("NA",) * 7

    # Calculate positions based on anchors
    guide_start = RE_site_pos - SEQ_CONST.guide_length
    donor_start = RE_site_pos + SEQ_CONST.RE_site_length
    bc0_end = msd_pos
    bc0_start = bc0_end - SEQ_CONST.bc0_length
    sps_end = bc0_start
    sps_start = sps_end - SEQ_CONST.SPS_length
    donor_end = sps_start

    # Validate positions
    if guide_start < 0 or donor_end <= donor_start:
        return ("NA",) * 7

    # Extract features
    leader = seq[:guide_start]
    guide = seq[guide_start:RE_site_pos]
    donor = seq[donor_start:donor_end]
    SPS = seq[sps_start:sps_end]
    bc0 = seq[bc0_start:bc0_end]
    msd = seq[msd_pos:msd_pos + SEQ_CONST.msd_length]

    # End-slice donor_bc0_fragment (consistent key for linking)
    donor_bc0_fragment = extract_donor_bc0_fragment(
        seq,
        msd_pos=msd_pos,
        msd_pattern=SEQ_CONST.msd,
        fragment_length=SEQ_CONST.donor_bc0_fragment_length,
        bc0_length=SEQ_CONST.bc0_length,
        sps_length=SEQ_CONST.SPS_length,
    )
    if donor_bc0_fragment is None:
        donor_bc0_fragment = "NA"

    return leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment


def process_collapsed_fastq(
    fastq_path: Path,
    inner_primers: str,
    sequencing_run_date: str,
    library_ID: str,
    bc0_counts: Dict[tuple, Dict[str, int]],
    all_dates: Set[str],
    min_quality: int = 0,
    quality_stats: Optional[Dict[str, int]] = None,
) -> Tuple[int, int]:
    """
    Process a single collapsed FASTQ file and accumulate BC0 counts.

    Args:
        fastq_path: Path to collapsed FASTQ file
        inner_primers: Primer pair ID (determines orientation)
        sequencing_run_date: Date string for count column naming
        library_ID: Library identifier
        bc0_counts: Accumulator dict to update (modified in place)
        all_dates: Set of all sequencing dates (for initializing count columns)
        min_quality: Minimum average Phred quality to keep (0 = no filtering)
        quality_stats: Optional dict to accumulate quality filtering stats

    Returns:
        Tuple of (unique_sequences_processed, total_counts_processed)
    """
    total_unique_seqs = 0
    total_counts = 0
    filtered_low_quality = 0
    filtered_counts = 0

    # Determine direction based on primer pair
    is_reverse = inner_primers in ("KR1884-KR1953", "KR1882-KR1953")
    seq_direction = "rev" if is_reverse else "fwd"
    counts_col = f"{seq_direction}_{sequencing_run_date[-4:]}"

    with pysam.FastxFile(str(fastq_path)) as fin:
        for entry in fin:
            total_unique_seqs += 1
            counts = int(entry.name)
            total_counts += counts
            seq = entry.sequence
            qual = entry.quality

            # Quality filtering (if enabled)
            if min_quality > 0 and qual:
                avg_qual = calculate_avg_phred(qual)
                if avg_qual < min_quality:
                    filtered_low_quality += 1
                    filtered_counts += counts
                    continue

            # Reverse complement reads from reverse primer pairs
            if is_reverse:
                seq = rev_comp(seq)
                qual = qual[::-1] if qual else ""
                # Add trailing bases for KR1882-KR1953
                if inner_primers == "KR1882-KR1953":
                    seq += "GAC"
                    qual += "EEE"

            # Parse sequence components
            leader, guide, donor, SPS, bc0, msd, donor_bc0_frag = parse_sequence(seq)
            bc0_info = (leader, guide, donor, SPS, bc0, msd, donor_bc0_frag, seq)

            # Initialize count dict for new bc0_info
            if bc0_info not in bc0_counts:
                bc0_counts[bc0_info] = {}
                # Initialize all possible count columns from all dates in keyfile
                for run_date in all_dates:
                    for direction in ["fwd", "rev"]:
                        col_name = f"{direction}_{run_date[-4:]}"
                        bc0_counts[bc0_info][col_name] = 0

            bc0_counts[bc0_info][counts_col] += counts

    # Update quality filtering stats if provided
    if quality_stats is not None:
        quality_stats["filtered_low_quality"] = (
            quality_stats.get("filtered_low_quality", 0) + filtered_low_quality
        )
        quality_stats["filtered_counts"] = (
            quality_stats.get("filtered_counts", 0) + filtered_counts
        )

    return total_unique_seqs, total_counts


def calculate_totals(bc0_counts: Dict[tuple, Dict[str, int]]) -> None:
    """
    Calculate total counts for each bc0_info entry.

    Adds total_counts, total_fwd_counts, total_rev_counts to each entry.

    Args:
        bc0_counts: Dict to update in place
    """
    for bc0_info, counts_dict in bc0_counts.items():
        total = sum(counts_dict.values())
        total_fwd = sum(v for k, v in counts_dict.items() if k.startswith("fwd_"))
        total_rev = sum(v for k, v in counts_dict.items() if k.startswith("rev_"))

        counts_dict["total_counts"] = total
        counts_dict["total_fwd_counts"] = total_fwd
        counts_dict["total_rev_counts"] = total_rev


def write_counts_file(
    output_path: Path,
    library_ID: str,
    bc0_counts: Dict[tuple, Dict[str, int]],
) -> None:
    """
    Write counts to TSV file.

    Args:
        output_path: Path to output file
        library_ID: Library identifier
        bc0_counts: Dict of bc0_info -> counts
    """
    # Determine all column names dynamically
    if not bc0_counts:
        return

    # Get all count column names from first entry
    sample_counts = next(iter(bc0_counts.values()))
    total_cols = ["total_counts", "total_fwd_counts", "total_rev_counts"]
    per_run_cols = sorted(k for k in sample_counts.keys() if k not in total_cols)
    seq_cols = ["leader", "guide", "donor", "SPS", "bc0", "msd", "donor_bc0_fragment", "seq"]

    header = ["library_ID"] + total_cols + per_run_cols + seq_cols

    with open(output_path, "w") as outfile:
        outfile.write("\t".join(header) + "\n")

        for bc0_info, counts_dict in bc0_counts.items():
            row = [library_ID]
            row.extend(str(counts_dict.get(col, 0)) for col in total_cols)
            row.extend(str(counts_dict.get(col, 0)) for col in per_run_cols)
            row.extend(str(val) for val in bc0_info)
            outfile.write("\t".join(row) + "\n")


# ============================================================================
# Main Processing
# ============================================================================

def process_library(
    library_ID: str,
    samples: List[SampleInfo],
    collapsed_fastq_dir: Path,
    counts_dir: Path,
    all_dates: Set[str],
    min_quality: int = 0,
) -> Dict[str, any]:
    """
    Process all samples for a single library.

    Args:
        library_ID: Library identifier
        samples: List of SampleInfo for this library
        collapsed_fastq_dir: Directory containing collapsed FASTQs
        counts_dir: Directory for output files
        all_dates: Set of all sequencing dates
        min_quality: Minimum average Phred quality (0 = no filtering)

    Returns:
        Dict with processing statistics
    """
    start_time = timeit.default_timer()
    bc0_counts: Dict[tuple, Dict[str, int]] = {}
    quality_stats: Dict[str, int] = {"filtered_low_quality": 0, "filtered_counts": 0}

    total_files = 0
    total_seqs = 0
    total_counts = 0
    missing_files = []

    for sample in samples:
        for seq_date in sample.sequencing_dates:
            # Build expected filename
            sample_name = f"{seq_date}_{sample.outer_primers}_{sample.sample_number}"
            fastq_path = collapsed_fastq_dir / f"{sample_name}_collapsed.fastq"

            if not fastq_path.exists():
                missing_files.append(str(fastq_path))
                continue

            total_files += 1
            n_seqs, n_counts = process_collapsed_fastq(
                fastq_path,
                sample.inner_primers,
                seq_date,
                library_ID,
                bc0_counts,
                all_dates,
                min_quality=min_quality,
                quality_stats=quality_stats,
            )
            total_seqs += n_seqs
            total_counts += n_counts

    # Calculate totals
    calculate_totals(bc0_counts)

    # Write output
    output_path = counts_dir / f"{library_ID}_guide_donor_bc0_counts.tsv"
    write_counts_file(output_path, library_ID, bc0_counts)

    elapsed = timeit.default_timer() - start_time

    return {
        "library_ID": library_ID,
        "files_processed": total_files,
        "unique_sequences": len(bc0_counts),
        "total_seqs_in_files": total_seqs,
        "total_counts": total_counts,
        "missing_files": missing_files,
        "elapsed_seconds": elapsed,
        "output_file": str(output_path),
        "filtered_low_quality": quality_stats["filtered_low_quality"],
        "filtered_counts": quality_stats["filtered_counts"],
    }


def main():
    """Main processing pipeline."""
    args = parse_args()

    # Get configuration
    config = get_config_from_env_or_defaults()

    # Override with command-line args if provided
    if args.project_dir:
        config["project_dir"] = args.project_dir
        config["processed_data_dir"] = args.project_dir / "processed_data" / args.datatype
        config["keyfile_dir"] = args.project_dir / "keyfiles" / args.datatype

    if args.keyfile:
        config["sample_keyfile"] = args.keyfile

    # Resolve keyfile path
    keyfile_name = config["sample_keyfile"]
    if Path(keyfile_name).is_absolute():
        sample_keyfile = Path(keyfile_name)
    else:
        sample_keyfile = config["keyfile_dir"] / keyfile_name

    # Setup derived paths
    collapsed_fastq_dir = config["processed_data_dir"] / "collapsed_fastq"
    counts_dir = config["processed_data_dir"] / "guide_donor_bc0"
    log_dir = config["processed_data_dir"] / "log"

    start_time = datetime.now()
    print("=" * 70)
    print("Step 02: Combine Collapsed FASTQs and Parse Guide-Donor-BC0")
    print(f"Started: {start_time}")
    print(f"Project: {config['project_dir']}")
    print(f"Keyfile: {sample_keyfile}")
    if args.min_quality > 0:
        print(f"Quality filter: min Phred Q{args.min_quality}")
    print("=" * 70)

    # Setup directories
    setup_directories([counts_dir, log_dir])

    # Validate paths
    if not sample_keyfile.exists():
        print(f"ERROR: Sample keyfile not found: {sample_keyfile}")
        return 1

    if not collapsed_fastq_dir.exists():
        print(f"ERROR: Collapsed FASTQ directory not found: {collapsed_fastq_dir}")
        return 1

    # Load sample keyfile
    print(f"\nLoading sample keyfile: {sample_keyfile}")
    library_to_samples, all_dates = load_sample_keyfile(sample_keyfile)
    print(f"Found sequencing dates: {', '.join(sorted(all_dates))}")

    # Filter to actual libraries (exclude 'empty')
    libraries = [lib for lib in library_to_samples.keys() if lib != "empty"]
    print(f"Found {len(libraries)} libraries to process: {', '.join(libraries)}")

    # Process each library
    print("\nProcessing libraries...")
    all_stats = []

    for library_ID in tqdm(libraries, desc="Libraries"):
        samples = library_to_samples[library_ID]
        stats = process_library(
            library_ID,
            samples,
            collapsed_fastq_dir,
            counts_dir,
            all_dates,
            min_quality=args.min_quality,
        )
        all_stats.append(stats)

        # Report
        print(f"\n  {library_ID}:")
        print(f"    Files processed: {stats['files_processed']}")
        print(f"    Unique BC0 combinations: {stats['unique_sequences']:,}")
        print(f"    Total counts: {stats['total_counts']:,}")
        if args.min_quality > 0 and stats.get("filtered_low_quality", 0) > 0:
            filtered_pct = 100 * stats["filtered_counts"] / stats["total_counts"] if stats["total_counts"] > 0 else 0
            print(f"    Filtered (low quality): {stats['filtered_low_quality']:,} seqs ({filtered_pct:.1f}% of counts)")
        print(f"    Time: {stats['elapsed_seconds']:.1f}s")

        if stats["missing_files"]:
            print(f"    WARNING: {len(stats['missing_files'])} missing files")
            for f in stats["missing_files"][:3]:
                print(f"      - {Path(f).name}")
            if len(stats["missing_files"]) > 3:
                print(f"      ... and {len(stats['missing_files']) - 3} more")

    # Summary
    total_elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print("Step B Complete")
    print(f"  Libraries processed: {len(libraries)}")
    print(f"  Output directory: {counts_dir}")
    print(f"  Total elapsed time: {total_elapsed}")
    print(f"{'=' * 70}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
