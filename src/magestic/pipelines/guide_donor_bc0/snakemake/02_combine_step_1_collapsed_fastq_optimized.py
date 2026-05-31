#!/usr/bin/env python3
"""
Step 02: Combine Collapsed FASTQs and Parse Guide-Donor-BC0 (OPTIMIZED)

OPTIMIZATIONS over original:
1. Parallel file processing - Process collapsed FASTQs concurrently using multiprocessing
2. Batch sequence parsing - Process sequences in batches for better memory efficiency
3. Pandas-based aggregation - More efficient than nested dicts for large datasets
4. Streaming file I/O - Process files without loading entirely into memory

Input: Collapsed FASTQ files from Step 01
Output: Per-library guide_donor_bc0_counts.tsv files

Author: Kevin R. Roy
"""

import os
import sys
import argparse
import timeit
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Set, Tuple, Optional
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import pandas as pd
import numpy as np
import pysam

# Import anchor-based extraction utilities
from magestic.utils.sequence_anchors import (
    find_anchor,
    find_anchor_from_end,
    extract_donor_bc0_fragment,
)


# ============================================================================
# Configuration
# ============================================================================

@dataclass(frozen=True)
class SequenceConstants:
    """Fixed sequence lengths and patterns for parsing reads."""
    guide_length: int = 20
    RE_site: str = "GTTTGAAGAGC"  # SpCas9/SpG
    RE_site_lbcas12a: str = "TTTCGAAGAGC"  # LbCas12a
    RE_site_length: int = 11
    bc0_length: int = 10
    SPS_length: int = 20
    msd: str = "AGGAAAACAGAC"  # Constant 3' anchor
    msd_length: int = 12
    donor_bc0_fragment_length: int = 60


SEQ_CONST = SequenceConstants()


@dataclass
class SampleInfo:
    """Information about a single sample from the keyfile."""
    sample_number: str
    library_ID: str
    inner_primers: str
    outer_primers: str
    sequencing_dates: List[str]
    polymerase: str


# ============================================================================
# Utility Functions
# ============================================================================

def rev_comp(seq: str) -> str:
    """Return reverse complement of DNA sequence."""
    return seq.translate(
        str.maketrans("ACGTacgtRYMKrymkVBHDvbhdNn ", "TGCAtgcaYRKMyrkmBVDHbvdhNn ")
    )[::-1]


def parse_sequence(seq: str) -> Tuple[str, str, str, str, str, str, str]:
    """
    Parse a sequence to extract guide, donor, BC0, and other components.
    Auto-detects SpCas9 vs LbCas12a RE sites.
    """
    # Find msd from end (most reliable anchor)
    msd_pos = find_anchor_from_end(seq, SEQ_CONST.msd, max_distance_from_end=30)
    if msd_pos is None:
        return ("NA",) * 7

    # Find RE site - try SpCas9/SpG first, then LbCas12a
    RE_site_pos = find_anchor(seq, SEQ_CONST.RE_site, end=min(150, len(seq)))
    if RE_site_pos is None:
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

    # Extract donor_bc0_fragment
    frag_end = msd_pos
    frag_start = frag_end - SEQ_CONST.donor_bc0_fragment_length
    if frag_start >= 0:
        donor_bc0_fragment = seq[frag_start:frag_end]
    else:
        donor_bc0_fragment = "NA"

    return leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment


def load_sample_keyfile(keyfile_path: Path) -> Tuple[Dict[str, List[SampleInfo]], Set[str]]:
    """Load sample keyfile and organize by library_ID."""
    library_to_samples: Dict[str, List[SampleInfo]] = {}
    all_dates: Set[str] = set()

    with open(keyfile_path) as infile:
        header = infile.readline().strip().split("\t")

        for line in infile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) < 5:
                continue

            sample_number = fields[0]
            library_ID = fields[1]
            inner_primers = fields[2]
            outer_primers = fields[3]
            dates_str = fields[4]
            polymerase = fields[5] if len(fields) > 5 else "Q5"

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
# Optimized Processing Functions
# ============================================================================

def process_single_fastq(args: Tuple) -> pd.DataFrame:
    """
    Process a single collapsed FASTQ file and return a DataFrame.

    This function is designed to be called in parallel.
    """
    fastq_path, inner_primers, seq_date, library_ID = args

    if not Path(fastq_path).exists():
        return pd.DataFrame()

    # Determine direction based on primer pair
    is_reverse = inner_primers in ("KR1884-KR1953", "KR1882-KR1953")
    seq_direction = "rev" if is_reverse else "fwd"
    counts_col = f"{seq_direction}_{seq_date[-4:]}"

    records = []

    with pysam.FastxFile(str(fastq_path)) as fin:
        for entry in fin:
            counts = int(entry.name)
            seq = entry.sequence

            # Reverse complement reads from reverse primer pairs
            if is_reverse:
                seq = rev_comp(seq)
                if inner_primers == "KR1882-KR1953":
                    seq += "GAC"

            # Parse sequence components
            leader, guide, donor, SPS, bc0, msd, donor_bc0_frag = parse_sequence(seq)

            records.append({
                'leader': leader,
                'guide': guide,
                'donor': donor,
                'SPS': SPS,
                'bc0': bc0,
                'msd': msd,
                'donor_bc0_fragment': donor_bc0_frag,
                'seq': seq,
                'counts': counts,
                'direction': seq_direction,
                'date': seq_date[-4:],
            })

    if not records:
        return pd.DataFrame()

    return pd.DataFrame(records)


def aggregate_counts(dfs: List[pd.DataFrame], all_dates: Set[str]) -> pd.DataFrame:
    """
    Aggregate counts from multiple DataFrames into final format.

    Groups by sequence components and pivots counts by direction_date.
    """
    if not dfs:
        return pd.DataFrame()

    # Combine all dataframes
    combined = pd.concat(dfs, ignore_index=True)

    if combined.empty:
        return pd.DataFrame()

    # Create grouping key from sequence components
    key_cols = ['leader', 'guide', 'donor', 'SPS', 'bc0', 'msd', 'donor_bc0_fragment', 'seq']

    # Create direction_date column for pivoting
    combined['dir_date'] = combined['direction'] + '_' + combined['date']

    # Group and aggregate
    grouped = combined.groupby(key_cols + ['dir_date'])['counts'].sum().reset_index()

    # Pivot to wide format
    pivoted = grouped.pivot_table(
        index=key_cols,
        columns='dir_date',
        values='counts',
        fill_value=0,
        aggfunc='sum'
    ).reset_index()

    # Ensure all date columns exist
    for date in all_dates:
        for direction in ['fwd', 'rev']:
            col = f"{direction}_{date[-4:]}"
            if col not in pivoted.columns:
                pivoted[col] = 0

    # Calculate totals
    date_cols = [c for c in pivoted.columns if c.startswith('fwd_') or c.startswith('rev_')]
    fwd_cols = [c for c in date_cols if c.startswith('fwd_')]
    rev_cols = [c for c in date_cols if c.startswith('rev_')]

    pivoted['total_fwd_counts'] = pivoted[fwd_cols].sum(axis=1) if fwd_cols else 0
    pivoted['total_rev_counts'] = pivoted[rev_cols].sum(axis=1) if rev_cols else 0
    pivoted['total_counts'] = pivoted['total_fwd_counts'] + pivoted['total_rev_counts']

    return pivoted


def process_library_parallel(
    library_ID: str,
    samples: List[SampleInfo],
    collapsed_fastq_dir: Path,
    counts_dir: Path,
    all_dates: Set[str],
    n_workers: int = 8,
) -> Dict:
    """
    Process all samples for a single library using parallel file processing.
    """
    start_time = timeit.default_timer()

    # Build list of (fastq_path, inner_primers, seq_date, library_ID) tuples
    file_args = []
    for sample in samples:
        for seq_date in sample.sequencing_dates:
            sample_name = f"{seq_date}_{sample.outer_primers}_{sample.sample_number}"
            fastq_path = collapsed_fastq_dir / f"{sample_name}_collapsed.fastq"
            file_args.append((str(fastq_path), sample.inner_primers, seq_date, library_ID))

    print(f"    Processing {len(file_args)} FASTQ files with {n_workers} workers...")

    # Process files in parallel
    dfs = []
    total_files = 0
    missing_files = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {executor.submit(process_single_fastq, args): args for args in file_args}

        for future in as_completed(futures):
            args = futures[future]
            try:
                df = future.result()
                if df.empty:
                    if not Path(args[0]).exists():
                        missing_files.append(args[0])
                else:
                    dfs.append(df)
                    total_files += 1
            except Exception as e:
                print(f"    Error processing {args[0]}: {e}")

    print(f"    Aggregating counts from {total_files} files...")

    # Aggregate all results
    result_df = aggregate_counts(dfs, all_dates)

    if result_df.empty:
        return {
            "library_ID": library_ID,
            "files_processed": 0,
            "unique_sequences": 0,
            "total_counts": 0,
            "missing_files": missing_files,
            "elapsed_seconds": timeit.default_timer() - start_time,
            "output_file": None,
        }

    # Add library_ID column
    result_df.insert(0, 'library_ID', library_ID)

    # Reorder columns
    key_cols = ['library_ID', 'total_counts', 'total_fwd_counts', 'total_rev_counts']
    date_cols = sorted([c for c in result_df.columns if c.startswith('fwd_') or c.startswith('rev_')])
    seq_cols = ['leader', 'guide', 'donor', 'SPS', 'bc0', 'msd', 'donor_bc0_fragment', 'seq']

    final_cols = key_cols + date_cols + seq_cols
    final_cols = [c for c in final_cols if c in result_df.columns]
    result_df = result_df[final_cols]

    # Sort by total_counts descending
    result_df = result_df.sort_values('total_counts', ascending=False)

    # Write output
    output_path = counts_dir / f"{library_ID}_guide_donor_bc0_counts.tsv"
    result_df.to_csv(output_path, sep='\t', index=False)

    elapsed = timeit.default_timer() - start_time

    return {
        "library_ID": library_ID,
        "files_processed": total_files,
        "unique_sequences": len(result_df),
        "total_counts": result_df['total_counts'].sum(),
        "missing_files": missing_files,
        "elapsed_seconds": elapsed,
        "output_file": str(output_path),
    }


# ============================================================================
# Main
# ============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Combine collapsed FASTQs and parse guide-donor-BC0 sequences (OPTIMIZED)"
    )
    parser.add_argument(
        "--project-dir", "-p",
        type=Path,
        help="Project directory path"
    )
    parser.add_argument(
        "--keyfile", "-k",
        type=str,
        help="Sample keyfile name or path"
    )
    parser.add_argument(
        "--datatype", "-d",
        type=str,
        default="guide_donor_bc0",
        help="Subproject/datatype name"
    )
    parser.add_argument(
        "--n-workers", "-j",
        type=int,
        default=8,
        help="Number of parallel workers (default: 8)"
    )
    return parser.parse_args()


def main():
    """Main processing pipeline."""
    args = parse_args()

    # Get project directory
    project_dir = args.project_dir or Path(os.environ.get(
        "PROJECT_DIR",
        "/path/to/projects/QTL/20210822_yL274-yL284_QTL_screen"
    ))
    datatype = args.datatype

    processed_data_dir = project_dir / "processed_data" / datatype
    keyfile_dir = project_dir / "keyfiles" / datatype

    print("=" * 70)
    print("Step 02: Combine Collapsed FASTQs and Parse Guide-Donor-BC0 (OPTIMIZED)")
    print(f"Started: {pd.Timestamp.now()}")
    print(f"Project: {project_dir}")
    print(f"Workers: {args.n_workers}")
    print("=" * 70)

    # Resolve keyfile path
    if args.keyfile:
        if Path(args.keyfile).is_absolute():
            sample_keyfile = Path(args.keyfile)
        else:
            sample_keyfile = keyfile_dir / args.keyfile
    else:
        sample_keyfile = keyfile_dir / os.environ.get("SAMPLE_KEYFILE", "step_1_sample_key_unified.tsv")

    print(f"\nLoading sample keyfile: {sample_keyfile}")
    library_to_samples, all_dates = load_sample_keyfile(sample_keyfile)

    print(f"Found sequencing dates: {', '.join(sorted(all_dates))}")
    print(f"Found {len(library_to_samples)} libraries to process: {', '.join(sorted(library_to_samples.keys()))}")

    # Setup directories
    collapsed_fastq_dir = processed_data_dir / "collapsed_fastq"
    counts_dir = processed_data_dir / "guide_donor_bc0"
    counts_dir.mkdir(parents=True, exist_ok=True)

    # Process libraries
    print("\nProcessing libraries...")
    total_start = timeit.default_timer()

    all_stats = []
    for library_ID, samples in sorted(library_to_samples.items()):
        print(f"\n  Library: {library_ID}")
        stats = process_library_parallel(
            library_ID,
            samples,
            collapsed_fastq_dir,
            counts_dir,
            all_dates,
            n_workers=args.n_workers,
        )
        all_stats.append(stats)

        print(f"    Files processed: {stats['files_processed']}")
        print(f"    Unique sequences: {stats['unique_sequences']:,}")
        print(f"    Total counts: {stats['total_counts']:,}")
        print(f"    Time: {stats['elapsed_seconds']:.1f}s")

    total_elapsed = timeit.default_timer() - total_start

    print("\n" + "=" * 70)
    print("Step 02 Complete (OPTIMIZED)")
    print(f"  Libraries processed: {len(all_stats)}")
    print(f"  Total elapsed time: {total_elapsed/60:.1f} minutes")
    print("=" * 70)


if __name__ == "__main__":
    main()
