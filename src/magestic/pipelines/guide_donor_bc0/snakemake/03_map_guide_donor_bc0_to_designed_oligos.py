#!/usr/bin/env python3
"""
Step 03: Map Guide-Donor-BC0 Sequences to Designed Oligos

Maps observed guide-donor-bc0 sequences to the designed oligo pool using
a multi-level matching strategy:

1. Exact match: guide + RE_site + donor matches designed sequence
2. Guide + donor intersection: both guide and donor match same oligo
3. Unambiguous donor: single oligo matches donor (perfect or middle region)
4. Unambiguous guide: single oligo matches guide

This handles sequencing errors and synthesis variations while maintaining
accurate oligo assignment.

Input: Per-library guide_donor_bc0_counts.tsv from Step 02
Output: Per-library matched_guide_donor_bc0_counts.tsv with oligo_name column

Performance optimizations (2026-02-19):
- Vectorized lookup table building (replaces iterrows)
- Parallel library processing via multiprocessing
- Chunked pandas I/O for large files

Usage:
    python 03_map_guide_donor_bc0_to_designed_oligos.py [--n-jobs N]
"""

import argparse
import os
import timeit
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime
from functools import partial
from multiprocessing import cpu_count
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from tqdm import tqdm

# ============================================================================
# Configuration - Read from environment variables set by SLURM script
# ============================================================================

def get_config_from_env():
    """Get project configuration from environment variables."""
    project_dir_str = os.environ.get("PROJECT_DIR")
    if not project_dir_str:
        raise ValueError(
            "PROJECT_DIR environment variable not set. "
            "Run via SLURM script or set: export PROJECT_DIR=/path/to/project"
        )

    project_dir = Path(project_dir_str)
    datatype = os.environ.get("DATATYPE", "guide_donor_bc0")

    processed_data_dir = project_dir / "processed_data" / datatype
    keyfile_dir = project_dir / "keyfiles" / datatype

    # Oligo design file from env or search keyfiles directory
    oligo_design_file = os.environ.get("OLIGO_DESIGN_FILE")
    if oligo_design_file:
        oligo_pool_files = [Path(oligo_design_file)]
    else:
        # Search for oligo design files in keyfiles directory
        oligo_pool_files = list(keyfile_dir.glob("*Twist*.tsv"))
        if not oligo_pool_files:
            raise ValueError(f"No oligo design files found in {keyfile_dir}")

    return {
        "project_dir": project_dir,
        "processed_data_dir": processed_data_dir,
        "keyfile_dir": keyfile_dir,
        "counts_dir": processed_data_dir / "guide_donor_bc0",
        "log_dir": processed_data_dir / "log",
        "oligo_pool_files": oligo_pool_files,
    }

# Defer config loading until main() to allow environment setup
CONFIG = None


# ============================================================================
# Sequence Constants
# ============================================================================

@dataclass(frozen=True)
class OligoConstants:
    """Constants for oligo sequence parsing and matching."""
    SPCAS9_CLONING_SITE: str = "GTTTGAAGAGC"
    PERFECT_LEADER: str = "GCAGGCGCGCC"
    # Slice positions within oligo_seq
    GUIDE_START: int = 20
    GUIDE_END: int = 40
    DONOR_START: int = 51
    DONOR_END: int = 180
    # Middle donor region for partial matching
    MIDDLE_DONOR_OFFSET: int = 30
    MIDDLE_DONOR_LENGTH: int = 60


OLIGO = OligoConstants()


# ============================================================================
# Lookup Table Builder
# ============================================================================

@dataclass
class OligoLookups:
    """Lookup dictionaries for oligo matching."""
    guide_donor_to_oligo: Dict[str, str]  # Exact match
    guide_to_oligos: Dict[str, List[str]]  # One guide -> many oligos
    donor_to_oligos: Dict[str, List[str]]  # Perfect donor -> many oligos
    middle_donor_to_oligos: Dict[str, List[str]]  # Partial donor -> many oligos


def load_oligo_pools(config: dict) -> pd.DataFrame:
    """
    Load and combine all oligo pool design files.

    Args:
        config: Configuration dictionary with oligo_pool_files

    Returns:
        DataFrame with oligo_name, oligo_seq, guide, donor, middle_donor columns
    """
    dfs = []

    for filepath in config["oligo_pool_files"]:
        if not filepath.exists():
            print(f"  WARNING: Oligo pool not found: {filepath}")
            continue

        print(f"  Loading: {filepath.name}")
        df = pd.read_csv(filepath, sep="\t")
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError("No oligo pool files found")

    combined = pd.concat(dfs, ignore_index=True)

    # Extract sequence components
    combined["guide"] = combined["oligo_seq"].str[OLIGO.GUIDE_START:OLIGO.GUIDE_END]
    combined["donor"] = combined["oligo_seq"].str[OLIGO.DONOR_START:OLIGO.DONOR_END].str.upper()
    combined["middle_donor"] = (
        combined["oligo_seq"]
        .str[OLIGO.DONOR_START + OLIGO.MIDDLE_DONOR_OFFSET:
             OLIGO.DONOR_START + OLIGO.MIDDLE_DONOR_OFFSET + OLIGO.MIDDLE_DONOR_LENGTH]
        .str.upper()
    )

    print(f"  Loaded {len(combined):,} oligos from {len(dfs)} files")
    return combined


def build_lookup_tables(oligo_df: pd.DataFrame) -> OligoLookups:
    """
    Build lookup dictionaries for efficient oligo matching.

    Uses vectorized pandas operations for ~10-50x speedup over iterrows.

    Args:
        oligo_df: DataFrame with oligo_name, guide, donor, middle_donor

    Returns:
        OligoLookups dataclass with all lookup dicts
    """
    # 1. Exact guide+donor match (vectorized)
    oligo_df["exact_key"] = (
        oligo_df["guide"] + OLIGO.SPCAS9_CLONING_SITE + oligo_df["donor"]
    )
    guide_donor_to_oligo = dict(zip(oligo_df["exact_key"], oligo_df["oligo_name"]))

    # 2. Guide -> oligos (vectorized groupby)
    guide_to_oligos = (
        oligo_df.groupby("guide")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # 3. Donor -> oligos (vectorized groupby)
    donor_to_oligos = (
        oligo_df.groupby("donor")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # 4. Middle donor -> oligos (vectorized groupby)
    # For simplicity, use the pre-computed middle_donor column
    middle_donor_to_oligos = (
        oligo_df.groupby("middle_donor")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    # Note: The original sliding window approach is unnecessary since
    # we only need to match against the middle region at a fixed position.
    # The middle_donor column already contains the correct 60bp region.

    return OligoLookups(
        guide_donor_to_oligo=guide_donor_to_oligo,
        guide_to_oligos=guide_to_oligos,
        donor_to_oligos=donor_to_oligos,
        middle_donor_to_oligos=middle_donor_to_oligos,
    )


# ============================================================================
# Oligo Matching
# ============================================================================

@dataclass
class MatchResult:
    """Result of oligo matching attempt."""
    oligo_name: str
    perfect_guide: bool
    perfect_donor: bool
    perfect_middle_donor: bool
    unambiguous_guide: bool
    unambiguous_donor: bool
    # Additional fields for ambiguous match logging
    n_guide_matches: int = 0
    n_donor_matches: int = 0
    ambiguous: bool = False


def match_to_oligo(guide: str, donor: str, lookups: OligoLookups) -> MatchResult:
    """
    Attempt to match guide and donor sequences to a designed oligo.

    Matching strategy (in order of preference):
    1. Exact guide+donor match
    2. Intersection of guide-matching and donor-matching oligos
    3. Unambiguous donor match (single oligo)
    4. Unambiguous guide match (single oligo)

    Args:
        guide: 20bp guide sequence
        donor: Donor sequence
        lookups: OligoLookups with all lookup dicts

    Returns:
        MatchResult with oligo_name and match quality flags
    """
    # Check for exact match first
    exact_key = guide + OLIGO.SPCAS9_CLONING_SITE + donor
    if exact_key in lookups.guide_donor_to_oligo:
        return MatchResult(
            oligo_name=lookups.guide_donor_to_oligo[exact_key],
            perfect_guide=True,
            perfect_donor=True,
            perfect_middle_donor=True,
            unambiguous_guide=True,
            unambiguous_donor=True,
            n_guide_matches=1,
            n_donor_matches=1,
            ambiguous=False,
        )

    # Check individual matches
    perfect_guide = guide in lookups.guide_to_oligos
    perfect_donor = donor in lookups.donor_to_oligos

    # Extract middle donor region
    middle_start = OLIGO.MIDDLE_DONOR_OFFSET
    middle_end = OLIGO.MIDDLE_DONOR_OFFSET + OLIGO.MIDDLE_DONOR_LENGTH
    middle_donor = donor[middle_start:middle_end] if len(donor) >= middle_end else ""
    perfect_middle_donor = middle_donor in lookups.middle_donor_to_oligos

    # Get candidate oligos
    guide_oligos = lookups.guide_to_oligos.get(guide, []) if perfect_guide else []

    if perfect_donor:
        donor_oligos = lookups.donor_to_oligos[donor]
    elif perfect_middle_donor:
        donor_oligos = lookups.middle_donor_to_oligos[middle_donor]
    else:
        donor_oligos = []

    n_guide_matches = len(guide_oligos)
    n_donor_matches = len(donor_oligos)
    unambiguous_guide = n_guide_matches == 1
    unambiguous_donor = n_donor_matches == 1

    # Try to find match
    oligo_name = "NA"
    ambiguous = False

    if guide_oligos and donor_oligos:
        # Find intersection
        intersection = [o for o in guide_oligos if o in donor_oligos]
        if len(intersection) == 1:
            oligo_name = intersection[0]
            unambiguous_guide = True
            unambiguous_donor = True
        elif len(intersection) > 1:
            # Multiple oligos match both guide AND donor - true ambiguity
            ambiguous = True
    elif unambiguous_donor:
        oligo_name = donor_oligos[0]
    elif unambiguous_guide:
        oligo_name = guide_oligos[0]
    elif n_guide_matches > 1 or n_donor_matches > 1:
        # Multiple matches but couldn't resolve
        ambiguous = True

    return MatchResult(
        oligo_name=oligo_name,
        perfect_guide=perfect_guide,
        perfect_donor=perfect_donor,
        perfect_middle_donor=perfect_middle_donor,
        unambiguous_guide=unambiguous_guide,
        unambiguous_donor=unambiguous_donor,
        n_guide_matches=n_guide_matches,
        n_donor_matches=n_donor_matches,
        ambiguous=ambiguous,
    )


# ============================================================================
# File Processing
# ============================================================================

def _match_row(row: pd.Series, lookups: OligoLookups) -> pd.Series:
    """
    Match a single row to designed oligo. Used for vectorized apply.

    Args:
        row: DataFrame row with guide, donor, leader columns
        lookups: OligoLookups for matching

    Returns:
        Series with match results
    """
    guide = row["guide"]
    donor = row["donor"]
    leader = row["leader"]

    if guide == "NA" or pd.isna(guide):
        return pd.Series({
            "scaffold": "NA",
            "perfect_leader": False,
            "perfect_guide": False,
            "perfect_donor": False,
            "perfect_middle_donor": False,
            "unambiguous_guide": False,
            "unambiguous_donor": False,
            "oligo_name": "NA",
            "ambiguous": False,
        })

    match_result = match_to_oligo(guide, donor, lookups)
    return pd.Series({
        "scaffold": "WT_SpCas9",
        "perfect_leader": leader == OLIGO.PERFECT_LEADER,
        "perfect_guide": match_result.perfect_guide,
        "perfect_donor": match_result.perfect_donor,
        "perfect_middle_donor": match_result.perfect_middle_donor,
        "unambiguous_guide": match_result.unambiguous_guide,
        "unambiguous_donor": match_result.unambiguous_donor,
        "oligo_name": match_result.oligo_name,
        "ambiguous": match_result.ambiguous,
    })


def process_counts_file(
    input_path: Path,
    output_path: Path,
    lookups: OligoLookups,
    chunk_size: int = 100000,
) -> Dict[str, int]:
    """
    Process a single counts file and match to designed oligos.

    Uses chunked pandas reading for memory efficiency and better I/O performance.

    Args:
        input_path: Path to input counts TSV
        output_path: Path to output matched TSV
        lookups: OligoLookups for matching
        chunk_size: Number of rows to process at a time (default 100k)

    Returns:
        Dict with processing statistics
    """
    start_time = timeit.default_timer()

    stats = {"total": 0, "matched": 0, "perfect_match": 0, "ambiguous": 0}
    first_chunk = True

    # Read in chunks for memory efficiency
    for chunk in pd.read_csv(input_path, sep="\t", chunksize=chunk_size):
        stats["total"] += len(chunk)

        # Rename columns to use consistent names
        # Input columns from Step 02 (at end): leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment, seq
        # Get column names to identify positions
        cols = chunk.columns.tolist()

        # Sequence fields are at the end (assuming fixed positions from end)
        chunk["leader"] = chunk.iloc[:, -8]
        chunk["guide"] = chunk.iloc[:, -7]
        chunk["donor"] = chunk.iloc[:, -6]
        chunk["bc0"] = chunk.iloc[:, -4]
        chunk["donor_bc0_fragment"] = chunk.iloc[:, -2]

        # Get library_ID and counts from beginning
        chunk["library_ID_col"] = chunk.iloc[:, 0]
        chunk["counts"] = chunk.iloc[:, 1]

        # Apply matching function to each row
        match_results = chunk.apply(lambda row: _match_row(row, lookups), axis=1)

        # Build output DataFrame
        output_df = pd.DataFrame({
            "library_ID": chunk["library_ID_col"],
            "scaffold": match_results["scaffold"],
            "bc0": chunk["bc0"],
            "counts": chunk["counts"],
            "perfect_leader": match_results["perfect_leader"],
            "leader": chunk["leader"],
            "perfect_guide": match_results["perfect_guide"],
            "guide": chunk["guide"],
            "perfect_donor": match_results["perfect_donor"],
            "perfect_middle_donor_region": match_results["perfect_middle_donor"],
            "donor": chunk["donor"],
            "donor_bc0_fragment": chunk["donor_bc0_fragment"],
            "unambiguous_guide": match_results["unambiguous_guide"],
            "unambiguous_donor": match_results["unambiguous_donor"],
            "oligo_name": match_results["oligo_name"],
        })

        # Update statistics
        matched_mask = output_df["oligo_name"] != "NA"
        stats["matched"] += matched_mask.sum()

        perfect_mask = (
            matched_mask &
            output_df["perfect_guide"].astype(bool) &
            output_df["perfect_donor"].astype(bool)
        )
        stats["perfect_match"] += perfect_mask.sum()

        # Track ambiguous matches (multiple candidates but couldn't resolve)
        ambiguous_mask = match_results["ambiguous"].astype(bool)
        stats["ambiguous"] += ambiguous_mask.sum()

        # Write output
        output_df.to_csv(
            output_path,
            sep="\t",
            index=False,
            header=first_chunk,
            mode="w" if first_chunk else "a",
        )
        first_chunk = False

        # Progress update
        elapsed = (timeit.default_timer() - start_time) / 60
        print(f"    Processed {stats['total']:,} sequences ({elapsed:.1f} min)")

    stats["elapsed_seconds"] = timeit.default_timer() - start_time
    return stats


def process_library_wrapper(
    library_args: Tuple[str, Path, Path, OligoLookups],
) -> Dict:
    """
    Wrapper for parallel library processing.

    Args:
        library_args: Tuple of (library_ID, input_path, output_path, lookups)

    Returns:
        Dict with processing statistics including library_ID
    """
    library_ID, input_path, output_path, lookups = library_args

    if not input_path.exists():
        return {"library_ID": library_ID, "error": f"Input file not found: {input_path.name}"}

    try:
        stats = process_counts_file(input_path, output_path, lookups)
        stats["library_ID"] = library_ID
        return stats
    except Exception as e:
        return {"library_ID": library_ID, "error": str(e)}


def get_library_list(config: dict) -> List[str]:
    """
    Get list of libraries to process from sample keyfile.

    Args:
        config: Configuration dictionary with keyfile_dir

    Returns:
        List of library IDs (excluding 'empty')
    """
    # Find sample keyfile - try common names
    keyfile_dir = config["keyfile_dir"]
    sample_keyfile = None
    for name in ["sample_key.tsv", "step_1_sample_key.tsv"]:
        candidate = keyfile_dir / name
        if candidate.exists():
            sample_keyfile = candidate
            break

    # Also check env variable
    sample_keyfile_env = os.environ.get("SAMPLE_KEYFILE")
    if sample_keyfile_env:
        sample_keyfile = keyfile_dir / sample_keyfile_env

    if not sample_keyfile or not sample_keyfile.exists():
        raise FileNotFoundError(f"No sample keyfile found in {keyfile_dir}")

    libraries = set()

    with open(sample_keyfile) as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 2:
                library_ID = fields[1]
                if library_ID != "empty":
                    libraries.add(library_ID)

    return sorted(libraries)


# ============================================================================
# Main Processing
# ============================================================================

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Map guide-donor-bc0 sequences to designed oligos"
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=1,
        help="Number of parallel jobs (default: 1, sequential processing)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=100000,
        help="Chunk size for reading input files (default: 100000)",
    )
    return parser.parse_args()


def main():
    """Main processing pipeline."""
    args = parse_args()

    # Load configuration from environment variables
    config = get_config_from_env()
    counts_dir = config["counts_dir"]

    # Set thread limits for nested parallelism prevention
    n_jobs = min(args.n_jobs, cpu_count())
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    start_time = datetime.now()
    print("=" * 70)
    print("Step 03: Map Guide-Donor-BC0 to Designed Oligos")
    print(f"Started: {start_time}")
    print(f"Project: {config['project_dir']}")
    print(f"Parallel jobs: {n_jobs}")
    print("=" * 70)

    # Validate paths
    if not counts_dir.exists():
        print(f"ERROR: Input directory not found: {counts_dir}")
        return

    # Load oligo pools and build lookups
    print("\nLoading oligo pool design files...")
    oligo_df = load_oligo_pools(config)

    print("\nBuilding lookup tables (vectorized)...")
    lookup_start = timeit.default_timer()
    lookups = build_lookup_tables(oligo_df)
    lookup_time = timeit.default_timer() - lookup_start
    print(f"  Exact matches: {len(lookups.guide_donor_to_oligo):,}")
    print(f"  Unique guides: {len(lookups.guide_to_oligos):,}")
    print(f"  Unique donors: {len(lookups.donor_to_oligos):,}")
    print(f"  Middle donors: {len(lookups.middle_donor_to_oligos):,}")
    print(f"  Build time: {lookup_time:.2f}s")

    # Get library list
    libraries = get_library_list(config)
    print(f"\nLibraries to process: {', '.join(libraries)}")

    # Prepare library arguments
    library_args = []
    for library_ID in libraries:
        input_path = counts_dir / f"{library_ID}_guide_donor_bc0_counts.tsv"
        output_path = counts_dir / f"{library_ID}_matched_guide_donor_bc0_counts.tsv"
        library_args.append((library_ID, input_path, output_path, lookups))

    # Process libraries (sequential or parallel)
    print(f"\nProcessing {len(libraries)} libraries...")
    all_stats = []

    if n_jobs == 1:
        # Sequential processing
        for args_tuple in tqdm(library_args, desc="Libraries"):
            library_ID = args_tuple[0]
            print(f"\n  Processing: {library_ID}")
            stats = process_library_wrapper(args_tuple)
            all_stats.append(stats)

            if "error" in stats:
                print(f"    ERROR: {stats['error']}")
            else:
                match_rate = 100 * stats["matched"] / stats["total"] if stats["total"] > 0 else 0
                perfect_rate = 100 * stats["perfect_match"] / stats["total"] if stats["total"] > 0 else 0
                print(f"    Total sequences: {stats['total']:,}")
                print(f"    Matched to oligo: {stats['matched']:,} ({match_rate:.1f}%)")
                print(f"    Perfect matches: {stats['perfect_match']:,} ({perfect_rate:.1f}%)")
                print(f"    Time: {stats['elapsed_seconds']:.1f}s")
    else:
        # Parallel processing
        print(f"  Using {n_jobs} parallel workers")
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = {
                executor.submit(process_library_wrapper, args_tuple): args_tuple[0]
                for args_tuple in library_args
            }

            for future in tqdm(as_completed(futures), total=len(futures), desc="Libraries"):
                library_ID = futures[future]
                try:
                    stats = future.result()
                    all_stats.append(stats)

                    if "error" in stats:
                        print(f"\n  {library_ID}: ERROR - {stats['error']}")
                    else:
                        match_rate = 100 * stats["matched"] / stats["total"] if stats["total"] > 0 else 0
                        print(f"\n  {library_ID}: {stats['total']:,} seqs, "
                              f"{match_rate:.1f}% matched, {stats['elapsed_seconds']:.1f}s")
                except Exception as e:
                    print(f"\n  {library_ID}: EXCEPTION - {e}")
                    all_stats.append({"library_ID": library_ID, "error": str(e)})

    # Summary
    total_elapsed = datetime.now() - start_time
    successful = [s for s in all_stats if "error" not in s]
    failed = [s for s in all_stats if "error" in s]

    print(f"\n{'=' * 70}")
    print("Step C Complete")
    print(f"  Libraries processed: {len(successful)}")
    if failed:
        print(f"  Libraries failed: {len(failed)}")
        for s in failed:
            print(f"    - {s['library_ID']}: {s['error']}")

    if successful:
        total_seqs = sum(s["total"] for s in successful)
        total_matched = sum(s["matched"] for s in successful)
        total_perfect = sum(s["perfect_match"] for s in successful)
        total_ambiguous = sum(s.get("ambiguous", 0) for s in successful)
        print(f"\n  Total sequences: {total_seqs:,}")
        print(f"  Total matched: {total_matched:,} ({100*total_matched/total_seqs:.1f}%)")
        print(f"  Perfect matches: {total_perfect:,} ({100*total_perfect/total_seqs:.1f}%)")
        if total_ambiguous > 0:
            print(f"  Ambiguous matches: {total_ambiguous:,} ({100*total_ambiguous/total_seqs:.2f}%)")

    print(f"\n  Output directory: {counts_dir}")
    print(f"  Total elapsed time: {total_elapsed}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
