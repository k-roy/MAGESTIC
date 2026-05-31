#!/usr/bin/env python3
"""
Step 03: Map Guide-Donor-BC0 Sequences to Designed Oligos (OPTIMIZED)

Maps observed guide-donor-bc0 sequences to the designed oligo pool using
a multi-level matching strategy:

1. Exact match: guide + RE_site + donor matches designed sequence
2. Guide + donor intersection: both guide and donor match same oligo
3. Unambiguous donor: single oligo matches donor (perfect or middle region)
4. Unambiguous guide: single oligo matches guide

OPTIMIZATION (2026-03-22):
- Two-pass approach: vectorized exact matching handles ~85% of rows
- Parallel fallback matching for remaining rows using joblib
- Increased chunk size (500k) for better I/O efficiency
- Progress tracking within library processing

Input: Per-library guide_donor_bc0_counts.tsv from Step 02
Output: Per-library matched_guide_donor_bc0_counts.tsv with oligo_name column

Usage:
    python 03_map_guide_donor_bc0_to_designed_oligos_optimized.py [--n-jobs N]

    # Test mode: process single library
    python 03_map_guide_donor_bc0_to_designed_oligos_optimized.py --test-library V546
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

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
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
    middle_donor_to_oligos = (
        oligo_df.groupby("middle_donor")["oligo_name"]
        .apply(lambda x: x.drop_duplicates().tolist())
        .to_dict()
    )

    return OligoLookups(
        guide_donor_to_oligo=guide_donor_to_oligo,
        guide_to_oligos=guide_to_oligos,
        donor_to_oligos=donor_to_oligos,
        middle_donor_to_oligos=middle_donor_to_oligos,
    )


# ============================================================================
# Oligo Matching (for fallback - non-exact matches)
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
    n_guide_matches: int = 0
    n_donor_matches: int = 0
    ambiguous: bool = False


def match_to_oligo(guide: str, donor: str, lookups: OligoLookups) -> MatchResult:
    """
    Attempt to match guide and donor sequences to a designed oligo.

    Matching strategy (in order of preference):
    1. Exact guide+donor match (handled by vectorized pass - should not reach here)
    2. Intersection of guide-matching and donor-matching oligos
    3. Unambiguous donor match (single oligo)
    4. Unambiguous guide match (single oligo)
    """
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
            ambiguous = True
    elif unambiguous_donor:
        oligo_name = donor_oligos[0]
    elif unambiguous_guide:
        oligo_name = guide_oligos[0]
    elif n_guide_matches > 1 or n_donor_matches > 1:
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
# Optimized File Processing
# ============================================================================

def match_fallback_batch(df: pd.DataFrame, lookups: OligoLookups) -> pd.DataFrame:
    """
    Apply multi-level fallback matching to a batch of unmatched rows.

    Args:
        df: DataFrame with guide, donor, leader columns
        lookups: OligoLookups for matching

    Returns:
        DataFrame with match results, indexed to match input
    """
    results = []
    for idx, row in df.iterrows():
        guide = row["guide"]
        donor = row["donor"]
        leader = row["leader"]

        if guide == "NA" or pd.isna(guide):
            results.append({
                "oligo_name": "NA",
                "perfect_guide": False,
                "perfect_donor": False,
                "perfect_middle_donor": False,
                "unambiguous_guide": False,
                "unambiguous_donor": False,
                "perfect_leader": False,
                "scaffold": "NA",
            })
        else:
            match = match_to_oligo(guide, donor, lookups)
            results.append({
                "oligo_name": match.oligo_name,
                "perfect_guide": match.perfect_guide,
                "perfect_donor": match.perfect_donor,
                "perfect_middle_donor": match.perfect_middle_donor,
                "unambiguous_guide": match.unambiguous_guide,
                "unambiguous_donor": match.unambiguous_donor,
                "perfect_leader": leader == OLIGO.PERFECT_LEADER,
                "scaffold": "WT_SpCas9" if match.oligo_name != "NA" else "NA",
            })

    return pd.DataFrame(results, index=df.index)


def process_counts_file_vectorized(
    input_path: Path,
    output_path: Path,
    lookups: OligoLookups,
    chunk_size: int = 500000,
    n_workers: int = 8,
) -> Dict[str, int]:
    """
    Process a single counts file with vectorized exact matching and parallel fallback.

    Two-pass approach:
    - Pass 1: Vectorized exact-key lookup (handles ~85% of rows)
    - Pass 2: Parallel fallback matching for remaining rows

    Args:
        input_path: Path to input counts TSV
        output_path: Path to output matched TSV
        lookups: OligoLookups for matching
        chunk_size: Number of rows to process at a time (default 500k)
        n_workers: Number of parallel workers for fallback matching

    Returns:
        Dict with processing statistics
    """
    start_time = timeit.default_timer()

    stats = {
        "total": 0,
        "matched": 0,
        "exact_match": 0,
        "fallback_match": 0,
        "unmatched": 0,
    }
    first_chunk = True
    chunk_num = 0

    # Read in chunks for memory efficiency
    for chunk in pd.read_csv(input_path, sep="\t", chunksize=chunk_size):
        chunk_num += 1
        chunk_start = timeit.default_timer()
        stats["total"] += len(chunk)

        # Parse column positions (same as original)
        # Input columns from Step 02 (at end): leader, guide, donor, SPS, bc0, msd, donor_bc0_fragment, seq
        chunk["leader"] = chunk.iloc[:, -8]
        chunk["guide"] = chunk.iloc[:, -7]
        chunk["donor"] = chunk.iloc[:, -6]
        chunk["bc0"] = chunk.iloc[:, -4]
        chunk["donor_bc0_fragment"] = chunk.iloc[:, -2]
        chunk["library_ID_col"] = chunk.iloc[:, 0]
        chunk["counts"] = chunk.iloc[:, 1]

        # ========================================
        # PASS 1: Vectorized exact-key matching
        # ========================================
        chunk["exact_key"] = chunk["guide"] + OLIGO.SPCAS9_CLONING_SITE + chunk["donor"]
        chunk["oligo_name"] = chunk["exact_key"].map(lookups.guide_donor_to_oligo)

        # Identify exact matches
        exact_matched = chunk["oligo_name"].notna()
        n_exact = exact_matched.sum()
        stats["exact_match"] += n_exact

        # Set flags for exact matches (all True for exact matches)
        chunk.loc[exact_matched, "perfect_guide"] = True
        chunk.loc[exact_matched, "perfect_donor"] = True
        chunk.loc[exact_matched, "perfect_middle_donor"] = True
        chunk.loc[exact_matched, "unambiguous_guide"] = True
        chunk.loc[exact_matched, "unambiguous_donor"] = True
        chunk.loc[exact_matched, "perfect_leader"] = chunk.loc[exact_matched, "leader"] == OLIGO.PERFECT_LEADER
        chunk.loc[exact_matched, "scaffold"] = "WT_SpCas9"

        # ========================================
        # PASS 2: Fallback matching for unmatched
        # ========================================
        unmatched_mask = ~exact_matched
        n_unmatched = unmatched_mask.sum()

        if n_unmatched > 0:
            unmatched_df = chunk.loc[unmatched_mask, ["guide", "donor", "leader"]].copy()

            # Determine parallelization strategy
            n_sub = min(n_workers, max(1, n_unmatched // 10000))

            if n_sub > 1 and n_unmatched > 10000:
                # Split into sub-chunks for parallel processing
                sub_indices = np.array_split(unmatched_df.index, n_sub)
                sub_dfs = [unmatched_df.loc[idx] for idx in sub_indices]

                results = Parallel(n_jobs=n_sub, backend="loky")(
                    delayed(match_fallback_batch)(sub_df, lookups)
                    for sub_df in sub_dfs
                )
                fallback_df = pd.concat(results)
            else:
                # Process sequentially for small batches
                fallback_df = match_fallback_batch(unmatched_df, lookups)

            # Merge fallback results back into main chunk
            for col in ["oligo_name", "perfect_guide", "perfect_donor",
                       "perfect_middle_donor", "unambiguous_guide",
                       "unambiguous_donor", "perfect_leader", "scaffold"]:
                chunk.loc[unmatched_mask, col] = fallback_df[col].values

            # Count fallback matches
            fallback_matched = chunk.loc[unmatched_mask, "oligo_name"] != "NA"
            stats["fallback_match"] += fallback_matched.sum()
            stats["unmatched"] += (~fallback_matched).sum()

        stats["matched"] = stats["exact_match"] + stats["fallback_match"]

        # Fill any remaining NA values (using infer_objects to avoid FutureWarning)
        chunk["oligo_name"] = chunk["oligo_name"].fillna("NA")
        chunk["scaffold"] = chunk["scaffold"].fillna("NA")
        for col in ["perfect_leader", "perfect_guide", "perfect_donor",
                    "perfect_middle_donor", "unambiguous_guide", "unambiguous_donor"]:
            chunk[col] = chunk[col].fillna(False).infer_objects(copy=False)

        # Build output DataFrame with correct column names
        output_df = pd.DataFrame({
            "library_ID": chunk["library_ID_col"],
            "scaffold": chunk["scaffold"],
            "bc0": chunk["bc0"],
            "counts": chunk["counts"],
            "perfect_leader": chunk["perfect_leader"],
            "leader": chunk["leader"],
            "perfect_guide": chunk["perfect_guide"],
            "guide": chunk["guide"],
            "perfect_donor": chunk["perfect_donor"],
            "perfect_middle_donor_region": chunk["perfect_middle_donor"],
            "donor": chunk["donor"],
            "donor_bc0_fragment": chunk["donor_bc0_fragment"],
            "unambiguous_guide": chunk["unambiguous_guide"],
            "unambiguous_donor": chunk["unambiguous_donor"],
            "oligo_name": chunk["oligo_name"],
        })

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
        chunk_time = timeit.default_timer() - chunk_start
        elapsed = timeit.default_timer() - start_time
        exact_pct = 100 * n_exact / len(chunk) if len(chunk) > 0 else 0
        print(f"    Chunk {chunk_num}: {len(chunk):,} rows, "
              f"{exact_pct:.1f}% exact match, {chunk_time:.1f}s "
              f"(total: {elapsed/60:.1f} min)")

    stats["elapsed_seconds"] = timeit.default_timer() - start_time
    return stats


def process_library_wrapper(
    library_args: Tuple[str, Path, Path, OligoLookups, int],
) -> Dict:
    """
    Wrapper for parallel library processing.

    Args:
        library_args: Tuple of (library_ID, input_path, output_path, lookups, n_workers)

    Returns:
        Dict with processing statistics including library_ID
    """
    library_ID, input_path, output_path, lookups, n_workers = library_args

    if not input_path.exists():
        return {"library_ID": library_ID, "error": f"Input file not found: {input_path.name}"}

    try:
        stats = process_counts_file_vectorized(
            input_path, output_path, lookups,
            chunk_size=500000, n_workers=n_workers
        )
        stats["library_ID"] = library_ID
        return stats
    except Exception as e:
        import traceback
        return {"library_ID": library_ID, "error": str(e), "traceback": traceback.format_exc()}


def get_library_list(config: dict) -> List[str]:
    """
    Get list of libraries to process from sample keyfile.
    """
    keyfile_dir = config["keyfile_dir"]
    sample_keyfile = None
    for name in ["sample_key.tsv", "step_1_sample_key.tsv", "step_1_sample_key_unified.tsv"]:
        candidate = keyfile_dir / name
        if candidate.exists():
            sample_keyfile = candidate
            break

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
        description="Map guide-donor-bc0 sequences to designed oligos (OPTIMIZED)"
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=8,
        help="Number of parallel workers for fallback matching (default: 8)",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500000,
        help="Chunk size for reading input files (default: 500000)",
    )
    parser.add_argument(
        "--test-library",
        type=str,
        default=None,
        help="Process only this library (for testing)",
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
    os.environ["LOKY_MAX_CPU_COUNT"] = str(n_jobs)

    start_time = datetime.now()
    print("=" * 70)
    print("Step 03: Map Guide-Donor-BC0 to Designed Oligos (OPTIMIZED)")
    print(f"Started: {start_time}")
    print(f"Project: {config['project_dir']}")
    print(f"Parallel workers: {n_jobs}")
    print(f"Chunk size: {args.chunk_size:,}")
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
    if args.test_library:
        libraries = [args.test_library]
        print(f"\nTEST MODE: Processing only {args.test_library}")
    else:
        libraries = get_library_list(config)
        print(f"\nLibraries to process: {', '.join(libraries)}")

    # Process libraries sequentially (parallelization is within each library)
    print(f"\nProcessing {len(libraries)} libraries...")
    all_stats = []

    for library_ID in libraries:
        input_path = counts_dir / f"{library_ID}_guide_donor_bc0_counts.tsv"
        output_path = counts_dir / f"{library_ID}_matched_guide_donor_bc0_counts.tsv"

        if args.test_library:
            # For test mode, use a different output file name
            output_path = counts_dir / f"{library_ID}_matched_guide_donor_bc0_counts_OPTIMIZED.tsv"

        print(f"\n  Processing: {library_ID}")
        print(f"    Input: {input_path.name}")
        print(f"    Output: {output_path.name}")

        stats = process_library_wrapper((library_ID, input_path, output_path, lookups, n_jobs))
        all_stats.append(stats)

        if "error" in stats:
            print(f"    ERROR: {stats['error']}")
            if "traceback" in stats:
                print(f"    {stats['traceback']}")
        else:
            total = stats["total"]
            exact_rate = 100 * stats["exact_match"] / total if total > 0 else 0
            fallback_rate = 100 * stats["fallback_match"] / total if total > 0 else 0
            match_rate = 100 * stats["matched"] / total if total > 0 else 0
            print(f"    Total sequences: {total:,}")
            print(f"    Exact matches: {stats['exact_match']:,} ({exact_rate:.1f}%)")
            print(f"    Fallback matches: {stats['fallback_match']:,} ({fallback_rate:.1f}%)")
            print(f"    Total matched: {stats['matched']:,} ({match_rate:.1f}%)")
            print(f"    Unmatched: {stats['unmatched']:,}")
            print(f"    Time: {stats['elapsed_seconds']:.1f}s ({stats['elapsed_seconds']/60:.1f} min)")

    # Summary
    total_elapsed = datetime.now() - start_time
    successful = [s for s in all_stats if "error" not in s]
    failed = [s for s in all_stats if "error" in s]

    print(f"\n{'=' * 70}")
    print("Step 03 Complete (OPTIMIZED)")
    print(f"  Libraries processed: {len(successful)}")
    if failed:
        print(f"  Libraries failed: {len(failed)}")
        for s in failed:
            print(f"    - {s['library_ID']}: {s['error']}")

    if successful:
        total_seqs = sum(s["total"] for s in successful)
        total_exact = sum(s["exact_match"] for s in successful)
        total_fallback = sum(s["fallback_match"] for s in successful)
        total_matched = sum(s["matched"] for s in successful)
        total_unmatched = sum(s["unmatched"] for s in successful)

        print(f"\n  Total sequences: {total_seqs:,}")
        print(f"  Exact matches: {total_exact:,} ({100*total_exact/total_seqs:.1f}%)")
        print(f"  Fallback matches: {total_fallback:,} ({100*total_fallback/total_seqs:.1f}%)")
        print(f"  Total matched: {total_matched:,} ({100*total_matched/total_seqs:.1f}%)")
        print(f"  Unmatched: {total_unmatched:,} ({100*total_unmatched/total_seqs:.1f}%)")

    print(f"\n  Output directory: {counts_dir}")
    print(f"  Total elapsed time: {total_elapsed}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
