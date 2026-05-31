#!/usr/bin/env python3
"""
Step 05: Assess Plasmid Library Representation

Analyzes the representation of designed oligos across plasmid libraries after
BC0 purity filtering. Aggregates counts at the oligo name level and identifies
missing or underrepresented oligos.

Key analyses:
- Oligo-level aggregation from guide-bc0 level data
- Missing oligo identification per library
- GC content bias assessment
- SPS-based filtering to correct top pool

Input: Purity-filtered guide-donor-bc0 counts from Step 04
Output: Oligo-level representation summaries per library

Usage:
    python 05_assess_plasmid_library_representation.py
"""

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import sys

import pandas as pd
from tqdm import tqdm

# ============================================================================
# Configuration - Read from environment variables set by SLURM script
# ============================================================================

import os

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

    purity_dir = processed_data_dir / "guide_donor_bc0_purity_filtered"

    return {
        "project_dir": project_dir,
        "processed_data_dir": processed_data_dir,
        "keyfile_dir": keyfile_dir,
        "purity_dir": purity_dir,
        "summary_dir": purity_dir / "summary_plots",
        "oligo_pool_files": oligo_pool_files,
    }

# File patterns
INPUT_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


@dataclass
class OligoConstants:
    """Constants for oligo sequence parsing."""
    DONOR_START: int = 51
    DONOR_END: int = 180
    SPS_START: int = 180
    SPS_LENGTH: int = 20


@dataclass
class LibraryConfig:
    """Configuration for library processing."""
    # Library to SPS (SubPool Sequence) mapping
    # SPS identifies which subpool of oligos is in each library
    library_sps_map: Dict[str, str] = field(default_factory=lambda: {
        "V629": "CTTGACTCGACAGTGAGAGC",  # 2024 Twist pool
        "V631": "CGTATGGAGATGACGTCAGG",  # 2024 Twist pool
        "V634": "CGATAGACGACTGGACAGCA",  # 2021 Twist pool
    })

    @property
    def expected_libraries(self) -> List[str]:
        return list(self.library_sps_map.keys())


# Global constants
OLIGO = OligoConstants()
LIB_CONFIG = LibraryConfig()


def setup_directories(env_config: dict) -> None:
    """Create necessary output directories."""
    env_config["summary_dir"].mkdir(parents=True, exist_ok=True)


def gc_content(seq: str) -> float:
    """
    Calculate the GC content of a DNA sequence.

    Args:
        seq: DNA sequence string

    Returns:
        GC content as fraction (0.0 to 1.0)
    """
    if not isinstance(seq, str) or len(seq) == 0:
        return 0.0
    seq_upper = seq.upper()
    gc_count = seq_upper.count("G") + seq_upper.count("C")
    return gc_count / len(seq_upper)


def load_purity_files(input_dir: Path, suffix: str) -> List[pd.DataFrame]:
    """
    Load all purity-filtered files matching the suffix.

    Args:
        input_dir: Directory containing purity files
        suffix: File suffix to match

    Returns:
        List of DataFrames with library column added
    """
    files = sorted(input_dir.glob(f"*{suffix}"))
    df_list = []

    for filepath in files:
        library = filepath.stem.replace(suffix.replace(".tsv", ""), "")
        print(f"  Loading: {filepath.name} (library: {library})")
        df = pd.read_csv(filepath, sep="\t")
        df["library"] = library
        df_list.append(df)

    return df_list


def load_designed_oligos(env_config: dict) -> pd.DataFrame:
    """
    Load and process all designed oligo pool files.

    Args:
        env_config: Configuration dictionary with oligo_pool_files

    Returns:
        Combined DataFrame with extracted sequence features
    """
    dfs = []

    for filepath in env_config["oligo_pool_files"]:
        if not filepath.exists():
            print(f"  WARNING: Oligo pool not found: {filepath}")
            continue

        print(f"  Loading: {filepath.name}")
        df = pd.read_csv(filepath, sep="\t")
        dfs.append(df)

    if not dfs:
        raise FileNotFoundError("No oligo pool files found")

    combined = pd.concat(dfs, ignore_index=True)

    # Extract sequence features
    combined["designed_donor"] = (
        combined["oligo_seq"]
        .str[OLIGO.DONOR_START:OLIGO.DONOR_END]
        .str.upper()
    )
    combined["SPS"] = (
        combined["oligo_seq"]
        .str[OLIGO.SPS_START:]
        .str.upper()
    )
    combined["oligo_seq"] = combined["oligo_seq"].str.upper()
    combined["gc_content"] = combined["designed_donor"].apply(gc_content)

    return combined


def process_donor_gc_content(df: pd.DataFrame, env_config: dict) -> pd.DataFrame:
    """
    Calculate GC content for donor sequences and report issues.

    Args:
        df: DataFrame with donor column
        env_config: Configuration dictionary with summary_dir

    Returns:
        DataFrame with gc_content column added
    """
    if "donor" not in df.columns:
        print("  Column 'donor' not found; skipping GC content calculation")
        return df

    df["gc_content"] = df["donor"].apply(gc_content)

    # Report non-string donors
    non_string_mask = ~df["donor"].apply(lambda x: isinstance(x, str))
    non_string_count = non_string_mask.sum()

    if non_string_count > 0:
        print(f"  WARNING: Found {non_string_count} non-string donor entries")
        non_string_df = df[non_string_mask]
        non_string_file = env_config["summary_dir"] / "non_string_donors.tsv"
        non_string_df.to_csv(non_string_file, sep="\t", index=False)

    return df


def aggregate_by_oligo_name(df: pd.DataFrame) -> pd.DataFrame:
    """
    Group counts at oligo_name level for each library.

    Args:
        df: DataFrame with library, oligo_name, grouped_guide_bc0_counts, bc0

    Returns:
        Aggregated DataFrame with total counts and unique BC0 count per oligo
    """
    required_cols = {"library", "oligo_name", "grouped_guide_bc0_counts", "bc0"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    return df.groupby(["library", "oligo_name"], as_index=False).agg(
        grouped_guide_bc0_counts=("grouped_guide_bc0_counts", "sum"),
        unique_bc0=("bc0", "nunique"),
    )


def add_missing_oligos(
    df: pd.DataFrame,
    designed_oligos: pd.DataFrame,
) -> pd.DataFrame:
    """
    Add missing oligos per library based on expected SPS with zero counts.

    Args:
        df: Aggregated oligo counts DataFrame
        designed_oligos: DataFrame with all designed oligos

    Returns:
        DataFrame with missing oligos added (zero counts)
    """
    df_parts = [df]

    for library, sps in LIB_CONFIG.library_sps_map.items():
        # Get oligos expected in this library (based on SPS)
        expected_oligos = designed_oligos[designed_oligos["SPS"] == sps]

        # Get oligos already present
        existing_oligos = df.loc[df["library"] == library, "oligo_name"]
        missing = expected_oligos[~expected_oligos["oligo_name"].isin(existing_oligos)]

        if missing.empty:
            print(f"  Library {library}: All expected oligos present")
            continue

        # Check for sequence overlap (same sequence, different name)
        if "oligo_seq" in df.columns:
            existing_seqs = df.loc[df["library"] == library, "oligo_seq"]
            overlap = missing[missing["oligo_seq"].isin(existing_seqs)]

            if not overlap.empty:
                print(f"  WARNING: Library {library} has {len(overlap)} missing oligos "
                      f"sharing oligo_seq with existing oligos")
                missing = missing[~missing.index.isin(overlap.index)]

        if not missing.empty:
            print(f"  Library {library}: Adding {len(missing)} missing oligos with zero counts")
            to_append = missing.assign(
                library=library,
                grouped_guide_bc0_counts=0,
                unique_bc0=0,
            )[["library", "oligo_name", "grouped_guide_bc0_counts", "unique_bc0",
               "oligo_seq", "designed_donor", "gc_content", "SPS"]]
            df_parts.append(to_append)

    return pd.concat(df_parts, ignore_index=True)


def filter_by_top_sps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter to only include the top SPS per library by summed counts.

    This handles cases where a library may have contamination from multiple
    SPS pools - we keep only the dominant one.

    Args:
        df: DataFrame with library, SPS, grouped_oligo_name_counts

    Returns:
        Filtered DataFrame with only top SPS per library
    """
    # Sum counts by library and SPS
    sps_totals = (
        df.groupby(["library", "SPS"], as_index=False)
        .agg(
            total_counts=("grouped_oligo_name_counts", "sum"),
            unique_bc0=("unique_bc0", "sum"),
        )
        .sort_values(["library", "total_counts"], ascending=[True, False])
    )

    # Keep only top SPS per library
    top_sps = sps_totals.drop_duplicates(subset=["library"], keep="first")

    print("\n  Top SPS per library:")
    for _, row in top_sps.iterrows():
        print(f"    {row['library']}: {row['SPS']} ({row['total_counts']:,} counts)")

    # Filter to top SPS
    return df.merge(
        top_sps[["library", "SPS"]],
        on=["library", "SPS"],
        how="inner",
    )


def save_outputs(df: pd.DataFrame, env_config: dict) -> None:
    """
    Save various output files.

    Args:
        df: Processed oligo-level DataFrame
        env_config: Configuration dictionary with summary_dir
    """
    summary_dir = env_config["summary_dir"]

    # All libraries combined
    all_file = summary_dir / "all_libraries_oligo_name_filtered.tsv"
    df.to_csv(all_file, sep="\t", index=False)
    print(f"\n  Saved: {all_file.name}")

    # Per library
    for library in df["library"].unique():
        lib_df = df[df["library"] == library]
        lib_file = summary_dir / f"{library}_oligo_name_filtered.tsv"
        lib_df.to_csv(lib_file, sep="\t", index=False)
        print(f"  Saved: {lib_file.name} ({len(lib_df):,} oligos)")

    # Missing oligos report
    missing = df[df["grouped_oligo_name_counts"] == 0]
    if not missing.empty:
        missing_file = summary_dir / "missing_oligos_with_zero_counts.tsv"
        missing.to_csv(missing_file, sep="\t", index=False)
        print(f"  Saved: {missing_file.name} ({len(missing):,} missing oligos)")


def main():
    """Main processing pipeline."""
    # Load configuration from environment
    env_config = get_config_from_env()
    purity_dir = env_config["purity_dir"]
    summary_dir = env_config["summary_dir"]

    start_time = datetime.now()
    print("=" * 70)
    print("Step 05: Assess Plasmid Library Representation")
    print(f"Started: {start_time}")
    print(f"Project: {env_config['project_dir']}")
    print("=" * 70)

    # Setup
    setup_directories(env_config)

    # Check input directory
    if not purity_dir.exists():
        print(f"ERROR: Input directory not found: {purity_dir}")
        sys.exit(1)

    # Load purity-filtered files
    print("\nLoading purity-filtered files...")
    df_list = load_purity_files(purity_dir, INPUT_SUFFIX)

    if not df_list:
        print("ERROR: No input files found")
        sys.exit(1)

    # Combine and process
    print(f"\nCombining {len(df_list)} files...")
    combined_df = pd.concat(df_list, ignore_index=True)
    print(f"  Total rows: {len(combined_df):,}")

    # Process GC content
    print("\nProcessing GC content...")
    combined_df = process_donor_gc_content(combined_df, env_config)

    # Load designed oligos
    print("\nLoading designed oligo pools...")
    designed_oligos = load_designed_oligos(env_config)
    print(f"  Total designed oligos: {len(designed_oligos):,}")

    # Aggregate by oligo name
    print("\nAggregating by oligo name...")
    oligo_df = aggregate_by_oligo_name(combined_df)
    print(f"  Unique library-oligo combinations: {len(oligo_df):,}")

    # Merge with designed oligos for annotations
    oligo_df = oligo_df.merge(
        designed_oligos[["oligo_name", "oligo_seq", "designed_donor", "gc_content", "SPS"]],
        on="oligo_name",
        how="left",
    )

    # Ensure uppercase sequences
    for col in ["SPS", "oligo_seq"]:
        if col in oligo_df.columns:
            oligo_df[col] = oligo_df[col].str.upper()

    # Add missing oligos
    print("\nChecking for missing oligos...")
    oligo_df = add_missing_oligos(oligo_df, designed_oligos)

    # Rename counts column for clarity
    oligo_df.rename(
        columns={"grouped_guide_bc0_counts": "grouped_oligo_name_counts"},
        inplace=True,
    )

    # Filter by top SPS per library
    print("\nFiltering to top SPS per library...")
    filtered_df = filter_by_top_sps(oligo_df)

    # Save outputs
    print("\nSaving outputs...")
    save_outputs(filtered_df, env_config)

    # Summary statistics
    print("\n" + "-" * 50)
    print("Library Summary:")
    for library in LIB_CONFIG.expected_libraries:
        lib_df = filtered_df[filtered_df["library"] == library]
        if lib_df.empty:
            print(f"  {library}: No data")
            continue

        n_oligos = len(lib_df)
        n_detected = (lib_df["grouped_oligo_name_counts"] > 0).sum()
        total_counts = lib_df["grouped_oligo_name_counts"].sum()
        mean_bc0 = lib_df["unique_bc0"].mean()

        print(f"  {library}:")
        print(f"    Oligos: {n_detected:,} / {n_oligos:,} detected")
        print(f"    Total counts: {total_counts:,}")
        print(f"    Mean BC0s per oligo: {mean_bc0:.1f}")

    # Final summary
    elapsed = datetime.now() - start_time
    print(f"\n{'=' * 70}")
    print("Step E Complete")
    print(f"  Output directory: {summary_dir}")
    print(f"  Elapsed time: {elapsed}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
