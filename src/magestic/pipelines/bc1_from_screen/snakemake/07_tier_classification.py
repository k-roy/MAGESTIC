#!/usr/bin/env python3
"""
Step 04: Abundance Tier Classification for BC1 Barcodes

Classifies BC1 barcodes into abundance tiers based on t0 read counts within each library.

Abundance Tier Definitions:
- abundance_high (≥100 reads): Individual BC1 analysis in DESeq2
- abundance_medium (20-99 reads): Aggregate by oligo before DESeq2
- abundance_low (<20 reads): Excluded from analysis (low confidence)

NOTE: These "abundance tiers" are distinct from the "confidence_tier" column
in bc1_donor_bc0_pipeline which classifies oligo attribution confidence.

Usage:
    python 04_tier_classification.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --bc1-reference bc1_reference_table.tsv \\
        --sample-key sample_key.tsv \\
        --library yL406 \\
        --output-tier-mapping tier_mapping.tsv \\
        --output-stats tier_stats.tsv

MEMORY (2026-08-06): the long-format counts file is read in CHUNKS and the t0
totals are accumulated incrementally, instead of materialising all ~122M rows
as one DataFrame. The BC1 reference table is read with `usecols` because only
`bc1` + `oligo_name` are ever used. Both changes are numerically inert -- see
planning/pipeline_provenance_rebuild/15_step07_memory_fix.md for the
byte-identity proof (md5 of the tier mapping + tier stats, before vs after, on
the full yL406 input).

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

# Resolve TierConfig without requiring the magestic package to be installed.
# Verified 2026-08-06: NONE of the conda envs the pipeline is driven with
# (snakemake / ihw_env / pydeseq2) has magestic installed -- only the base
# anaconda3 python does -- so the bare `from magestic...` import made this step
# die with ModuleNotFoundError under exactly the interpreters Snakemake uses.
# Fall back to the sibling config.py, which imports nothing but stdlib + yaml.
try:
    from magestic.pipelines.bc1_from_screen.config import TierConfig
    from magestic.pipelines.bc1_from_screen.library_scope import resolve_library_samples
    from magestic.pipelines.bc1_from_screen.reference_projection import project_bc1_reference
except ImportError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from config import TierConfig
    from library_scope import resolve_library_samples
    from reference_projection import project_bc1_reference

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Default rows per read_csv chunk when streaming the long-format counts file.
# MEASURED on the full yL406 file (122,647,244 rows, SLURM 37782771):
#   chunksize 10,000,000 -> peak RSS 1.257 GiB, wall 0:58.01
#   chunksize  2,000,000 -> peak RSS 0.434 GiB, wall 0:58.66
# i.e. 2M is 2.9x lighter for no measurable time cost, and the output is
# byte-identical at both (and at 7 and 50,000) -- chunksize is a memory knob
# only.
DEFAULT_CHUNKSIZE = 2_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Classify BC1s into abundance tiers based on t0 counts"
    )
    parser.add_argument(
        "--bc1-counts", type=Path, required=True,
        help="Input BC1 counts long-format TSV"
    )
    parser.add_argument(
        "--bc1-reference", type=Path, required=True,
        help="BC1 reference table with oligo assignments"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition, generations, sample_name columns"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name (e.g., yL437, yL442, yL406)"
    )
    parser.add_argument(
        "--output-tier-mapping", type=Path, required=True,
        help="Output BC1 tier mapping TSV"
    )
    parser.add_argument(
        "--output-stats", type=Path, required=True,
        help="Output tier statistics TSV"
    )
    parser.add_argument(
        "--tier-1-min", type=int, default=100,
        help="Minimum t0 reads for abundance_high (default: 100)"
    )
    parser.add_argument(
        "--tier-2-min", type=int, default=20,
        help="Minimum t0 reads for abundance_medium (default: 20)"
    )
    parser.add_argument(
        "--t0-condition", type=str, default="YPD",
        help="Condition to identify t0 samples (default: YPD)"
    )
    parser.add_argument(
        "--screen-dates", type=str, nargs="+",
        help="Screen dates to filter for this library (optional)"
    )
    parser.add_argument(
        "--chunksize", type=int, default=DEFAULT_CHUNKSIZE,
        help=f"Rows per chunk when streaming the counts file "
             f"(default: {DEFAULT_CHUNKSIZE:,}). Memory only -- the result is "
             f"independent of this value."
    )
    return parser.parse_args()


def identify_t0_samples(
    sample_key: pd.DataFrame,
    t0_condition: str = "YPD",
    screen_dates: list = None,
    sample_col: str = "sample_name",
) -> list:
    """
    Identify t0 sample names from the sample key.

    t0 samples are those with:
    - condition matching t0_condition (default: YPD)
    - generations == 0
    - Optionally filtered by screen_dates

    Extracted verbatim from calculate_t0_counts() (2026-08-06) so the streaming
    and in-memory paths share ONE definition of "which samples are t0".
    """
    # Compare generations numerically without truncation: 0.5 is not a baseline, and a
    # non-numeric value (NaN after coercion) never equals 0.
    generations_numeric = pd.to_numeric(sample_key["generations"], errors="coerce")
    is_t0_condition = sample_key["condition"] == t0_condition
    unparsed = is_t0_condition & generations_numeric.isna()
    if unparsed.any():
        logger.warning(
            f"  {int(unparsed.sum())} '{t0_condition}' row(s) have a non-numeric generations "
            f"value and are not treated as t0 (e.g. "
            f"{sample_key.loc[unparsed, sample_col].astype(str).tolist()[:5]})")
    t0_mask = is_t0_condition & (generations_numeric == 0)

    # Optionally filter by screen dates
    if screen_dates:
        t0_mask = t0_mask & (sample_key["screen_date"].isin(screen_dates))

    t0_samples = sample_key.loc[t0_mask, sample_col].unique().tolist()

    if not t0_samples:
        raise ValueError(
            f"No t0 samples found for condition='{t0_condition}', generations=0"
            + (f", screen_dates={screen_dates}" if screen_dates else "")
        )

    return t0_samples


def calculate_t0_counts(
    bc1_counts: pd.DataFrame,
    sample_key: pd.DataFrame,
    t0_condition: str = "YPD",
    screen_dates: list = None,
    bc1_col: str = "bc1",
    sample_col: str = "sample_name",
    count_col: str = "counts",
) -> pd.DataFrame:
    """
    Calculate total t0 counts per BC1 from an already-loaded counts DataFrame.

    Retained unchanged as the reference (in-memory) implementation and for any
    caller that already holds the frame; `main()` uses the streaming variant
    below, which is proven to produce a byte-identical result.

    Args:
        bc1_counts: Long-format BC1 counts (bc1, sample_name, counts)
        sample_key: Sample metadata with condition, generations, sample_name
        t0_condition: Condition value identifying t0 samples
        screen_dates: Optional list of screen dates to filter
        bc1_col: Column name for BC1 barcodes
        sample_col: Column name for sample identifiers
        count_col: Column name for read counts

    Returns:
        DataFrame with bc1 and t0_total columns
    """
    t0_samples = identify_t0_samples(sample_key, t0_condition, screen_dates, sample_col)

    logger.info(f"  Found {len(t0_samples)} t0 samples")

    # Filter BC1 counts to t0 samples
    t0_counts = bc1_counts[bc1_counts[sample_col].isin(t0_samples)]

    # Sum counts per BC1 across all t0 samples
    bc1_t0_totals = (
        t0_counts
        .groupby(bc1_col)[count_col]
        .sum()
        .reset_index()
        .rename(columns={count_col: "t0_total"})
    )

    return bc1_t0_totals


def calculate_t0_counts_streaming(
    bc1_counts_path: Path,
    sample_key: pd.DataFrame,
    t0_condition: str = "YPD",
    screen_dates: list = None,
    bc1_col: str = "bc1",
    sample_col: str = "sample_name",
    count_col: str = "counts",
    chunksize: int = DEFAULT_CHUNKSIZE,
) -> pd.DataFrame:
    """
    Streaming equivalent of calculate_t0_counts() that never holds the whole
    long-format counts file in memory.

    Identical by construction:
      * the t0 sample set comes from the same identify_t0_samples()
      * each chunk is filtered with the same `.isin(t0_samples)` predicate
      * per-chunk `groupby(bc1)[counts].sum()` partial sums are added exactly
        (integer addition; no float accumulation is introduced)
      * the final frame is produced by the SAME
        `groupby(bc1)[counts].sum().reset_index().rename(...)` call, so the row
        ORDER (groupby sort=True) and the column dtype are whatever pandas
        would have produced from the one-shot path.

    Returns:
        DataFrame with bc1 and t0_total columns
    """
    t0_samples = identify_t0_samples(sample_key, t0_condition, screen_dates, sample_col)

    logger.info(f"  Found {len(t0_samples)} t0 samples")

    totals = {}
    n_rows = 0
    n_t0_rows = 0
    n_chunks = 0

    reader = pd.read_csv(
        bc1_counts_path,
        sep="\t",
        usecols=[bc1_col, sample_col, count_col],
        chunksize=chunksize,
    )
    with reader:
        for chunk in reader:
            n_chunks += 1
            n_rows += len(chunk)
            sub = chunk[chunk[sample_col].isin(t0_samples)]
            del chunk
            if len(sub):
                n_t0_rows += len(sub)
                partial = sub.groupby(bc1_col)[count_col].sum()
                for key, value in partial.items():
                    prev = totals.get(key)
                    totals[key] = value if prev is None else prev + value
                del partial
            del sub

    logger.info(
        f"  Streamed {n_rows:,} rows in {n_chunks} chunk(s) of {chunksize:,}; "
        f"{n_t0_rows:,} rows from t0 samples"
    )

    if totals:
        accumulated = pd.DataFrame(
            {bc1_col: list(totals.keys()), count_col: list(totals.values())}
        )
    else:
        accumulated = pd.DataFrame(
            {bc1_col: pd.Series(dtype=object), count_col: pd.Series(dtype="int64")}
        )
    del totals

    bc1_t0_totals = (
        accumulated
        .groupby(bc1_col)[count_col]
        .sum()
        .reset_index()
        .rename(columns={count_col: "t0_total"})
    )

    return bc1_t0_totals


def classify_bc1_tiers(
    bc1_t0_counts: pd.DataFrame,
    tier_config: TierConfig,
    bc1_col: str = "bc1",
    t0_col: str = "t0_total",
) -> pd.DataFrame:
    """
    Classify BC1s into tiers based on t0 read counts.

    Args:
        bc1_t0_counts: DataFrame with bc1 and t0_total columns
        tier_config: TierConfig with tier thresholds
        bc1_col: Column name for BC1 barcodes
        t0_col: Column name for t0 counts

    Returns:
        DataFrame with bc1, t0_total, and tier columns
    """
    result = bc1_t0_counts.copy()
    result["tier"] = result[t0_col].apply(tier_config.get_tier)
    return result


def create_tier_summary(tier_mapping: pd.DataFrame, library: str) -> pd.DataFrame:
    """
    Create summary statistics for abundance tier classification.

    Args:
        tier_mapping: DataFrame with tier assignments
        library: Library name

    Returns:
        DataFrame with tier summary statistics
    """
    summary_data = []

    for tier in ["abundance_high", "abundance_medium", "abundance_low"]:
        tier_data = tier_mapping[tier_mapping["tier"] == tier]
        summary_data.append({
            "library": library,
            "tier": tier,
            "n_bc1s": len(tier_data),
            "total_t0_reads": tier_data["t0_total"].sum() if len(tier_data) > 0 else 0,
            "mean_t0_reads": tier_data["t0_total"].mean() if len(tier_data) > 0 else 0,
            "median_t0_reads": tier_data["t0_total"].median() if len(tier_data) > 0 else 0,
        })

    return pd.DataFrame(summary_data)


def main():
    args = parse_args()

    logger.info(f"Step 04: Abundance tier classification for library {args.library}")
    logger.info(f"  Input BC1 counts: {args.bc1_counts}")
    logger.info(f"  BC1 reference: {args.bc1_reference}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  t0 condition: {args.t0_condition}")
    logger.info(f"  Tier thresholds: abundance_high >= {args.tier_1_min}, abundance_medium >= {args.tier_2_min}")

    # Detect count column name (handles both 'counts' and 'sample_counts')
    # from the HEADER ONLY -- the file itself is streamed below.
    counts_columns = list(pd.read_csv(args.bc1_counts, sep="\t", nrows=0).columns)
    if "counts" in counts_columns:
        count_col = "counts"
    elif "sample_counts" in counts_columns:
        count_col = "sample_counts"
    else:
        raise ValueError("BC1 counts file must have 'counts' or 'sample_counts' column")
    logger.info(f"  Using count column: {count_col}")

    # Load BC1 reference table.
    # Only `bc1` and `oligo_name` are ever used (see the merge below), so read
    # just those two columns -- the full table carries the ~235 bp `donor`
    # sequence for every BC1 and is ~2 orders of magnitude larger in memory.
    logger.info("Loading BC1 reference table...")
    ref_columns = list(pd.read_csv(args.bc1_reference, sep="\t", nrows=0).columns)
    if "bc1" in ref_columns and "oligo_name" in ref_columns:
        ref_table = pd.read_csv(args.bc1_reference, sep="\t", usecols=["bc1", "oligo_name"])
    else:
        ref_table = pd.read_csv(args.bc1_reference, sep="\t")
    logger.info(f"  Loaded {len(ref_table):,} BC1s")

    # Load sample key
    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} samples")

    # Verify required columns
    required_cols = ["sample_name", "condition", "generations"]
    missing = [c for c in required_cols if c not in sample_key.columns]
    if missing:
        raise ValueError(f"Sample key missing required columns: {missing}")

    # Scope the baseline to --library with the same rule step 08 uses. Without this a
    # multi-library key summed every library's t0 reads into this library's tiers.
    try:
        keep_samples = resolve_library_samples(sample_key, args.library)
    except ValueError as error:
        logger.error(str(error))
        return 1
    if keep_samples is not None:
        sample_key = sample_key[sample_key["sample_name"].astype(str).isin(set(keep_samples))]
        logger.info(f"  Sample key scoped to library {args.library}: {len(sample_key):,} samples")

    # Configure tier thresholds
    tier_config = TierConfig(
        tier_1_min_reads=args.tier_1_min,
        tier_2_min_reads=args.tier_2_min,
    )

    # Calculate t0 counts by STREAMING the long-format counts file.
    logger.info(f"Calculating t0 counts (condition={args.t0_condition}, generations=0)...")
    t0_counts = calculate_t0_counts_streaming(
        args.bc1_counts,
        sample_key,
        t0_condition=args.t0_condition,
        screen_dates=args.screen_dates,
        count_col=count_col,
        chunksize=args.chunksize,
    )
    logger.info(f"  Found {len(t0_counts):,} BC1s with t0 counts")

    # Classify into tiers
    logger.info("Classifying BC1s into tiers...")
    tier_mapping = classify_bc1_tiers(t0_counts, tier_config)

    # Add oligo information from the reference table, one row per BC1. A BC1 with
    # several reference oligos abstains (oligo_name empty, candidates kept) instead of
    # multiplying its tier row, which double-counted its reads in the tier summary.
    if "bc1" in ref_table.columns and "oligo_name" in ref_table.columns:
        projection = project_bc1_reference(ref_table)
        tier_mapping = tier_mapping.merge(
            projection[["bc1", "oligo_name", "bc1_reference_n_oligos", "bc1_reference_oligo_candidates"]],
            on="bc1", how="left", validate="many_to_one")
        tier_mapping["bc1_reference_n_oligos"] = tier_mapping["bc1_reference_n_oligos"].fillna(0).astype(int)
        n_with_oligo = tier_mapping["oligo_name"].notna().sum()
        n_ambiguous = int((tier_mapping["bc1_reference_n_oligos"] > 1).sum())
        logger.info(f"  {n_with_oligo:,} BC1s mapped to one oligo; {n_ambiguous:,} have several "
                    f"reference oligos and are left unassigned")

    # Create summary statistics
    tier_stats = create_tier_summary(tier_mapping, args.library)

    # Log summary
    for tier in ["abundance_high", "abundance_medium", "abundance_low"]:
        n_bc1s = len(tier_mapping[tier_mapping["tier"] == tier])
        logger.info(f"  {tier}: {n_bc1s:,} BC1s")

    # Ensure output directories exist
    args.output_tier_mapping.parent.mkdir(parents=True, exist_ok=True)
    args.output_stats.parent.mkdir(parents=True, exist_ok=True)

    # Save outputs
    logger.info(f"Saving tier mapping to {args.output_tier_mapping}")
    tier_mapping.to_csv(args.output_tier_mapping, sep="\t", index=False)

    logger.info(f"Saving tier statistics to {args.output_stats}")
    tier_stats.to_csv(args.output_stats, sep="\t", index=False)

    logger.info("Step 04: Abundance tier classification complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
