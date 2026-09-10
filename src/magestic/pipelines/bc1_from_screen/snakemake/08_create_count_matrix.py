#!/usr/bin/env python3
"""
Step 08: Create DESeq2 Count Matrix with T0 Consolidation

Creates a wide-format count matrix combining the two analysis abundance tiers
(high + medium) for a library, with configurable t0 sample consolidation.

NOTE: Abundance tiers are distinct from confidence tiers in bc1_donor_bc0_pipeline.

T0 Consolidation Modes:
- screen_date: Pool all t0 samples from the same screen_date (yL406, yL437)
- pseudo_replicates: Create N pseudo-replicates from t0 samples (yL442)
- none: No consolidation

-------------------------------------------------------------------------------
2026-08-06 FIXES (planning/pipeline_provenance_rebuild/10_step08_and_detector_fixes.md)
-------------------------------------------------------------------------------
D2  Tier filter is now VOCABULARY-AGNOSTIC. It previously hardcoded
    ["abundance_high", "abundance_medium"], which matches nothing in the yL406 /
    yL437 / yL442 tier mappings on Oak (they carry tier_1 / tier_2) -- the step
    then produced `Matrix shape: 0 BC1s x 0 samples` and **exited 0**. Hardcoding
    ["tier_1","tier_2"] instead would have broken the OTHER half of the corpus
    (20240806_yL410-yL433 carries abundance_*), so both spellings are accepted via
    magestic.pipelines.bc1_from_screen.config.ANALYSIS_TIERS.

D2b FAIL-LOUD GUARDS. Four separate empty-result conditions now exit non-zero with
    a specific message instead of writing a header-only matrix at exit 0
    (docs/ARCHITECTURE.md: no silent success).

D3  LIBRARY-SAMPLE FILTER RESTORED. The Feb-2026 production run filtered the long
    counts to the target library's samples before pivoting
    (yL437_count_matrix_v2.log: 69,854,629 -> 13,072,159 rows over 192 samples);
    the tracked script had lost it and would pivot across other libraries' samples.
    Screens whose sample key carries no library column (e.g. yL406, single-library)
    log a loud warning and are left unfiltered -- never silently, never fatally.

MEM Streaming pivot. The long file is no longer materialised: it is read with
    `chunksize`, each chunk is reduced to the keep-set via Index.get_indexer, and
    counts accumulate into preallocated int32 columns with np.bincount. Peak RSS
    on yL406 goes 13.629 -> 1.545 GiB with byte-identical output (audit
    planning/pipeline_provenance_rebuild/05_deseq2_memory_audit.md).

Usage:
    # Screen-date consolidation (yL437)
    python 05_create_count_matrix.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --tier-mapping tier_mapping.tsv \\
        --sample-key sample_key.tsv \\
        --library yL437 \\
        --t0-consolidation-mode screen_date \\
        --output count_matrix.tsv

    # Pseudo-replicate consolidation (yL442)
    python 05_create_count_matrix.py \\
        --bc1-counts bc1_counts_long.tsv \\
        --tier-mapping tier_mapping.tsv \\
        --sample-key sample_key.tsv \\
        --library yL442 \\
        --t0-consolidation-mode pseudo_replicates \\
        --t0-pseudo-rep-count 5 \\
        --output count_matrix.tsv

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

# Single source of truth for the tier vocabulary (see D2 above).
#
# This file is deployed in TWO trees that spell the same module differently: the
# installed package (`magestic.pipelines.bc1_from_screen.config`) and the
# standalone common/scripts tree (`bc1_from_screen_pipeline.config`). On top of
# that, the conda env the Snakemake `PYTHON` points at does NOT necessarily have
# magestic installed (verified 2026-08-06: the `snakemake` env does not, so the
# package copy of step 07 dies with ModuleNotFoundError). Resolve in three steps
# so both deployed copies stay byte-identical AND neither depends on an install:
#   1. installed package;
#   2. standalone `bc1_from_screen_pipeline` package on the parent-of-parent dir;
#   3. the sibling `config.py` that sits one level up from this snakemake/ dir --
#      true in BOTH trees, and config.py imports nothing but the stdlib + yaml.
try:
    from magestic.pipelines.bc1_from_screen.config import ANALYSIS_TIERS
except ImportError:
    _here = Path(__file__).resolve()
    try:
        sys.path.insert(0, str(_here.parents[2]))
        from bc1_from_screen_pipeline.config import ANALYSIS_TIERS
    except ImportError:
        sys.path.insert(0, str(_here.parents[1]))
        from config import ANALYSIS_TIERS

INT32_MAX = np.iinfo(np.int32).max

#: Column names, in priority order, that may carry the library assignment in a
#: sample key. `library_ID` is what the 20250912 combined_key.tsv uses (and what
#: the Feb-2026 production log filtered on).
LIBRARY_COL_CANDIDATES = ("library_ID", "library_id", "library")

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def consolidate_t0_by_screen_date(
    wide_matrix: pd.DataFrame,
    sample_key: pd.DataFrame,
    library: str,
    sample_name_col: str = "sample_name",
    generation_col: str = "generations",
    screen_date_col: str = "screen_date",
    library_col: str = "library_ID"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Consolidate t0 samples by screen_date.

    All t0 samples from the same screen_date are pooled together.
    This is appropriate when all conditions on the same date share
    the same starting population (e.g., yL406, yL437).

    Returns:
        (consolidated_matrix, metadata_df)
    """
    # Filter sample_key to library
    if library_col in sample_key.columns:
        lib_key = sample_key[sample_key[library_col] == library].copy()
    else:
        lib_key = sample_key.copy()

    # Identify t0 samples (generation == 0)
    t0_mask = lib_key[generation_col].astype(str).str.strip() == "0"
    t0_key = lib_key[t0_mask]

    if len(t0_key) == 0:
        logger.warning("No t0 samples found in sample_key")
        return wide_matrix, pd.DataFrame()

    # Get unique screen dates
    screen_dates = t0_key[screen_date_col].unique()
    logger.info(f"  Found {len(t0_key)} t0 samples across {len(screen_dates)} screen dates")

    # Create consolidated t0 columns (one per screen_date)
    consolidated_cols = {}
    metadata_rows = []

    for screen_date in screen_dates:
        date_mask = t0_key[screen_date_col] == screen_date
        date_samples = t0_key[date_mask][sample_name_col].tolist()

        # Filter to samples that exist in matrix
        date_samples = [s for s in date_samples if s in wide_matrix.columns]

        if len(date_samples) == 0:
            logger.warning(f"  No t0 samples found in matrix for screen_date {screen_date}")
            continue

        consolidated_name = f"t0_{screen_date}_consolidated"
        consolidated_cols[consolidated_name] = wide_matrix[date_samples].sum(axis=1)

        metadata_rows.append({
            "sample_name": consolidated_name,
            "screen_date": screen_date,
            "generation": 0,
            "n_samples_consolidated": len(date_samples),
            "source_samples": ",".join(date_samples)
        })

        logger.info(f"  Screen date {screen_date}: {len(date_samples)} t0 samples -> {consolidated_name}")

    # Remove original t0 samples, add consolidated
    all_t0_samples = t0_key[sample_name_col].tolist()
    non_t0_cols = [c for c in wide_matrix.columns if c not in all_t0_samples]

    result = pd.concat([
        wide_matrix[non_t0_cols],
        pd.DataFrame(consolidated_cols, index=wide_matrix.index)
    ], axis=1)

    return result, pd.DataFrame(metadata_rows)


def consolidate_t0_as_pseudo_replicates(
    wide_matrix: pd.DataFrame,
    sample_key: pd.DataFrame,
    library: str,
    n_pseudo_reps: int,
    seed: int = 42,
    sample_name_col: str = "sample_name",
    generation_col: str = "generations",
    library_col: str = "library_ID"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create pseudo-replicates from t0 samples to match treatment replicate counts.

    This is appropriate when all t0 samples are from the same screen_date
    but there are more t0 samples than condition replicates (e.g., yL442).

    Args:
        wide_matrix: BC1 x samples count matrix
        sample_key: Sample metadata
        library: Library name
        n_pseudo_reps: Number of pseudo-replicates to create
        seed: Random seed for reproducible sample assignment

    Returns:
        (consolidated_matrix, metadata_df)
    """
    # Filter sample_key to library
    if library_col in sample_key.columns:
        lib_key = sample_key[sample_key[library_col] == library].copy()
    else:
        lib_key = sample_key.copy()

    # Identify t0 samples (generation == 0)
    t0_mask = lib_key[generation_col].astype(str).str.strip() == "0"
    t0_key = lib_key[t0_mask]

    t0_samples = t0_key[sample_name_col].tolist()
    # Filter to samples that exist in matrix
    t0_samples = [s for s in t0_samples if s in wide_matrix.columns]

    if len(t0_samples) == 0:
        logger.warning("No t0 samples found")
        return wide_matrix, pd.DataFrame()

    logger.info(f"  Creating {n_pseudo_reps} pseudo-replicates from {len(t0_samples)} t0 samples")

    # Shuffle and partition
    np.random.seed(seed)
    shuffled = np.random.permutation(t0_samples)

    base_size = len(t0_samples) // n_pseudo_reps
    remainder = len(t0_samples) % n_pseudo_reps

    groups = []
    start = 0
    for i in range(n_pseudo_reps):
        # First 'remainder' groups get one extra sample
        size = base_size + (1 if i < remainder else 0)
        groups.append(shuffled[start:start+size].tolist())
        start += size

    # Create consolidated columns
    consolidated_cols = {}
    metadata_rows = []

    for i, group_samples in enumerate(groups, 1):
        pseudo_name = f"t0_pseudorep_{i}"
        consolidated_cols[pseudo_name] = wide_matrix[group_samples].sum(axis=1)

        metadata_rows.append({
            "sample_name": pseudo_name,
            "pseudo_rep_id": i,
            "generation": 0,
            "n_samples_consolidated": len(group_samples),
            "source_samples": ",".join(group_samples)
        })

        logger.info(f"  Pseudo-replicate {i}: {len(group_samples)} samples -> {pseudo_name}")

    # Remove original t0 samples, add pseudo-replicates
    non_t0_cols = [c for c in wide_matrix.columns if c not in t0_samples]

    result = pd.concat([
        wide_matrix[non_t0_cols],
        pd.DataFrame(consolidated_cols, index=wide_matrix.index)
    ], axis=1)

    return result, pd.DataFrame(metadata_rows)


def resolve_library_samples(
    sample_key: pd.DataFrame,
    library: str,
    sample_name_col: str = "sample_name",
) -> Optional[List[str]]:
    """
    Resolve the set of sample_names belonging to ``library`` (defect D3).

    Tries each of LIBRARY_COL_CANDIDATES in order and uses the first column that
    actually CONTAINS ``library``. A column that exists but never mentions the
    library is not a match -- e.g. yL406's key has a `yL` column holding trafo
    IDs (`trafo_1-1`), which would otherwise resolve to an empty sample set and
    wrongly abort a perfectly good single-library screen.

    Returns:
        List of sample names, or None if no candidate column carries the library
        (caller must warn loudly and skip filtering, NOT filter to nothing).
    """
    if sample_name_col not in sample_key.columns:
        logger.warning(
            f"  Sample key has no '{sample_name_col}' column "
            f"(has {list(sample_key.columns)[:12]}...); cannot apply the library filter."
        )
        return None

    for col in LIBRARY_COL_CANDIDATES:
        if col not in sample_key.columns:
            continue
        mask = sample_key[col].astype(str).str.strip() == library
        if not mask.any():
            logger.info(
                f"  Column '{col}' exists but contains no '{library}' "
                f"(values seen: {sorted(sample_key[col].astype(str).unique())[:8]}); "
                f"trying next candidate."
            )
            continue
        samples = sample_key.loc[mask, sample_name_col].dropna().astype(str).unique().tolist()
        logger.info(f"  Samples for library {library}: {len(samples):,} (from column '{col}')")
        return samples

    logger.warning(
        f"  NO library column found for '{library}' among {list(LIBRARY_COL_CANDIDATES)} "
        f"in the sample key -- the library filter is NOT being applied and the matrix "
        f"will span every sample present in the counts file. This is correct for a "
        f"single-library screen (e.g. yL406) and WRONG for a multi-library one; if this "
        f"screen has more than one library, add a 'library_ID' column to the sample key."
    )
    return None


def build_wide_matrix_streaming(
    bc1_counts_path: Path,
    keep_bc1: np.ndarray,
    keep_samples: Optional[List[str]],
    chunksize: int = 5_000_000,
) -> Tuple[pd.DataFrame, int, int]:
    """
    Stream the long-format counts into a dense (bc1 x sample) int32 matrix.

    Equivalent to (and byte-identical with) the previous
    ``read_csv(whole file).pivot_table(index='bc1', columns='sample_name',
    values='counts', aggfunc='sum', fill_value=0)`` at ~1/8.8th the peak RSS
    (13.629 -> 1.545 GiB on yL406). Ordering is constructed to match pivot_table:
    index = sorted observed bc1, columns = sorted observed sample_name.

    Both filters are applied INSIDE the chunk loop, so the unfiltered frame never
    exists.

    Args:
        bc1_counts_path: long-format TSV with bc1 / sample_name / counts
        keep_bc1: sorted unique bc1 keep-set (the analysis tiers)
        keep_samples: sample_names to retain, or None for "all samples"
        chunksize: rows per read_csv chunk

    Returns:
        (wide DataFrame indexed by bc1, rows_read, rows_kept)
    """
    keep_index = pd.Index(keep_bc1)          # hashed once, reused for every chunk
    n_keep = len(keep_index)
    sample_filter = set(keep_samples) if keep_samples is not None else None

    sample_to_col: Dict[str, int] = {}
    cols: List[Optional[np.ndarray]] = []
    row_seen = np.zeros(n_keep, dtype=bool)
    n_rows_read = 0
    n_rows_kept = 0

    reader = pd.read_csv(
        bc1_counts_path,
        sep="\t",
        usecols=["bc1", "sample_name", "counts"],
        dtype={"bc1": "str", "sample_name": "str", "counts": "int32"},
        chunksize=chunksize,
    )
    for chunk in reader:
        n_rows_read += len(chunk)

        pos = keep_index.get_indexer(chunk["bc1"].to_numpy())   # -1 where absent
        m = pos >= 0
        if sample_filter is not None:
            m &= chunk["sample_name"].isin(sample_filter).to_numpy()
        if not m.any():
            del pos, m, chunk
            continue

        pos = pos[m]
        cnt = chunk["counts"].to_numpy()[m]
        smp = chunk["sample_name"].to_numpy()[m]
        del chunk, m
        n_rows_kept += pos.size
        row_seen[pos] = True

        # Group surviving rows by sample: factorize (C hash) then sort int codes.
        codes, uniques = pd.factorize(smp)
        del smp
        order = np.argsort(codes, kind="stable")
        codes_s, pos_s, cnt_s = codes[order], pos[order], cnt[order]
        del codes, order, pos, cnt
        bounds = np.flatnonzero(codes_s[1:] != codes_s[:-1]) + 1
        starts = np.concatenate(([0], bounds))
        ends = np.concatenate((bounds, [codes_s.size]))
        for s, e in zip(starts, ends):
            name = uniques[codes_s[s]]
            col = sample_to_col.get(name)
            if col is None:
                col = len(cols)
                sample_to_col[name] = col
                cols.append(np.zeros(n_keep, dtype=np.int32))
            # bincount sums duplicate (bc1, sample) rows exactly as aggfunc="sum"
            add = np.bincount(pos_s[s:e], weights=cnt_s[s:e], minlength=n_keep)
            if add.max() > INT32_MAX:
                raise OverflowError(
                    "per-(bc1,sample) count exceeds int32 -- widen the accumulator"
                )
            cols[col] += add.astype(np.int32)
        del codes_s, pos_s, cnt_s, uniques

    # Assemble in pivot_table's ordering.
    sample_names = sorted(sample_to_col)
    n_obs = int(row_seen.sum())
    mat = np.empty((n_obs, len(sample_names)), dtype=np.int32)
    for j, name in enumerate(sample_names):
        k = sample_to_col[name]
        mat[:, j] = cols[k][row_seen]
        cols[k] = None                       # free as we go
    del cols, sample_to_col

    wide = pd.DataFrame(mat, columns=sample_names, copy=False)
    wide.index = pd.Index(keep_index[row_seen].to_numpy(), name="bc1")
    wide.columns.name = "sample_name"        # matches pivot_table's column axis name
    return wide, n_rows_read, n_rows_kept


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create wide-format count matrix for DESeq2 with T0 consolidation"
    )
    parser.add_argument(
        "--bc1-counts", type=Path, required=True,
        help="Input BC1 counts long-format TSV"
    )
    parser.add_argument(
        "--tier-mapping", type=Path, required=True,
        help="BC1 tier mapping TSV"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition/generation info"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name to process"
    )
    parser.add_argument(
        "--output", type=Path, required=True,
        help="Output wide-format count matrix TSV"
    )
    parser.add_argument(
        "--metadata-output", type=Path,
        help="Output path for consolidation metadata (optional)"
    )
    parser.add_argument(
        "--t0-consolidation-mode", type=str,
        choices=["screen_date", "pseudo_replicates", "none"],
        default="screen_date",
        help="T0 consolidation strategy (default: screen_date)"
    )
    parser.add_argument(
        "--t0-pseudo-rep-count", type=int,
        help="Number of t0 pseudo-replicates to create (for pseudo_replicates mode)"
    )
    parser.add_argument(
        "--screen-dates", type=str, nargs="+",
        help="Screen dates for this library (optional filter)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for pseudo-replicate assignment (default: 42)"
    )
    parser.add_argument(
        "--chunksize", type=int, default=5_000_000,
        help="Rows per read_csv chunk for the streaming pivot (default: 5,000,000). "
             "Lower it to trade wall time for peak memory."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    logger.info(f"Step 08: Creating count matrix for library {args.library}")
    logger.info(f"  BC1 counts: {args.bc1_counts}")
    logger.info(f"  Tier mapping: {args.tier_mapping}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  T0 consolidation mode: {args.t0_consolidation_mode}")

    # Validate args
    if args.t0_consolidation_mode == "pseudo_replicates" and not args.t0_pseudo_rep_count:
        logger.error("--t0-pseudo-rep-count is required for pseudo_replicates mode")
        return 1
    if args.screen_dates:
        logger.error(
            "--screen-dates is no longer supported: it filtered by substring match on "
            "sample_name, which is not a library scoping rule. Scope the run with the "
            "sample key's library column instead (defect D3)."
        )
        return 1

    # ---- Load the two small inputs; the BIG long file is streamed later ------
    logger.info("Loading tier mapping...")
    tier_mapping = pd.read_csv(
        args.tier_mapping, sep="\t", usecols=["bc1", "tier"],
        dtype={"bc1": "str", "tier": "str"},
    )
    logger.info(f"  Loaded {len(tier_mapping):,} BC1 tier assignments")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} sample entries")

    # ---- D2: tier keep-set, accepting BOTH tier vocabularies ----------------
    tier_bc1s = np.sort(
        tier_mapping.loc[tier_mapping["tier"].isin(ANALYSIS_TIERS), "bc1"].unique()
    )
    observed_tiers = sorted(tier_mapping["tier"].dropna().unique().tolist())
    logger.info(
        f"  BC1s in analysis tiers {list(ANALYSIS_TIERS)}: {len(tier_bc1s):,} "
        f"(tier values present in the mapping: {observed_tiers})"
    )

    # GUARD 1 -- empty tier keep-set
    if len(tier_bc1s) == 0:
        logger.error(
            "EMPTY TIER KEEP-SET: no BC1 in the tier mapping carries any of the accepted "
            f"analysis tiers {list(ANALYSIS_TIERS)}. The mapping contains {observed_tiers}. "
            "Refusing to write an empty count matrix (docs/ARCHITECTURE.md: no silent success)."
        )
        return 1

    # ---- D3: library-sample keep-set ---------------------------------------
    keep_samples = resolve_library_samples(sample_key, args.library)

    # ---- Streaming pivot (see build_wide_matrix_streaming) ------------------
    logger.info("Creating wide-format count matrix (streaming)...")
    wide_matrix, n_rows_read, n_rows_kept = build_wide_matrix_streaming(
        args.bc1_counts, tier_bc1s, keep_samples, chunksize=args.chunksize,
    )
    logger.info(f"  Loaded {n_rows_read:,} rows")
    logger.info(f"  After tier + library filters: {n_rows_kept:,} rows")

    # GUARD 2/3 -- everything filtered away
    if n_rows_kept == 0:
        if keep_samples is None:
            logger.error(
                f"EMPTY AFTER TIER FILTER: none of the {n_rows_read:,} count rows carry a "
                f"BC1 from the {len(tier_bc1s):,}-BC1 analysis-tier keep-set. The tier "
                "mapping and the counts file do not describe the same screen."
            )
        else:
            logger.error(
                f"EMPTY AFTER LIBRARY FILTER: none of the {n_rows_read:,} count rows match "
                f"both the analysis-tier keep-set ({len(tier_bc1s):,} BC1s) and the "
                f"{len(keep_samples):,} sample(s) assigned to library '{args.library}'. "
                "Check that the sample key's library column and the counts file's "
                "sample_name values come from the same screen."
            )
        return 1

    logger.info(
        f"  Matrix shape: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples"
    )

    # Apply T0 consolidation
    metadata_df = pd.DataFrame()

    if args.t0_consolidation_mode == "screen_date":
        logger.info("Consolidating t0 samples by screen_date...")
        wide_matrix, metadata_df = consolidate_t0_by_screen_date(
            wide_matrix, sample_key, args.library
        )
        logger.info(f"  Consolidated matrix: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")

    elif args.t0_consolidation_mode == "pseudo_replicates":
        logger.info(f"Creating {args.t0_pseudo_rep_count} t0 pseudo-replicates...")
        wide_matrix, metadata_df = consolidate_t0_as_pseudo_replicates(
            wide_matrix, sample_key, args.library,
            n_pseudo_reps=args.t0_pseudo_rep_count,
            seed=args.seed
        )
        logger.info(f"  Consolidated matrix: {wide_matrix.shape[0]:,} BC1s x {wide_matrix.shape[1]:,} samples")

    else:
        logger.info("No t0 consolidation applied")

    # Add tier information to output
    tier_info = tier_mapping[tier_mapping["bc1"].isin(wide_matrix.index)][["bc1", "tier"]].drop_duplicates()

    # Reset index for output
    wide_matrix = wide_matrix.reset_index()

    # Merge tier info
    wide_matrix = wide_matrix.merge(tier_info, on="bc1", how="left")

    # Reorder columns: bc1, tier, then samples
    cols = ["bc1", "tier"] + [c for c in wide_matrix.columns if c not in ["bc1", "tier"]]
    wide_matrix = wide_matrix[cols]

    # GUARD 4 -- an empty (or sample-less) final matrix must never be written as
    # a "successful" header-only file for step 10 to consume.
    n_sample_cols = wide_matrix.shape[1] - 2
    if len(wide_matrix) == 0 or n_sample_cols == 0:
        logger.error(
            f"EMPTY COUNT MATRIX: {len(wide_matrix):,} BC1s x {n_sample_cols:,} samples "
            "after t0 consolidation. Nothing was written. Refusing to hand a header-only "
            "matrix to step 10 (docs/ARCHITECTURE.md: no silent success)."
        )
        return 1

    # Ensure output directory exists
    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Save output
    logger.info(f"Saving count matrix to {args.output}")
    wide_matrix.to_csv(args.output, sep="\t", index=False)

    # Save metadata if requested
    if args.metadata_output and len(metadata_df) > 0:
        logger.info(f"Saving consolidation metadata to {args.metadata_output}")
        args.metadata_output.parent.mkdir(parents=True, exist_ok=True)
        metadata_df.to_csv(args.metadata_output, sep="\t", index=False)

    # Log summary by abundance tier -- iterate over the tiers actually present so
    # this works under either vocabulary (D2).
    for tier in [t for t in ANALYSIS_TIERS if (wide_matrix["tier"] == t).any()]:
        n_tier = (wide_matrix["tier"] == tier).sum()
        logger.info(f"  {tier}: {n_tier:,} BC1s")
    n_untiered = wide_matrix["tier"].isna().sum()
    if n_untiered:
        logger.warning(
            f"  {n_untiered:,} BC1s in the matrix have NO tier after the merge -- "
            "the tier mapping has duplicate or missing bc1 entries."
        )

    logger.info("Step 08: Count matrix creation complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
