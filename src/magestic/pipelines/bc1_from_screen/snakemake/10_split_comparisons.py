#!/usr/bin/env python3
"""
Step 10: Split Count Tables by Comparison

Splits the wide-format count matrix into per-comparison reference/treatment
pairs for DESeq2 analysis.

IMPORTANT: Uses screen_date-based t0 consolidation.
All gen=0 samples from the same screen_date are from the same glycerol aliquot
(technical replicates), so they are consolidated as reference for all
comparisons on that screen_date.

Creates:
- comparison_manifest.tsv: List of all comparisons with metadata
- reference/{comparison}_reference.tsv: Reference (t0) counts
- treatment/{comparison}_treatment.tsv: Treatment counts

Usage:
    python 10_split_comparisons.py \
        --count-matrix count_matrix.tsv \
        --sample-key sample_key.tsv \
        --library yL406 \
        --output-dir counts_for_deseq/yL406

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split count matrix into per-comparison tables"
    )
    parser.add_argument(
        "--count-matrix", type=Path, required=True,
        help="Input wide-format count matrix TSV"
    )
    parser.add_argument(
        "--sample-key", type=Path, required=True,
        help="Sample key TSV with condition/generation/screen_date info"
    )
    parser.add_argument(
        "--library", type=str, required=True,
        help="Library name"
    )
    parser.add_argument(
        "--output-dir", type=Path, required=True,
        help="Output directory for split tables"
    )
    parser.add_argument(
        "--screen-date-column", type=str, default="screen_date",
        help="Column name for screen date (default: screen_date)"
    )
    parser.add_argument(
        "--condition-column", type=str, default="condition",
        help="Column name for condition (default: condition)"
    )
    parser.add_argument(
        "--generation-column", type=str, default="generations",
        help="Column name for generation (default: generations)"
    )
    parser.add_argument(
        "--sample-name-column", type=str, default="sample_name",
        help="Column name for sample name (default: sample_name)"
    )
    parser.add_argument(
        "--consolidate-by-screen-date", action="store_true", default=True,
        help="Consolidate all gen=0 samples from same screen_date as reference (default: True)"
    )
    parser.add_argument(
        "--no-consolidate", dest="consolidate_by_screen_date", action="store_false",
        help="Use condition-specific t0 only (not recommended)"
    )
    parser.add_argument(
        "--split-by-yeast-sublibrary", action="store_true", default=False,
        help="Split comparisons by yeast_sublibrary (requires --bc1-reference-table)"
    )
    parser.add_argument(
        "--bc1-reference-table", type=Path,
        help="BC1 reference table with yeast_sublibrary annotations (required if --split-by-yeast-sublibrary)"
    )
    parser.add_argument(
        "--yeast-sublibrary-column", type=str, default="oligo_name",
        help="Column to extract yeast_sublibrary from (default: oligo_name, first underscore-delimited field)"
    )
    return parser.parse_args()


def build_sample_info(sample_key: pd.DataFrame,
                      sample_columns: List[str],
                      sample_name_col: str,
                      condition_col: str,
                      generation_col: str,
                      screen_date_col: str) -> Dict:
    """
    Build sample info dictionary from sample key.

    Returns:
        {sample_name: {
            "condition": str,
            "generation": int,
            "screen_date": str
        }}
    """
    sample_info = {}

    # Build lookup from sample key
    key_lookup = {}
    if sample_name_col in sample_key.columns:
        for _, row in sample_key.iterrows():
            sample_name = row.get(sample_name_col)
            if pd.isna(sample_name):
                continue

            # Get generation, handling 'none' or other non-numeric values
            gen_value = row.get(generation_col, row.get("generation", 0))
            try:
                if pd.isna(gen_value) or str(gen_value).lower() == "none":
                    gen_value = 0
                else:
                    gen_value = int(gen_value)
            except (ValueError, TypeError):
                gen_value = 0

            # Get screen_date
            screen_date = row.get(screen_date_col, "unknown")
            if pd.isna(screen_date):
                screen_date = "unknown"
            else:
                screen_date = str(screen_date)

            # Get condition
            condition = row.get(condition_col, "unknown")
            if pd.isna(condition):
                condition = "unknown"

            key_lookup[sample_name] = {
                "condition": condition,
                "generation": gen_value,
                "screen_date": screen_date
            }

    # Map sample columns to info
    for col in sample_columns:
        if col in key_lookup:
            sample_info[col] = key_lookup[col]
        else:
            # Not found in sample key - use defaults
            sample_info[col] = {
                "condition": "unknown",
                "generation": 0,
                "screen_date": "unknown"
            }

    return sample_info


def extract_yeast_sublibrary_mapping(bc1_reference: pd.DataFrame, yeast_sublibrary_col: str) -> Dict[str, str]:
    """
    Extract BC1 -> yeast_sublibrary mapping from reference table.

    For oligo_name column, extracts prefix before first underscore (e.g., "SpG", "LbCas12a").

    Returns:
        {bc1: yeast_sublibrary}
    """
    if yeast_sublibrary_col not in bc1_reference.columns:
        logger.error(f"Column '{yeast_sublibrary_col}' not found in BC1 reference table")
        return {}

    # Vectorized extraction - much faster than iterrows
    df = bc1_reference[["bc1", yeast_sublibrary_col]].dropna(subset=["bc1"]).copy()

    if yeast_sublibrary_col == "oligo_name":
        # Extract prefix before first underscore
        df["yeast_sublib"] = df[yeast_sublibrary_col].fillna("unknown").astype(str).str.split("_").str[0]
    else:
        df["yeast_sublib"] = df[yeast_sublibrary_col].fillna("unknown").astype(str)

    return dict(zip(df["bc1"], df["yeast_sublib"]))


def group_samples_by_screen_date(sample_info: Dict) -> Dict[str, Dict]:
    """
    Group samples by screen_date for t0 consolidation.

    All gen=0 samples from the same screen_date are from the same
    glycerol aliquot (technical replicates), so we consolidate them.

    Returns:
        {screen_date: {
            "ref": [list of gen=0 sample names],
            "treat": {comparison_id: [list of treatment sample names]}
        }}
    """
    by_screen_date = {}

    for sample, info in sample_info.items():
        screen_date = info["screen_date"]

        if screen_date not in by_screen_date:
            by_screen_date[screen_date] = {"ref": [], "treat": {}}

        if info["generation"] == 0:
            by_screen_date[screen_date]["ref"].append(sample)
        else:
            condition = info["condition"]
            gen = info["generation"]
            key = f"{condition}_{gen}"

            if key not in by_screen_date[screen_date]["treat"]:
                by_screen_date[screen_date]["treat"][key] = []
            by_screen_date[screen_date]["treat"][key].append(sample)

    return by_screen_date


def create_comparisons_by_screen_date(
    count_matrix_filtered: pd.DataFrame,
    sample_info: Dict,
    args,
    reference_dir: Path,
    treatment_dir: Path,
    sublib_suffix: str,
    yeast_sublib
) -> List[Dict]:
    """
    Create comparisons using screen_date-based t0 consolidation.

    Returns list of manifest rows.
    """
    manifest_rows = []
    bc1_col = "bc1"

    by_screen_date = group_samples_by_screen_date(sample_info)

    logger.info("Screen-date summary:")
    for sd, data in sorted(by_screen_date.items()):
        n_ref = len(data["ref"])
        n_treat_groups = len(data["treat"])
        logger.info(f"  {sd}: {n_ref} reference samples, {n_treat_groups} treatment groups")

    # Create comparisons grouped by screen_date
    for screen_date, data in by_screen_date.items():
        reference_samples = data["ref"]

        if not reference_samples:
            logger.warning(f"  No reference (gen=0) samples for screen_date {screen_date}, skipping")
            continue

        for comparison_id, treatment_samples in data["treat"].items():
            if not treatment_samples:
                continue

            # Sanitize for filename - include screen_date to make unique
            safe_id = comparison_id.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")
            full_id = f"{screen_date}_{safe_id}"  # Include screen_date for uniqueness

            logger.info(f"  Processing: {comparison_id} (screen_date={screen_date})")
            logger.info(f"    Reference: {len(reference_samples)} samples (all gen=0 from {screen_date})")
            logger.info(f"    Treatment: {len(treatment_samples)} samples")

            # Reference counts - ALL gen=0 from this screen_date
            ref_cols = [bc1_col] + reference_samples
            ref_df = count_matrix_filtered[ref_cols].copy()
            ref_df = ref_df.set_index(bc1_col)

            ref_file = reference_dir / f"{args.library}{sublib_suffix}_{full_id}_reference.tsv"
            ref_df.to_csv(ref_file, sep="\t")

            # Treatment counts
            treat_cols = [bc1_col] + treatment_samples
            treat_df = count_matrix_filtered[treat_cols].copy()
            treat_df = treat_df.set_index(bc1_col)

            treat_file = treatment_dir / f"{args.library}{sublib_suffix}_{full_id}_treatment.tsv"
            treat_df.to_csv(treat_file, sep="\t")

            # Parse condition and generation
            parts = comparison_id.rsplit("_", 1)
            condition = parts[0] if len(parts) > 1 else comparison_id
            generation = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

            manifest_rows.append({
                "library": args.library,
                "yeast_sublibrary": yeast_sublib if yeast_sublib else "all",
                "comparison_id": f"{args.library}{sublib_suffix}_{full_id}",
                "condition": condition,
                "generation": generation,
                "screen_date": screen_date,
                "reference_file": ref_file.name,
                "treatment_file": treat_file.name,
                "n_reference_samples": len(reference_samples),
                "n_treatment_samples": len(treatment_samples),
                "n_bc1s": len(ref_df),
            })

    return manifest_rows


def create_comparisons_no_consolidation(
    count_matrix_filtered: pd.DataFrame,
    sample_info: Dict,
    args,
    reference_dir: Path,
    treatment_dir: Path,
    sublib_suffix: str,
    yeast_sublib
) -> List[Dict]:
    """
    Create comparisons using ALL gen=0 samples as reference (no screen_date grouping).

    Returns list of manifest rows.
    """
    manifest_rows = []
    bc1_col = "bc1"

    logger.info("Using all gen=0 samples as pooled reference (no screen_date consolidation)")

    # Identify all reference and treatment samples
    reference_samples = [s for s, info in sample_info.items() if info["generation"] == 0]
    treatment_samples = [s for s, info in sample_info.items() if info["generation"] > 0]

    logger.info(f"  Reference samples (all gen=0): {len(reference_samples)}")
    logger.info(f"  Treatment samples: {len(treatment_samples)}")

    if not reference_samples:
        logger.warning("  No reference (gen=0) samples found!")
        return manifest_rows

    # Group treatment by condition and generation
    comparisons = {}
    for sample in treatment_samples:
        info = sample_info[sample]
        key = f"{info['condition']}_{info['generation']}"
        if key not in comparisons:
            comparisons[key] = []
        comparisons[key].append(sample)

    logger.info(f"  Creating {len(comparisons)} comparisons")

    for comparison_id, treatment_list in comparisons.items():
        safe_id = comparison_id.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")

        # Reference counts (ALL gen=0)
        ref_cols = [bc1_col] + reference_samples
        ref_df = count_matrix_filtered[ref_cols].copy()
        ref_df = ref_df.set_index(bc1_col)

        ref_file = reference_dir / f"{args.library}{sublib_suffix}_{safe_id}_reference.tsv"
        ref_df.to_csv(ref_file, sep="\t")

        # Treatment counts
        treat_cols = [bc1_col] + treatment_list
        treat_df = count_matrix_filtered[treat_cols].copy()
        treat_df = treat_df.set_index(bc1_col)

        treat_file = treatment_dir / f"{args.library}{sublib_suffix}_{safe_id}_treatment.tsv"
        treat_df.to_csv(treat_file, sep="\t")

        parts = comparison_id.rsplit("_", 1)
        condition = parts[0] if len(parts) > 1 else comparison_id
        generation = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

        manifest_rows.append({
            "library": args.library,
            "yeast_sublibrary": yeast_sublib if yeast_sublib else "all",
            "comparison_id": f"{args.library}{sublib_suffix}_{safe_id}",
            "condition": condition,
            "generation": generation,
            "screen_date": "pooled",
            "reference_file": ref_file.name,
            "treatment_file": treat_file.name,
            "n_reference_samples": len(reference_samples),
            "n_treatment_samples": len(treatment_list),
            "n_bc1s": len(ref_df),
        })

        logger.info(f"  {comparison_id}: ref={len(reference_samples)}, treat={len(treatment_list)}")

    return manifest_rows


def main():
    args = parse_args()

    logger.info(f"Step 10: Splitting comparisons for library {args.library}")
    logger.info(f"  Count matrix: {args.count_matrix}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  Output directory: {args.output_dir}")
    logger.info(f"  Screen-date consolidation: {args.consolidate_by_screen_date}")
    logger.info(f"  Screen-date column: {args.screen_date_column}")

    # Load data
    logger.info("Loading count matrix...")
    count_matrix = pd.read_csv(args.count_matrix, sep="\t")
    logger.info(f"  Loaded {len(count_matrix):,} BC1s x {len(count_matrix.columns)} columns")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} sample entries")
    logger.info(f"  Sample key columns: {list(sample_key.columns)}")

    # Handle yeast_sublibrary splitting
    yeast_sublibrary_bc1_mapping = None
    yeast_sublibraries_to_process = [None]  # Default: process all BC1s together

    if args.split_by_yeast_sublibrary:
        if not args.bc1_reference_table:
            logger.error("--bc1-reference-table required when using --split-by-yeast-sublibrary")
            return 1

        logger.info(f"Loading BC1 reference table for yeast_sublibrary annotations...")
        bc1_reference = pd.read_csv(args.bc1_reference_table, sep="\t")
        logger.info(f"  Loaded {len(bc1_reference):,} BC1s")

        yeast_sublibrary_bc1_mapping = extract_yeast_sublibrary_mapping(bc1_reference, args.yeast_sublibrary_column)
        logger.info(f"  Extracted yeast_sublibrary annotations for {len(yeast_sublibrary_bc1_mapping):,} BC1s")

        # Determine unique sublibraries
        yeast_sublibraries_to_process = sorted(set(yeast_sublibrary_bc1_mapping.values()))
        logger.info(f"  Sublibraries to process: {yeast_sublibraries_to_process}")

        # Show BC1 counts per yeast_sublibrary
        from collections import Counter
        sublibrary_counts = Counter(yeast_sublibrary_bc1_mapping.values())
        for sublib in sorted(sublibrary_counts.keys()):
            logger.info(f"    {sublib}: {sublibrary_counts[sublib]:,} BC1s")

    # Identify sample columns (excluding metadata)
    meta_cols = ["bc1", "tier", "oligo_name"]
    sample_columns = [c for c in count_matrix.columns if c not in meta_cols]
    logger.info(f"  Sample columns: {len(sample_columns)}")

    # Build sample info from sample key
    sample_info = build_sample_info(
        sample_key,
        sample_columns,
        args.sample_name_column,
        args.condition_column,
        args.generation_column,
        args.screen_date_column
    )

    # Check for screen_date availability
    screen_dates = set(info["screen_date"] for info in sample_info.values())
    logger.info(f"  Unique screen_dates: {screen_dates}")

    if "unknown" in screen_dates and len(screen_dates) == 1:
        logger.warning("  No screen_date information found in sample key!")
        logger.warning("  Falling back to no consolidation (using all gen=0 as pooled reference)")
        args.consolidate_by_screen_date = False

    # Create output directories
    reference_dir = args.output_dir / "reference"
    treatment_dir = args.output_dir / "treatment"
    reference_dir.mkdir(parents=True, exist_ok=True)
    treatment_dir.mkdir(parents=True, exist_ok=True)

    # Manifest rows
    manifest_rows = []
    bc1_col = "bc1"

    # Loop over sublibraries (or process all BC1s together if not splitting)
    for yeast_sublib in yeast_sublibraries_to_process:
        # Filter count matrix by yeast_sublibrary if splitting
        if yeast_sublib is None:
            # Process all BC1s together
            count_matrix_filtered = count_matrix.copy()
            sublib_suffix = ""
            logger.info(f"\nProcessing all BC1s together...")
        else:
            # Filter to BC1s from this yeast_sublibrary
            bc1s_in_yeast_sublibrary = [bc1 for bc1, sl in yeast_sublibrary_bc1_mapping.items() if sl == yeast_sublib]
            count_matrix_filtered = count_matrix[count_matrix[bc1_col].isin(bc1s_in_yeast_sublibrary)].copy()
            sublib_suffix = f"_{yeast_sublib}"
            logger.info(f"\nProcessing yeast_sublibrary: {yeast_sublib}")
            logger.info(f"  Filtered to {len(count_matrix_filtered):,} BC1s")

            if len(count_matrix_filtered) == 0:
                logger.warning(f"  No BC1s found for yeast_sublibrary {yeast_sublib}, skipping")
                continue

        # Create comparisons using appropriate method
        if args.consolidate_by_screen_date:
            rows = create_comparisons_by_screen_date(
                count_matrix_filtered, sample_info, args,
                reference_dir, treatment_dir, sublib_suffix, yeast_sublib
            )
        else:
            rows = create_comparisons_no_consolidation(
                count_matrix_filtered, sample_info, args,
                reference_dir, treatment_dir, sublib_suffix, yeast_sublib
            )

        manifest_rows.extend(rows)

    # Save manifest
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_file = args.output_dir / "comparison_manifest.tsv"
    manifest_df.to_csv(manifest_file, sep="\t", index=False)
    logger.info(f"\nSaved manifest with {len(manifest_df)} comparisons to {manifest_file}")

    # Summary
    if len(manifest_df) > 0:
        ref_counts = manifest_df["n_reference_samples"]
        logger.info(f"Reference sample counts: min={ref_counts.min()}, max={ref_counts.max()}, "
                   f"median={ref_counts.median():.0f}")

    logger.info("Step 10: Split complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
