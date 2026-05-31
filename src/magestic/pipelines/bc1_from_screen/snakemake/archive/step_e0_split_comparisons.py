#!/usr/bin/env python3
"""
Step E0: Split Count Tables by Comparison

Splits the wide-format count matrix into per-comparison reference/treatment
pairs for DESeq2 analysis.

Creates:
- comparison_manifest.tsv: List of all comparisons with metadata
- reference/{comparison}_reference.tsv: Reference (t0) counts
- treatment/{comparison}_treatment.tsv: Treatment counts

Usage:
    python step_e0_split_comparisons.py \\
        --count-matrix count_matrix.tsv \\
        --sample-key sample_key.tsv \\
        --library yL437 \\
        --output-dir counts_for_deseq/yL437

Author: Kevin R. Roy
"""

import argparse
import logging
import sys
from pathlib import Path

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
        help="Sample key TSV with condition/generation info"
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
        "--reference-condition", type=str, default="t0",
        help="Reference condition name pattern (default: t0)"
    )
    return parser.parse_args()


def identify_sample_conditions(sample_key: pd.DataFrame, columns: list) -> dict:
    """
    Map sample columns to their conditions and generations.

    Returns:
        {sample_name: {"condition": str, "generation": int}}
    """
    sample_info = {}

    for col in columns:
        # Skip non-sample columns
        if col in ["bc1", "tier", "oligo_name"]:
            continue

        # Try to match with sample key
        if "sample_name" in sample_key.columns:
            match = sample_key[sample_key["sample_name"] == col]
            if len(match) > 0:
                row = match.iloc[0]
                sample_info[col] = {
                    "condition": row.get("condition", "unknown"),
                    "generation": row.get("generation", 0),
                }
                continue

        # Parse from column name if not in sample key
        # Common patterns: {date}_{condition}_{generation} or {condition}_gen{N}
        parts = col.replace("_", "-").split("-")
        condition = "unknown"
        generation = 0

        for i, part in enumerate(parts):
            # Look for generation markers
            if part.lower().startswith("gen"):
                try:
                    generation = int(part[3:])
                except ValueError:
                    pass
            elif part.lower() == "t0" or part == "0":
                generation = 0
            elif part.isdigit() and i == len(parts) - 1:
                generation = int(part)
            # Look for condition (non-numeric, non-date)
            elif not part.isdigit() and len(part) > 2 and not part.startswith("20"):
                condition = part

        sample_info[col] = {"condition": condition, "generation": generation}

    return sample_info


def main():
    args = parse_args()

    logger.info(f"Step E0: Splitting comparisons for library {args.library}")
    logger.info(f"  Count matrix: {args.count_matrix}")
    logger.info(f"  Sample key: {args.sample_key}")
    logger.info(f"  Output directory: {args.output_dir}")

    # Load data
    logger.info("Loading count matrix...")
    count_matrix = pd.read_csv(args.count_matrix, sep="\t")
    logger.info(f"  Loaded {len(count_matrix):,} BC1s x {len(count_matrix.columns)} columns")

    logger.info("Loading sample key...")
    sample_key = pd.read_csv(args.sample_key, sep="\t")
    logger.info(f"  Loaded {len(sample_key):,} sample entries")

    # Identify sample columns
    sample_columns = [c for c in count_matrix.columns if c not in ["bc1", "tier", "oligo_name"]]
    logger.info(f"  Sample columns: {len(sample_columns)}")

    # Map samples to conditions/generations
    sample_info = identify_sample_conditions(sample_key, sample_columns)

    # Identify reference (t0) and treatment samples
    reference_samples = [s for s, info in sample_info.items()
                        if info["generation"] == 0 or args.reference_condition in s.lower()]
    treatment_samples = [s for s, info in sample_info.items()
                        if s not in reference_samples and info["generation"] > 0]

    logger.info(f"  Reference samples: {len(reference_samples)}")
    logger.info(f"  Treatment samples: {len(treatment_samples)}")

    # Group treatment samples by condition and generation
    comparisons = {}
    for sample in treatment_samples:
        info = sample_info[sample]
        key = f"{info['condition']}_{info['generation']}"
        if key not in comparisons:
            comparisons[key] = []
        comparisons[key].append(sample)

    logger.info(f"  Unique comparisons: {len(comparisons)}")

    # Create output directories
    reference_dir = args.output_dir / "reference"
    treatment_dir = args.output_dir / "treatment"
    reference_dir.mkdir(parents=True, exist_ok=True)
    treatment_dir.mkdir(parents=True, exist_ok=True)

    # Create manifest
    manifest_rows = []

    # Split and save for each comparison
    bc1_col = "bc1"
    meta_cols = [c for c in ["bc1", "tier", "oligo_name"] if c in count_matrix.columns]

    for comparison_id, treatment_samples_list in comparisons.items():
        logger.info(f"  Processing comparison: {comparison_id}")

        # Reference counts
        ref_cols = meta_cols + reference_samples
        ref_df = count_matrix[ref_cols].copy()
        ref_df = ref_df.set_index(bc1_col)
        ref_df = ref_df[[c for c in ref_df.columns if c not in ["tier", "oligo_name"]]]

        ref_file = reference_dir / f"{args.library}_{comparison_id}_reference.tsv"
        ref_df.to_csv(ref_file, sep="\t")

        # Treatment counts
        treat_cols = meta_cols + treatment_samples_list
        treat_df = count_matrix[treat_cols].copy()
        treat_df = treat_df.set_index(bc1_col)
        treat_df = treat_df[[c for c in treat_df.columns if c not in ["tier", "oligo_name"]]]

        treat_file = treatment_dir / f"{args.library}_{comparison_id}_treatment.tsv"
        treat_df.to_csv(treat_file, sep="\t")

        # Parse condition and generation from comparison_id
        parts = comparison_id.rsplit("_", 1)
        condition = parts[0] if len(parts) > 1 else comparison_id
        generation = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0

        manifest_rows.append({
            "library": args.library,
            "comparison_id": f"{args.library}_{comparison_id}",
            "condition": condition,
            "generation": generation,
            "reference_file": ref_file.name,
            "treatment_file": treat_file.name,
            "n_reference_samples": len(reference_samples),
            "n_treatment_samples": len(treatment_samples_list),
        })

        logger.info(f"    Reference: {len(reference_samples)} samples -> {ref_file.name}")
        logger.info(f"    Treatment: {len(treatment_samples_list)} samples -> {treat_file.name}")

    # Save manifest
    manifest_df = pd.DataFrame(manifest_rows)
    manifest_file = args.output_dir / "comparison_manifest.tsv"
    manifest_df.to_csv(manifest_file, sep="\t", index=False)
    logger.info(f"Saved manifest with {len(manifest_df)} comparisons to {manifest_file}")

    logger.info("Step E0: Split complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
