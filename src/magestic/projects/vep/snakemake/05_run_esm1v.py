#!/usr/bin/env python3
"""
Run ESM1-v model on protein variants.

This script is called by Snakemake or can be run standalone.
ESM1-v scores missense variants using masked language modeling
to compute delta log-likelihood (mutant - WT).

Usage:
    python 05_run_esm1v.py --input variants.tsv --output scores.tsv

Author: Kevin R. Roy
Date: 2026-02-20
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports

from magestic.projects.vep.models.esm1v import ESM1vModel, run_esm1v_batch, DELETERIOUS_THRESHOLD
from magestic.projects.vep.core import load_metadata, write_results


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ESM1-v protein variant effect prediction"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input TSV file with columns: protein_id, sequence, position, wt_aa, mut_aa"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output TSV file for scores"
    )
    parser.add_argument(
        "--metadata", "-m",
        type=Path,
        help="Metadata TSV file (optional)"
    )
    parser.add_argument(
        "--device",
        default="cuda:0",
        help="PyTorch device (default: cuda:0)"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=DELETERIOUS_THRESHOLD,
        help=f"Score threshold for deleterious prediction (default: {DELETERIOUS_THRESHOLD})"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=== ESM1-v Protein Variant Scoring ===")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Deleterious threshold: {args.threshold}")

    # Load variants
    import pandas as pd
    variants_df = pd.read_csv(args.input, sep='\t')

    required_cols = ['protein_id', 'sequence', 'position', 'wt_aa', 'mut_aa']
    missing = set(required_cols) - set(variants_df.columns)
    if missing:
        print(f"ERROR: Missing required columns: {missing}")
        print(f"Available columns: {list(variants_df.columns)}")
        sys.exit(1)

    print(f"Loaded {len(variants_df)} variants")

    # Load metadata if provided
    metadata = None
    if args.metadata and args.metadata.exists():
        metadata = load_metadata(args.metadata)
        print(f"Loaded metadata with {len(metadata)} entries")

    # Run ESM1-v
    results = run_esm1v_batch(
        variants_df=variants_df,
        device=args.device,
        threshold=args.threshold
    )

    # Write results
    write_results(results, args.output, "ESM1-v")

    # Summary statistics
    if 'prediction' in results.columns:
        n_deleterious = (results['prediction'] == 'deleterious').sum()
        n_neutral = (results['prediction'] == 'neutral').sum()
        print(f"\nPrediction summary:")
        print(f"  Deleterious (< {args.threshold}): {n_deleterious} ({100*n_deleterious/len(results):.1f}%)")
        print(f"  Neutral (>= {args.threshold}): {n_neutral} ({100*n_neutral/len(results):.1f}%)")


if __name__ == "__main__":
    main()
