#!/usr/bin/env python3
"""
Run Shorkie model on DNA sequences.

This script is called by Snakemake or can be run standalone.
Shorkie is an 8-fold ensemble model for RNA expression prediction
using 16kb input sequences.

Usage:
    python 04_run_shorkie.py --input sequences.fasta --output predictions.tsv

Author: Kevin R. Roy
Date: 2026-02-20
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports

from magestic.projects.vep.models.shorkie import ShorkieModel, run_shorkie_batch, SEQ_LENGTH
from magestic.projects.vep.core import load_fasta_directory, load_metadata, write_results


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Shorkie 8-fold ensemble RNA expression prediction"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input FASTA file or directory containing 16kb FASTA files"
    )
    parser.add_argument(
        "--output", "-o",
        type=Path,
        required=True,
        help="Output TSV file for predictions"
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
        "--folds",
        type=int,
        nargs="+",
        default=list(range(8)),
        help="Which folds to use (default: all 8)"
    )
    parser.add_argument(
        "--wt-cache",
        type=Path,
        help="Path to cache WT predictions for efficiency"
    )
    parser.add_argument(
        "--incremental",
        action="store_true",
        help="Skip test IDs already in output file"
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to existing output file instead of overwriting"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=== Shorkie 8-Fold Ensemble Prediction ===")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Expected sequence length: {SEQ_LENGTH} bp")

    # Load sequences
    if args.input.is_dir():
        sequences = load_fasta_directory(args.input, glob_pattern="*_tests_16kb.fasta")
    else:
        from core import read_fasta
        sequences = read_fasta(args.input)

    print(f"Loaded {len(sequences)} sequences")

    # Handle incremental processing
    existing_ids = set()
    if args.incremental and args.output.exists():
        import pandas as pd
        existing = pd.read_csv(args.output, sep='\t')
        existing_ids = set(existing['test_id'])
        sequences = {k: v for k, v in sequences.items() if k not in existing_ids}
        print(f"Incremental mode: {len(sequences)} new sequences to process")

    if len(sequences) == 0:
        print("No sequences to process")
        return

    # Load metadata if provided
    metadata = None
    if args.metadata and args.metadata.exists():
        metadata = load_metadata(args.metadata)
        print(f"Loaded metadata with {len(metadata)} entries")

    # Run Shorkie
    results = run_shorkie_batch(
        sequences=sequences,
        metadata=metadata,
        device=args.device,
        folds=args.folds,
        wt_cache_path=args.wt_cache
    )

    # Handle append mode
    if args.append and args.output.exists():
        import pandas as pd
        existing_df = pd.read_csv(args.output, sep='\t')
        results = pd.concat([existing_df, results], ignore_index=True)

    # Write results
    write_results(results, args.output, "Shorkie")

    print(f"\nCompleted: {len(results)} total predictions")


if __name__ == "__main__":
    main()
