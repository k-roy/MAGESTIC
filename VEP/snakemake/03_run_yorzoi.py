#!/usr/bin/env python3
"""
Run Yorzoi model on DNA sequences.

This script is called by Snakemake or can be run standalone.
Yorzoi predicts RNA expression using a Borzoi-like architecture
trained on yeast RNA-seq data.

Usage:
    python 03_run_yorzoi.py --input sequences.fasta --output predictions.tsv

Author: Kevin R. Roy
Date: 2026-02-20
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from models.yorzoi import YorzoiModel, run_yorzoi_batch
from core import load_fasta_directory, load_metadata, write_results


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Yorzoi RNA expression prediction"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Input FASTA file or directory containing FASTA files"
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
        "--wt-cache",
        type=Path,
        help="Path to cache WT predictions for efficiency"
    )
    parser.add_argument(
        "--incremental",
        action="store_true",
        help="Skip test IDs already in output file"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=== Yorzoi RNA Expression Prediction ===")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")

    # Load sequences
    if args.input.is_dir():
        sequences = load_fasta_directory(args.input, glob_pattern="*_tests_*.fasta")
    else:
        from core import read_fasta
        sequences = read_fasta(args.input)

    print(f"Loaded {len(sequences)} sequences")

    # Handle incremental processing
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

    # Run Yorzoi
    results = run_yorzoi_batch(
        sequences=sequences,
        metadata=metadata,
        device=args.device,
        wt_cache_path=args.wt_cache
    )

    # Write results
    write_results(results, args.output, "Yorzoi")

    print(f"\nCompleted: {len(results)} sequences processed")


if __name__ == "__main__":
    main()
