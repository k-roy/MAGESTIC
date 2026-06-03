#!/usr/bin/env python3
"""
Run Evo2 model on DNA sequences.

This script is called by Snakemake or can be run standalone.
Evo2 scores sequences using DNA language model perplexity.

Usage:
    python 02_run_evo2.py --input sequences.fasta --output scores.tsv

Author: Kevin R. Roy
Date: 2026-02-20
"""

import argparse
import sys
from pathlib import Path

# Add parent directory to path for imports

from magestic.projects.vep.models.evo2 import run_evo2_predictions, generate_slurm_script
from magestic.projects.vep.core import load_fasta_directory, load_metadata, write_results


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Evo2 DNA language model scoring"
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
        help="Output TSV file for scores"
    )
    parser.add_argument(
        "--metadata", "-m",
        type=Path,
        help="Metadata TSV file (optional)"
    )
    parser.add_argument(
        "--model",
        choices=["evo2_7b", "evo2_40b"],
        default="evo2_7b",
        help="Evo2 model size (default: evo2_7b)"
    )
    parser.add_argument(
        "--device",
        default="cuda:0",
        help="PyTorch device (default: cuda:0)"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1,
        help="Batch size for inference (default: 1)"
    )
    parser.add_argument(
        "--save-logprobs",
        action="store_true",
        help="Save per-position log probabilities"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print(f"=== Evo2 {args.model} Scoring ===")
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")

    # Load sequences
    if args.input.is_dir():
        sequences = load_fasta_directory(args.input, glob_pattern="*_tests_*.fasta")
    else:
        from core import read_fasta
        sequences = read_fasta(args.input)

    print(f"Loaded {len(sequences)} sequences")

    # Load metadata if provided
    metadata = None
    if args.metadata and args.metadata.exists():
        metadata = load_metadata(args.metadata)
        print(f"Loaded metadata with {len(metadata)} entries")

    # Run Evo2
    # NOTE: Evo2 model runs via Singularity container with file-based I/O
    # Write sequences to temp file, run predictions, parse results
    import tempfile
    import os

    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
        temp_input = f.name
        for seq_id, seq in sequences.items():
            f.write(f"{seq}\n")

    temp_output = temp_input.replace('.txt', '_scores.tsv')

    try:
        # Determine FP8 usage from model name
        use_fp8 = "H100" in args.device or "h100" in args.device.lower()
        model_size = "7b" if "7b" in args.model else "40b"

        results_df = run_evo2_predictions(
            input_file=temp_input,
            output_file=temp_output,
            model_size=model_size,
            use_fp8=use_fp8,
            batch_size=args.batch_size,
        )

        # Map results back to sequence IDs
        import pandas as pd
        seq_ids = list(sequences.keys())
        results_df['test_id'] = [seq_ids[i] for i in results_df['sequence_index']]
        results = results_df[['test_id', 'score']].copy()
        results.columns = ['test_id', 'evo2_score']

    finally:
        # Cleanup temp files
        if os.path.exists(temp_input):
            os.remove(temp_input)
        if os.path.exists(temp_output):
            os.remove(temp_output)

    # Write results
    write_results(results, args.output, f"Evo2 {args.model}")

    print(f"\nCompleted: {len(results)} sequences scored")


if __name__ == "__main__":
    main()
