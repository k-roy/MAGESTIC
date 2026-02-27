#!/usr/bin/env python3
"""
Backfill VEP Cache from Existing Results

Import existing VEP prediction results into the central cache so that
future runs can skip already-computed sequences.

Usage:
    python backfill_cache.py \
        --results shorkie_pairwise_16kb_predictions.tsv \
        --fasta-dir pairwise_tests_50kb/fasta_16kb_shorkie/ \
        --model shorkie \
        --project bloom_qtl_pairwise

Author: Kevin R. Roy
Date: 2026-02-22
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, Optional
from tqdm import tqdm

import pandas as pd

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
from vep_cache import VEPCache, create_model_params_file


# Default cache location
DEFAULT_CACHE_DIR = Path("/path/to/common/vep_cache")

# Model-specific configurations
MODEL_CONFIGS = {
    "shorkie": {
        "window_size": 16384,
        "model_version": "shorkie_v1.0",
        "score_column": "logSED",
        "test_id_column": "test_id",
        "fasta_pattern": "*_tests_16kb.fasta",
        "params": {
            "window_size": 16384,
            "model_type": "shorkie_8fold_ensemble",
            "scoring": "logSED"
        }
    },
    "yorzoi": {
        "window_size": 50000,
        "model_version": "yorzoi_v1.0",
        "score_column": "log2fc",
        "test_id_column": "test_id",
        "fasta_pattern": "*_tests_50kb.fasta",
        "params": {
            "window_size": 50000,
            "model_type": "yorzoi_expression",
            "scoring": "log2_fold_change"
        }
    },
    "evo2_7b": {
        "window_size": 50000,
        "model_version": "evo2_7b_v1.0",
        "score_column": "delta_score",
        "test_id_column": "test_id",
        "fasta_pattern": "*_tests_50kb.fasta",
        "params": {
            "window_size": 50000,
            "model_type": "evo2_7b",
            "scoring": "delta_log_likelihood"
        }
    },
    "evo2_40b": {
        "window_size": 50000,
        "model_version": "evo2_40b_v1.0",
        "score_column": "delta_score",
        "test_id_column": "test_id",
        "fasta_pattern": "*_tests_50kb.fasta",
        "params": {
            "window_size": 50000,
            "model_type": "evo2_40b",
            "scoring": "delta_log_likelihood"
        }
    },
    "esm1v": {
        "window_size": None,  # Protein-based, no window
        "model_version": "esm1v_v1.0",
        "score_column": "delta_logits",
        "test_id_column": "test_id",
        "fasta_pattern": None,  # Uses protein sequences
        "params": {
            "model_type": "esm1v_protein",
            "scoring": "delta_logits"
        }
    }
}


def load_fasta_sequences(fasta_dir: Path, pattern: str = "*.fasta") -> Dict[str, str]:
    """Load all sequences from FASTA files in a directory.

    Args:
        fasta_dir: Directory containing FASTA files
        pattern: Glob pattern for FASTA files

    Returns:
        Dictionary mapping sequence ID -> sequence
    """
    sequences = {}

    for fasta_file in fasta_dir.glob(pattern):
        current_id = None
        current_seq = []

        with open(fasta_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        sequences[current_id] = "".join(current_seq)
                    # Extract ID from header (first word after >)
                    current_id = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)

            if current_id:
                sequences[current_id] = "".join(current_seq)

    return sequences


def backfill_from_results(
    results_tsv: Path,
    fasta_dir: Optional[Path],
    model: str,
    project_name: str,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    protein_column: Optional[str] = None,
    batch_size: int = 1000,
    dry_run: bool = False
) -> int:
    """Backfill cache from existing results file.

    Args:
        results_tsv: Path to TSV file with prediction results
        fasta_dir: Directory containing FASTA files with sequences
        model: Model name (shorkie, yorzoi, evo2_7b, etc.)
        project_name: Name of the project for metadata
        cache_dir: Base cache directory
        protein_column: For ESM1-v, column containing protein sequence
        batch_size: Number of entries to insert per batch
        dry_run: If True, don't actually insert (just show what would be done)

    Returns:
        Number of entries backfilled
    """
    # Get model configuration
    if model not in MODEL_CONFIGS:
        raise ValueError(f"Unknown model: {model}. Choose from: {list(MODEL_CONFIGS.keys())}")

    config = MODEL_CONFIGS[model]

    print(f"\n{'=' * 60}")
    print(f"Backfilling {model} cache from {results_tsv.name}")
    print(f"{'=' * 60}")

    # Load results
    print(f"\nLoading results from {results_tsv}...")
    results_df = pd.read_csv(results_tsv, sep="\t")
    print(f"  Loaded {len(results_df)} results")

    # Check required columns
    test_id_col = config["test_id_column"]
    score_col = config["score_column"]

    if test_id_col not in results_df.columns:
        raise ValueError(f"Results missing {test_id_col} column")

    if score_col not in results_df.columns:
        # Try alternative score columns
        alt_cols = ["score", "logSED", "delta_score", "delta_logits", "log2fc"]
        for alt in alt_cols:
            if alt in results_df.columns:
                score_col = alt
                break
        else:
            raise ValueError(f"Results missing score column ({score_col})")

    print(f"  Using score column: {score_col}")

    # Filter to valid results (non-null scores)
    valid_results = results_df[results_df[score_col].notna()]
    print(f"  Valid results (non-null score): {len(valid_results)}")

    # Load sequences
    if model == "esm1v" and protein_column:
        # For ESM1-v, sequences are in the results file
        print(f"\nUsing protein sequences from column: {protein_column}")
        sequences = dict(zip(valid_results[test_id_col], valid_results[protein_column]))
    elif fasta_dir:
        print(f"\nLoading sequences from {fasta_dir}...")
        pattern = config["fasta_pattern"] or "*.fasta"
        sequences = load_fasta_sequences(fasta_dir, pattern)
        print(f"  Loaded {len(sequences)} sequences")
    else:
        raise ValueError("Must provide either fasta_dir or protein_column")

    # Initialize cache
    cache = VEPCache(cache_dir, model)
    params_hash = cache.hash_params(config["params"])

    print(f"\nModel params hash: {params_hash}")
    print(f"Existing cache entries: {cache.count(params_hash)}")

    if dry_run:
        print("\n[DRY RUN] Would backfill the following:")

    # Process in batches
    entries = []
    skipped_no_seq = 0
    skipped_existing = 0

    # Check existing entries first
    print("\nChecking for existing cache entries...")
    test_ids = valid_results[test_id_col].tolist()
    seq_hashes = {}
    for test_id in tqdm(test_ids, desc="Hashing sequences"):
        if test_id in sequences:
            seq_hashes[test_id] = cache.hash_sequence(sequences[test_id])

    # Batch check existing
    existing = cache.get_batch(list(seq_hashes.values()), params_hash)
    existing_hashes = {h for h, v in existing.items() if v is not None}

    print(f"\nProcessing {len(valid_results)} results...")
    for _, row in tqdm(valid_results.iterrows(), total=len(valid_results), desc="Backfilling"):
        test_id = row[test_id_col]

        # Check if we have the sequence
        if test_id not in sequences:
            skipped_no_seq += 1
            continue

        sequence = sequences[test_id]
        seq_hash = seq_hashes.get(test_id) or cache.hash_sequence(sequence)

        # Skip if already in cache
        if seq_hash in existing_hashes:
            skipped_existing += 1
            continue

        # Build result dict from row
        result = row.to_dict()
        result["score"] = row[score_col]
        result["seq_length"] = len(sequence)

        entries.append((seq_hash, result, test_id))

        # Insert batch
        if len(entries) >= batch_size and not dry_run:
            cache.put_batch(
                entries, params_hash, project_name,
                model_version=config["model_version"],
                window_size=config["window_size"]
            )
            entries = []

    # Insert remaining
    if entries and not dry_run:
        cache.put_batch(
            entries, params_hash, project_name,
            model_version=config["model_version"],
            window_size=config["window_size"]
        )

    # Final stats
    n_backfilled = len(valid_results) - skipped_no_seq - skipped_existing

    print(f"\n{'=' * 60}")
    print("Backfill Summary")
    print(f"{'=' * 60}")
    print(f"Total results:        {len(valid_results)}")
    print(f"Skipped (no seq):     {skipped_no_seq}")
    print(f"Skipped (existing):   {skipped_existing}")
    print(f"Backfilled:           {n_backfilled}")
    print(f"Final cache size:     {cache.count(params_hash)}")

    if dry_run:
        print("\n[DRY RUN] No changes made")

    return n_backfilled


def main():
    parser = argparse.ArgumentParser(
        description="Backfill VEP cache from existing results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Backfill Shorkie predictions
  python backfill_cache.py \\
      --results shorkie_pairwise_16kb_predictions.tsv \\
      --fasta-dir fasta_16kb_shorkie/ \\
      --model shorkie \\
      --project bloom_qtl_pairwise

  # Backfill ESM1-v predictions (protein sequences in results)
  python backfill_cache.py \\
      --results esm1v_predictions.tsv \\
      --model esm1v \\
      --protein-column wt_protein \\
      --project bloom_qtl_magestic
        """
    )

    parser.add_argument(
        "--results", "-r",
        type=Path,
        required=True,
        help="Path to TSV file with prediction results"
    )
    parser.add_argument(
        "--fasta-dir", "-f",
        type=Path,
        help="Directory containing FASTA files with sequences"
    )
    parser.add_argument(
        "--model", "-m",
        type=str,
        required=True,
        choices=list(MODEL_CONFIGS.keys()),
        help="Model name"
    )
    parser.add_argument(
        "--project", "-p",
        type=str,
        required=True,
        help="Project name for metadata"
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help=f"Cache directory (default: {DEFAULT_CACHE_DIR})"
    )
    parser.add_argument(
        "--protein-column",
        type=str,
        help="For ESM1-v: column containing protein sequence"
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Batch size for insertions (default: 1000)"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without making changes"
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.results.exists():
        print(f"Error: Results file not found: {args.results}")
        sys.exit(1)

    if args.fasta_dir and not args.fasta_dir.exists():
        print(f"Error: FASTA directory not found: {args.fasta_dir}")
        sys.exit(1)

    if args.model == "esm1v" and not args.protein_column:
        print("Error: ESM1-v requires --protein-column")
        sys.exit(1)

    if args.model != "esm1v" and not args.fasta_dir:
        print(f"Error: {args.model} requires --fasta-dir")
        sys.exit(1)

    # Run backfill
    n_backfilled = backfill_from_results(
        results_tsv=args.results,
        fasta_dir=args.fasta_dir,
        model=args.model,
        project_name=args.project,
        cache_dir=args.cache_dir,
        protein_column=args.protein_column,
        batch_size=args.batch_size,
        dry_run=args.dry_run
    )

    print(f"\nDone! Backfilled {n_backfilled} predictions.")


if __name__ == "__main__":
    main()
