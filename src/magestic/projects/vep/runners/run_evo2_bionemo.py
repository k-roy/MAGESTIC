#!/usr/bin/env python3
"""
Unified Evo2 scoring using BioNeMo with VEP cache integration.

Supports both Evo2-7B (single GPU) and Evo2-40B (4x GPU tensor-parallel).
Uses sequence hash-based caching to avoid redundant computation.

Usage:
    # Score with 7B model (single GPU)
    python run_evo2_bionemo.py \
        --input sequences.fasta \
        --output scores.tsv \
        --model 7b

    # Score with 40B model (4x GPU)
    python run_evo2_bionemo.py \
        --input sequences.fasta \
        --output scores.tsv \
        --model 40b \
        --tensor-parallel-size 4

Author: Kevin R. Roy
Date: 2026-03-15
"""
import argparse
import json
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from magestic.projects.vep.core.vep_cache import VEPCache

# Configuration
BASE_DIR = Path("/path/to")
BIONEMO_DIR = BASE_DIR / "software/bionemo"
CACHE_DIR = BASE_DIR / "common/vep_cache"

CHECKPOINTS = {
    "7b": BIONEMO_DIR / "checkpoints/evo2_7b_nemo2",
    "40b": BIONEMO_DIR / "checkpoints/evo2_40b_nemo2",
}

# Model parameters for cache key
MODEL_PARAMS = {
    "7b": {
        "model_size": "7b",
        "seq_length": 50000,
        "reduce_method": "mean",
        "framework": "bionemo",
        "version": "2.7.1"
    },
    "40b": {
        "model_size": "40b",
        "seq_length": 50000,
        "reduce_method": "mean",
        "framework": "bionemo",
        "version": "2.7.1"
    }
}


def parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """Parse FASTA file into list of (seq_id, sequence) tuples."""
    sequences = []
    current_id = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append((current_id, ''.join(current_seq)))
                # Extract first word as seq_id
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_id is not None:
            sequences.append((current_id, ''.join(current_seq)))

    return sequences


def check_cache(
    sequences: List[Tuple[str, str]],
    cache: VEPCache,
    params_hash: str
) -> Tuple[Dict[str, float], List[Tuple[str, str]]]:
    """Check cache for existing predictions.

    Returns:
        cached_scores: Dict mapping seq_id -> score for cached sequences
        uncached: List of (seq_id, sequence) tuples that need scoring
    """
    # Compute sequence hashes
    seq_hashes = [(seq_id, seq, cache.hash_sequence(seq)) for seq_id, seq in sequences]

    # Batch lookup
    all_hashes = [h for _, _, h in seq_hashes]
    cache_results = cache.get_batch(all_hashes, params_hash)

    cached_scores = {}
    uncached = []

    for seq_id, seq, seq_hash in seq_hashes:
        result = cache_results.get(seq_hash)
        if result is not None:
            cached_scores[seq_id] = result["score"]
        else:
            uncached.append((seq_id, seq, seq_hash))

    return cached_scores, uncached


def run_bionemo_inference(
    sequences: List[Tuple[str, str, str]],  # (seq_id, sequence, seq_hash)
    model_size: str,
    tensor_parallel_size: int,
    cache: VEPCache,
    params_hash: str,
    source_project: str
) -> Dict[str, float]:
    """Run BioNeMo inference on uncached sequences.

    Note: This function is meant to be called from within a SLURM job
    with the BioNeMo container. The actual inference is done by
    calling bionemo.evo2.run.predict.
    """
    if not sequences:
        return {}

    # Create temporary FASTA for inference
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        temp_fasta = Path(f.name)
        for seq_id, seq, _ in sequences:
            f.write(f">{seq_id}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

    # Create temporary output directory
    with tempfile.TemporaryDirectory() as temp_dir:
        output_dir = Path(temp_dir)

        # Import and run BioNeMo inference
        # This assumes we're running inside the BioNeMo container
        try:
            from bionemo.evo2.run.predict import predict

            ckpt_dir = CHECKPOINTS[model_size]

            predict(
                fasta_path=str(temp_fasta),
                ckpt_dir=str(ckpt_dir),
                output_dir=str(output_dir),
                tensor_parallel_size=tensor_parallel_size,
                model_size=model_size,
                output_log_prob_seqs=True,
                log_prob_collapse_option="mean",
                micro_batch_size=1
            )

            # Read results
            import torch

            pred_file = output_dir / "predictions__rank_0__dp_rank_0.pt"
            seq_map_file = output_dir / "seq_idx_map.json"

            if not pred_file.exists():
                raise RuntimeError(f"Predictions file not found: {pred_file}")

            preds = torch.load(pred_file, map_location='cpu')
            with open(seq_map_file) as f:
                seq_map = json.load(f)

            # Extract scores
            log_probs = preds['log_probs_seqs'].numpy()

            scores = {}
            cache_entries = []

            for seq_id, seq, seq_hash in sequences:
                if seq_id in seq_map:
                    idx = seq_map[seq_id]
                    score = float(log_probs[idx])
                    scores[seq_id] = score

                    # Add to cache
                    cache_entries.append((
                        seq_hash,
                        {"score": score, "seq_length": len(seq)},
                        seq_id
                    ))

            # Batch store to cache
            if cache_entries:
                cache.put_batch(
                    cache_entries,
                    params_hash,
                    source_project,
                    model_version=f"evo2_{model_size}",
                    window_size=50000
                )

            return scores

        finally:
            # Cleanup temp FASTA
            temp_fasta.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(
        description="Score sequences using Evo2 via BioNeMo with VEP cache"
    )
    parser.add_argument(
        '--input', '-i', required=True,
        help='Input FASTA file'
    )
    parser.add_argument(
        '--output', '-o', required=True,
        help='Output TSV file with seq_id and score columns'
    )
    parser.add_argument(
        '--model', '-m', required=True, choices=['7b', '40b'],
        help='Model size: 7b (single GPU) or 40b (4x GPU)'
    )
    parser.add_argument(
        '--tensor-parallel-size', '-t', type=int, default=None,
        help='Tensor parallel size (default: 1 for 7b, 4 for 40b)'
    )
    parser.add_argument(
        '--project', '-p', default='vep_scoring',
        help='Project name for cache tracking (default: vep_scoring)'
    )
    parser.add_argument(
        '--no-cache', action='store_true',
        help='Disable caching (run all sequences)'
    )
    parser.add_argument(
        '--cache-only', action='store_true',
        help='Only return cached results, skip inference'
    )
    args = parser.parse_args()

    # Set default tensor parallel size
    if args.tensor_parallel_size is None:
        args.tensor_parallel_size = 4 if args.model == '40b' else 1

    # Initialize
    input_path = Path(args.input)
    output_path = Path(args.output)

    print(f"=== Evo2-{args.model.upper()} Scoring with VEP Cache ===", file=sys.stderr)
    print(f"Input: {input_path}", file=sys.stderr)
    print(f"Output: {output_path}", file=sys.stderr)
    print(f"Model: evo2_{args.model}", file=sys.stderr)
    print(f"Tensor parallel: {args.tensor_parallel_size}", file=sys.stderr)
    print(f"Project: {args.project}", file=sys.stderr)
    print("", file=sys.stderr)

    # Parse input
    sequences = parse_fasta(input_path)
    print(f"Loaded {len(sequences)} sequences", file=sys.stderr)

    # Initialize cache
    cache = VEPCache(CACHE_DIR, f"evo2_{args.model}")
    params = MODEL_PARAMS[args.model]
    params_hash = cache.hash_params(params)

    # Check cache
    if args.no_cache:
        cached_scores = {}
        uncached = [(seq_id, seq, cache.hash_sequence(seq)) for seq_id, seq in sequences]
    else:
        cached_scores, uncached = check_cache(sequences, cache, params_hash)

    print(f"Cache hits: {len(cached_scores)}", file=sys.stderr)
    print(f"Need scoring: {len(uncached)}", file=sys.stderr)

    # Run inference on uncached sequences
    if uncached and not args.cache_only:
        print(f"\nRunning inference on {len(uncached)} sequences...", file=sys.stderr)
        t0 = time.time()

        new_scores = run_bionemo_inference(
            uncached,
            args.model,
            args.tensor_parallel_size,
            cache,
            params_hash,
            args.project
        )

        elapsed = time.time() - t0
        print(f"Inference completed in {elapsed:.1f}s ({len(new_scores)/elapsed:.2f} seq/s)", file=sys.stderr)

        # Merge with cached scores
        all_scores = {**cached_scores, **new_scores}
    else:
        all_scores = cached_scores

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write("seq_id\tscore\n")
        for seq_id, seq in sequences:
            score = all_scores.get(seq_id, float('nan'))
            f.write(f"{seq_id}\t{score}\n")

    print(f"\nOutput written to: {output_path}", file=sys.stderr)
    print(f"Total sequences: {len(sequences)}", file=sys.stderr)
    print(f"Scored: {len(all_scores)}", file=sys.stderr)


if __name__ == '__main__':
    main()
