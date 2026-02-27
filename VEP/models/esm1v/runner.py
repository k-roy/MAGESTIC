"""
ESM1-v Model Runner

ESM1-v is a protein language model for scoring missense variant effects.
It uses masked language modeling to compute the effect of amino acid
substitutions on protein fitness.

Requirements:
    - conda environment: esm
    - PyTorch with CUDA support
    - fair-esm package

Author: Kevin R. Roy
Date: 2026-02-20
"""

import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from ..base import VEPModel, SLURMConfig, PredictionResult


# Default paths
BASE_DIR = Path("/path/to")

# ESM1-v thresholds for pathogenicity
DELETERIOUS_THRESHOLD = -7.5  # Below this = deleterious
# Note: ESM1-v scores are negative log-odds, more negative = more deleterious


class ESM1vModel(VEPModel):
    """
    ESM1-v protein variant effect prediction model.

    ESM1-v uses masked language modeling to score the effect of
    amino acid substitutions. The delta log-likelihood between
    wild-type and mutant amino acids indicates pathogenicity.

    Usage:
        model = ESM1vModel()
        model.load_model()
        result = model.predict_variant(wt_protein, position, wt_aa, mut_aa, "test_1")
    """

    def __init__(
        self,
        model_name: str = "esm1v_t33_650M_UR90S_1",
        device: str = "cuda:0",
        deleterious_threshold: float = DELETERIOUS_THRESHOLD,
    ):
        super().__init__(model_name="esm1v", device=device)
        self.esm_model_name = model_name
        self.deleterious_threshold = deleterious_threshold
        self.alphabet = None
        self.batch_converter = None

    def load_model(self) -> None:
        """Load ESM1-v model and alphabet."""
        import torch
        import esm

        print(f"Loading ESM1-v model: {self.esm_model_name}...")
        start = time.time()

        self.model, self.alphabet = esm.pretrained.load_model_and_alphabet(
            self.esm_model_name
        )
        self.model = self.model.to(self.device)
        self.model.eval()
        self.batch_converter = self.alphabet.get_batch_converter()

        print(f"  Model loaded in {time.time() - start:.1f}s")
        self._is_loaded = True

    def predict_variant(
        self,
        protein_sequence: str,
        position: int,
        wt_aa: str,
        mut_aa: str,
        test_id: str,
        use_sliding_window: bool = False,
        window_size: int = 1022,
        **kwargs
    ) -> PredictionResult:
        """
        Score a single missense variant.

        Args:
            protein_sequence: Full protein sequence
            position: 1-based position of the mutation
            wt_aa: Wild-type amino acid
            mut_aa: Mutant amino acid
            test_id: Identifier for this test
            use_sliding_window: Use sliding window for long proteins
            window_size: Size of sliding window (ESM1-v max is 1022)

        Returns:
            PredictionResult with delta logits score
        """
        import torch

        if not self._is_loaded:
            self.load_model()

        # Validate input
        if position < 1 or position > len(protein_sequence):
            raise ValueError(f"Position {position} out of range [1, {len(protein_sequence)}]")

        # 0-indexed position
        pos_idx = position - 1

        # Verify WT amino acid matches
        if protein_sequence[pos_idx] != wt_aa:
            raise ValueError(
                f"Expected {wt_aa} at position {position}, "
                f"found {protein_sequence[pos_idx]}"
            )

        # Determine if we need sliding window
        used_window = False
        if len(protein_sequence) > window_size:
            if use_sliding_window:
                # Extract window centered on mutation
                window_start = max(0, pos_idx - window_size // 2)
                window_end = min(len(protein_sequence), window_start + window_size)
                window_start = max(0, window_end - window_size)

                window_seq = protein_sequence[window_start:window_end]
                window_pos = pos_idx - window_start
                used_window = True
            else:
                # Truncate to first window_size residues
                window_seq = protein_sequence[:window_size]
                window_pos = pos_idx
                if pos_idx >= window_size:
                    raise ValueError(
                        f"Position {position} outside window. "
                        f"Use use_sliding_window=True for long proteins."
                    )
        else:
            window_seq = protein_sequence
            window_pos = pos_idx

        # Prepare sequence for ESM
        data = [("protein", window_seq)]
        _, _, batch_tokens = self.batch_converter(data)
        batch_tokens = batch_tokens.to(self.device)

        # Mask the target position (+1 for BOS token)
        masked_tokens = batch_tokens.clone()
        masked_tokens[0, window_pos + 1] = self.alphabet.mask_idx

        # Get model predictions
        with torch.no_grad():
            results = self.model(masked_tokens, repr_layers=[])
            logits = results["logits"]

        # Get logits at masked position
        position_logits = logits[0, window_pos + 1]

        # Get indices for WT and mutant amino acids
        wt_idx = self.alphabet.get_idx(wt_aa)
        mut_idx = self.alphabet.get_idx(mut_aa)

        # Delta logits: positive = neutral, negative = deleterious
        wt_logit = position_logits[wt_idx].item()
        mut_logit = position_logits[mut_idx].item()
        delta_logits = mut_logit - wt_logit

        # Classify based on threshold
        prediction = "deleterious" if delta_logits < self.deleterious_threshold else "neutral"

        return PredictionResult(
            test_id=test_id,
            score=delta_logits,
            prediction=prediction,
            extra_metrics={
                "wt_logit": wt_logit,
                "mut_logit": mut_logit,
                "protein_length": len(protein_sequence),
                "used_sliding_window": used_window,
                "variant_desc": f"{wt_aa}{position}{mut_aa}",
            }
        )

    def predict_sequence(
        self,
        sequence: str,
        test_id: str,
        **kwargs
    ) -> PredictionResult:
        """
        Not applicable for ESM1-v - use predict_variant instead.

        Raises:
            NotImplementedError
        """
        raise NotImplementedError(
            "ESM1-v scores variants, not sequences. "
            "Use predict_variant(protein_seq, position, wt_aa, mut_aa, test_id) instead."
        )

    def predict_variants_batch(
        self,
        variants: List[Dict[str, Any]],
        batch_size: int = 16,
        use_sliding_window: bool = True,
        window_size: int = 1022,
    ) -> List[PredictionResult]:
        """
        Score multiple missense variants in batches for ~4x throughput.

        FIX (2026-02-23): Added batch processing to improve GPU utilization.
        Previous implementation processed variants one at a time.

        Args:
            variants: List of dicts with keys:
                - protein_sequence: Full protein sequence
                - position: 1-based position of mutation
                - wt_aa: Wild-type amino acid
                - mut_aa: Mutant amino acid
                - test_id: Identifier for this test
            batch_size: Number of variants to process in parallel (default 16)
            use_sliding_window: Use sliding window for long proteins
            window_size: Size of sliding window (ESM1-v max is 1022)

        Returns:
            List of PredictionResults in same order as input
        """
        import torch

        if not self._is_loaded:
            self.load_model()

        results = []

        # Process in batches
        for batch_start in range(0, len(variants), batch_size):
            batch_end = min(batch_start + batch_size, len(variants))
            batch_variants = variants[batch_start:batch_end]

            # Prepare batch data
            batch_data = []
            batch_meta = []  # Store metadata for each variant

            for var in batch_variants:
                protein_seq = var['protein_sequence']
                position = var['position']
                pos_idx = position - 1
                wt_aa = var['wt_aa']
                mut_aa = var['mut_aa']
                test_id = var['test_id']

                # Validate
                if pos_idx < 0 or pos_idx >= len(protein_seq):
                    batch_meta.append({
                        'error': f"Position {position} out of range [1, {len(protein_seq)}]",
                        'test_id': test_id,
                    })
                    continue

                if protein_seq[pos_idx] != wt_aa:
                    batch_meta.append({
                        'error': f"Expected {wt_aa} at position {position}, found {protein_seq[pos_idx]}",
                        'test_id': test_id,
                    })
                    continue

                # Handle long proteins with sliding window
                used_window = False
                if len(protein_seq) > window_size:
                    if use_sliding_window:
                        window_start = max(0, pos_idx - window_size // 2)
                        window_end = min(len(protein_seq), window_start + window_size)
                        window_start = max(0, window_end - window_size)
                        window_seq = protein_seq[window_start:window_end]
                        window_pos = pos_idx - window_start
                        used_window = True
                    else:
                        window_seq = protein_seq[:window_size]
                        window_pos = pos_idx
                        if pos_idx >= window_size:
                            batch_meta.append({
                                'error': f"Position {position} outside window",
                                'test_id': test_id,
                            })
                            continue
                else:
                    window_seq = protein_seq
                    window_pos = pos_idx

                batch_data.append((f"protein_{len(batch_data)}", window_seq))
                batch_meta.append({
                    'test_id': test_id,
                    'window_pos': window_pos,
                    'wt_aa': wt_aa,
                    'mut_aa': mut_aa,
                    'protein_length': len(protein_seq),
                    'used_window': used_window,
                    'position': position,
                })

            # Skip if no valid variants in batch
            if not batch_data:
                for meta in batch_meta:
                    if 'error' in meta:
                        results.append(PredictionResult(
                            test_id=meta['test_id'],
                            score=float('nan'),
                            prediction='error',
                            extra_metrics={'error': meta['error']},
                        ))
                continue

            # Convert batch to tokens
            _, _, batch_tokens = self.batch_converter(batch_data)
            batch_tokens = batch_tokens.to(self.device)

            # Create masked versions for each variant
            masked_tokens = batch_tokens.clone()
            valid_indices = []
            for i, meta in enumerate(batch_meta):
                if 'error' not in meta:
                    masked_tokens[len(valid_indices), meta['window_pos'] + 1] = self.alphabet.mask_idx
                    valid_indices.append(i)

            # Get predictions for entire batch
            with torch.no_grad():
                logits = self.model(masked_tokens, repr_layers=[])["logits"]

            # Extract results for each variant
            result_idx = 0
            for i, meta in enumerate(batch_meta):
                if 'error' in meta:
                    results.append(PredictionResult(
                        test_id=meta['test_id'],
                        score=float('nan'),
                        prediction='error',
                        extra_metrics={'error': meta['error']},
                    ))
                else:
                    position_logits = logits[result_idx, meta['window_pos'] + 1]

                    wt_idx = self.alphabet.get_idx(meta['wt_aa'])
                    mut_idx = self.alphabet.get_idx(meta['mut_aa'])

                    wt_logit = position_logits[wt_idx].item()
                    mut_logit = position_logits[mut_idx].item()
                    delta_logits = mut_logit - wt_logit

                    prediction = "deleterious" if delta_logits < self.deleterious_threshold else "neutral"

                    results.append(PredictionResult(
                        test_id=meta['test_id'],
                        score=delta_logits,
                        prediction=prediction,
                        extra_metrics={
                            'wt_logit': wt_logit,
                            'mut_logit': mut_logit,
                            'protein_length': meta['protein_length'],
                            'used_sliding_window': meta['used_window'],
                            'variant_desc': f"{meta['wt_aa']}{meta['position']}{meta['mut_aa']}",
                        }
                    ))
                    result_idx += 1

            # Clear GPU memory between batches
            torch.cuda.empty_cache()

        return results

    def get_slurm_config(self) -> SLURMConfig:
        """Return SLURM resource requirements for ESM1-v."""
        return SLURMConfig(
            partition="gpu",
            gres="gpu:1",
            constraint="GPU_MEM:32GB",
            mem_gb=32,
            cpus=8,
            time_hours=8,
            use_container=False,
            conda_env="esm",
        )

    def _generate_run_command(
        self,
        input_file: Path,
        output_file: Path,
        **kwargs
    ) -> str:
        """Generate ESM1-v run command for SLURM script."""
        scripts_dir = BASE_DIR / "common/scripts/variant_effect_prediction/snakemake"

        return f'''# Run ESM1-v variant scoring
python {scripts_dir}/05_run_esm1v.py \\
    --input {input_file} \\
    --output {output_file}
'''


def run_esm1v_batch(
    variants_df: pd.DataFrame,
    output_file: Path,
    wt_protein_col: str = "wt_protein",
    mut_protein_col: str = "mut_protein",
    position_col: str = "aa_position",
    wt_aa_col: str = "wt_aa",
    mut_aa_col: str = "mut_aa",
    test_id_col: str = "test_id",
    use_sliding_window: bool = True,
    batch_size: int = 16,
    device: str = "cuda:0",
) -> pd.DataFrame:
    """
    Run ESM1-v predictions on a batch of protein variants.

    FIX (2026-02-23): Uses batched prediction for ~4x throughput improvement.
    Previously processed variants one at a time.

    Args:
        variants_df: DataFrame with variant information
        output_file: Output TSV path
        wt_protein_col: Column name for WT protein sequence
        mut_protein_col: Column name for mutant protein sequence (optional)
        position_col: Column name for amino acid position (1-based)
        wt_aa_col: Column name for wild-type amino acid
        mut_aa_col: Column name for mutant amino acid
        test_id_col: Column name for test identifier
        use_sliding_window: Use sliding window for long proteins
        batch_size: Number of variants to process in parallel (default 16)
        device: CUDA device

    Returns:
        DataFrame with prediction results
    """
    model = ESM1vModel(device=device)
    model.load_model()

    # Prepare variants list for batched processing
    variants_list = []
    row_metadata = []  # Store extra columns from input rows

    for idx, row in variants_df.iterrows():
        variants_list.append({
            'protein_sequence': row[wt_protein_col],
            'position': int(row[position_col]),
            'wt_aa': row[wt_aa_col],
            'mut_aa': row[mut_aa_col],
            'test_id': row[test_id_col],
        })

        # Store extra metadata to add back later
        extra = {}
        for col in ['target_locus', 'orf', 'gene', 'ref_strain', 'alt_strain']:
            if col in row:
                extra[col] = row[col]
        row_metadata.append(extra)

    # Run batched predictions with progress bar
    print(f"Processing {len(variants_list)} variants in batches of {batch_size}...")
    from tqdm import tqdm

    all_results = []
    n_batches = (len(variants_list) + batch_size - 1) // batch_size

    for batch_idx in tqdm(range(n_batches), desc="ESM1-v batches"):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, len(variants_list))
        batch_variants = variants_list[batch_start:batch_end]

        batch_results = model.predict_variants_batch(
            variants=batch_variants,
            batch_size=batch_size,
            use_sliding_window=use_sliding_window,
        )

        # Add extra metadata from input rows
        for i, result in enumerate(batch_results):
            result_dict = result.to_dict()
            result_dict.update(row_metadata[batch_start + i])
            all_results.append(result_dict)

    # Save results
    df = pd.DataFrame(all_results)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)

    # Summary
    if 'prediction' in df.columns:
        print(f"\nESM1-v Results Summary:")
        print(f"  Total: {len(df)}")
        print(f"  Deleterious: {(df['prediction'] == 'deleterious').sum()}")
        print(f"  Neutral: {(df['prediction'] == 'neutral').sum()}")
        if (df['prediction'] == 'error').any():
            print(f"  Errors: {(df['prediction'] == 'error').sum()}")

    return df


if __name__ == "__main__":
    print("ESM1-v Model Test")
    print("=" * 50)

    # Test SLURM config
    model = ESM1vModel()
    config = model.get_slurm_config()
    print(f"SLURM config:")
    print(f"  Partition: {config.partition}")
    print(f"  GPU: {config.gres}")
    print(f"  Memory: {config.mem_gb}G")
    print(f"  Conda env: {config.conda_env}")

    print(f"\nDeleterious threshold: {model.deleterious_threshold}")
    print("(Variants with delta_logits < threshold are classified as deleterious)")
