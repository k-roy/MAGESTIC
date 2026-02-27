"""
Shorkie Model Runner

Shorkie is an 8-fold ensemble model for predicting RNA expression effects
in yeast. It uses a Borzoi-like architecture with 16kb input sequences.

Requirements:
    - conda environment: shorkie_env
    - Model weights in: common/models/shorkie/model_f{0-7}_best.h5
    - baskerville-yeast package

Author: Kevin R. Roy
Date: 2026-02-20
"""

import json
import os
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from ..base import VEPModel, SLURMConfig, PredictionResult


# Suppress TensorFlow warnings
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"


# Default paths
BASE_DIR = Path("/path/to")
SHORKIE_DIR = BASE_DIR / "software/shorkie-paper"
MODEL_DIR = BASE_DIR / "common/models/shorkie"
PARAMS_FILE = SHORKIE_DIR / "model/shorkie/params.json"

# Model parameters
SEQ_LENGTH = 16384  # Required input length
NUM_FOLDS = 8
NUM_FEATURES = 170  # 4 nucleotides + 165 species + 1 special token

# Edge exclusion for analysis (center 10kb)
EDGE_EXCLUSION = 3192  # (16384 - 10000) / 2
ANALYSIS_WINDOW = 10000


def one_hot_encode(sequence: str, num_features: int = NUM_FEATURES) -> np.ndarray:
    """
    Convert DNA sequence to one-hot encoding for Shorkie (vectorized for speed).

    Shorkie expects 170 features: 4 nucleotides + 165 species + 1 special token.
    Only the first 4 features (nucleotides) are used for DNA encoding.

    Args:
        sequence: DNA sequence string
        num_features: Number of features (default 170)

    Returns:
        One-hot encoded array of shape (seq_len, num_features)
    """
    seq_len = len(sequence)
    one_hot = np.zeros((seq_len, num_features), dtype=np.float32)

    # Convert sequence to numpy array of bytes
    seq_array = np.frombuffer(sequence.upper().encode(), dtype='S1')

    # Create mapping array: A=0, C=1, G=2, T=3, others=4
    char_to_idx = np.full(256, 4, dtype=np.int8)
    char_to_idx[65] = 0  # A
    char_to_idx[67] = 1  # C
    char_to_idx[71] = 2  # G
    char_to_idx[84] = 3  # T

    # Convert characters to indices
    indices = char_to_idx[seq_array.view(np.uint8)]

    # Set valid bases
    valid_mask = indices < 4
    valid_positions = np.where(valid_mask)[0]
    one_hot[valid_positions, indices[valid_mask]] = 1.0

    # Set N/unknown to 0.25 for first 4 features
    invalid_positions = np.where(~valid_mask)[0]
    one_hot[invalid_positions, :4] = 0.25

    return one_hot


class ShorkieModel(VEPModel):
    """
    Shorkie RNA expression prediction model.

    Shorkie is an 8-fold ensemble that predicts RNA coverage tracks.
    It takes 16kb sequences and outputs binned predictions (16bp bins).

    Usage:
        model = ShorkieModel()
        model.load_model()
        result = model.predict_sequence(sequence, "test_1")
    """

    def __init__(
        self,
        model_dir: Path = MODEL_DIR,
        params_file: Path = PARAMS_FILE,
        num_folds: int = NUM_FOLDS,
        device: str = "cuda:0"
    ):
        super().__init__(model_name="shorkie", device=device)
        self.model_dir = Path(model_dir)
        self.params_file = Path(params_file)
        self.num_folds = num_folds
        self.models = []

    def load_model(self) -> None:
        """Load Shorkie model ensemble."""
        try:
            from baskerville import seqnn
        except ImportError:
            raise ImportError(
                "baskerville-yeast not installed.\n"
                "Activate shorkie_env: conda activate shorkie_env"
            )

        # Load model parameters
        with open(self.params_file) as f:
            params = json.load(f)

        params_model = params["model"]
        params_model["num_features"] = 165 + 5  # Species features

        print(f"Loading {self.num_folds} Shorkie model folds...")
        self.models = []

        for fold in range(self.num_folds):
            model_file = self.model_dir / f"model_f{fold}_best.h5"

            if not model_file.exists():
                raise FileNotFoundError(
                    f"Model weights not found: {model_file}\n"
                    f"Model weights should be in: {self.model_dir}"
                )

            print(f"  Loading fold {fold}...", end=" ", flush=True)
            seqnn_model = seqnn.SeqNN(params_model)
            seqnn_model.restore(str(model_file), trunk=False, by_name=True)
            seqnn_model.build_ensemble(True, [0])
            self.models.append(seqnn_model)
            print("done")

        print(f"Loaded {len(self.models)} models")
        self._is_loaded = True

    def _predict_single(self, sequence_one_hot: np.ndarray) -> np.ndarray:
        """
        Predict coverage tracks for a sequence using model ensemble.

        FIX (2026-02-23): Added memory-efficient streaming mode that loads
        one model at a time to reduce GPU memory usage from ~800MB to ~100MB.

        Args:
            sequence_one_hot: One-hot encoded sequence (L, num_features)

        Returns:
            Predicted coverage tracks, averaged across folds (L_out, num_tracks)
        """
        # Add batch dimension
        seq_batch = np.expand_dims(sequence_one_hot, 0)

        if self.models:
            # All models loaded - use fast path
            predictions = []
            for model in self.models:
                pred = model.predict(seq_batch, verbose=0)
                predictions.append(pred)

            # Average across folds
            avg_pred = np.mean(predictions, axis=0)
            return avg_pred[0]
        else:
            # Memory-efficient streaming mode
            return self._predict_streaming(seq_batch)

    def _predict_streaming(self, seq_batch: np.ndarray) -> np.ndarray:
        """
        Memory-efficient prediction that loads one model at a time.

        This reduces peak GPU memory from ~800MB to ~100MB by:
        1. Loading one model fold
        2. Running inference
        3. Releasing the model before loading the next

        Args:
            seq_batch: Batched input (1, L, num_features)

        Returns:
            Averaged predictions across all folds
        """
        import gc
        try:
            import tensorflow as tf
        except ImportError:
            raise ImportError("TensorFlow not available for streaming mode")

        from baskerville import seqnn

        # Load model parameters once
        with open(self.params_file) as f:
            params = json.load(f)
        params_model = params["model"]
        params_model["num_features"] = 165 + 5

        # Accumulate predictions
        accumulated_pred = None
        n_folds = 0

        for fold in range(self.num_folds):
            model_file = self.model_dir / f"model_f{fold}_best.h5"
            if not model_file.exists():
                continue

            # Load single model
            seqnn_model = seqnn.SeqNN(params_model)
            seqnn_model.restore(str(model_file), trunk=False, by_name=True)
            seqnn_model.build_ensemble(True, [0])

            # Run prediction
            pred = seqnn_model.predict(seq_batch, verbose=0)

            # Accumulate
            if accumulated_pred is None:
                accumulated_pred = pred.astype(np.float64)
            else:
                accumulated_pred += pred

            n_folds += 1

            # Release model to free GPU memory
            del seqnn_model
            tf.keras.backend.clear_session()
            gc.collect()

        if accumulated_pred is None or n_folds == 0:
            raise RuntimeError("No model folds loaded successfully")

        # Average across folds
        avg_pred = accumulated_pred / n_folds
        return avg_pred[0].astype(np.float32)

    def predict_sequence(
        self,
        sequence: str,
        test_id: str,
        wt_prediction: Optional[np.ndarray] = None,
        **kwargs
    ) -> PredictionResult:
        """
        Run Shorkie prediction on a single sequence.

        Args:
            sequence: DNA sequence (must be 16384 bp)
            test_id: Identifier for this test
            wt_prediction: Optional cached WT prediction for comparison

        Returns:
            PredictionResult with logSED score and coverage metrics
        """
        if not self._is_loaded:
            self.load_model()

        # Validate sequence length
        if len(sequence) != SEQ_LENGTH:
            raise ValueError(
                f"Sequence length {len(sequence)} != required {SEQ_LENGTH}"
            )

        # Encode and predict
        one_hot = one_hot_encode(sequence)
        prediction = self._predict_single(one_hot)

        # Calculate metrics
        if wt_prediction is not None:
            metrics = self._calculate_metrics(wt_prediction, prediction)
        else:
            # Just return total coverage if no WT comparison
            metrics = {
                "sum_coverage": float(prediction.sum()),
                "logSED": 0.0,
                "coverage_ratio": 1.0,
            }

        return PredictionResult(
            test_id=test_id,
            score=metrics.get("logSED", 0.0),
            prediction=None,  # Shorkie doesn't have discrete categories
            extra_metrics=metrics
        )

    def _calculate_metrics(
        self,
        wt_pred: np.ndarray,
        test_pred: np.ndarray,
        edge_exclusion: int = EDGE_EXCLUSION
    ) -> Dict[str, float]:
        """
        Calculate expression effect metrics from coverage predictions.

        Args:
            wt_pred: WT coverage predictions (L_out, num_tracks)
            test_pred: Test coverage predictions
            edge_exclusion: Positions to exclude from edges (in bp)

        Returns:
            Dict with metric values
        """
        # Shorkie output is binned (16bp bins)
        bin_size = 16
        start_bin = edge_exclusion // bin_size
        end_bin = wt_pred.shape[0] - start_bin

        # Extract center window
        wt_center = wt_pred[start_bin:end_bin, :]
        test_center = test_pred[start_bin:end_bin, :]

        # Sum coverage over gene region and all tracks
        wt_sum = float(wt_center.sum())
        test_sum = float(test_center.sum())

        # logSED = log2(sum_test + 1) - log2(sum_wt + 1)
        logSED = np.log2(test_sum + 1) - np.log2(wt_sum + 1)
        coverage_ratio = test_sum / (wt_sum + 1)

        # Per-position log2FC (sum across tracks first)
        wt_pos = wt_center.sum(axis=1)
        test_pos = test_center.sum(axis=1)
        pseudo = 1.0
        per_pos_log2fc = np.log2((test_pos + pseudo) / (wt_pos + pseudo))
        max_abs_log2fc = float(np.max(np.abs(per_pos_log2fc)))
        mean_abs_log2fc = float(np.mean(np.abs(per_pos_log2fc)))

        return {
            "wt_sum_coverage": wt_sum,
            "test_sum_coverage": test_sum,
            "logSED": float(logSED),
            "coverage_ratio": float(coverage_ratio),
            "max_abs_log2fc": max_abs_log2fc,
            "mean_abs_log2fc": mean_abs_log2fc,
        }

    def get_slurm_config(self) -> SLURMConfig:
        """Return SLURM resource requirements for Shorkie."""
        return SLURMConfig(
            partition="gpu",
            gres="gpu:1",
            constraint="GPU_MEM:32GB",
            mem_gb=32,
            cpus=8,
            time_hours=12,
            use_container=False,
            conda_env="shorkie_env",
        )

    def _generate_run_command(
        self,
        input_file: Path,
        output_file: Path,
        **kwargs
    ) -> str:
        """Generate Shorkie run command for SLURM script."""
        scripts_dir = BASE_DIR / "common/scripts/variant_effect_prediction/snakemake"

        return f'''# Run Shorkie predictions
export TF_FORCE_GPU_ALLOW_GROWTH=true

python {scripts_dir}/04_run_shorkie.py \\
    --input {input_file} \\
    --output {output_file}
'''


def run_shorkie_batch(
    sequences: Dict[str, str],
    metadata: pd.DataFrame,
    output_file: Path,
    model_dir: Path = MODEL_DIR,
    params_file: Path = PARAMS_FILE,
) -> pd.DataFrame:
    """
    Run Shorkie predictions on a batch of sequences.

    This is the main entry point for batch processing with WT caching.

    Args:
        sequences: {test_id: sequence} mapping
        metadata: DataFrame with test_id, orf, background_strain columns
        output_file: Output TSV path
        model_dir: Directory containing model weights
        params_file: Path to params.json

    Returns:
        DataFrame with prediction results
    """
    model = ShorkieModel(model_dir=model_dir, params_file=params_file)
    model.load_model()

    # Cache WT predictions
    wt_cache = {}  # (orf, strain) -> prediction

    results = []
    from tqdm import tqdm

    for _, row in tqdm(metadata.iterrows(), total=len(metadata), desc="Shorkie"):
        test_id = row['test_id']
        orf = row['orf']
        background_strain = row['background_strain']
        test_type = row.get('test_type', 'unknown')

        if test_id not in sequences:
            print(f"  Warning: {test_id} not found in sequences")
            continue

        test_seq = sequences[test_id]

        # Get or compute WT prediction
        wt_key = (orf, background_strain)
        if wt_key not in wt_cache:
            wt_id = f"{orf}_{background_strain}_wt"
            if wt_id in sequences:
                wt_seq = sequences[wt_id]
                wt_one_hot = one_hot_encode(wt_seq)
                wt_pred = model._predict_single(wt_one_hot)
                wt_cache[wt_key] = wt_pred
            else:
                print(f"  Warning: WT {wt_id} not found")
                continue

        wt_pred = wt_cache[wt_key]

        # Predict test sequence (skip if it's WT itself)
        if test_type == 'wt':
            test_pred = wt_pred
        else:
            test_one_hot = one_hot_encode(test_seq)
            test_pred = model._predict_single(test_one_hot)

        # Calculate metrics
        metrics = model._calculate_metrics(wt_pred, test_pred)

        result = {
            'test_id': test_id,
            'orf': orf,
            'gene': row.get('gene', orf),
            'background_strain': background_strain,
            'test_type': test_type,
            'n_variants': row.get('n_variants', 0),
            **metrics
        }
        results.append(result)

    # Save results
    df = pd.DataFrame(results)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep='\t', index=False)

    return df


if __name__ == "__main__":
    print("Shorkie Model Test")
    print("=" * 50)

    # Check if model weights exist
    if not MODEL_DIR.exists():
        print(f"Model directory not found: {MODEL_DIR}")
    else:
        print(f"Model directory: {MODEL_DIR}")
        weights = list(MODEL_DIR.glob("model_f*_best.h5"))
        print(f"Found {len(weights)} model weight files")

    if not PARAMS_FILE.exists():
        print(f"Params file not found: {PARAMS_FILE}")
    else:
        print(f"Params file: {PARAMS_FILE}")

    # Test SLURM config
    model = ShorkieModel()
    config = model.get_slurm_config()
    print(f"\nSLURM config:")
    print(f"  Partition: {config.partition}")
    print(f"  GPU: {config.gres}")
    print(f"  Memory: {config.mem_gb}G")
    print(f"  Conda env: {config.conda_env}")
