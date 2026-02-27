"""
Yorzoi Model Runner

Yorzoi is a Borzoi-like model for yeast RNA expression prediction.
It uses a tiled prediction approach for sequences longer than its
input window (4992 bp).

Requirements:
    - conda environment: yorzoi_local
    - HuggingFace model: tom-ellis-lab/yorzoi

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
HF_HOME = BASE_DIR / "common/software/huggingface"
TORCH_HOME = BASE_DIR / "common/software/torch_hub"

# Model parameters
MODEL_INPUT_LENGTH = 4992
MODEL_OUTPUT_LENGTH = 3000
MODEL_RESOLUTION = 10
FLANKING = (MODEL_INPUT_LENGTH - MODEL_OUTPUT_LENGTH) // 2  # 996 bp
DEFAULT_STEP_SIZE = 1500

# Edge exclusion for analysis (center 45kb for 50kb sequences)
EDGE_EXCLUSION = 2500
ANALYSIS_WINDOW = 45000


def one_hot_encode(sequence: str) -> "torch.Tensor":
    """One-hot encode DNA sequence for Yorzoi (vectorized for speed)."""
    import torch

    # Convert sequence to numpy array of characters
    seq_array = np.frombuffer(sequence.upper().encode(), dtype='S1')

    # Create mapping array: A=0, C=1, G=2, T=3, others=4
    # Using ord values: A=65, C=67, G=71, T=84
    char_to_idx = np.full(256, 4, dtype=np.int8)
    char_to_idx[65] = 0  # A
    char_to_idx[67] = 1  # C
    char_to_idx[71] = 2  # G
    char_to_idx[84] = 3  # T

    # Convert characters to indices
    indices = char_to_idx[seq_array.view(np.uint8)]

    # Create one-hot encoding
    seq_len = len(sequence)
    encoded = np.zeros((seq_len, 4), dtype=np.float32)

    # Set valid bases
    valid_mask = indices < 4
    valid_positions = np.where(valid_mask)[0]
    encoded[valid_positions, indices[valid_mask]] = 1.0

    # Set N/unknown to 0.25
    invalid_positions = np.where(~valid_mask)[0]
    encoded[invalid_positions, :] = 0.25

    return torch.from_numpy(encoded)


def borzoi_transform_inv(y: "torch.Tensor") -> "torch.Tensor":
    """Inverse of Borzoi training transform."""
    import torch
    expd = torch.where(y <= 384, y, 384 + torch.pow(y - 384, 2))
    x = torch.pow(expd, 1 / 0.75)
    return x


def unbin_tensor(y_binned: "torch.Tensor", resolution: int = MODEL_RESOLUTION) -> "torch.Tensor":
    """Upsample binned predictions to base resolution."""
    y_unbinned = y_binned.repeat_interleave(resolution, dim=-1)
    y_unbinned = y_unbinned / resolution
    return y_unbinned


class YorzoiModel(VEPModel):
    """
    Yorzoi RNA expression prediction model.

    Yorzoi predicts RNA coverage tracks using a tiled approach for
    sequences longer than 4992 bp.

    Usage:
        model = YorzoiModel()
        model.load_model()
        result = model.predict_sequence(sequence, "test_1")
    """

    def __init__(
        self,
        step_size: int = DEFAULT_STEP_SIZE,
        device: str = "cuda:0"
    ):
        super().__init__(model_name="yorzoi", device=device)
        self.step_size = step_size

    def load_model(self) -> None:
        """Load Yorzoi model from HuggingFace."""
        import torch

        print("Loading Yorzoi model from HuggingFace...")
        start = time.time()

        from yorzoi.model.borzoi import Borzoi
        self.model = Borzoi.from_pretrained("tom-ellis-lab/yorzoi")
        self.model = self.model.to(self.device)
        self.model.eval()

        print(f"  Model loaded in {time.time() - start:.1f}s")
        self._is_loaded = True

    def _calculate_tile_positions(
        self,
        seq_length: int,
    ) -> List[Tuple[int, int, int, int]]:
        """Calculate tile positions for full sequence coverage."""
        tiles = []
        output_pos = 0

        while output_pos + MODEL_OUTPUT_LENGTH <= seq_length:
            window_start = output_pos - FLANKING

            if window_start < 0:
                window_start = 0
                output_start = FLANKING
            else:
                output_start = output_pos

            window_end = window_start + MODEL_INPUT_LENGTH
            output_end = output_start + MODEL_OUTPUT_LENGTH

            if window_end > seq_length:
                window_end = seq_length
                window_start = seq_length - MODEL_INPUT_LENGTH
                output_start = window_start + FLANKING
                output_end = min(output_start + MODEL_OUTPUT_LENGTH, seq_length)

            tiles.append((window_start, window_end, output_start, output_end))
            output_pos += self.step_size

            if output_end >= seq_length:
                break

        if tiles and tiles[-1][3] < seq_length - FLANKING:
            window_end = seq_length
            window_start = seq_length - MODEL_INPUT_LENGTH
            output_start = window_start + FLANKING
            output_end = output_start + MODEL_OUTPUT_LENGTH
            tiles.append((window_start, window_end, output_start, output_end))

        return tiles

    def _predict_single_tile(
        self,
        sequence: str,
        window_start: int,
        window_end: int
    ) -> np.ndarray:
        """Run model inference on a single tile."""
        import torch

        window_seq = sequence[window_start:window_end]

        if len(window_seq) < MODEL_INPUT_LENGTH:
            pad_needed = MODEL_INPUT_LENGTH - len(window_seq)
            window_seq = 'N' * (pad_needed // 2) + window_seq + 'N' * (pad_needed - pad_needed // 2)

        encoded = one_hot_encode(window_seq).unsqueeze(0).to(self.device)

        with torch.no_grad():
            with torch.autocast(device_type="cuda"):
                raw_output = self.model(encoded)

        processed = unbin_tensor(borzoi_transform_inv(raw_output), resolution=MODEL_RESOLUTION)
        return processed[0].cpu().numpy()

    def predict_tiled(self, sequence: str) -> Tuple[np.ndarray, np.ndarray]:
        """
        Predict full sequence using tiling approach.

        FIX (2026-02-23): Added GPU memory clearing between tiles to prevent
        memory fragmentation during long predictions (50kb sequences = ~30 tiles).

        Returns:
            Tuple of (plus_strand, minus_strand) predictions
        """
        import torch
        import gc

        seq_length = len(sequence)
        tiles = self._calculate_tile_positions(seq_length)

        plus_accum = np.zeros(seq_length, dtype=np.float64)
        minus_accum = np.zeros(seq_length, dtype=np.float64)
        coverage = np.zeros(seq_length, dtype=np.float64)

        # Handle edge regions
        padded_start = 'N' * FLANKING + sequence[:MODEL_INPUT_LENGTH - FLANKING]
        tile_pred = self._predict_single_tile(padded_start, 0, MODEL_INPUT_LENGTH)
        plus_tile = np.mean(tile_pred[:81], axis=0)
        minus_tile = np.mean(tile_pred[81:], axis=0)
        plus_accum[0:FLANKING] += plus_tile[0:FLANKING]
        minus_accum[0:FLANKING] += minus_tile[0:FLANKING]
        coverage[0:FLANKING] += 1

        padded_end = sequence[-(MODEL_INPUT_LENGTH - FLANKING):] + 'N' * FLANKING
        tile_pred = self._predict_single_tile(padded_end, 0, MODEL_INPUT_LENGTH)
        plus_tile = np.mean(tile_pred[:81], axis=0)
        minus_tile = np.mean(tile_pred[81:], axis=0)
        plus_accum[-FLANKING:] += plus_tile[-FLANKING:]
        minus_accum[-FLANKING:] += minus_tile[-FLANKING:]
        coverage[-FLANKING:] += 1

        for i, (window_start, window_end, output_start, output_end) in enumerate(tiles):
            tile_pred = self._predict_single_tile(sequence, window_start, window_end)

            plus_tile = np.mean(tile_pred[:81], axis=0)
            minus_tile = np.mean(tile_pred[81:], axis=0)

            tile_output_start = output_start - (window_start + FLANKING)
            tile_output_end = tile_output_start + (output_end - output_start)

            plus_accum[output_start:output_end] += plus_tile[tile_output_start:tile_output_end]
            minus_accum[output_start:output_end] += minus_tile[tile_output_start:tile_output_end]
            coverage[output_start:output_end] += 1

            # Clear GPU cache every 10 tiles to prevent memory fragmentation
            if (i + 1) % 10 == 0:
                torch.cuda.empty_cache()
                gc.collect()

        mask = coverage > 0
        plus_avg = np.zeros(seq_length)
        minus_avg = np.zeros(seq_length)
        plus_avg[mask] = plus_accum[mask] / coverage[mask]
        minus_avg[mask] = minus_accum[mask] / coverage[mask]

        # Final cleanup after all tiles
        torch.cuda.empty_cache()
        gc.collect()

        return plus_avg, minus_avg

    def predict_sequence(
        self,
        sequence: str,
        test_id: str,
        wt_prediction: Optional[Tuple[np.ndarray, np.ndarray]] = None,
        **kwargs
    ) -> PredictionResult:
        """
        Run Yorzoi prediction on a single sequence.

        Args:
            sequence: DNA sequence (any length >= MODEL_INPUT_LENGTH)
            test_id: Identifier for this test
            wt_prediction: Optional cached WT prediction (plus, minus) for comparison

        Returns:
            PredictionResult with log2FC score and coverage metrics
        """
        if not self._is_loaded:
            self.load_model()

        # Predict test sequence
        test_plus, test_minus = self.predict_tiled(sequence)

        # Calculate metrics
        if wt_prediction is not None:
            wt_plus, wt_minus = wt_prediction
            metrics = self._calculate_metrics(wt_plus, wt_minus, test_plus, test_minus)
        else:
            metrics = {
                "wt_mean_coverage": 0.0,
                "test_mean_coverage": float((test_plus + test_minus).mean()),
                "log2fc": 0.0,
                "coverage_ratio": 1.0,
            }

        return PredictionResult(
            test_id=test_id,
            score=metrics.get("log2fc", 0.0),
            prediction=None,
            extra_metrics=metrics
        )

    def _calculate_metrics(
        self,
        wt_plus: np.ndarray,
        wt_minus: np.ndarray,
        test_plus: np.ndarray,
        test_minus: np.ndarray,
        edge_exclusion: int = EDGE_EXCLUSION
    ) -> Dict[str, float]:
        """Calculate expression effect metrics from coverage predictions."""
        # Extract center window
        start = edge_exclusion
        end = len(wt_plus) - edge_exclusion

        wt_plus_center = wt_plus[start:end]
        wt_minus_center = wt_minus[start:end]
        test_plus_center = test_plus[start:end]
        test_minus_center = test_minus[start:end]

        # Total coverage (both strands)
        wt_total = wt_plus_center + wt_minus_center
        test_total = test_plus_center + test_minus_center

        # Mean coverage
        wt_mean = np.mean(wt_total)
        test_mean = np.mean(test_total)

        # Log2 fold change (with pseudocount)
        pseudo = 1.0
        log2fc = np.log2((test_mean + pseudo) / (wt_mean + pseudo))
        coverage_ratio = test_mean / (wt_mean + pseudo)

        # Per-position log2FC
        per_pos_log2fc = np.log2((test_total + pseudo) / (wt_total + pseudo))
        max_abs_log2fc = float(np.max(np.abs(per_pos_log2fc)))
        mean_abs_log2fc = float(np.mean(np.abs(per_pos_log2fc)))

        # Strand-specific
        plus_log2fc = np.log2((np.mean(test_plus_center) + pseudo) / (np.mean(wt_plus_center) + pseudo))
        minus_log2fc = np.log2((np.mean(test_minus_center) + pseudo) / (np.mean(wt_minus_center) + pseudo))

        return {
            'wt_mean_coverage': float(wt_mean),
            'test_mean_coverage': float(test_mean),
            'log2fc': float(log2fc),
            'coverage_ratio': float(coverage_ratio),
            'max_abs_log2fc': max_abs_log2fc,
            'mean_abs_log2fc': mean_abs_log2fc,
            'plus_log2fc': float(plus_log2fc),
            'minus_log2fc': float(minus_log2fc),
        }

    def get_slurm_config(self) -> SLURMConfig:
        """Return SLURM resource requirements for Yorzoi."""
        return SLURMConfig(
            partition="gpu",
            gres="gpu:1",
            constraint="GPU_SKU:H100_SXM5&GPU_MEM:80GB",
            mem_gb=64,
            cpus=8,
            time_hours=24,
            use_container=False,
            conda_env="yorzoi_local",
        )

    def _generate_run_command(
        self,
        input_file: Path,
        output_file: Path,
        **kwargs
    ) -> str:
        """Generate Yorzoi run command for SLURM script."""
        scripts_dir = BASE_DIR / "common/scripts/variant_effect_prediction/snakemake"

        return f'''# Set HuggingFace cache directories
export HF_HOME={HF_HOME}
export TRANSFORMERS_CACHE=$HF_HOME/hub
export TORCH_HOME={TORCH_HOME}

# Run Yorzoi predictions
python {scripts_dir}/03_run_yorzoi.py \\
    --input {input_file} \\
    --output {output_file}
'''


def run_yorzoi_batch(
    sequences: Dict[str, str],
    metadata: pd.DataFrame,
    output_file: Path,
    device: str = "cuda:0",
) -> pd.DataFrame:
    """
    Run Yorzoi predictions on a batch of sequences.

    This is the main entry point for batch processing with WT caching.

    Args:
        sequences: {test_id: sequence} mapping
        metadata: DataFrame with test_id, orf, background_strain columns
        output_file: Output TSV path
        device: CUDA device

    Returns:
        DataFrame with prediction results
    """
    model = YorzoiModel(device=device)
    model.load_model()

    # Cache WT predictions
    wt_cache = {}  # (orf, strain) -> (plus, minus)

    results = []
    from tqdm import tqdm

    for _, row in tqdm(metadata.iterrows(), total=len(metadata), desc="Yorzoi"):
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
                wt_plus, wt_minus = model.predict_tiled(wt_seq)
                wt_cache[wt_key] = (wt_plus, wt_minus)
            else:
                print(f"  Warning: WT {wt_id} not found")
                continue

        wt_pred = wt_cache[wt_key]

        # Predict test sequence (skip if it's WT itself)
        if test_type == 'wt':
            test_plus, test_minus = wt_pred
        else:
            test_plus, test_minus = model.predict_tiled(test_seq)

        # Calculate metrics
        metrics = model._calculate_metrics(wt_pred[0], wt_pred[1], test_plus, test_minus)

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
    print("Yorzoi Model Test")
    print("=" * 50)

    # Test SLURM config
    model = YorzoiModel()
    config = model.get_slurm_config()
    print(f"SLURM config:")
    print(f"  Partition: {config.partition}")
    print(f"  GPU: {config.gres}")
    print(f"  Constraint: {config.constraint}")
    print(f"  Memory: {config.mem_gb}G")
    print(f"  Conda env: {config.conda_env}")
