"""
Abstract Base Class for VEP Models

All VEP model runners should inherit from VEPModel and implement
the abstract methods.

Includes sequence-hash based caching to avoid recomputing predictions
for identical sequences with different IDs.

Author: Kevin R. Roy
Date: 2026-02-20
"""

from abc import ABC, abstractmethod
from collections import defaultdict
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import pandas as pd


def hash_sequence(seq: str) -> str:
    """Create MD5 hash of sequence content (case-insensitive)."""
    return hashlib.md5(seq.upper().strip().encode()).hexdigest()


@dataclass
class SequenceCache:
    """
    Cache for sequence predictions based on sequence content hash.

    This prevents recomputing predictions for identical sequences that
    have different IDs (e.g., same variant in different genomic windows).

    Usage:
        cache = SequenceCache()

        # Process sequences with deduplication
        for test_id, sequence in sequences.items():
            seq_hash = cache.register(test_id, sequence)

            if cache.has_result(seq_hash):
                # Use cached result
                result = cache.get_result(seq_hash, test_id)
            else:
                # Compute and cache
                result = model.predict(sequence)
                cache.set_result(seq_hash, result)
    """

    # hash -> canonical test_id (first occurrence)
    _hash_to_canonical: Dict[str, str] = field(default_factory=dict)

    # hash -> [all test_ids with this sequence]
    _hash_to_ids: Dict[str, List[str]] = field(default_factory=lambda: defaultdict(list))

    # hash -> sequence (for verification)
    _hash_to_seq: Dict[str, str] = field(default_factory=dict)

    # hash -> cached PredictionResult
    _hash_to_result: Dict[str, 'PredictionResult'] = field(default_factory=dict)

    def register(self, test_id: str, sequence: str) -> str:
        """
        Register a sequence and return its hash.

        Tracks all test_ids that share the same sequence.
        """
        seq_hash = hash_sequence(sequence)
        self._hash_to_ids[seq_hash].append(test_id)

        if seq_hash not in self._hash_to_canonical:
            self._hash_to_canonical[seq_hash] = test_id
            self._hash_to_seq[seq_hash] = sequence

        return seq_hash

    def has_result(self, seq_hash: str) -> bool:
        """Check if result is cached for this sequence hash."""
        return seq_hash in self._hash_to_result

    def get_result(self, seq_hash: str, test_id: str) -> 'PredictionResult':
        """Get cached result, updating the test_id."""
        cached = self._hash_to_result[seq_hash]
        # Return a copy with the requested test_id
        return PredictionResult(
            test_id=test_id,
            score=cached.score,
            prediction=cached.prediction,
            extra_metrics=cached.extra_metrics.copy() if cached.extra_metrics else None
        )

    def set_result(self, seq_hash: str, result: 'PredictionResult') -> None:
        """Cache a result for this sequence hash."""
        self._hash_to_result[seq_hash] = result

    def get_unique_sequences(self) -> Dict[str, str]:
        """Return {canonical_id: sequence} for unique sequences only."""
        return {
            self._hash_to_canonical[h]: self._hash_to_seq[h]
            for h in self._hash_to_canonical
        }

    def expand_results(self, results: List['PredictionResult']) -> List['PredictionResult']:
        """
        Expand results from unique sequences back to all original IDs.

        Takes predictions indexed by canonical_id and duplicates them
        for all other IDs that share the same sequence.
        """
        # Build lookup by canonical_id
        result_by_canonical = {r.test_id: r for r in results}

        expanded = []
        for seq_hash, ids in self._hash_to_ids.items():
            canonical = self._hash_to_canonical[seq_hash]
            if canonical not in result_by_canonical:
                continue

            base_result = result_by_canonical[canonical]
            for test_id in ids:
                expanded.append(PredictionResult(
                    test_id=test_id,
                    score=base_result.score,
                    prediction=base_result.prediction,
                    extra_metrics=base_result.extra_metrics.copy() if base_result.extra_metrics else None
                ))

        return expanded

    def get_stats(self) -> Dict[str, int]:
        """Return deduplication statistics."""
        total_ids = sum(len(ids) for ids in self._hash_to_ids.values())
        unique_seqs = len(self._hash_to_canonical)
        return {
            "total_ids": total_ids,
            "unique_sequences": unique_seqs,
            "duplicates": total_ids - unique_seqs,
            "cached_results": len(self._hash_to_result)
        }

    def save_mapping(self, path: Path) -> None:
        """Save ID mapping to TSV file for result expansion."""
        with open(path, 'w') as f:
            f.write("seq_hash\tcanonical_id\tn_duplicates\tall_ids\n")
            for seq_hash, ids in sorted(self._hash_to_ids.items()):
                canonical = self._hash_to_canonical[seq_hash]
                f.write(f"{seq_hash}\t{canonical}\t{len(ids)}\t{json.dumps(ids)}\n")


@dataclass
class SLURMConfig:
    """SLURM resource configuration for a VEP model."""

    partition: str = "gpu"
    gres: str = "gpu:1"
    constraint: Optional[str] = None
    mem_gb: int = 64
    cpus: int = 8
    time_hours: int = 24

    # Container settings (for Singularity-based models like Evo2)
    use_container: bool = False
    container_path: Optional[Path] = None

    # Conda environment (for conda-based models)
    conda_env: Optional[str] = None

    def to_slurm_header(self, job_name: str, log_dir: Path) -> str:
        """Generate SLURM header directives."""
        lines = [
            "#!/bin/bash",
            f"#SBATCH --job-name={job_name}",
            f"#SBATCH --partition={self.partition}",
            f"#SBATCH --gres={self.gres}",
        ]

        if self.constraint:
            lines.append(f'#SBATCH --constraint="{self.constraint}"')

        lines.extend([
            f"#SBATCH --mem={self.mem_gb}G",
            f"#SBATCH --cpus-per-task={self.cpus}",
            f"#SBATCH --time={self.time_hours}:00:00",
            f"#SBATCH --output={log_dir}/{job_name}_%j.log",
            f"#SBATCH --error={log_dir}/{job_name}_%j.err",
        ])

        return "\n".join(lines)


@dataclass
class PredictionResult:
    """Standardized prediction result from any VEP model."""

    test_id: str
    score: float
    prediction: Optional[str] = None  # "deleterious", "neutral", etc.

    # Model-specific metrics (optional)
    extra_metrics: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for DataFrame construction."""
        result = {
            "test_id": self.test_id,
            "score": self.score,
        }
        if self.prediction:
            result["prediction"] = self.prediction
        if self.extra_metrics:
            result.update(self.extra_metrics)
        return result


class VEPModel(ABC):
    """
    Abstract base class for Variant Effect Prediction models.

    Subclasses must implement:
        - load_model(): Load model weights
        - predict_sequence(): Run prediction on a single sequence
        - get_slurm_config(): Return SLURM resource requirements

    Optional overrides:
        - predict_batch(): Batch prediction (default: iterate over predict_sequence)
        - cleanup(): Release model resources
    """

    def __init__(self, model_name: str, device: str = "cuda:0", use_cache: bool = True):
        """
        Initialize VEP model.

        Args:
            model_name: Human-readable model name (e.g., "evo2_7b", "yorzoi")
            device: PyTorch device string
            use_cache: If True, deduplicate sequences by content hash to avoid
                       recomputing predictions for identical sequences with different IDs
        """
        self.model_name = model_name
        self.device = device
        self.model = None
        self._is_loaded = False
        self.use_cache = use_cache
        self._cache: Optional[SequenceCache] = None

    @abstractmethod
    def load_model(self) -> None:
        """
        Load model weights.

        Should set self.model and self._is_loaded = True.
        """
        pass

    @abstractmethod
    def predict_sequence(
        self,
        sequence: str,
        test_id: str,
        **kwargs
    ) -> PredictionResult:
        """
        Run prediction on a single sequence.

        Args:
            sequence: DNA or protein sequence
            test_id: Identifier for this test
            **kwargs: Model-specific parameters

        Returns:
            PredictionResult with score and optional prediction label
        """
        pass

    @abstractmethod
    def get_slurm_config(self) -> SLURMConfig:
        """
        Return SLURM resource requirements for this model.

        Returns:
            SLURMConfig with appropriate resources
        """
        pass

    def predict_batch(
        self,
        sequences: Dict[str, str],
        expand_duplicates: bool = True,
        **kwargs
    ) -> List[PredictionResult]:
        """
        Run prediction on a batch of sequences with automatic deduplication.

        Uses sequence-hash caching to avoid recomputing predictions for
        identical sequences with different IDs. This can dramatically reduce
        compute time when the same variant appears in multiple genomic windows.

        Args:
            sequences: {test_id: sequence} mapping
            expand_duplicates: If True, expand results to include all original IDs.
                              If False, return only unique sequence results.
            **kwargs: Model-specific parameters

        Returns:
            List of PredictionResult objects (one per input ID if expand_duplicates=True)
        """
        if not self._is_loaded:
            self.load_model()

        if not self.use_cache:
            # No caching - process all sequences
            results = []
            for test_id, sequence in sequences.items():
                result = self.predict_sequence(sequence, test_id, **kwargs)
                results.append(result)
            return results

        # Initialize cache and register all sequences
        self._cache = SequenceCache()
        for test_id, sequence in sequences.items():
            self._cache.register(test_id, sequence)

        # Get statistics
        stats = self._cache.get_stats()
        if stats["duplicates"] > 0:
            print(f"Deduplication: {stats['total_ids']} sequences -> {stats['unique_sequences']} unique "
                  f"({stats['duplicates']} duplicates, {100*stats['duplicates']/stats['total_ids']:.0f}% reduction)")

        # Process only unique sequences
        unique_results = []
        unique_sequences = self._cache.get_unique_sequences()

        for test_id, sequence in unique_sequences.items():
            result = self.predict_sequence(sequence, test_id, **kwargs)
            unique_results.append(result)

        # Expand back to all IDs if requested
        if expand_duplicates and stats["duplicates"] > 0:
            return self._cache.expand_results(unique_results)
        else:
            return unique_results

    def save_dedup_mapping(self, path: Path) -> None:
        """Save the deduplication mapping to a file for later result expansion."""
        if self._cache:
            self._cache.save_mapping(path)
        else:
            raise ValueError("No cache available. Run predict_batch first with use_cache=True.")

    def predict_to_dataframe(
        self,
        sequences: Dict[str, str],
        **kwargs
    ) -> pd.DataFrame:
        """
        Run predictions and return as DataFrame.

        Args:
            sequences: {test_id: sequence} mapping
            **kwargs: Model-specific parameters

        Returns:
            DataFrame with columns: test_id, score, [prediction], [extra_metrics...]
        """
        results = self.predict_batch(sequences, **kwargs)
        return pd.DataFrame([r.to_dict() for r in results])

    def cleanup(self) -> None:
        """
        Release model resources.

        Override for models that need explicit cleanup (e.g., GPU memory).
        """
        self.model = None
        self._is_loaded = False

    def __enter__(self):
        """Context manager entry."""
        self.load_model()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.cleanup()
        return False

    def generate_slurm_script(
        self,
        input_file: Path,
        output_file: Path,
        job_name: Optional[str] = None,
        log_dir: Optional[Path] = None,
        **extra_params
    ) -> str:
        """
        Generate SLURM submission script for this model.

        Args:
            input_file: Path to input sequences
            output_file: Path for output predictions
            job_name: SLURM job name (default: model_name)
            log_dir: Directory for logs (default: output_file.parent/logs)
            **extra_params: Additional parameters for the script

        Returns:
            Complete SLURM script as string
        """
        slurm_config = self.get_slurm_config()
        job_name = job_name or f"{self.model_name}_pred"
        log_dir = log_dir or output_file.parent / "logs"

        # Generate header
        script = slurm_config.to_slurm_header(job_name, log_dir)
        script += "\n\n"

        # Add environment setup
        if slurm_config.conda_env:
            script += f'''# Environment setup (edit conda path for your HPC)
eval "$(conda shell.bash hook)"
conda activate {slurm_config.conda_env}

echo "Python: $(which python)"
echo "Python version: $(python --version)"
'''

        script += f'''
echo "=== {self.model_name} Prediction ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
echo "Start time: $(date)"
echo ""
echo "Input: {input_file}"
echo "Output: {output_file}"
echo ""

mkdir -p {log_dir}
mkdir -p {output_file.parent}

'''

        # Model-specific run command (should be overridden in subclasses)
        script += self._generate_run_command(input_file, output_file, **extra_params)

        script += f'''
echo ""
echo "End time: $(date)"
'''

        return script

    def _generate_run_command(
        self,
        input_file: Path,
        output_file: Path,
        **kwargs
    ) -> str:
        """
        Generate the model-specific run command.

        Override in subclasses to provide the actual command.
        """
        return f"# TODO: Add {self.model_name} run command\n"
