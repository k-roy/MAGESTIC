"""
Abstract Base Class for VEP Models

All VEP model runners should inherit from VEPModel and implement
the abstract methods.

Author: Kevin R. Roy
Date: 2026-02-20
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import pandas as pd


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

    def __init__(self, model_name: str, device: str = "cuda:0"):
        """
        Initialize VEP model.

        Args:
            model_name: Human-readable model name (e.g., "evo2_7b", "yorzoi")
            device: PyTorch device string
        """
        self.model_name = model_name
        self.device = device
        self.model = None
        self._is_loaded = False

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
        **kwargs
    ) -> List[PredictionResult]:
        """
        Run prediction on a batch of sequences.

        Default implementation iterates over predict_sequence().
        Override for models with native batching support.

        Args:
            sequences: {test_id: sequence} mapping
            **kwargs: Model-specific parameters

        Returns:
            List of PredictionResult objects
        """
        if not self._is_loaded:
            self.load_model()

        results = []
        for test_id, sequence in sequences.items():
            result = self.predict_sequence(sequence, test_id, **kwargs)
            results.append(result)

        return results

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
            script += f'''# Environment setup
eval "$(/path/to/anaconda3/bin/conda shell.bash hook)"
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
