"""
VEP Pipeline Configuration

Provides VEPConfig dataclass for standardized project configuration.

Usage:
    from config import VEPConfig

    # From YAML
    config = VEPConfig.from_yaml("project_config.yaml")

    # Programmatic
    config = VEPConfig(
        project_name="my_project",
        project_dir=Path("/path/to/project"),
        fasta_50kb_dir=Path("/path/to/fasta_50kb"),
        output_dir=Path("/path/to/output"),
    )

Author: Kevin R. Roy
Date: 2026-02-20
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any
import yaml


# Default paths
DEFAULT_CONDA_BASE = Path("/path/to/anaconda3")
DEFAULT_EVO2_CONTAINER = Path("/path/to/common/software/evo2_with_dependencies.sif")


@dataclass
class ModelConfig:
    """Configuration for a single VEP model."""

    enabled: bool = True
    conda_env: Optional[str] = None
    container_path: Optional[Path] = None

    # SLURM resources
    partition: str = "gpu"
    gres: str = "gpu:1"
    constraint: Optional[str] = None
    mem_gb: int = 64
    cpus: int = 8
    time_hours: int = 24

    # Model-specific options
    extra_options: Dict[str, Any] = field(default_factory=dict)


@dataclass
class VEPConfig:
    """
    Configuration for a VEP pipeline run.

    Required fields:
        project_name: Human-readable project name
        project_dir: Base directory for the project

    Input paths (at least one required):
        fasta_50kb_dir: Directory with 50kb FASTA files (for Evo2, Yorzoi)
        fasta_16kb_dir: Directory with 16kb FASTA files (for Shorkie)
        protein_variants_file: TSV with protein variants (for ESM1-v)

    Output:
        output_dir: Base output directory (default: project_dir/vep_results)
    """

    # Required
    project_name: str
    project_dir: Path

    # Input paths
    fasta_50kb_dir: Optional[Path] = None
    fasta_16kb_dir: Optional[Path] = None
    protein_variants_file: Optional[Path] = None
    metadata_file: Optional[Path] = None

    # Output
    output_dir: Optional[Path] = None

    # Models to run (default: all available based on inputs)
    models_to_run: List[str] = field(default_factory=lambda: ["evo2", "yorzoi", "shorkie", "esm1v"])

    # Model-specific configuration
    evo2: ModelConfig = field(default_factory=lambda: ModelConfig(
        conda_env=None,  # Uses container
        container_path=DEFAULT_EVO2_CONTAINER,
        constraint="GPU_SKU:H100_SXM5&GPU_MEM:80GB",
        mem_gb=64,
        time_hours=12,
    ))

    yorzoi: ModelConfig = field(default_factory=lambda: ModelConfig(
        conda_env="yorzoi_local",
        constraint="GPU_SKU:A100_SXM4|GPU_SKU:H100_SXM5",  # Requires Ampere+ for flash-attention
        mem_gb=64,
        time_hours=12,
    ))

    shorkie: ModelConfig = field(default_factory=lambda: ModelConfig(
        conda_env="shorkie",
        mem_gb=64,
        time_hours=12,
    ))

    esm1v: ModelConfig = field(default_factory=lambda: ModelConfig(
        conda_env="esm1v",
        mem_gb=32,
        time_hours=2,
    ))

    # Evo2-specific options
    evo2_model: str = "evo2_7b"  # or "evo2_40b"
    save_position_logprobs: bool = True

    # Processing options
    incremental: bool = False  # Skip already-processed sequences
    batch_size: int = 1

    # SLURM defaults
    default_partition: str = "gpu"

    def __post_init__(self):
        """Validate and set defaults."""
        self.project_dir = Path(self.project_dir)

        # Set default output dir
        if self.output_dir is None:
            self.output_dir = self.project_dir / "vep_results"
        else:
            self.output_dir = Path(self.output_dir)

        # Convert paths
        if self.fasta_50kb_dir:
            self.fasta_50kb_dir = Path(self.fasta_50kb_dir)
        if self.fasta_16kb_dir:
            self.fasta_16kb_dir = Path(self.fasta_16kb_dir)
        if self.protein_variants_file:
            self.protein_variants_file = Path(self.protein_variants_file)
        if self.metadata_file:
            self.metadata_file = Path(self.metadata_file)

        # Validate at least one input
        if not any([self.fasta_50kb_dir, self.fasta_16kb_dir, self.protein_variants_file]):
            raise ValueError("At least one input path required: fasta_50kb_dir, fasta_16kb_dir, or protein_variants_file")

        # Auto-detect which models can run based on available inputs
        self._validate_models()

    def _validate_models(self):
        """Check which models can run based on available inputs."""
        available_models = []

        if self.fasta_50kb_dir and self.fasta_50kb_dir.exists():
            available_models.extend(["evo2", "yorzoi"])

        if self.fasta_16kb_dir and self.fasta_16kb_dir.exists():
            available_models.append("shorkie")

        if self.protein_variants_file and self.protein_variants_file.exists():
            available_models.append("esm1v")

        # Filter requested models to available ones
        self.models_to_run = [m for m in self.models_to_run if m in available_models]

        if not self.models_to_run:
            raise ValueError(
                f"No models can run with provided inputs.\n"
                f"Available inputs: "
                f"50kb={self.fasta_50kb_dir}, "
                f"16kb={self.fasta_16kb_dir}, "
                f"protein={self.protein_variants_file}"
            )

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "VEPConfig":
        """Load configuration from YAML file."""
        yaml_path = Path(yaml_path)

        if not yaml_path.exists():
            raise FileNotFoundError(f"Config file not found: {yaml_path}")

        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        # Handle nested model configs
        for model_name in ["evo2", "yorzoi", "shorkie", "esm1v"]:
            if model_name in data and isinstance(data[model_name], dict):
                data[model_name] = ModelConfig(**data[model_name])

        return cls(**data)

    def to_yaml(self, yaml_path: Path) -> None:
        """Save configuration to YAML file."""
        yaml_path = Path(yaml_path)
        yaml_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to dict, handling nested dataclasses
        data = self._to_dict()

        with open(yaml_path, 'w') as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    def _to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for YAML serialization."""
        result = {}

        for key, value in self.__dict__.items():
            if isinstance(value, Path):
                result[key] = str(value)
            elif isinstance(value, ModelConfig):
                result[key] = {
                    k: str(v) if isinstance(v, Path) else v
                    for k, v in value.__dict__.items()
                }
            else:
                result[key] = value

        return result

    def get_output_path(self, model: str, filename: str) -> Path:
        """Get output path for a model's results."""
        return self.output_dir / f"{model}_results" / filename

    def get_snakemake_config(self) -> Dict[str, Any]:
        """Get configuration dict for Snakemake."""
        return {
            "VEP_DIR": str(Path(__file__).parent),
            "FASTA_50KB_DIR": str(self.fasta_50kb_dir) if self.fasta_50kb_dir else None,
            "FASTA_16KB_DIR": str(self.fasta_16kb_dir) if self.fasta_16kb_dir else None,
            "PROTEIN_VARIANTS_FILE": str(self.protein_variants_file) if self.protein_variants_file else None,
            "METADATA_FILE": str(self.metadata_file) if self.metadata_file else None,
            "OUTPUT_DIR": str(self.output_dir),
            "EVO2_MODEL": self.evo2_model,
            "MODELS_TO_RUN": self.models_to_run,
        }
