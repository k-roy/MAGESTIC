"""
Core VEP Pipeline Utilities

This module provides shared functionality used across all VEP models:

- fasta: FASTA file I/O (via sequence_utils)
- metadata: Test metadata loading and handling
- slurm: SLURM script generation
- results: Result writing and formatting

Author: Kevin R. Roy
Date: 2026-02-20
"""

from pathlib import Path
from typing import Dict, List, Optional, Union
import pandas as pd
import sys

# Handle imports whether running as package or via sys.path
# Import directly from fasta_io submodule (avoids mappy dependency from alignment_variants)
try:
    # Try relative import first (when installed as package)
    from ..sequence_utils.fasta_io import read_fasta, write_fasta, iter_fasta
except ImportError:
    # Fall back to direct file import (when running via sys.path.insert)
    import importlib.util
    _fasta_io_path = Path(__file__).parent.parent / "sequence_utils" / "fasta_io.py"
    _spec = importlib.util.spec_from_file_location("fasta_io", _fasta_io_path)
    _fasta_io = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_fasta_io)
    read_fasta = _fasta_io.read_fasta
    write_fasta = _fasta_io.write_fasta
    iter_fasta = _fasta_io.iter_fasta


def load_fasta_directory(
    fasta_dir: Path,
    glob_pattern: str = "*_tests_*.fasta"
) -> Dict[str, str]:
    """
    Load all sequences from FASTA files in a directory.

    Commonly used patterns:
        - "*_tests_50kb.fasta" for Evo2/Yorzoi
        - "*_tests_16kb.fasta" for Shorkie

    Args:
        fasta_dir: Directory containing FASTA files
        glob_pattern: Glob pattern to match files

    Returns:
        {test_id: sequence} for all sequences across all files
    """
    fasta_dir = Path(fasta_dir)
    all_sequences = {}

    for fasta_file in sorted(fasta_dir.glob(glob_pattern)):
        sequences = read_fasta(fasta_file)
        all_sequences.update(sequences)

    return all_sequences


def load_metadata(
    metadata_file: Path,
    required_columns: Optional[List[str]] = None
) -> pd.DataFrame:
    """
    Load test metadata TSV file.

    Standard columns:
        - test_id: Unique identifier for each test
        - orf: Target ORF (e.g., YNL085W)
        - gene: Gene name (e.g., MKT1)
        - background_strain: Reference strain (e.g., BYa, RMx)
        - test_type: Type of test (wt, single_variant, all_variants)
        - n_variants: Number of variants in this test

    Args:
        metadata_file: Path to metadata TSV
        required_columns: Columns that must be present (raises ValueError if missing)

    Returns:
        DataFrame with test metadata
    """
    metadata_file = Path(metadata_file)

    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    df = pd.read_csv(metadata_file, sep='\t')

    if required_columns:
        missing = set(required_columns) - set(df.columns)
        if missing:
            raise ValueError(
                f"Metadata missing required columns: {missing}\n"
                f"Available columns: {list(df.columns)}"
            )

    return df


def write_results(
    results: pd.DataFrame,
    output_file: Path,
    model_name: str,
    include_summary: bool = True
) -> None:
    """
    Write VEP results to TSV file with optional summary.

    Args:
        results: DataFrame with prediction results
        output_file: Output TSV path
        model_name: Model name for logging
        include_summary: Print summary statistics
    """
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    results.to_csv(output_file, sep='\t', index=False)

    if include_summary:
        print(f"\n{'=' * 60}")
        print(f"{model_name} Results Summary")
        print(f"{'=' * 60}")
        print(f"Output: {output_file}")
        print(f"Total predictions: {len(results)}")

        if 'prediction' in results.columns:
            print(f"\nPrediction distribution:")
            for pred, count in results['prediction'].value_counts().items():
                print(f"  {pred}: {count} ({100*count/len(results):.1f}%)")

        if 'score' in results.columns:
            print(f"\nScore statistics:")
            print(f"  Mean: {results['score'].mean():.4f}")
            print(f"  Std:  {results['score'].std():.4f}")
            print(f"  Min:  {results['score'].min():.4f}")
            print(f"  Max:  {results['score'].max():.4f}")


def generate_slurm_header(
    job_name: str,
    partition: str = "gpu",
    gres: str = "gpu:1",
    constraint: Optional[str] = None,
    mem_gb: int = 64,
    cpus: int = 8,
    time_hours: int = 24,
    log_dir: Optional[Path] = None
) -> str:
    """
    Generate SLURM script header.

    Args:
        job_name: SLURM job name
        partition: SLURM partition
        gres: GPU resource request
        constraint: GPU constraint (e.g., "GPU_SKU:H100_SXM5&GPU_MEM:80GB")
        mem_gb: Memory in GB
        cpus: Number of CPUs
        time_hours: Wall time in hours
        log_dir: Directory for log files

    Returns:
        SLURM header as string
    """
    log_dir = log_dir or Path("logs")

    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --partition={partition}",
        f"#SBATCH --gres={gres}",
    ]

    if constraint:
        lines.append(f'#SBATCH --constraint="{constraint}"')

    lines.extend([
        f"#SBATCH --mem={mem_gb}G",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --time={time_hours}:00:00",
        f"#SBATCH --output={log_dir}/{job_name}_%j.log",
        f"#SBATCH --error={log_dir}/{job_name}_%j.err",
    ])

    return "\n".join(lines)


def generate_conda_init(conda_env: str) -> str:
    """
    Generate proper conda initialization for SLURM scripts.

    Uses eval + conda shell hook instead of source ~/.bashrc
    to ensure conda PATH is properly set.

    Args:
        conda_env: Name of conda environment to activate

    Returns:
        Shell commands for conda initialization
    """
    return f'''# Environment setup - proper conda initialization for batch jobs
# Set CONDA_EXE or edit this path for your HPC installation
eval "$(conda shell.bash hook)"
conda activate {conda_env}

# Verify environment
echo "Python: $(which python)"
echo "Python version: $(python --version)"
echo "Conda env: $CONDA_DEFAULT_ENV"
'''
