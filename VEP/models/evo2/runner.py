"""
evo2 Runner Module

Run evo2 DNA language model predictions via Singularity container.

Author: Kevin R. Roy
Last Updated: 2026-02-10

CRITICAL: Uses --cleanenv flag to prevent conda PATH pollution
that causes Triton compiler failures.

Working configurations validated 2026-02-10:
- H100 with FP8: scores -0.97111857, -0.97105616
- A40 without FP8: scores -0.9703312, -0.9702675
"""

import subprocess
from pathlib import Path
from typing import Optional, List
import pandas as pd

from .config import (
    EVO2_DIR,
    CONTAINER,
    MODEL_CONFIGS,
    SLURM_CONSTRAINTS,
    DEFAULT_INFERENCE_PARAMS,
    get_config_path,
    get_checkpoint_path,
    check_gpu_compatibility,
)


def run_evo2_predictions(
    input_file: Path,
    output_file: Path,
    model_size: str = "7b",
    use_fp8: bool = True,
    batch_size: int = 1,
    debug: bool = True,
    extra_bind_paths: Optional[List[Path]] = None,
) -> pd.DataFrame:
    """
    Run evo2 predictions using Singularity container.

    IMPORTANT: This runs synchronously and requires GPU access.
    For cluster jobs, use generate_slurm_script() instead.

    Args:
        input_file: Plain text file with one DNA sequence per line
        output_file: Output TSV file path
        model_size: "7b" or "40b"
        use_fp8: True for H100 (CC 9.0+), False for A100/A40
        batch_size: Batch size (default 1 for long sequences)
        debug: Enable verbose logging
        extra_bind_paths: Additional paths to bind into container

    Returns:
        DataFrame with columns: sequence_index, score

    Raises:
        RuntimeError: If evo2 inference fails
        FileNotFoundError: If input file or container doesn't exist
    """
    # Validate inputs
    input_file = Path(input_file)
    output_file = Path(output_file)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if not CONTAINER.exists():
        raise FileNotFoundError(f"Singularity container not found: {CONTAINER}")

    # Ensure output directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Get config and checkpoint paths
    config_path = get_config_path(model_size, use_fp8)
    checkpoint_path = get_checkpoint_path(model_size)

    # Build bind paths
    bind_paths = [
        f"{EVO2_DIR}:{EVO2_DIR}",
        f"{input_file.parent}:{input_file.parent}",
        f"{output_file.parent}:{output_file.parent}",
        f"/tmp:{EVO2_DIR}/.triton",  # Triton cache
    ]

    if extra_bind_paths:
        for path in extra_bind_paths:
            bind_paths.append(f"{path}:{path}")

    # Build singularity command
    # CRITICAL: --cleanenv prevents conda PATH pollution that causes
    # Triton compiler to find wrong gcc
    cmd = [
        "singularity", "exec",
        "--cleanenv",  # CRITICAL: Prevents conda PATH pollution
        "--nv",        # Enable NVIDIA GPU support
    ]

    # Add bind paths
    for bind_path in bind_paths:
        cmd.extend(["--bind", bind_path])

    # Set working directory and container
    cmd.extend([
        "--pwd", str(EVO2_DIR),
        str(CONTAINER),
    ])

    # Add inference command
    cmd.extend([
        "python3", "inference.py",
        "--config_path", config_path,
        "--checkpoint_path", checkpoint_path,
        "--input_file", str(input_file),
        "--output_file", str(output_file),
        "--batch_size", str(batch_size),
    ])

    if debug:
        cmd.append("--debug")

    # Run inference
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=str(EVO2_DIR),
    )

    if result.returncode != 0:
        error_msg = f"evo2 inference failed with exit code {result.returncode}\n"
        error_msg += f"STDERR:\n{result.stderr}\n"
        error_msg += f"STDOUT:\n{result.stdout}"
        raise RuntimeError(error_msg)

    # Read and return results
    if not output_file.exists():
        raise RuntimeError(f"Output file not created: {output_file}")

    return pd.read_csv(output_file, sep='\t')


def generate_slurm_script(
    input_file: Path,
    output_file: Path,
    job_name: str = "evo2_pred",
    model_size: str = "7b",
    gpu_type: str = "H100",
    time_hours: int = 8,
    partition: str = "owners",
    extra_bind_paths: Optional[List[Path]] = None,
) -> str:
    """
    Generate SLURM script for evo2 job.

    Args:
        input_file: Input sequences file
        output_file: Output scores file
        job_name: SLURM job name
        model_size: "7b" or "40b"
        gpu_type: "H100", "A100", or "A40"
        time_hours: Wall time in hours
        partition: SLURM partition
        extra_bind_paths: Additional paths to bind

    Returns:
        Complete SLURM script as string
    """
    input_file = Path(input_file)
    output_file = Path(output_file)

    # Determine config based on GPU type
    use_fp8 = gpu_type in ["H100", "H200", "RTX_6000_Ada", "L40S"]
    config_path = get_config_path(model_size, use_fp8)
    checkpoint_path = get_checkpoint_path(model_size)

    # Get SLURM constraint
    constraint = SLURM_CONSTRAINTS.get(gpu_type, f"GPU_SKU:{gpu_type}")

    # Memory requirement
    mem_gb = MODEL_CONFIGS[model_size]["min_gpu_mem_gb"] + 32  # Extra buffer

    # Build bind paths
    bind_paths = [
        '"$PWD:$PWD"',
        f'"{input_file.parent}:{input_file.parent}"',
        f'"{output_file.parent}:{output_file.parent}"',
        '"/tmp:$PWD/.triton"',
    ]
    if extra_bind_paths:
        for path in extra_bind_paths:
            bind_paths.append(f'"{path}:{path}"')

    bind_args = " \\\n        --bind ".join(bind_paths)

    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --constraint="{constraint}"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem={mem_gb}G
#SBATCH --time={time_hours}:00:00
#SBATCH --output=logs/{job_name}_%j.out
#SBATCH --error=logs/{job_name}_%j.err

# evo2 Inference Script
# Generated by variant_effect_prediction.models.evo2.runner
# GPU Type: {gpu_type} ({"FP8" if use_fp8 else "non-FP8"} mode)

set -e
mkdir -p logs

echo "=========================================="
echo "evo2 {model_size.upper()} Inference"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "GPU Type: {gpu_type}"
echo "Config: {config_path} ({"FP8" if use_fp8 else "non-FP8"})"
date
echo ""

# Show GPU info
nvidia-smi -L
nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader
echo ""

cd {EVO2_DIR}

echo "Input: {input_file}"
echo "Output: {output_file}"
echo ""

# CRITICAL: --cleanenv prevents conda PATH pollution
# This fixes Triton compiler failures where it finds the wrong gcc
singularity exec \\
    --cleanenv \\
    --nv \\
    --bind {bind_args} \\
    --pwd "$PWD" \\
    {CONTAINER} \\
    python3 inference.py \\
        --config_path {config_path} \\
        --checkpoint_path {checkpoint_path} \\
        --input_file {input_file} \\
        --output_file {output_file} \\
        --debug

echo ""
if [ -f "{output_file}" ]; then
    echo "SUCCESS!"
    echo "Output lines: $(wc -l < {output_file})"
    head -5 {output_file}
else
    echo "ERROR: Output file not created"
    exit 1
fi
echo "=========================================="
echo "Completed at: $(date)"
'''
    return script


def generate_array_slurm_script(
    input_dir: Path,
    output_dir: Path,
    batch_pattern: str = "batch_*.seq",
    job_name: str = "evo2_array",
    model_size: str = "7b",
    gpu_type: str = "H100",
    time_hours: int = 2,
    partition: str = "owners",
    array_range: str = "1-100",
) -> str:
    """
    Generate SLURM array job script for processing multiple batch files.

    Args:
        input_dir: Directory containing batch files
        output_dir: Directory for output files
        batch_pattern: Glob pattern for batch files
        job_name: SLURM job name
        model_size: "7b" or "40b"
        gpu_type: "H100", "A100", or "A40"
        time_hours: Wall time per task
        partition: SLURM partition
        array_range: SLURM array range (e.g., "1-100")

    Returns:
        Complete SLURM array script as string
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    # Determine config based on GPU type
    use_fp8 = gpu_type in ["H100", "H200", "RTX_6000_Ada", "L40S"]
    config_path = get_config_path(model_size, use_fp8)
    checkpoint_path = get_checkpoint_path(model_size)

    # Get SLURM constraint
    constraint = SLURM_CONSTRAINTS.get(gpu_type, f"GPU_SKU:{gpu_type}")
    mem_gb = MODEL_CONFIGS[model_size]["min_gpu_mem_gb"] + 32

    script = f'''#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --constraint="{constraint}"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem={mem_gb}G
#SBATCH --time={time_hours}:00:00
#SBATCH --array={array_range}
#SBATCH --output=logs/{job_name}_%A_%a.out
#SBATCH --error=logs/{job_name}_%A_%a.err

# evo2 Array Job Script
# Generated by variant_effect_prediction.models.evo2.runner

set -e
mkdir -p logs

INPUT_DIR="{input_dir}"
OUTPUT_DIR="{output_dir}"
mkdir -p "$OUTPUT_DIR"

# Get input file for this task
INPUT_FILE="$INPUT_DIR/{batch_pattern.replace('*', '${SLURM_ARRAY_TASK_ID}')}"
OUTPUT_FILE="$OUTPUT_DIR/scores_${{SLURM_ARRAY_TASK_ID}}.tsv"

echo "=========================================="
echo "evo2 Array Task $SLURM_ARRAY_TASK_ID"
echo "=========================================="
echo "Input: $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
date

if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

cd {EVO2_DIR}

singularity exec \\
    --cleanenv \\
    --nv \\
    --bind "$PWD:$PWD" \\
    --bind "$INPUT_DIR:$INPUT_DIR" \\
    --bind "$OUTPUT_DIR:$OUTPUT_DIR" \\
    --bind "/tmp:$PWD/.triton" \\
    --pwd "$PWD" \\
    {CONTAINER} \\
    python3 inference.py \\
        --config_path {config_path} \\
        --checkpoint_path {checkpoint_path} \\
        --input_file "$INPUT_FILE" \\
        --output_file "$OUTPUT_FILE" \\
        --debug

echo "Task $SLURM_ARRAY_TASK_ID completed at: $(date)"
'''
    return script


if __name__ == "__main__":
    from pathlib import Path

    print("evo2 Runner Module Test")
    print("=" * 50)

    # Create test paths
    test_input = Path("/tmp/evo2_test_input.txt")
    test_output = Path("/tmp/evo2_test_output.tsv")

    # Create a small test input
    with open(test_input, 'w') as f:
        f.write("ATGCATGCATGC" * 100 + "\n")  # 1200 bp

    print(f"Test input: {test_input}")
    print(f"Test output: {test_output}")
    print()

    # Generate SLURM script for H100
    print("Generated SLURM script for H100:")
    print("-" * 40)
    script = generate_slurm_script(
        test_input, test_output,
        job_name="evo2_test",
        gpu_type="H100",
        time_hours=1,
    )
    print(script[:1000] + "...")
    print()

    # Generate SLURM script for A100 (non-FP8)
    print("Generated SLURM script for A100 (non-FP8):")
    print("-" * 40)
    script_a100 = generate_slurm_script(
        test_input, test_output,
        job_name="evo2_test_a100",
        gpu_type="A100",
        time_hours=1,
    )
    print(script_a100[:1000] + "...")
    print()

    print("NOTE: To run, save script to .slurm file and submit with sbatch")
    print("H100 needed for FP8, A100/A40 work with non-FP8 config")
