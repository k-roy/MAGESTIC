"""
evo2 Model Module

DNA language model for genomic fitness prediction.
- Uses Singularity container (vortex.v3.sif)
- Works on H100 (FP8) and A100/A40 (non-FP8)
- Models: 7B (14GB checkpoint), 40B (82GB checkpoint)

Author: Kevin R. Roy
Last Updated: 2026-02-10

Usage:
    from variant_effect_prediction.models.evo2 import (
        generate_slurm_script,
        check_gpu_compatibility,
        EVO2_DIR,
    )

    # Generate SLURM script for H100
    script = generate_slurm_script(
        input_file="/path/to/sequences.txt",
        output_file="/path/to/scores.tsv",
        gpu_type="H100",
    )

    # Or for A100 (uses non-FP8 config)
    script = generate_slurm_script(
        input_file="/path/to/sequences.txt",
        output_file="/path/to/scores.tsv",
        gpu_type="A100",
    )
"""

from .config import (
    EVO2_DIR,
    CONTAINER,
    CHECKPOINTS_DIR,
    MODEL_CONFIGS,
    GPU_MATRIX,
    SLURM_CONSTRAINTS,
    get_config_path,
    get_checkpoint_path,
    check_gpu_compatibility,
)

from .runner import (
    run_evo2_predictions,
    generate_slurm_script,
    generate_array_slurm_script,
)

__all__ = [
    # Config
    "EVO2_DIR",
    "CONTAINER",
    "CHECKPOINTS_DIR",
    "MODEL_CONFIGS",
    "GPU_MATRIX",
    "SLURM_CONSTRAINTS",
    "get_config_path",
    "get_checkpoint_path",
    "check_gpu_compatibility",
    # Runner
    "run_evo2_predictions",
    "generate_slurm_script",
    "generate_array_slurm_script",
]
