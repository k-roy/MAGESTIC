"""
evo2 Configuration

Paths and settings for running evo2 DNA language model.

Author: Kevin R. Roy
Last Updated: 2026-02-10

Working Configurations (validated 2026-02-10):
- H100 with FP8: scores -0.97111857, -0.97105616
- A40 without FP8: scores -0.9703312, -0.9702675
- Scores are within ~0.001, confirming non-FP8 mode works well
"""

from pathlib import Path
from typing import Dict, Any

# Base paths
# Set VEP_BASE_DIR environment variable to override
import os
BASE_DIR = Path(os.environ.get("VEP_BASE_DIR", Path(__file__).parent))
EVO2_DIR = BASE_DIR / "common/software/evo2"
CONTAINER = EVO2_DIR / "vortex.v3.sif"
CHECKPOINTS_DIR = EVO2_DIR / "checkpoints"
CONFIGS_DIR = EVO2_DIR / "configs"

# Model configurations
# Critical: inner_mlp_size must be 11264 to match evo2_7b.pt checkpoint
MODEL_CONFIGS: Dict[str, Dict[str, Any]] = {
    "7b": {
        "checkpoint": "checkpoints/evo2_7b.pt",
        "config_fp8": "configs/evo2-7b-1m.yml",        # For H100 (CC 9.0+)
        "config_no_fp8": "configs/evo2-7b-1m-no-fp8.yml",  # For A100/A40 (CC 8.0-8.6)
        "min_gpu_mem_gb": 32,  # Safe operation needs 32GB+
        "checkpoint_size_gb": 14,
    },
    "40b": {
        "checkpoint": "checkpoints/evo2_40b.pt",
        "config_fp8": "configs/evo2-40b-1m.yml",
        "config_no_fp8": None,  # Not available - 40B requires H100+
        "min_gpu_mem_gb": 80,
        "checkpoint_size_gb": 82,
    },
}

# GPU Compatibility Matrix (Compute Capability determines FP8 support)
# FP8 requires CC >= 8.9 (Ada Lovelace or Hopper architecture)
GPU_MATRIX = {
    "H200": {"cc": 9.0, "vram_gb": 141, "fp8": True, "7b": True, "40b": True},
    "H100_SXM5": {"cc": 9.0, "vram_gb": 80, "fp8": True, "7b": True, "40b": True},
    "H100_PCIe": {"cc": 9.0, "vram_gb": 80, "fp8": True, "7b": True, "40b": True},
    "RTX_6000_Ada": {"cc": 8.9, "vram_gb": 48, "fp8": True, "7b": True, "40b": False},
    "L40S": {"cc": 8.9, "vram_gb": 48, "fp8": True, "7b": True, "40b": False},
    "A100_SXM4": {"cc": 8.0, "vram_gb": 80, "fp8": False, "7b": True, "40b": False},
    "A100_PCIe": {"cc": 8.0, "vram_gb": 40, "fp8": False, "7b": True, "40b": False},
    "A40": {"cc": 8.6, "vram_gb": 48, "fp8": False, "7b": True, "40b": False},
    "A30": {"cc": 8.0, "vram_gb": 24, "fp8": False, "7b": False, "40b": False},
    "V100_32GB": {"cc": 7.0, "vram_gb": 32, "fp8": False, "7b": False, "40b": False},
    "V100_16GB": {"cc": 7.0, "vram_gb": 16, "fp8": False, "7b": False, "40b": False},
    "TITAN_Xp": {"cc": 6.1, "vram_gb": 12, "fp8": False, "7b": False, "40b": False},
}

# SLURM partition constraints for targeting specific GPUs
SLURM_CONSTRAINTS = {
    "H100": "GPU_SKU:H100_SXM5",
    "A100": "GPU_SKU:A100_SXM4",
    "A40": "GPU_SKU:A40",
    "80GB": "GPU_MEM:80GB",  # Gets H100 or A100-80GB
    "48GB": "GPU_MEM:48GB",  # Gets A40 or L40S
}

# Inference parameters
DEFAULT_INFERENCE_PARAMS = {
    "batch_size": 1,
    "prepend_bos": False,
    "reduced_method": "mean",  # 'mean' or 'sum' for per-position scores
    "average_reverse_complement": False,
    "debug": True,
}


def get_config_path(model_size: str = "7b", use_fp8: bool = True) -> str:
    """
    Get the appropriate config file path.

    Args:
        model_size: "7b" or "40b"
        use_fp8: True for H100 (CC 9.0+), False for A100/A40 (CC 8.0-8.6)

    Returns:
        Relative path to config file (relative to EVO2_DIR)
    """
    config = MODEL_CONFIGS[model_size]
    if use_fp8:
        return config["config_fp8"]
    else:
        if config["config_no_fp8"] is None:
            raise ValueError(f"Model {model_size} requires FP8 (H100 GPU)")
        return config["config_no_fp8"]


def get_checkpoint_path(model_size: str = "7b") -> str:
    """Get checkpoint file path (relative to EVO2_DIR)."""
    return MODEL_CONFIGS[model_size]["checkpoint"]


def check_gpu_compatibility(gpu_name: str, model_size: str = "7b") -> tuple[bool, str]:
    """
    Check if a GPU is compatible with the specified model.

    Args:
        gpu_name: GPU name (e.g., "H100_SXM5", "A100_SXM4")
        model_size: "7b" or "40b"

    Returns:
        Tuple of (is_compatible, reason)
    """
    if gpu_name not in GPU_MATRIX:
        return False, f"Unknown GPU: {gpu_name}"

    gpu = GPU_MATRIX[gpu_name]
    model = MODEL_CONFIGS[model_size]

    # Check VRAM
    if gpu["vram_gb"] < model["min_gpu_mem_gb"]:
        return False, f"Insufficient VRAM: {gpu['vram_gb']}GB < {model['min_gpu_mem_gb']}GB required"

    # Check model-specific compatibility
    if model_size == "7b" and not gpu["7b"]:
        return False, f"GPU {gpu_name} not compatible with 7B model"
    if model_size == "40b" and not gpu["40b"]:
        return False, f"GPU {gpu_name} not compatible with 40B model"

    # Determine if FP8 is available
    use_fp8 = gpu["fp8"]
    config_type = "FP8" if use_fp8 else "non-FP8"

    return True, f"Compatible using {config_type} config"


if __name__ == "__main__":
    print("evo2 Configuration Test")
    print("=" * 50)
    print(f"EVO2_DIR: {EVO2_DIR}")
    print(f"Container exists: {CONTAINER.exists()}")
    print(f"Checkpoints dir exists: {CHECKPOINTS_DIR.exists()}")
    print()

    print("Model Configurations:")
    for model, config in MODEL_CONFIGS.items():
        print(f"  {model}:")
        print(f"    Checkpoint: {config['checkpoint']}")
        print(f"    FP8 config: {config['config_fp8']}")
        print(f"    Non-FP8 config: {config['config_no_fp8']}")
        print(f"    Min GPU mem: {config['min_gpu_mem_gb']}GB")
    print()

    print("GPU Compatibility Check (7B model):")
    for gpu_name in ["H100_SXM5", "A100_SXM4", "A40", "V100_32GB"]:
        compat, reason = check_gpu_compatibility(gpu_name, "7b")
        status = "OK" if compat else "FAIL"
        print(f"  {gpu_name}: [{status}] {reason}")
