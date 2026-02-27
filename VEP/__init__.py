"""
Variant Effect Prediction (VEP) Pipeline
=========================================

A modular framework for predicting variant effects using multiple models:
- evo2: DNA language model (genomic fitness)
- esm1v: Protein language model (missense effects)
- shorkie: RNA expression predictor (8-fold ensemble)
- yorzoi: RNA expression predictor (Borzoi-like)

Modules:
    config: VEPConfig dataclass for project configuration
    core: Shared utilities (FASTA I/O, metadata, results, SLURM)
    locus_preparation: Build genomic loci and insert variants
    sequence_utils: Common sequence operations
    models: VEP model runners (abstract base + implementations)
    snakemake: Include-able Snakemake rules and step scripts

Usage:
    # Snakemake
    snakemake --configfile project_config.yaml

    # Programmatic
    from variant_effect_prediction import VEPConfig
    from variant_effect_prediction.models.yorzoi import YorzoiModel

Author: Kevin R. Roy
"""

from . import locus_preparation
from . import sequence_utils
from . import models
from . import core
from .config import VEPConfig, ModelConfig

__all__ = [
    "VEPConfig",
    "ModelConfig",
    "locus_preparation",
    "sequence_utils",
    "models",
    "core",
]

__version__ = "1.0.0"
