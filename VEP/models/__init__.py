"""
VEP Models Module

Model runners for variant effect prediction:
- evo2: DNA language model (7B/40B, requires H100 GPU)
- esm1v: Protein language model (650M params)
- shorkie: RNA expression predictor (TensorFlow, 8-fold ensemble)
- yorzoi: RNA expression predictor (HuggingFace transformer)

Each model submodule provides:
- runner.py: Main prediction runner
- config.py: Model-specific configuration

Author: Kevin R. Roy
"""

from . import evo2
from . import esm1v
from . import shorkie
from . import yorzoi
