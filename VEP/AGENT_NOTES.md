# Variant Effect Prediction Pipeline - Agent Notes

**Date:** 2026-02-10
**Author:** Kevin R. Roy
**Purpose:** Guide for agents working on VEP pipeline modules

---

## Overview

This is a modular pipeline for predicting variant effects using multiple models. Each module can be developed independently.

### Module Structure

```
variant_effect_prediction/
├── __init__.py
├── AGENT_NOTES.md              # This file
│
├── locus_preparation/          # Module 1: Build loci & insert variants
│   ├── __init__.py
│   └── AGENT_NOTES.md
│
├── sequence_utils/             # Module 2: Shared sequence operations
│   ├── __init__.py
│   └── AGENT_NOTES.md
│
└── models/                     # Module 3: VEP model runners
    ├── __init__.py
    ├── evo2/
    │   ├── __init__.py
    │   └── AGENT_NOTES.md
    ├── esm1v/
    │   ├── __init__.py
    │   └── AGENT_NOTES.md
    ├── shorkie/
    │   ├── __init__.py
    │   └── AGENT_NOTES.md
    └── yorzoi/
        ├── __init__.py
        └── AGENT_NOTES.md
```

---

## Source Scripts to Extract From

All existing scripts are in:
```
/path/to/projects/QTL/variant_effect_prediction/scripts/20260205_in_silico_QTN_dissection/
```

Key scripts:
| Script | Functions to Extract |
|--------|---------------------|
| `12_build_all_bloom_loci_s288c.py` | locus_builder, coordinate_mapping |
| `14_generate_variant_dissection_tests.py` | variant_insertion |
| `18_generate_trio_dissection_tests.py` | trio-based variant insertion |
| `19_run_trio_dissection_models.py` | Model runners (evo2, shorkie, yorzoi) |
| `20_run_trio_esm1v.py` | ESM1-v runner, translation |

---

## Data Locations

### Reference Data
```
/path/to/common/
├── reference_genomes/
│   └── MAGESTIC_background_strain.fasta
├── annotation_files/
│   └── MAGESTIC_background_strain_annotations.gff
└── published_data/
    └── bloom_et_al_2019/
        └── Bloom_QTL_gene_level.tsv
```

### Processed Data
```
/path/to/projects/QTL/variant_effect_prediction/processed_data/20260205_in_silico_QTN_dissection/
├── bloom_all_loci_s288c/           # Built loci
├── trio_variant_dissection/pilot/  # Trio test sequences
└── evo2_variant_dissection/        # evo2 batch files
```

---

## Model Requirements Summary

| Model | GPU Required | Conda Env | Key Dependency |
|-------|--------------|-----------|----------------|
| evo2 | H100 only (FP8) | N/A (Singularity) | vortex.v3.sif |
| esm1v | Any 16GB+ | `esm` | fair-esm |
| shorkie | Any GPU | `shorkie_env` | TensorFlow |
| yorzoi | A100/H100 | `yorzoi_local` | flash-attn |

---

## Agent Assignment

Each agent should focus on ONE module:
1. **Agent A**: `locus_preparation/` - Extract locus building functions
2. **Agent B**: `sequence_utils/` - Extract common utilities
3. **Agent C**: `models/evo2/` - Fix Singularity issues, create runner
4. **Agent D**: `models/esm1v/` - Create protein prediction runner
5. **Agent E**: `models/shorkie/` - Create TensorFlow runner
6. **Agent F**: `models/yorzoi/` - Create HuggingFace runner

---

## Testing Pattern

Each module should have a simple test:
```python
# At bottom of each module file
if __name__ == "__main__":
    # Quick test with minimal data
    print("Testing module...")
    # Test code here
```

---

## Import Pattern for Project Scripts

```python
import sys
from pathlib import Path
BASE_DIR = Path("/path/to")
sys.path.insert(0, str(BASE_DIR / "common/scripts"))

from variant_effect_prediction.locus_preparation import build_locus
from variant_effect_prediction.sequence_utils import translate, one_hot_encode
from variant_effect_prediction.models.evo2 import run_evo2_predictions
```

---

*Last updated: 2026-02-10*
