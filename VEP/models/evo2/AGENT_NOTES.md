# evo2 Model Module - Agent Notes

**Date:** 2026-02-10
**Purpose:** Run evo2 DNA language model predictions
**Status:** WORKING - 7B model validated on H100/A100/A40; 40B requires multi-GPU

---

## Quick Reference

| Model | Min VRAM | GPU Requirement | Single GPU Options | Config |
|-------|----------|-----------------|-------------------|--------|
| **7B** | 32GB | Ampere+ (CC ≥8.0) | H100, A100, A40, L40S | `evo2-7b-1m.yml` or `evo2-7b-1m-no-fp8.yml` |
| **40B** | >80GB | Ampere+ + Multi-GPU | ❌ Needs 2×H100 | `evo2-40b-1m.yml` (not tested) |

**⚠️ V100 (CC 7.0) does NOT work** - FlashAttention requires Ampere architecture or newer!

---

## Working Configurations (Validated 2026-02-10)

| Job ID | GPU | Config | Scores | Key Fix |
|--------|-----|--------|--------|---------|
| 15463182 | H100 80GB | evo2-7b-1m.yml | ✅ 719 scores | `--cleanenv` |
| 15462819 | A40 48GB | evo2-7b-1m-no-fp8.yml | ✅ 719 scores | `--cleanenv` + non-FP8 |
| 15484440 | A40 48GB | auto-detect | ✅ 719 scores | Auto GPU detection |
| 15483880 | H100 80GB | evo2-7b-1m.yml | ✅ 719 scores | Production pilot |

**Scores are within ~0.001**, confirming non-FP8 mode works well for inference.

---

## GPU Compatibility Matrix

| GPU | VRAM | Compute Cap | FP8 Support | 7B Model | 40B Model | Config |
|-----|------|-------------|-------------|----------|-----------|--------|
| **H100** | 80GB | 9.0 | ✅ Yes | ✅ Works | ❌ OOM | `evo2-7b-1m.yml` |
| **A100-80GB** | 80GB | 8.0 | ❌ No | ✅ Works | ❌ OOM | `evo2-7b-1m-no-fp8.yml` |
| **A100-40GB** | 40GB | 8.0 | ❌ No | ✅ Works | ❌ OOM | `evo2-7b-1m-no-fp8.yml` |
| **A40** | 48GB | 8.6 | ❌ No | ✅ Works | ❌ OOM | `evo2-7b-1m-no-fp8.yml` |
| **L40S** | 48GB | 8.9 | ✅ Yes | ✅ Works | ❌ OOM | `evo2-7b-1m.yml` |
| **V100-32GB** | 32GB | 7.0 | ❌ No | ❌ FlashAttn! | ❌ OOM | FlashAttention requires CC 8.0+ |
| **RTX 3090** | 24GB | 8.6 | ❌ No | ❌ OOM | ❌ OOM | - |
| **TITAN Xp** | 12GB | 6.1 | ❌ No | ❌ OOM | ❌ OOM | - |

**Key insights:**
- **FlashAttention requires Ampere+ (CC ≥ 8.0)** - V100 and older won't work!
- FP8 requires Compute Capability >= 8.9 (Ada Lovelace, Hopper, or newer)
- 7B model needs ~32GB VRAM minimum + Ampere+ GPU
- 40B model checkpoint is 82GB → needs multi-GPU setup
- Use `GPU_MEM:48GB` constraint to ensure L40S/H100 (avoids V100 assignment)

---

## Critical Fixes Discovered

### Fix 1: `--cleanenv` Flag (MOST IMPORTANT)

**Problem:** Triton compiler finds conda's gcc in PATH instead of container's gcc.
```
subprocess.CalledProcessError: Command '['/path/to/anaconda3/bin/x86_64-conda-linux-gnu-cc'...]' returned non-zero exit status 1.
```

**Solution:** Use `--cleanenv` flag with singularity:
```bash
singularity exec --cleanenv --nv ...
```

This prevents the host environment from polluting the container's PATH.

### Fix 2: FP8 Config for Older GPUs

**Problem:** FP8 assertion error on A40/A100:
```
AssertionError: Device compute capability 8.9 or higher required for FP8
```

**Solution:** Use non-FP8 config:
```bash
--config_path configs/evo2-7b-1m-no-fp8.yml
```

### Fix 3: Correct Model Config

**Problem:** Size mismatch error:
```
size mismatch for layers.0.in_linear.weight
```

**Solution:** Use the `1m` configs (not `8k`):
- `evo2-7b-1m.yml` has `inner_mlp_size: 11264` (correct)
- `evo2-7b-8k.yml` has wrong size (don't use)

### Fix 4: Triton Cache Binding

**Problem:** Triton can't write cache files.

**Solution:** Bind `/tmp` to `.triton`:
```bash
--bind "/tmp:$PWD/.triton"
```

---

## Module Files

```
common/scripts/variant_effect_prediction/models/evo2/
├── __init__.py       # Module exports
├── AGENT_NOTES.md    # This file
├── config.py         # Paths, GPU matrix, model configs
└── runner.py         # Singularity-based runner
```

---

## Auto-Detecting GPU Script (RECOMMENDED)

For maximum flexibility across partitions, use auto-detection instead of hardcoded GPU constraints:

```bash
#!/bin/bash
#SBATCH --job-name=evo2_7b
#SBATCH --partition=gpu  # or owners
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00

set -e
mkdir -p logs

# Auto-detect GPU compute capability
GPU_CC=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1)
CC_MAJOR=$(echo $GPU_CC | cut -d. -f1)

echo "GPU Compute Capability: $GPU_CC"

if [ "$CC_MAJOR" -ge 9 ]; then
    CONFIG="configs/evo2-7b-1m.yml"        # FP8 for H100 (CC 9.0+)
    echo "Using FP8 config"
else
    CONFIG="configs/evo2-7b-1m-no-fp8.yml" # Non-FP8 for A40/A100 (CC 8.x)
    echo "Using non-FP8 config"
fi

cd /path/to/common/software/evo2

singularity exec --cleanenv --nv \
    --bind "$PWD:$PWD" \
    --bind "$INPUT_DIR:$INPUT_DIR" \
    --bind "/tmp:$PWD/.triton" \
    --pwd "$PWD" \
    vortex.v3.sif \
    python3 inference.py \
        --config_path $CONFIG \
        --checkpoint_path checkpoints/evo2_7b.pt \
        --input_file "$INPUT_FILE" \
        --output_file "$OUTPUT_FILE"
```

**Benefits:**
- Avoids failed submissions when specific GPU types aren't available
- Works on both `gpu` and `owners` partitions
- Automatically picks optimal config based on actual GPU assigned

---

## Quick Usage

### Generate SLURM script programmatically

```python
import sys
from pathlib import Path

BASE_DIR = Path("/path/to")
sys.path.insert(0, str(BASE_DIR / "common/scripts"))

from variant_effect_prediction.models.evo2 import generate_slurm_script

# For H100 (FP8 mode)
script = generate_slurm_script(
    input_file="/path/to/sequences.txt",
    output_file="/path/to/scores.tsv",
    gpu_type="H100",
)

# For A100/A40 (non-FP8 mode)
script = generate_slurm_script(
    input_file="/path/to/sequences.txt",
    output_file="/path/to/scores.tsv",
    gpu_type="A100",  # Uses non-FP8 config automatically
)
```

### Check GPU compatibility

```python
from variant_effect_prediction.models.evo2.config import check_gpu_compatibility

compatible, reason = check_gpu_compatibility("A100_SXM4", "7b")
print(f"A100: {reason}")  # "Compatible using non-FP8 config"

compatible, reason = check_gpu_compatibility("H100_SXM5", "40b")
print(f"H100 40B: {reason}")  # "40B model requires >80GB VRAM"
```

---

## Input/Output Format

**Input:** Plain text file, one DNA sequence per line
```
ATGCATGCATGC...
GCTAGCTAGCTA...
```

**Output:** TSV file with scores
```
score
-0.97111857
-0.97105616
```

Score is mean log-likelihood per position (negative, higher = better).

---

## Known Issues & Solutions

### Issue 1: Triton Compiler Error
- **Symptom:** `subprocess.CalledProcessError` with conda compiler path
- **Solution:** Use `--cleanenv` flag (already in runner.py)

### Issue 2: FP8 Assertion Error
- **Symptom:** `AssertionError: Device compute capability 8.9 or higher required`
- **Solution:** Use `evo2-7b-1m-no-fp8.yml` config

### Issue 3: Size Mismatch Error
- **Symptom:** `size mismatch for layers.0.in_linear.weight`
- **Solution:** Use config with `inner_mlp_size: 11264` (1m configs have this)

### Issue 4: OOM Error on 7B
- **Symptom:** `torch.OutOfMemoryError: CUDA out of memory`
- **Solution:** Need GPU with >= 32GB VRAM for 7B model

### Issue 5: OOM Error on 40B
- **Symptom:** `torch.OutOfMemoryError` on H100 80GB
- **Cause:** 40B checkpoint is 82GB, exceeds single GPU memory
- **Solution:** Requires multi-GPU setup (not yet implemented)

### Issue 6: Git Submodule Error
- **Symptom:** `git submodule update` fails
- **Solution:** Use direct `singularity exec` command, not `run_singularity` wrapper

### Issue 7: SLURM Constraint Unavailable
- **Symptom:** `sbatch: error: Requested node configuration is not available`
- **Solution:** Use auto-detecting GPU script without hardcoded constraints

### Issue 8: FlashAttention Requires Ampere+ GPUs (CRITICAL!)
- **Symptom:** `RuntimeError: FlashAttention only supports Ampere GPUs or newer.`
- **Cause:** V100 (CC 7.0) and older GPUs don't support FlashAttention
- **Why it matters:** Even with `GPU_MEM:32GB` constraint, V100-32GB gets assigned but crashes at runtime
- **Solution:** Use `GPU_MEM:48GB` constraint to ensure L40S (CC 8.9) or H100 (CC 9.0)
- **Alternative:** Use `--constraint="GPU_GEN:AMPERE"` or `--constraint="GPU_GEN:ADA"` on clusters that support it
- **Example error (Job 15544813):**
  ```
  Node: cluster-node (Tesla V100-SXM2-32GB)
  GPU Compute Capability: 7.0
  RuntimeError: FlashAttention only supports Ampere GPUs or newer.
  ```
- **Key insight:** evo2 7B requires BOTH: (1) ≥32GB VRAM AND (2) Ampere+ architecture (CC ≥8.0)

---

## Container Details

**Location:** `/path/to/common/software/evo2/vortex.v3.sif`

**Checkpoints:**
- `checkpoints/evo2_7b.pt` (14GB) - Works on single GPU
- `checkpoints/evo2_40b.pt` (82GB) - Needs multi-GPU

**Configs:**
- `configs/evo2-7b-1m.yml` - FP8 mode (H100, L40S)
- `configs/evo2-7b-1m-no-fp8.yml` - Non-FP8 mode (A100, A40)
- `configs/evo2-40b-1m.yml` - 40B model (multi-GPU only)
- `configs/evo2-7b-8k.yml` - WRONG `inner_mlp_size`, don't use

---

## SLURM Constraints for Targeting GPUs

```bash
# H100 (FP8 mode)
#SBATCH --constraint="GPU_SKU:H100_SXM5"

# A100 (non-FP8 mode)
#SBATCH --constraint="GPU_SKU:A100_SXM4"

# A40 (non-FP8 mode) - often unavailable
#SBATCH --constraint="GPU_SKU:A40"

# Any 80GB GPU (H100 or A100-80GB)
#SBATCH --constraint="GPU_MEM:80GB"
```

**Note:** Hardcoded constraints often fail due to GPU unavailability. Use auto-detection instead.

---

## 40B Model Status

**Current Status:** 🔄 Multi-GPU setup ready for testing (2026-02-11)

**Single GPU:** ❌ OOM on H100 80GB (checkpoint is 77GB)

**Multi-GPU Solution:** The vortex code has **built-in pipeline parallelism**!

The model auto-detects `torch.cuda.device_count()` and distributes layers across GPUs:
```python
# From vortex/model/model.py lines 637-645
num_gpus = torch.cuda.device_count()
layers_per_gpu = math.ceil(config.num_layers / num_gpus)
# 40B has 50 layers → 25 layers per GPU with 2 GPUs
```

### Multi-GPU SLURM Configuration

```bash
#!/bin/bash
#SBATCH --job-name=evo2_40b
#SBATCH --partition=owners
#SBATCH --gres=gpu:2              # Request 2 GPUs
#SBATCH --constraint="GPU_MEM:80GB"  # Ensures H100 or A100-80GB
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G                 # Need extra RAM for checkpoint loading
#SBATCH --time=4:00:00

# The model will auto-detect 2 GPUs and distribute 50 layers (25 per GPU)
# Each GPU will use ~38.5GB for model weights + activations
```

### GPU Requirements for 40B

| Setup | GPUs | VRAM per GPU | Total | Status |
|-------|------|--------------|-------|--------|
| 2× H100 | 2 | 80GB | 160GB | ✅ Should work |
| 2× A100-80GB | 2 | 80GB | 160GB | ✅ Should work |
| 4× A40 | 4 | 48GB | 192GB | ⚠️ Untested |

### Config Files

- `evo2-40b-1m.yml` - FP8 mode for H100 (CC 9.0+)
- `evo2-40b-1m-no-fp8.yml` - Non-FP8 mode for A100 (CC 8.x)

### Test Scripts

Located in `evo2_batch_pairwise/`:
- `run_evo2_40b_quick_test.slurm` - 10-sequence validation test
- `run_evo2_40b_test.slurm` - Full batch test (500 sequences)

### Memory Considerations

- **Checkpoint:** 77GB (loaded to CPU first, then distributed to GPUs)
- **System RAM:** Request 256GB to safely load checkpoint
- **Per-GPU activation memory:** ~10-15GB per 15kb sequence

---

## Pilot Job Results (2026-02-10)

### Completed Jobs

| Model | Job ID | GPU | Time | Sequences | Status |
|-------|--------|-----|------|-----------|--------|
| evo2 7B | 15483880 | H100 | 13:32 | 719 | ✅ Complete |
| evo2 7B auto | 15484440 | A40 | 34:52 | 719 | ✅ Complete |
| Yorzoi | 15484327 | H100 | 4:17 | 719 | ✅ Complete |
| ESM1-v | 15486355 | - | - | - | 🔄 Running |
| Shorkie | 15486523 | - | - | - | ⏳ Pending |

### Pilot Data
- **Location:** `projects/QTL/20260210_variant_effect_prediction/processed_data/20260205_in_silico_QTN_dissection/trio_variant_dissection/pilot/`
- **Sequences:** 719 (8 genes × ~90 variant tests each)
- **Format:** One DNA sequence per line (15,000-15,209 bp)

---

## Performance Benchmarks

| GPU | Config | 719 Sequences | Per-Sequence |
|-----|--------|---------------|--------------|
| H100 | FP8 | ~13 min | ~1.1 sec |
| A40 | non-FP8 | ~35 min | ~2.9 sec |

---

## Debugging Checklist

When evo2 fails, check in this order:

1. **GPU Memory:** Is there enough VRAM? (32GB+ for 7B)
2. **Compute Capability:** Is it >= 8.9 for FP8? If not, use non-FP8 config
3. **`--cleanenv` flag:** Is it present? (fixes Triton/gcc issues)
4. **Config file:** Using `1m` variant? (not `8k`)
5. **Triton cache:** Is `/tmp` bound to `.triton`?
6. **File permissions:** Can output directory be written?

---

## References

- Full debugging history: `/path/to/common/software/evo2/EVO2_PARAMETER_ANALYSIS.md`
- ArcInstitute/evo2: https://github.com/ArcInstitute/evo2
- NVIDIA NIM Prerequisites: https://docs.nvidia.com/nim/bionemo/evo2/latest/prerequisites.html

---

*Last updated: 2026-02-11 - Added FlashAttention CC 8.0+ requirement (Issue 8)*
