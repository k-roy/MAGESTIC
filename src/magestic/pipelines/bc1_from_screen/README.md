# BC1-From-Screen Pipeline

A portable pipeline for processing BC1 barcode counts from pooled CRISPR screens with library-specific tiered analysis.

## Overview

This pipeline processes BC1 barcode sequencing data from MAGESTIC screens to identify significant phenotype-associated variants. It supports:

- **Library-specific processing**: Different libraries (e.g., yL437, yL442) processed separately
- **Tiered BC1 analysis**: BC1s stratified by t0 read counts into tiers
- **Multi-timepoint screens**: Time series analysis with fitness slopes
- **Sherlock best practices**: Scratch-based I/O for SLURM jobs

## Features

- **Tier-aware hit calling**: Different thresholds for tier_1 vs tier_2
- **T0 consolidation**: Combines technical replicates for robust normalization
- **Adaptive thresholds**: Data-driven baseMean and concordance thresholds
- **Library-specific DESeq2**: Proper "median of ratios" normalization per library

## Installation

The pipeline is located at:
```
/path/to/common/scripts/bc1_from_screen_pipeline/
```

Add to your Python path:
```python
import sys
sys.path.insert(0, '/path/to/common/scripts')

from bc1_from_screen_pipeline.config import BC1FromScreenConfig, create_spg_qtl_screen_config
from bc1_from_screen_pipeline.pipeline import run_pipeline
```

## Usage

### Command Line

```bash
# Basic usage
python common/scripts/bc1_from_screen_pipeline/pipeline/run_pipeline.py \
    --project-dir /path/to/project

# With scratch for SLURM jobs (recommended)
python common/scripts/bc1_from_screen_pipeline/pipeline/run_pipeline.py \
    --project-dir /path/to/project \
    --use-scratch

# Run specific steps
python common/scripts/bc1_from_screen_pipeline/pipeline/run_pipeline.py \
    --project-dir /path/to/project \
    --steps tier,matrix,deseq2

# Use pre-defined SpG QTL screen config
python common/scripts/bc1_from_screen_pipeline/pipeline/run_pipeline.py \
    --use-spg-config
```

### Python API

```python
from pathlib import Path
from bc1_from_screen_pipeline.config import (
    BC1FromScreenConfig,
    LibraryConfig,
    TierConfig,
    create_spg_qtl_screen_config,
)
from bc1_from_screen_pipeline.pipeline import run_pipeline

# Use pre-defined config
config = create_spg_qtl_screen_config()

# Or create custom config
config = BC1FromScreenConfig(
    project_name="QTL",
    screen_name="my_screen",
    project_dir=Path("/path/to/project"),
    libraries=[
        LibraryConfig(
            name="myLib",
            screen_dates=["20251001", "20251015"],
            bottlenecked=False,
            timepoints=[0, 21],
            t0_consolidation_groups=4,
        ),
    ],
    use_scratch=True,  # For SLURM jobs
)

# Run pipeline
success = run_pipeline(config)
```

### SLURM Job Script

Use the template at `common/scripts/templates/sherlock_sbatch_template.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=bc1_screen
#SBATCH --partition=mylab
#SBATCH --account=mylab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00

source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate base

python /path/to/common/scripts/bc1_from_screen_pipeline/pipeline/run_pipeline.py \
    --use-spg-config \
    --use-scratch
```

## Module Structure

```
bc1_from_screen_pipeline/
├── __init__.py
├── config.py              # Configuration management + Sherlock support
├── README.md              # This file
├── core/
│   ├── __init__.py
│   ├── tier_classification.py    # BC1 tier assignment by t0 counts
│   ├── count_matrix.py           # DESeq2 matrix creation, t0 consolidation
│   ├── hit_calling.py            # Tier-aware hit calling
│   ├── variance_calibration.py   # Heteroscedasticity correction (synonymous variants)
│   ├── library_calibration.py    # Per-library bias calculation and correction
│   ├── empirical_thresholds.py   # Data-driven threshold calibration per condition
│   ├── noise_assessment.py       # Noisy condition/loci identification, penalties
│   ├── slope_fitness.py          # Multi-timepoint fitness via linear regression
│   ├── time_series.py            # Time series analysis framework
│   ├── causality.py              # Causality evaluation with evidence flags
│   ├── annotation.py             # Variant annotation using harmonized oligo table
│   ├── sublibrary.py             # Yeast sublibrary classification, nuclease metadata
│   ├── aggregation.py            # Hit aggregation by category, oligo, amino acid
│   └── ihw.py                    # Independent hypothesis weighting (IHW)
├── snakemake/             # Pipeline step scripts (01-18)
├── pipeline/
│   ├── __init__.py
│   └── run_pipeline.py    # Main pipeline runner
└── utils/
    └── __init__.py
```

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 01 | `validate_collapsed.py` | Validate collapsed FASTQ files |
| 02 | `parse_bc1_reads.py` | Parse BC1 barcodes from reads |
| 03 | `sample_qc.py` | Sample QC filtering |
| 04 | `combine_counts.py` | Combine per-sample counts |
| 05 | `error_detection.py` | Detect and flag sequencing errors |
| 06 | `absorb_error_counts.py` | Absorb error counts into parent BC1s |
| 07 | `tier_classification.py` | Classify BC1s into tiers by t0 counts |
| 08 | `create_count_matrix.py` | Create wide-format count matrix with t0 consolidation |
| 09 | `combine_technical_replicates.py` | Combine technical replicates |
| 10 | `split_comparisons.py` | Split count tables by comparison |
| 11 | `deseq2_single.py` | Run DESeq2 per library per comparison |
| 12 | `slope_fitness.py` | Multi-timepoint fitness estimation (linear regression) |
| 13 | `annotate_results.py` | Merge DESeq2 results with variant annotations |
| 14 | `hit_calling.py` | Tier-aware hit calling per comparison |
| 15 | `aggregate_hits.py` | Aggregate hits across comparisons |
| 16 | `rarefaction.py` | Rarefaction analysis for saturation curves |
| 17 | `deseq_rarefaction.py` | DESeq2 on rarefied data |
| 18 | `validate_sublibrary_composition.py` | Validate sublibrary composition |

## Tier Classification

BC1s are classified based on t0 read counts **within each library**:

| Tier | Read Threshold | Treatment |
|------|---------------|-----------|
| **tier_1** (abundance_high) | ≥100 reads | Individual BC1 DESeq2 |
| **tier_2** (abundance_medium) | 20-99 reads | Aggregated by oligo (pseudo-BC1s) |
| **tier_3** (abundance_low) | <20 reads | Flagged as low-confidence only |

**Note:** The code internally uses `abundance_high`, `abundance_medium`, `abundance_low` as tier identifiers. The tier_1/tier_2/tier_3 naming is used for configuration and output compatibility.

### Tiered BC1 Strategy

**Why tiers matter:**
- High-abundance BC1s (tier_1) have enough reads for individual statistical analysis
- Medium-abundance BC1s (tier_2) lack individual power but can be pooled by oligo
- Low-abundance BC1s (tier_3) are too sparse for reliable analysis

**Attribution for Tier 2:**
- Only BC1s with `closest_middle_donor_edit_distance == 0` are attributable
- These can be aggregated by `closest_oligo_name` (creates "pseudo-BC1s")
- Pseudo-BC1s are then analyzed alongside tier_1 BC1s in DESeq2

## Multi-Timepoint Analysis

For screens with multiple timepoints (e.g., 0, 6, 12, 21 generations):

### T0 Consolidation

The `consolidate_t0_samples()` function combines t0 technical replicates:
- Groups t0 samples into N consolidated replicates (default: 4)
- Uses summation to combine counts
- Tracks metadata: `n_samples_consolidated`, `source_samples`

### Slope Fitness (Step 12)

For multi-timepoint data, fitness is estimated via linear regression:

```python
from bc1_from_screen_pipeline.core.slope_fitness import calculate_slope_fitness

result = calculate_slope_fitness(
    timepoints=[0, 6, 12, 21],
    log2fc_values=[0, -0.5, -1.2, -2.1],
)
# Returns: slope, slope_se, r_squared, z_score, p_value, linearity_score, monotonicity
```

### Time Series Configuration

```yaml
libraries:
  - name: yL437
    timepoints: [0, 6, 12, 21]  # Multi-timepoint screen
    t0_consolidation_groups: 4
```

## Hit Calling Thresholds

| Parameter | Tier 1 | Tier 2 |
|-----------|--------|--------|
| Z-score threshold | 3.0 | 3.75 |
| Concordance threshold | 0.75 | 0.85 |
| Min BC1s aggregated | N/A | 2 |

### Adaptive Thresholds

The pipeline supports adaptive thresholds that adjust based on data characteristics:

| Config Parameter | Default | Description |
|------------------|---------|-------------|
| `use_adaptive_base_mean` | True | Adjust baseMean threshold per condition |
| `base_mean_percentile` | 5.0 | Percentile for baseMean cutoff |
| `strong_z_for_relaxed_bc1` | 5.0 | Z-score threshold for relaxed BC1 requirement |
| `low_cv_for_relaxed_bc1` | 0.3 | CV threshold for relaxed BC1 requirement |
| `high_concordance_for_relaxed_bc1` | 0.9 | Concordance for relaxed BC1 requirement |

## Advanced Analysis Modules

### Library Calibration

Corrects for per-library systematic biases:

```python
from bc1_from_screen_pipeline.core.library_calibration import (
    calculate_library_bias,
    center_by_library,
    run_bias_correction_pipeline,
)

# Calculate and apply bias correction
bias_stats = calculate_library_bias(df, library_col="yeast_sublibrary")
corrected_df = center_by_library(df, bias_stats)
```

### Variance Calibration

Uses synonymous variants as negative controls to calibrate expected variance:

```python
from bc1_from_screen_pipeline.core.variance_calibration import (
    calibrate_variance_from_synonymous,
)
```

### Noise Assessment

Identifies and penalizes noisy conditions and loci:

```python
from bc1_from_screen_pipeline.core.noise_assessment import (
    identify_noisy_conditions,
    identify_noisy_loci,
    apply_noise_penalties,
    run_noise_assessment_pipeline,
)
```

### Causality Framework

Evaluates causal evidence for variant effects:

```python
from bc1_from_screen_pipeline.core.causality import (
    CausalityEvidence,
    CausalityEvaluator,
    aggregate_to_aa_level,
)
```

## Sherlock Best Practices

This pipeline follows Stanford Sherlock HPC best practices:

1. **Use environment variables** for storage paths (`$OAK`, `$SCRATCH`)
2. **Stage data to Scratch** for high-bandwidth I/O during jobs
3. **Copy results back to Oak** for archival storage
4. **Initialize conda in sbatch scripts**, not `~/.bashrc`

When using `--use-scratch`:
- Input files are copied from Oak to `$SCRATCH`
- All processing happens on Scratch (75GB/s bandwidth)
- Final results are copied back to Oak

## Configuration

### YAML Configuration

Save configuration to a file:

```python
config.to_yaml(Path("my_config.yaml"))
```

Example `config.yaml`:
```yaml
project_name: QTL
screen_name: 20250912_SpG_QTL_screen
project_dir: /path/to/projects/QTL/20250912_SpG_QTL_screen

libraries:
  - name: yL437
    screen_dates: ["20251014", "20251018"]
    bottlenecked: true
    timepoints: [0, 6, 12, 21]
    t0_consolidation_groups: 4
  - name: yL442
    screen_dates: ["20251212", "20251216"]
    bottlenecked: false
    timepoints: [0, 21]
    t0_consolidation_groups: 8

tier_config:
  tier_1_min_reads: 100
  tier_2_min_reads: 20

hit_calling_config:
  tier_1_z_threshold: 3.0
  tier_1_concordance_threshold: 0.75
  tier_2_z_threshold: 3.75
  tier_2_concordance_threshold: 0.85
```

## Output Files

| File | Location | Description |
|------|----------|-------------|
| `bc1_tier_mapping_{library}.tsv` | bc1_counts/ | BC1 → tier assignment |
| `comparison_manifest.tsv` | counts_for_deseq/ | All comparisons |
| `{comparison}_DESeq2_results.tsv` | DESeq2_consolidated_t0/{library}/ | Raw DESeq2 |
| `{comparison}_annotated.tsv` | per_comparison/{library}/ | With tier + variant info |
| `hits_summary_by_tier.tsv` | hit_calling_tiered/ | Hit summary |

## Requirements

- Python 3.8+
- pandas
- numpy
- scipy
- pydeseq2 (for DESeq2 step)
- tqdm
- pyyaml

## Author

Kevin Roy Lab, Stanford University
