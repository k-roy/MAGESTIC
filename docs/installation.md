# Installation

MAGESTIC requires **Python ≥ 3.10**.

## Core install

```bash
pip install magestic
```

The core install pulls in the design and genomics stack (`numpy`, `pandas`,
`scipy`, `biopython`, `pyyaml`, `tqdm`, `matplotlib`, `seaborn`, `statsmodels`)
— enough for [library design](quickstart/design.md) and the shared genomics /
plotting utilities.

## Optional extras

Each analysis pipeline declares its own dependency group so you only install
what you use:

| Extra | Installs | Needed for |
|---|---|---|
| `align` | `rapidfuzz`, `python-levenshtein`, `pysam` | `guide_donor_bc0`, `bc1_donor_bc0` (read mapping & fuzzy donor matching) |
| `screen` | `pydeseq2`, `joblib`, `openpyxl` | `bc1_from_screen` (DESeq2 screen analysis) |
| `redi` | `openpyxl`, `pyarrow`, `python-levenshtein` | REDI clone-array engine |
| `snakemake` | `snakemake` | running any pipeline as a Snakemake workflow |
| `ihw` | `rpy2` (needs **R**) | Independent Hypothesis Weighting in the screen pipeline |
| `full` | `align` + `screen` + `snakemake` + `redi` | everything except `ihw` |
| `dev` | `pytest`, `pytest-cov`, `black`, `ruff`, `mypy` | development & tests |

```bash
pip install "magestic[align]"
pip install "magestic[screen]"
pip install "magestic[full]"
```

!!! note "IHW is optional"
    The `[ihw]` extra wraps the R/Bioconductor `IHW` package via `rpy2`. If R or
    `rpy2` is unavailable, `bc1_from_screen` automatically falls back to
    Benjamini–Hochberg multiple-testing correction — the pipeline still runs.

## Install from source

```bash
git clone https://github.com/k-roy/MAGESTIC.git
cd MAGESTIC
pip install -e ".[full,dev]"
```

The package uses the `src/` layout with a `hatchling` build backend; the
editable install exposes all the `magestic-*` console scripts.

## Verify

```bash
python -c "import magestic; print(magestic.__version__)"
magestic-pam-search --help
```

## HPC / cluster

The pipelines are built to run as **Snakemake** workflows and SLURM job arrays
for large screens. See [HPC / SLURM](user_guide/hpc_slurm.md) for cluster
submission patterns. On clusters where you cannot `pip install` system-wide, a
per-project conda environment works well:

```bash
conda create -n magestic python=3.11
conda activate magestic
pip install -e ".[full]"
```
