"""
Step A — trim + merge + collapse + counts-table.

Generates the per-sample `*_R1_R2_counts.tsv` table consumed by step_c.

This is a thin Python orchestrator around the shell pipeline (bbduk → bbmerge
→ collapse_identical_merged_reads_min_counts.py → collapsed_reads_to_counts_table.py).
The orchestrator is intended to be called from Snakemake (per-sample rule) or
SLURM (per-sample array task); it does NOT replace `step_a_trim_merge_collapse_counts_tables.sh`
which is preserved as the canonical reference for the exact bbduk/bbmerge
arguments.

Legacy reference:
    /oak/.../REDI/20250312_NNS_complex_satmut/scripts/step_a_trim_merge_collapse_counts_tables.sh
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from ..core.config import REDIConfig

# bbduk trimming literal — TruSeq adapter (5'→3')
TRUSEQ_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA"

# Path to the general-scripts dir (shared with NRD1 / bc1_donor_bc0). May be
# overridden via env var REDI_GENERAL_SCRIPTS_DIR.
DEFAULT_GENERAL_SCRIPTS_DIR = (
    "/path/to/"
    "scripts_and_keyfiles/general_scripts/"
)


def run_step_a_one_sample(
    prefix: str,
    config: REDIConfig,
    *,
    general_scripts_dir: str = DEFAULT_GENERAL_SCRIPTS_DIR,
    min_counts: int = 1,
    bbduk_ftl: int = 20,
    overwrite: bool = True,
) -> None:
    """Run trim + merge + collapse + counts-table for one prefix.

    Idempotent (skips if `counts_tables / <prefix>_R1_R2_counts.tsv` exists).
    """
    paths = config.paths
    paths.trimmed_fastq.mkdir(parents=True, exist_ok=True)
    paths.merged_fastq.mkdir(parents=True, exist_ok=True)
    paths.collapsed_fastq.mkdir(parents=True, exist_ok=True)
    paths.counts_tables.mkdir(parents=True, exist_ok=True)

    final_out = paths.counts_tables / f"{prefix}_R1_R2_counts.tsv"
    if final_out.exists():
        return

    fq1 = paths.raw_fastq / f"{prefix}_R1.fastq.gz"
    fq2 = paths.raw_fastq / f"{prefix}_R2.fastq.gz"
    trim1 = paths.trimmed_fastq / f"{prefix}_trimmed_R1.fastq"
    trim2 = paths.trimmed_fastq / f"{prefix}_trimmed_R2.fastq"
    merged = paths.merged_fastq / f"{prefix}_trimmed_merged.fastq"
    collapsed = paths.collapsed_fastq / f"{prefix}_collapsed.fastq"

    _bbduk(fq1, fq2, trim1, trim2, ftl=bbduk_ftl, overwrite=overwrite)
    _bbmerge(trim1, trim2, merged, overwrite=overwrite)
    _collapse(general_scripts_dir, merged, collapsed, min_counts=min_counts)
    _counts_table(general_scripts_dir, collapsed, prefix, final_out)


def _bbduk(
    in1: Path,
    in2: Path,
    out1: Path,
    out2: Path,
    *,
    ftl: int = 20,
    overwrite: bool = True,
) -> None:
    """Adapter + quality trim. Mirrors legacy step_a sh args exactly."""
    cmd = [
        "bbduk.sh",
        f"in1={in1}",
        f"in2={in2}",
        f"out1={out1}",
        f"out2={out2}",
        f"literal={TRUSEQ_ADAPTER}",
        f"overwrite={'true' if overwrite else 'false'}",
        f"ftl={ftl}",
        "ktrim=r",
        "k=23",
        "mink=11",
        "hdist=1",
        "tpe",
        "tbo",
    ]
    subprocess.run(cmd, check=True)


def _bbmerge(in1: Path, in2: Path, out: Path, *, overwrite: bool = True) -> None:
    cmd = [
        "bbmerge.sh",
        f"in1={in1}",
        f"in2={in2}",
        f"out={out}",
        f"overwrite={'true' if overwrite else 'false'}",
    ]
    subprocess.run(cmd, check=True)


def _collapse(
    general_scripts_dir: str,
    merged: Path,
    collapsed: Path,
    *,
    min_counts: int = 1,
) -> None:
    script = Path(general_scripts_dir) / "collapse_identical_merged_reads_min_counts.py"
    subprocess.run(
        ["python", str(script), str(merged), str(collapsed), str(min_counts)],
        check=True,
    )


def _counts_table(
    general_scripts_dir: str,
    collapsed: Path,
    sample_name: str,
    out: Path,
) -> None:
    script = Path(general_scripts_dir) / "collapsed_reads_to_counts_table.py"
    subprocess.run(
        ["python", str(script), str(collapsed), sample_name, str(out)],
        check=True,
    )
