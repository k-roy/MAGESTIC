"""Select one SLURM array task's sample through the shared WGS manifest reader.

The array template used to read physical manifest line N+1 with awk. That diverged
from `core.manifest.load_manifest`, which the direct and Snakemake paths use: a
space truncated a sample ID, comment and blank lines consumed task IDs, reordered
columns selected a BAM path, and CONFIG and MANIFEST could name different cohorts.
This module is now the only place the shell route enumerates samples, so array task
N is the Nth sample that `run_step_1` and the Snakefile would run.

    python -m magestic.pipelines.wgs.slurm.task_select --config C --manifest M --task N
    python -m magestic.pipelines.wgs.slurm.task_select --config C --manifest M --count

`--count` sizes the array (`sbatch --array=1-<count>`). Every refusal exits 2 with a
message on stderr and prints nothing on stdout, so the template never runs a task
for a row it could not identify.

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

from ..config import load_config
from ..pipeline._common import manifest_samples

_TASK_ID = re.compile(r"[1-9][0-9]*")


class TaskSelectionError(ValueError):
    """The array task cannot be mapped to exactly one manifest sample."""


def cohort_samples(config_path: str | Path, manifest_path: str | Path) -> list[str]:
    """Manifest samples in fan-out order, after proving MANIFEST is CONFIG's manifest.

    The stage CLI reads `paths.manifest` from CONFIG, so a different MANIFEST would
    size and index the array from one cohort while running samples from another.
    """
    config = load_config(config_path)
    declared = Path(config.paths.manifest)
    given = Path(manifest_path)
    try:
        same = os.path.samefile(declared, given)
    except OSError as exc:
        raise TaskSelectionError(f"cannot compare MANIFEST {given} with CONFIG's paths.manifest {declared}: {exc}") from exc
    if not same:
        raise TaskSelectionError(
            f"MANIFEST {given} is not the manifest CONFIG declares ({declared}); "
            "the array would be indexed from one cohort and run from another"
        )
    return [record.sample for record in manifest_samples(config)]


def select_sample(samples: list[str], task_id: str) -> str:
    """The sample for a 1-based SLURM array task ID; anything else is refused."""
    if not _TASK_ID.fullmatch(task_id or ""):
        raise TaskSelectionError(f"SLURM_ARRAY_TASK_ID must be a positive integer, got {task_id!r}")
    index = int(task_id)
    if index > len(samples):
        raise TaskSelectionError(
            f"array task {index} exceeds the {len(samples)} manifest sample(s); "
            "size the array with --count"
        )
    return samples[index - 1]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    parser.add_argument("--config", required=True)
    parser.add_argument("--manifest", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--task", help="1-based SLURM_ARRAY_TASK_ID")
    mode.add_argument("--count", action="store_true", help="print the number of manifest samples")
    args = parser.parse_args(argv)
    try:
        samples = cohort_samples(args.config, args.manifest)
        print(len(samples) if args.count else select_sample(samples, args.task))
    except (TaskSelectionError, ValueError, OSError) as exc:
        print(f"task_select: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
