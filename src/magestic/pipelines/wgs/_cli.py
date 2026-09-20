"""
magestic.pipelines.wgs._cli — CLI entry points.

Registered in pyproject.toml `[project.scripts]`:

    magestic-wgs               umbrella → umbrella_main
    magestic-wgs-bc1           step_1   → bc1_main
    magestic-wgs-donor-precall step_2   → donor_precall_main   (GATED — see below)
    magestic-wgs-sv-plasmid    step_3   → sv_plasmid_main
    magestic-wgs-reconcile     step_4   → reconcile_main
    magestic-wgs-join          step_5   → join_outcomes_main

The Snakemake path (B) invokes these same console scripts from its `shell:`
blocks, so there is exactly one entry point per stage and the byte-identical
parity contract (ARCHITECTURE_POLICY §4) is structural rather than a convention
someone has to maintain.

Nothing in this file touches TRACE. `trace_crispr.wgs` is called from
`pipeline/step_*.py` and `core/*.py` only; this layer parses arguments, calls
one step function, prints the path(s) it wrote and returns an exit code.

Step-function contract — the import each `*_main` lazy-imports in its body:

    from .pipeline.step_1_bc1_call      import run_step_1   # (config, sample=None) -> Path
    from .pipeline.step_2_donor_precall import run_step_2   # (config, sample=None) -> Path
    from .pipeline.step_3_sv_plasmid    import run_step_3   # (config, sample=None) -> Iterable[Path]
    from .pipeline.step_4_reconcile     import run_step_4   # (config, sample=None) -> Path
    from .pipeline.step_5_join_outcomes import run_step_5   # (config) -> Iterable[Path]

`--sample` selects a single sample (one Snakemake wildcard / one SLURM array
task); omitting it runs every sample in the manifest and writes the
cohort-level table. **The fan-out loop and the per-task skip-check
(ARCHITECTURE_POLICY §6) live in `pipeline/`, never here** — a loop in this file
would give paths (A) and (B) two different orchestrations and §4 parity would
stop being structural.

Per ARCHITECTURE_POLICY §7 (no silent success) no step call is wrapped in
`try`/`except`: a failure propagates as a traceback and a non-zero exit status.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .config import load_config


def _common_parser(prog: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog=prog)
    p.add_argument("--config", required=True, type=Path)
    return p


def _fanout_parser(prog: str) -> argparse.ArgumentParser:
    """`_common_parser` plus the per-sample fan-out selector.

    The fan-out unit is the sample for stages 1-4 (planning/166 §3).
    """
    p = _common_parser(prog)
    p.add_argument(
        "--sample",
        default=None,
        help=(
            "Sample name as it appears in the manifest — runs that one sample "
            "(one Snakemake wildcard / one SLURM array task). Omit to run every "
            "sample in the manifest and write the cohort-level table."
        ),
    )
    return p


def umbrella_main(argv: list[str] | None = None) -> int:
    """`magestic-wgs` umbrella — runs the full pipeline."""
    from .pipeline.run_pipeline import main as run_main

    return run_main(argv)


def bc1_main(argv: list[str] | None = None) -> int:
    """`magestic-wgs-bc1` — stage 1 de novo bc1 calling from the aligned BAM."""
    from .pipeline.step_1_bc1_call import run_step_1

    args = _fanout_parser("magestic-wgs-bc1").parse_args(argv)
    print(run_step_1(load_config(args.config), sample=args.sample))
    return 0


def donor_precall_main(argv: list[str] | None = None) -> int:
    """`magestic-wgs-donor-precall` — stage 2 donor-window editing outcome.

    GATED (planning/166 §8): the SYNTHETIC_ERROR axis ships disabled for the
    2021 pool, no gVCF exists to separate `designed_variant_unedited` from
    `unknown`, and TRACE's underscore-contig crash makes the stage incompatible
    with the mandated composite reference. The refusal lives in `step_2` /
    `qc.pipeline_qc`, not here.
    """
    from .pipeline.step_2_donor_precall import run_step_2

    args = _fanout_parser("magestic-wgs-donor-precall").parse_args(argv)
    print(run_step_2(load_config(args.config), sample=args.sample))
    return 0


def sv_plasmid_main(argv: list[str] | None = None) -> int:
    """`magestic-wgs-sv-plasmid` — stage 3 plasmid integration + donor-read cleaning."""
    from .pipeline.step_3_sv_plasmid import run_step_3

    args = _fanout_parser("magestic-wgs-sv-plasmid").parse_args(argv)
    for path in run_step_3(load_config(args.config), sample=args.sample):
        print(path)
    return 0


def reconcile_main(argv: list[str] | None = None) -> int:
    """`magestic-wgs-reconcile` — stage 4 per-cell edit reconciliation."""
    from .pipeline.step_4_reconcile import run_step_4

    args = _fanout_parser("magestic-wgs-reconcile").parse_args(argv)
    print(run_step_4(load_config(args.config), sample=args.sample))
    return 0


def join_outcomes_main(argv: list[str] | None = None) -> int:
    """`magestic-wgs-join` — stage 5 join of stage 2 x stage 4.

    Cohort-level, so it takes no `--sample`: the agreement matrix and the
    disagreement list are the whole point (planning/166 D2).
    """
    from .pipeline.step_5_join_outcomes import run_step_5

    args = _common_parser("magestic-wgs-join").parse_args(argv)
    for path in run_step_5(load_config(args.config)):
        print(path)
    return 0
