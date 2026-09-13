"""
Top-level MAGESTIC WGS pipeline orchestrator — the DIRECT path (A).

Runs steps 1-5 end-to-end on one config. The SNAKEMAKE path (B) invokes the same
per-step console scripts from its `shell:` blocks, so each stage has exactly one
entry point and dual-path byte-identity (`ARCHITECTURE_POLICY` §4) is structural
rather than a convention someone has to maintain.

Usage:
    magestic-wgs --config path/to/<cohort>_wgs.yaml [--steps 1,3,4]
                 [--from-step 1 --end-step 5] [--use-scratch] [--dry-run] [--debug]

**The stages are a diamond, not a chain** (spec §3, D2): 1 -> 4 is a hard
dependency (`reconcile_cell` scores a cell against its bc1 multiplicity, and its
`0/0/"high"` defaults are wrong rather than neutral); Stage 3's Half A is
independent while Half B needs an oligo assignment; Stage 2 is independent of both
and feeds only the Stage 5 join. `STEPS` is a valid TOPOLOGICAL order, not a claim
that each step consumes its predecessor.

**Step 2 is GATED, and a gated run exits non-zero.** Three measured gates (no
synthetic-donor input, no gVCF anywhere, and a crash on any contig name containing
an underscore — which every mandated composite reference has) stand between it and
a correct answer. A gated step is recorded as a `"GATED: ..."` string in the
artifacts dict, never as a path, and the run exits **2** unless `--allow-gated` is
given. Exit 0 from this pipeline means every requested step produced counted rows.

Author: Kevin R. Roy
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from ..config import WGSConfig, load_config
from ._common import GateError

logger = logging.getLogger(__name__)

STEPS = ("1", "2", "3", "4", "5")

STEP_TITLES = {
    "1": "bc1 calling (TRACE bc1_consensus)",
    "2": "donor-outcome pre-calling (TRACE donor_outcome) — GATED",
    "3": "plasmid integration + donor-read cleaning (TRACE sv_plasmid)",
    "4": "cell reconciliation (TRACE reconcile)",
    "5": "Stage 2 x Stage 4 outcome join (orchestration-owned)",
}

ARTIFACT_KEYS = {
    "1": "1_bc1_calls",
    "2": "2_donor_outcome",
    "3": "3_sv_plasmid",
    "4": "4_reconcile",
    "5": "5_outcome_join",
}

EXIT_OK = 0
EXIT_FAILED = 1
EXIT_GATED = 2


def select_steps(
    steps: str | None = None,
    from_step: str = "1",
    end_step: str = "5",
) -> list[str]:
    """Resolve `--steps` / `--from-step` / `--end-step` into an ordered step list.

    `--steps` wins when given; the result is always in canonical `STEPS` order, so
    a subset selected in any order runs in dependency order.
    """
    if steps:
        wanted = [s.strip() for s in steps.split(",") if s.strip()]
        unknown = [s for s in wanted if s not in STEPS]
        if unknown:
            raise ValueError(f"unknown step(s) {unknown}; valid steps are {list(STEPS)}")
        return [s for s in STEPS if s in set(wanted)]
    if from_step not in STEPS or end_step not in STEPS:
        raise ValueError(f"steps must be one of {STEPS}")
    start, end = STEPS.index(from_step), STEPS.index(end_step) + 1
    if start >= end:
        raise ValueError(f"--from-step {from_step} is after --end-step {end_step}")
    return list(STEPS[start:end])


def planned_outputs(config: WGSConfig, steps: list[str]) -> dict[str, Path]:
    """Resolved cohort output path per requested step, touching no filesystem.

    What `--dry-run` prints. Importing the step modules is what pulls in TRACE and
    pysam, so a dry run without the `magestic[wgs]` extra installed fails at import
    — which is the honest answer to "would this run?".
    """
    from . import (
        step_1_bc1_call,
        step_2_donor_precall,
        step_3_sv_plasmid,
        step_4_reconcile,
        step_5_join_outcomes,
    )

    resolvers = {
        "1": step_1_bc1_call.table_path,
        "2": step_2_donor_precall.table_path,
        "3": step_3_sv_plasmid.table_path,
        "4": step_4_reconcile.table_path,
        "5": step_5_join_outcomes.table_path,
    }
    return {ARTIFACT_KEYS[s]: resolvers[s](config) for s in steps}


def run(
    config: WGSConfig,
    from_step: str = "1",
    to_step: str = "5",
    *,
    steps: str | None = None,
) -> dict[str, object]:
    """Run the requested steps in order. Returns `{artifact key: Path | "GATED: ..."}`.

    A `GateError` is RECORDED, not raised, so the remaining steps still run and the
    gate is visible in the artifacts dict as a string rather than a path. Every
    other exception propagates: a failed step must reach the caller, and `main()`
    turns it into a non-zero exit (§7 — no silent success).
    """
    from . import (
        step_1_bc1_call,
        step_2_donor_precall,
        step_3_sv_plasmid,
        step_4_reconcile,
        step_5_join_outcomes,
    )

    runners = {
        "1": step_1_bc1_call.run_step_1,
        "2": step_2_donor_precall.run_step_2,
        "3": step_3_sv_plasmid.run_step_3,
        "4": step_4_reconcile.run_step_4,
        "5": step_5_join_outcomes.run_step_5,
    }
    artifacts: dict[str, object] = {}
    for step in select_steps(steps, from_step, to_step):
        logger.info("=" * 70)
        logger.info("step %s — %s", step, STEP_TITLES[step])
        logger.info("=" * 70)
        key = ARTIFACT_KEYS[step]
        try:
            artifacts[key] = runners[step](config)
        except GateError as exc:
            logger.warning("step %s is GATED and did not run: %s", step, exc)
            artifacts[key] = f"GATED: {exc}"
    return artifacts


def main(argv: list[str] | None = None) -> int:
    """Entry point for the `magestic-wgs` umbrella CLI."""
    p = argparse.ArgumentParser(
        prog="magestic-wgs",
        description="MAGESTIC WGS pipeline (orchestration over TRACE).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Stage 1 only, over every sample in the manifest
  magestic-wgs --config workflow/config/screens/<cohort>_wgs.yaml --steps 1

  # The runnable stages (step 2 is GATED; see the module docstring)
  magestic-wgs --config <cohort>_wgs.yaml --steps 1,3,4

  # Show what would be written, touch nothing
  magestic-wgs --config <cohort>_wgs.yaml --dry-run

Exit codes: 0 = every requested step produced counted rows,
            1 = a step failed, 2 = a requested step is GATED.
        """,
    )
    p.add_argument("--config", required=True, type=Path, help="Cohort config YAML")
    p.add_argument(
        "--project-dir",
        type=Path,
        default=None,
        help="Override the project root the config resolves paths against",
    )
    p.add_argument(
        "--steps",
        default=None,
        help=f"Comma-separated subset, e.g. 1,3,4. Available: {','.join(STEPS)}",
    )
    p.add_argument("--from-step", default="1", choices=STEPS)
    p.add_argument("--end-step", default="5", choices=STEPS)
    p.add_argument("--use-scratch", action="store_true", help="Stage I/O through $SCRATCH")
    p.add_argument("--allow-gated", action="store_true", help="Exit 0 even if a step is GATED")
    p.add_argument("--dry-run", action="store_true", help="Print planned outputs; write nothing")
    p.add_argument("--debug", action="store_true", help="Verbose logging")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    try:
        requested = select_steps(args.steps, args.from_step, args.end_step)
    except ValueError as exc:
        logger.error("%s", exc)
        return EXIT_FAILED

    config = load_config(
        args.config, project_dir=args.project_dir, use_scratch=args.use_scratch
    )

    if args.dry_run:
        print(f"config: {args.config}")
        for step in requested:
            print(f"  step {step} — {STEP_TITLES[step]}")
        for key, path in planned_outputs(config, requested).items():
            print(f"    {key}: {path}")
        print("dry run: nothing was written")
        return EXIT_OK

    # Scoped to the REQUESTED steps: a stage nobody asked to run must not fail the
    # run for a missing input. The pre-registered acceptance run is `--steps 1` on a
    # cohort whose Stage-3/4 inputs are gated and absent (spec §7.3, §8).
    errors = config.validate(steps=requested)
    if errors:
        for error in errors:
            logger.error("config error: %s", error)
        return EXIT_FAILED
    config.ensure_directories()

    artifacts = run(
        config, from_step=args.from_step, to_step=args.end_step, steps=args.steps
    )
    for key, val in artifacts.items():
        if isinstance(val, (list, tuple)):
            for path in val:
                print(f"  {key}: {path}")
        else:
            print(f"  {key}: {val}")

    gated = [k for k, v in artifacts.items() if isinstance(v, str)]
    if gated and not args.allow_gated:
        logger.error(
            "%d of %d requested step(s) were GATED and produced no output: %s. "
            "Exiting %d so this is not read as a successful run.",
            len(gated),
            len(artifacts),
            ", ".join(gated),
            EXIT_GATED,
        )
        return EXIT_GATED
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
