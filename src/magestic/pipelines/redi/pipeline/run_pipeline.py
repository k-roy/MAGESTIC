"""
Top-level REDI pipeline orchestrator.

Runs step_b through step_h end-to-end on a single config (step_a is per-sample
and typically driven by Snakemake / SLURM array).

Usage from CLI: `magestic-redi --config path/to/screen.yaml [--from STEP] [--to STEP]`.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..core.config import REDIConfig, load_config
from .step_b_combine_keys import run_step_b
from .step_c_process_bcs import run_step_c
from .step_d_purity import run_step_d
from .step_e_annotate import run_step_e
from .step_f_cross_array import run_step_f
from .step_g_wgs_check import run_step_g
from .step_h_plot import run_step_h

STEPS = ("b", "c", "d", "e", "f", "g", "h")


def run(config: REDIConfig, from_step: str = "b", to_step: str = "h") -> dict[str, object]:
    """Run the pipeline from `from_step` through `to_step`. Returns dict of
    artifacts produced."""
    if from_step not in STEPS or to_step not in STEPS:
        raise ValueError(f"steps must be one of {STEPS}")
    start = STEPS.index(from_step)
    end = STEPS.index(to_step) + 1

    artifacts: dict[str, object] = {}
    combined_key_path = config.paths.keyfiles / "combined_key.tsv"

    for step in STEPS[start:end]:
        if step == "b":
            artifacts["b_combined_key"] = run_step_b(config)
            combined_key_path = artifacts["b_combined_key"]
        elif step == "c":
            artifacts["c_combined"], artifacts["c_magestic_only"] = run_step_c(
                config, combined_key_path
            )
        elif step == "d":
            artifacts["d_top_clone"] = run_step_d(config)
        elif step == "e":
            artifacts["e_annotated"] = run_step_e(config)
        elif step == "f":
            artifacts["f_cross_array"] = run_step_f(config)
        elif step == "g":
            artifacts["g_wgs_check"] = run_step_g(config)
        elif step == "h":
            artifacts["h_plots"] = run_step_h(config)
    return artifacts


def main(argv: list[str] | None = None) -> int:
    """Entry point for `magestic-redi` umbrella CLI."""
    p = argparse.ArgumentParser(prog="magestic-redi")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--from", dest="from_step", default="b", choices=STEPS)
    p.add_argument("--to", dest="to_step", default="h", choices=STEPS)
    args = p.parse_args(argv)
    config = load_config(args.config)
    artifacts = run(config, from_step=args.from_step, to_step=args.to_step)
    for key, val in artifacts.items():
        print(f"  {key}: {val}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
