"""
magestic.pipelines.redi._cli — CLI entry points.

Registered in pyproject.toml `[project.scripts]`:

    magestic-redi              umbrella → umbrella_main
    magestic-redi-trim         step_a   → trim_main
    magestic-redi-combine-keys step_b   → combine_keys_main
    magestic-redi-extract-bc   step_c   → extract_bc_main
    magestic-redi-purity       step_d   → purity_main
    magestic-redi-annotate     step_e   → annotate_main
    magestic-redi-cross-array  step_f   → cross_array_main
    magestic-redi-wgs-check    step_g   → wgs_check_main
    magestic-redi-rearray      step_h_i → rearray_main
"""

from __future__ import annotations

import argparse
from pathlib import Path

from .core.config import load_config


def _common_parser(prog: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog=prog)
    p.add_argument("--config", required=True, type=Path)
    return p


def umbrella_main(argv: list[str] | None = None) -> int:
    """`magestic-redi` umbrella — runs the full pipeline."""
    from .pipeline.run_pipeline import main as run_main

    return run_main(argv)


def trim_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-trim` — step_a per-sample trim/merge/collapse/counts."""
    from .pipeline.step_a_trim_counts import run_step_a_one_sample

    p = _common_parser("magestic-redi-trim")
    p.add_argument("--sample", required=True, help="FASTQ prefix (no _R1.fastq.gz suffix)")
    args = p.parse_args(argv)
    config = load_config(args.config)
    run_step_a_one_sample(args.sample, config)
    return 0


def combine_keys_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-combine-keys` — step_b combined_key.tsv."""
    from .pipeline.step_b_combine_keys import run_step_b

    args = _common_parser("magestic-redi-combine-keys").parse_args(argv)
    print(run_step_b(load_config(args.config)))
    return 0


def extract_bc_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-extract-bc` — step_c bc1/REDI_bc extraction."""
    from .pipeline.step_c_process_bcs import run_step_c

    args = _common_parser("magestic-redi-extract-bc").parse_args(argv)
    config = load_config(args.config)
    out, mago_only = run_step_c(
        config, config.paths.keyfiles / "combined_key.tsv"
    )
    print(out)
    print(mago_only)
    return 0


def purity_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-purity` — step_d colony purity."""
    from .pipeline.step_d_purity import run_step_d

    args = _common_parser("magestic-redi-purity").parse_args(argv)
    print(run_step_d(load_config(args.config)))
    return 0


def annotate_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-annotate` — step_e annotate with bc1 reference."""
    from .pipeline.step_e_annotate import run_step_e

    args = _common_parser("magestic-redi-annotate").parse_args(argv)
    print(run_step_e(load_config(args.config)))
    return 0


def cross_array_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-cross-array` — step_f cross-array agreement."""
    from .pipeline.step_f_cross_array import run_step_f

    args = _common_parser("magestic-redi-cross-array").parse_args(argv)
    for path in run_step_f(load_config(args.config)):
        print(path)
    return 0


def wgs_check_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-wgs-check` — step_g WGS comparison."""
    from .pipeline.step_g_wgs_check import run_step_g

    args = _common_parser("magestic-redi-wgs-check").parse_args(argv)
    print(run_step_g(load_config(args.config)))
    return 0


def rearray_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-rearray` — step_h/i clone selection + rearray.

    Project-specific selection logic lives in `REDI/screens/<screen>/driver.py`;
    this CLI runs the engine plotting (step_h) and points the user at the
    screen driver for selection.
    """
    from .pipeline.step_h_plot import run_step_h

    args = _common_parser("magestic-redi-rearray").parse_args(argv)
    config = load_config(args.config)
    for k, v in run_step_h(config).items():
        print(f"  {k}: {v}")
    print(
        "Project-specific clone selection: see "
        f"REDI/screens/{config.screen_name}/driver.py"
    )
    return 0
