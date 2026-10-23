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
    magestic-redi-rearray      step_h (legacy command name) → rearray_main
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
    """`magestic-redi-rearray` — registered-position plots (legacy name).

    Collection-specific selection and robot worklists are not implemented by
    this command. Stage I exposes a generic ranking helper only.
    """
    from .pipeline.step_h_plot import run_step_h

    args = _common_parser("magestic-redi-rearray").parse_args(argv)
    config = load_config(args.config)
    for k, v in run_step_h(config).items():
        print(f"  {k}: {v}")
    print(
        "Plots complete. Collection-specific clone selection and robot worklists "
        "require a separate implementation; this command does not generate them."
    )
    return 0


def _plate_list(text: str) -> list[str]:
    """'5169-5172,5180' -> ['5169', '5170', '5171', '5172', '5180']."""
    out = []
    for part in text.split(','):
        part = part.strip()
        if '-' in part:
            a, b = part.split('-', 1)
            if not (a.isdigit() and b.isdigit()) or int(b) < int(a):
                raise SystemExit(f'bad plate range {part!r}')
            out += [str(n) for n in range(int(a), int(b) + 1)]
        elif part:
            out.append(part)
    return out


def worklist_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-worklist` — Stage I: explicit picks -> PIXL run files (new generation directory)."""
    import csv
    from .core.worklist import WorklistPlan, edge_wells
    from .pipeline.step_i_worklist import run_step_i_worklist

    p = argparse.ArgumentParser(prog='magestic-redi-worklist')
    p.add_argument('--picks', required=True, type=Path, help='TSV: source_plate, source_row, source_col[, source_density]')
    p.add_argument('--target-plates', required=True, help="pre-printed target barcodes, e.g. '5169-5186'")
    p.add_argument('--outdir', required=True, type=Path, help='new, empty generation directory')
    p.add_argument('--prefix', required=True)
    p.add_argument('--target-density', type=int, default=384, choices=(96, 384, 1536))
    p.add_argument('--fill-order', default='row', choices=('row', 'column'))
    p.add_argument('--exclude-edge', action='store_true', help='leave outer rows and columns of each target empty')
    p.add_argument('--exclude-wells', default='', help="extra empty target wells, e.g. 'A1,P24'")
    p.add_argument('--first-target-well-index', type=int, default=0)
    p.add_argument('--max-sources-per-run', type=int, default=2)
    p.add_argument('--source-rename', type=Path, help='TSV old_source_plate, new_source_plate (restamped plates)')
    p.add_argument('--keep-input-order', action='store_true')
    a = p.parse_args(argv)
    exclude = edge_wells(a.target_density) if a.exclude_edge else set()
    for w in filter(None, (x.strip() for x in a.exclude_wells.split(','))):
        exclude.add((w.rstrip('0123456789'), int(w[len(w.rstrip('0123456789')):])))
    rename = {}
    if a.source_rename:
        with open(a.source_rename, newline='') as fh:
            rename = {r['old_source_plate']: r['new_source_plate'] for r in csv.DictReader(fh, delimiter='\t')}
    plan = WorklistPlan(target_plates=_plate_list(a.target_plates), target_density=a.target_density,
                        fill_order=a.fill_order, exclude_target_wells=exclude,
                        first_target_well_index=a.first_target_well_index,
                        max_sources_per_run=a.max_sources_per_run, source_rename=rename,
                        keep_input_order=a.keep_input_order)
    print(run_step_i_worklist(a.picks, plan, a.outdir, a.prefix))
    return 0


def worklist_validate_main(argv: list[str] | None = None) -> int:
    """`magestic-redi-worklist-validate DIR` — re-check a written generation; exit 1 on any mismatch."""
    import json
    from .pipeline.step_i_worklist import validate_generation

    p = argparse.ArgumentParser(prog='magestic-redi-worklist-validate')
    p.add_argument('outdir', type=Path)
    report = validate_generation(p.parse_args(argv).outdir)
    print(json.dumps(report, indent=2))
    return 0 if report['ok'] else 1
