"""
Top-level REDI pipeline orchestrator.

Runs step_b through step_h end-to-end on a single config (step_a is per-sample
and typically driven by Snakemake / SLURM array).

Usage from CLI: `magestic-redi --config path/to/screen.yaml [--from STEP] [--to STEP]`.
"""

from __future__ import annotations

import argparse
import json
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


def registered_sample_names(config: REDIConfig) -> list[str]:
    """Derive file prefixes from explicit dated sample/assay registration.

    The sample key does not need a precomputed sample_name. Empty or conflicting
    registrations fail rather than creating a vacuous workflow fan-out.
    """
    from ..core.identity import keyed_frame
    from ..readers.keyfiles import combine_keyfiles
    keys = config.paths.keyfiles
    metadata = combine_keyfiles(keys/config.sample_key,keys/config.gdna_plate_key,
                                keys/config.redi_plate_barcode_key)
    date = str(config.sequencing_date)
    if 'sequencing_date' in metadata and not metadata.sequencing_date.eq(date).all():
        raise ValueError('declared sequencing_date conflicts with configured dated assay registration')
    names = date+'_'+metadata.outer_primers+'_'+metadata.sample_number
    if 'sample_name' in metadata and not metadata.sample_name.equals(names):
        raise ValueError('declared sample_name conflicts with dated assay/sample registration')
    metadata = keyed_frame(metadata.assign(sample_name=names),['sample_name'],label='workflow sample registration')
    if metadata.empty:
        raise ValueError('workflow requires at least one registered sample')
    result = metadata.sample_name.tolist()
    if any(name in ('.','..') or any(char in name for char in '/\\\x00\t\r\n') for name in result):
        raise ValueError('registered sample prefix must be a safe file basename')
    return result


def verify_completed_artifacts(primary_outputs) -> None:
    """Fail closed on incomplete, stale or modified workflow output sets.

    Snakemake schedules by declared dependencies, while this final gate checks
    source/input/output content hashes, including optional stage companions.
    A failed gate requires rebuilding the named affected stage, not accepting
    its exit code or an existing filename as completion evidence.
    """
    from ..core.io import source_digest
    from magestic.utils.provenance import md5_file
    current_source = source_digest()
    digests = {}
    # Stage A's receipt records the BBTools class/JNI tree and Java it ran with. A
    # warm workflow never re-enters Stage A, so the gate compares that identity with
    # the installation on PATH now. Computed once per call, not per sample.
    tool_identity = None
    for primary in primary_outputs:
        primary = Path(primary)
        proof_path = Path(str(primary)+'.provenance.json')
        try:
            proof = json.loads(proof_path.read_text())
            if (proof['pipeline'] != 'redi' or Path(proof['output']).resolve() != primary.resolve()
                    or proof['code']['redi_source_sha256'] != current_source):
                raise ValueError('source or output identity mismatch')
            if proof.get('step') == 'a_trim_counts':
                if tool_identity is None:
                    from .step_a_trim_counts import _tool_identity
                    tool_identity = _tool_identity()
                if proof['parameters'].get('tool_implementation') != tool_identity:
                    raise ValueError('Stage A tool implementation (BBTools classes/JNI or Java) '
                                     'differs from the installation on PATH; rebuild Stage A')
            outputs = proof['outputs']
            paths = [Path(item['path']).resolve() for item in outputs]
            if not outputs or primary.resolve() not in paths or len(paths) != len(set(paths)):
                raise ValueError('missing or duplicated output attestation')
            for item in [*proof['inputs'],*outputs]:
                path = Path(item['path']).resolve()
                if path not in digests:digests[path] = md5_file(path)
                if item['md5'] != digests[path]:
                    raise ValueError(f'content changed: {path}')
        except (OSError,ValueError,KeyError,TypeError) as exc:
            raise ValueError(f'incomplete or stale REDI stage attestation {proof_path}: {exc}') from exc


def run(config: REDIConfig, from_step: str = "b", to_step: str = "h") -> dict[str, object]:
    """Run the pipeline from `from_step` through `to_step`. Returns dict of
    artifacts produced."""
    if from_step not in STEPS or to_step not in STEPS:
        raise ValueError(f"steps must be one of {STEPS}")
    start = STEPS.index(from_step)
    end = STEPS.index(to_step) + 1
    if start >= end:
        raise ValueError("from_step must be before or equal to to_step; reversed pipeline range")

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
