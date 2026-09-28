"""Stage I worklist writer: explicit picks -> PIXL run files, master transfer table, picking summary.

    magestic-redi-worklist --picks picks.tsv --target-plates 5169-5186 --outdir DIR --prefix NAME

The pick table needs source_plate, source_row, source_col (and source_density when not 384); every other
column is carried into the master table unchanged. The output directory is a single worklist generation:
it must not exist yet, so a regenerated worklist can never sit beside, or silently replace, the one the
robot executed (the 2022 QTL defect, planning/184_). Publication happens only after every run file is
re-parsed and reproduces the master table exactly; the provenance sidecar of the master table is the
completion evidence.
"""
from __future__ import annotations

from pathlib import Path
import json
import tempfile

import pandas as pd

from ..core.io import provenance, publish_tables
from ..core.worklist import WorklistPlan, assign, run_file_lines, run_table, validate


def run_step_i_worklist(picks_path: str | Path, plan: WorklistPlan, outdir: str | Path, prefix: str) -> Path:
    outdir = Path(outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise FileExistsError(f'{outdir} already holds a worklist generation; write a new generation directory')
    if not prefix or '/' in prefix:
        raise ValueError('prefix must be a nonempty file-name stem')
    picks = pd.read_csv(picks_path, sep='\t', dtype=str, keep_default_na=False)
    transfers = assign(picks, plan)
    master = outdir / f'{prefix}_transfers.tsv'
    prov = provenance('i_worklist', master)
    prov.add_input(picks_path, role='explicit_pick_table')
    summary = run_table(transfers, plan.max_sources_per_run)
    with tempfile.TemporaryDirectory(prefix='redi_worklist_') as tmp:
        staged, run_files = [], {}
        for run in sorted(transfers.run.unique()):
            name = f'{prefix}_batch_1_{run}.csv'
            path = Path(tmp) / name
            lines = run_file_lines(transfers, run, plan.max_sources_per_run)
            path.write_text('\n'.join(lines) + '\n')
            run_files[run] = path
            staged.append((path, outdir / name, len(lines) - 1 - plan.max_sources_per_run))
        report = validate(transfers, run_files, plan.max_sources_per_run)
        if not report['ok']:
            raise AssertionError('worklist failed its round-trip check: ' + '; '.join(report['problems'][:5]))
        needed = Path(tmp) / f'{prefix}_plates_needed.txt'
        lines = ['targets: ' + ' '.join(plan.target_plates[:transfers.target_plate.nunique()]),
                 'sources, in order of first use: ' + ' '.join(dict.fromkeys(transfers.source_plate))]
        if plan.source_rename:
            lines.append('restamped sources: ' + ' '.join(f'{a}->{b}' for a, b in plan.source_rename.items()))
        needed.write_text('\n'.join(lines) + '\n')
        staged.append((needed, outdir / needed.name, None))
        prov.set_parameters(target_plates=plan.target_plates, target_density=plan.target_density,
                            fill_order=plan.fill_order, excluded_target_wells=sorted(f'{r}{c}' for r, c in plan.exclude_target_wells),
                            first_target_well_index=plan.first_target_well_index,
                            max_sources_per_run=plan.max_sources_per_run, source_rename=plan.source_rename,
                            keep_input_order=plan.keep_input_order, prefix=prefix)
        prov.set_counts(picks=len(transfers), runs=int(transfers.run.nunique()),
                        target_plates=int(transfers.target_plate.nunique()),
                        source_plates=int(transfers.source_plate.nunique()))
        prov.note(json.dumps({'round_trip_validation': report}, sort_keys=True))
        publish_tables([(master, transfers), (outdir / f'{prefix}_picking_summary.tsv', summary)], prov,
                       files=staged)
    return master


def validate_generation(outdir: str | Path) -> dict:
    """Re-check a published generation from disk: run files against the master table."""
    outdir = Path(outdir)
    masters = sorted(outdir.glob('*_transfers.tsv'))
    if len(masters) != 1:
        raise ValueError(f'{outdir}: expected one *_transfers.tsv, found {len(masters)}')
    transfers = pd.read_csv(masters[0], sep='\t', dtype=str, keep_default_na=False)
    prefix = masters[0].name[:-len('_transfers.tsv')]
    run_files = {int(p.stem.rsplit('_', 1)[1]): p for p in outdir.glob(f'{prefix}_batch_1_*.csv')}
    summary = pd.read_csv(outdir / f'{prefix}_picking_summary.tsv', sep='\t', dtype=str)
    max_sources = sum(c.endswith('_source_plate') and not c.startswith('total_') for c in summary.columns)
    transfers['run'] = transfers.run.astype(int)
    return validate(transfers, run_files, max_sources)
