"""PIXL rearray worklists: deterministic source-to-target assignment, legacy CSV writing, and a
round-trip validator.

This fills the Stage I robot-output gap (docs: REDI/docs/pipeline.md §I). It does not decide which clones
to pick; the caller supplies an explicit pick table. It decides where each pick goes, splits the transfers
into PIXL runs, writes the run files in the format the lab's PIXL accepted in 2024-2025, and proves the
written files reproduce the intended transfers exactly.

Why the validator exists: the 2022 QTL collection annotation is wrong for plates 5033-5054 because a
worklist folder was regenerated mid-campaign and the annotation read back a different generation than the
one executed (planning/184_, planning/220_ §5.6). Here the master transfer table and the run files are
written together, a generation directory is never rewritten, and every run file is re-parsed and compared
with the master table before publication.

Run-file format (one file per PIXL run; legacy 20240225 SaCas9-X1 and 20250620 NNS rearrays):

    5102,SBS,384,Target            one target plate
    4860,SBS,1536,Source           first source plate and its density
    4861,SBS,1536,Source           second source plate, or 0000,SBS,384,Source when the nest is empty
    4860,A,12,5102,A,1             source,row,col,target,row,col  in target fill order
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import csv
import re

import pandas as pd

# Row labels are lists, never strings: `'AB' in 'ABC...'` is a substring match.
ROWS_384 = list('ABCDEFGHIJKLMNOP')
ROWS_1536 = [*'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'AA', 'AB', 'AC', 'AD', 'AE', 'AF']
GEOMETRY = {96: (list('ABCDEFGH'), 12), 384: (ROWS_384, 24), 1536: (ROWS_1536, 48)}
EMPTY_NEST = '0000'
_PLATE = re.compile(r'[A-Za-z0-9_-]+')


def wells(density: int, order: str = 'row') -> list[tuple[str, int]]:
    """All wells of a plate in fill order: 'row' = A1, A2 ... A24, B1; 'column' = A1, B1 ... P1, A2."""
    if density not in GEOMETRY:
        raise ValueError(f'unsupported plate density {density}')
    rows, ncol = GEOMETRY[density]
    if order == 'row':
        return [(r, c) for r in rows for c in range(1, ncol + 1)]
    if order == 'column':
        return [(r, c) for c in range(1, ncol + 1) for r in rows]
    raise ValueError("fill order must be 'row' or 'column'")


def edge_wells(density: int) -> set[tuple[str, int]]:
    rows, ncol = GEOMETRY[density]
    return {(r, c) for r in rows for c in range(1, ncol + 1)
            if r in (rows[0], rows[-1]) or c in (1, ncol)}


def _natural(plate: str):
    return (0, int(plate), '') if plate.isdigit() else (1, 0, plate)


@dataclass
class WorklistPlan:
    """Explicit parameters; nothing is inferred from the data."""

    target_plates: list[str]                     # target barcodes in fill order (pre-printed labels)
    target_density: int = 384
    fill_order: str = 'row'
    exclude_target_wells: set[tuple[str, int]] = field(default_factory=set)
    first_target_well_index: int = 0             # resume a partly filled first target plate
    max_sources_per_run: int = 2                 # the PIXL has two source nests
    source_rename: dict[str, str] = field(default_factory=dict)   # old source barcode -> restamped barcode
    keep_input_order: bool = False               # default: sort picks by source plate, row, column


def _validate_picks(picks: pd.DataFrame) -> pd.DataFrame:
    need = {'source_plate', 'source_row', 'source_col'}
    missing = need - set(picks)
    if missing:
        raise ValueError(f'pick table lacks {sorted(missing)}')
    p = picks.copy()
    p['source_plate'] = p.source_plate.astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
    p['source_row'] = p.source_row.astype(str).str.strip().str.upper()
    p['source_col'] = p.source_col.astype(str).str.strip().str.replace(r'\.0$', '', regex=True)
    p['source_density'] = (p.source_density.astype(str).str.strip() if 'source_density' in p else '384')
    p.loc[p.source_density.eq(''), 'source_density'] = '384'
    if not p.source_plate.str.fullmatch(_PLATE).all() or p.source_plate.eq(EMPTY_NEST).any():
        raise ValueError('source plate barcodes must be nonempty identifiers other than 0000')
    for dens, rows in p.groupby('source_density'):
        if not dens.isdigit() or int(dens) not in GEOMETRY:
            raise ValueError(f'unsupported source density {dens!r}')
        valid_rows, ncol = GEOMETRY[int(dens)]
        bad = ~rows.source_row.isin(valid_rows) | ~rows.source_col.str.fullmatch(r'\d+') | \
              ~rows.source_col.where(rows.source_col.str.fullmatch(r'\d+'), '0').astype(int).between(1, ncol)
        if bad.any():
            raise ValueError(f'invalid {dens}-well source coordinates, e.g. {rows[bad].iloc[0][["source_plate", "source_row", "source_col"]].tolist()}')
    if p.groupby('source_plate').source_density.nunique().gt(1).any():
        raise ValueError('a source plate cannot have two densities')
    dup = p.duplicated(['source_plate', 'source_row', 'source_col'], keep=False)
    if dup.any():
        raise ValueError(f'{int(dup.sum())} pick rows name the same source well more than once')
    return p


def assign(picks: pd.DataFrame, plan: WorklistPlan) -> pd.DataFrame:
    """Return the master transfer table: one row per pick with its target well and PIXL run."""
    if plan.max_sources_per_run < 1:
        raise ValueError('max_sources_per_run must be >= 1')
    if len(set(plan.target_plates)) != len(plan.target_plates) or not plan.target_plates:
        raise ValueError('target plate barcodes must be nonempty and unique')
    p = _validate_picks(picks)
    p['original_source_plate'] = p.source_plate
    p['source_plate'] = p.source_plate.map(lambda s: plan.source_rename.get(s, s))
    if p.duplicated(['source_plate', 'source_row', 'source_col']).any():
        raise ValueError('source renaming merges two picks into one well')
    if set(p.source_plate) & set(plan.target_plates):
        raise ValueError('a plate cannot be both a source and a target')
    if not plan.keep_input_order:
        rows_order = {r: i for i, r in enumerate(ROWS_1536)}
        p = p.assign(_k1=p.source_plate.map(_natural), _k2=p.source_row.map(rows_order),
                     _k3=p.source_col.astype(int)).sort_values(['_k1', '_k2', '_k3'], kind='stable')
        p = p.drop(columns=['_k1', '_k2', '_k3'])
    slots = [(t, r, c) for t in plan.target_plates for (r, c) in wells(plan.target_density, plan.fill_order)
             if (r, c) not in plan.exclude_target_wells]
    slots = slots[plan.first_target_well_index:]
    if len(p) > len(slots):
        raise ValueError(f'{len(p)} picks exceed {len(slots)} available target wells; add target plates')
    p = p.reset_index(drop=True)
    p['target_plate'] = [s[0] for s in slots[:len(p)]]
    p['target_row'] = [s[1] for s in slots[:len(p)]]
    p['target_col'] = [s[2] for s in slots[:len(p)]]
    p['target_density'] = plan.target_density
    run_ids, n = [], 0
    for target, block in p.groupby('target_plate', sort=False):
        sources = list(dict.fromkeys(block.source_plate))
        groups = [sources[i:i + plan.max_sources_per_run] for i in range(0, len(sources), plan.max_sources_per_run)]
        which = {}
        for g in groups:
            n += 1
            for s in g:
                which[s] = n
        run_ids += [which[s] for s in block.source_plate]
    p['run'] = run_ids
    return p


def run_table(transfers: pd.DataFrame, max_sources: int) -> pd.DataFrame:
    """PIXL picking summary, legacy column names (processed_data/rearray/*_picking_summary.tsv)."""
    rows = []
    for run, g in transfers.groupby('run'):
        sources = list(dict.fromkeys(g.source_plate))
        counts = g.source_plate.value_counts()
        rec = {'batch_number': 1, 'rearray_file_number_suffix': run, 'Target_plate_ID': g.target_plate.iloc[0]}
        names = ['first', 'second', 'third', 'fourth'][:max_sources]
        plates = [sources[i] if i < len(sources) else EMPTY_NEST for i in range(max_sources)]
        rec.update({f'{n}_source_plate': s for n, s in zip(names, plates)})         # legacy column order:
        rec.update({f'total_colonies_from_{n}_source_plate': int(counts.get(s, 0))  # plates, then counts
                    for n, s in zip(names, plates)})
        rows.append(rec)
    return pd.DataFrame(rows)


def run_file_lines(transfers: pd.DataFrame, run: int, max_sources: int) -> list[str]:
    g = transfers[transfers.run.eq(run)]
    targets = g.target_plate.unique()
    if len(targets) != 1:
        raise ValueError(f'run {run} spans {len(targets)} target plates')
    lines = [f'{targets[0]},SBS,{int(g.target_density.iloc[0])},Target']
    sources = list(dict.fromkeys(g.source_plate))
    density = dict(zip(g.source_plate, g.source_density))
    for i in range(max_sources):
        lines.append(f'{sources[i]},SBS,{density[sources[i]]},Source' if i < len(sources)
                     else f'{EMPTY_NEST},SBS,384,Source')
    for r in g.itertuples():
        lines.append(f'{r.source_plate},{r.source_row},{r.source_col},{r.target_plate},{r.target_row},{r.target_col}')
    return lines


def parse_run_file(path: str | Path) -> tuple[list[tuple[str, str, str]], pd.DataFrame]:
    """Read a PIXL run file back: header (barcode, density, role) and body transfers."""
    header, body = [], []
    with open(path, newline='') as fh:
        for n, row in enumerate(csv.reader(fh), 1):
            if len(row) == 4 and row[2] and row[3] in ('Target', 'Source'):
                if body:
                    raise ValueError(f'{path}:{n}: header line after transfer lines')
                header.append((row[0], row[2], row[3]))
            elif len(row) == 6:
                body.append(row)
            else:
                raise ValueError(f'{path}:{n}: not a header or transfer line')
    frame = pd.DataFrame(body, columns=['source_plate', 'source_row', 'source_col',
                                        'target_plate', 'target_row', 'target_col'])
    return header, frame


def validate(transfers: pd.DataFrame, run_files: dict[int, Path], max_sources: int) -> dict:
    """Re-parse every run file and require it to reproduce the master table exactly."""
    problems, seen = [], []
    for run, path in sorted(run_files.items()):
        header, body = parse_run_file(path)
        targets = [h for h in header if h[2] == 'Target']
        sources = [h for h in header if h[2] == 'Source']
        if len(targets) != 1:
            problems.append(f'run {run}: {len(targets)} target header lines')
        if len(sources) != max_sources:
            problems.append(f'run {run}: {len(sources)} source header lines, expected {max_sources}')
        declared = {h[0] for h in sources} - {EMPTY_NEST}
        used = set(body.source_plate)
        if used != declared:
            problems.append(f'run {run}: body sources {sorted(used)} differ from header {sorted(declared)}')
        if targets and set(body.target_plate) != {targets[0][0]}:
            problems.append(f'run {run}: body targets differ from header')
        density = {h[0]: int(h[1]) for h in header}
        for r in body.itertuples():
            for plate, row, col in ((r.source_plate, r.source_row, r.source_col),
                                    (r.target_plate, r.target_row, r.target_col)):
                if density.get(plate) not in GEOMETRY:
                    problems.append(f'run {run}: plate {plate} has no header line with a valid density')
                    continue
                rows, ncol = GEOMETRY[density[plate]]
                if row not in rows or not col.isdigit() or not 1 <= int(col) <= ncol:
                    problems.append(f'run {run}: {plate} {row}{col} outside a {density[plate]}-well plate')
        seen.append(body.assign(run=run))
    written = pd.concat(seen, ignore_index=True) if seen else pd.DataFrame()
    key = ['source_plate', 'source_row', 'source_col', 'target_plate', 'target_row', 'target_col']
    master = transfers[key + ['run']].astype(str)
    written = written[key + ['run']].astype(str) if len(written) else written
    if written.duplicated(['target_plate', 'target_row', 'target_col']).any():
        problems.append('a target well is written twice')
    if written.duplicated(['source_plate', 'source_row', 'source_col']).any():
        problems.append('a source well is picked twice')
    a = set(map(tuple, master.values)); b = set(map(tuple, written.values))
    if a != b:
        problems.append(f'run files and master table differ: {len(a - b)} master-only, {len(b - a)} file-only')
    return {'ok': not problems, 'problems': problems, 'transfers_in_files': len(written),
            'transfers_in_master': len(master), 'runs': len(run_files)}
