"""Plot measured registered-position counts without inventing clone independence.

The legacy filenames retain ``clones`` for compatibility. Quantities represent
registered array positions, which may be repeat quadrants of the same isolate.
An independent-replicate gate is applied only when requested and documented in
an explicit input column; matching barcodes or positions do not supply that gate.
"""
from __future__ import annotations

from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from ..core.annotate import annotate_with_bc1_reference
from ..core.config import REDIConfig
from ..core.identity import POSITION_COLUMNS, count_integer, keyed_frame, valid_barcode
from ..core.io import provenance, publish_tables
from ..readers.keyfiles import _read_key


def _numeric(frame, column):
    if column not in frame:
        raise ValueError(f"plot filtering requires {column}")
    raw = frame[column]
    missing = raw.isna() | raw.astype('string').str.strip().eq('').fillna(False)
    values = pd.to_numeric(raw.where(~missing), errors='coerce')
    if values[~missing].isna().any() or not np.isfinite(values[~missing].astype(float)).all():
        raise ValueError(f"{column}: expected finite numeric values or explicit missing data")
    return values


def _counts(frame, column):
    if column not in frame:
        raise ValueError(f"plot filtering requires {column}")
    values = []
    for value in frame[column]:
        values.append(np.nan if pd.isna(value) or str(value).strip() == '' else count_integer(value, column))
    return pd.Series(values, index=frame.index)


def filter_high_quality_clones(annotated_cross_array_df: pd.DataFrame, array_a_name: str,
                               array_b_name: str, min_counts: int = 20,
                               min_purity: float = 0.95, min_replicates: int = 0,
                               bc1_length: int = 20) -> pd.DataFrame:
    """Filter by measured barcode counts/read fractions and an explicit gate.

    ``min_replicates=0`` disables the independent-replicate gate. A positive gate
    requires ``num_independent_replicates``; absence never means it passed.
    Array names are literal column labels, not expressions evaluated by query().
    """
    if type(min_counts) is not int or min_counts < 0:
        raise ValueError('min_counts must be a nonnegative integer')
    if not isinstance(min_purity, (int, float)) or isinstance(min_purity, bool) or not 0 <= min_purity <= 1:
        raise ValueError('min_purity must be a finite read fraction in [0,1]')
    if type(min_replicates) is not int or min_replicates < 0:
        raise ValueError('min_replicates must be a nonnegative integer')
    df = annotated_cross_array_df
    passes = []
    for name in (array_a_name, array_b_name):
        counts = _counts(df, name+'_REDI_bc_bc1_counts')
        fractions = _numeric(df, name+'_colony_purity')
        if not fractions.dropna().between(0, 1).all():
            raise ValueError('observed barcode read fractions must be in [0,1]')
        observation = name+'_bc1'
        if observation not in df:
            raise ValueError(f'plot filtering requires barcode observation {observation}')
        passes.append((counts >= min_counts) & (fractions >= min_purity)
                      & valid_barcode(df[observation],bc1_length))
    mask = passes[0] | passes[1]
    if min_replicates:
        if 'num_independent_replicates' not in df:
            raise ValueError('requested independent-replicate gate requires num_independent_replicates')
        mask &= _counts(df, 'num_independent_replicates') >= min_replicates
    return df.loc[mask].copy()


def _observed_identity(df, array_a_name, array_b_name, *, bc1_length=20):
    if df.columns.duplicated().any():
        raise ValueError('plot input has duplicate column names')
    out = df.copy()
    a, b = array_a_name+'_bc1', array_b_name+'_bc1'
    if a not in out or b not in out:
        raise ValueError('both array barcode observation columns are required')
    valid_a, valid_b = valid_barcode(out[a], bc1_length), valid_barcode(out[b], bc1_length)
    conflict = valid_a & valid_b & out[a].astype('string').ne(out[b].astype('string')).fillna(False)
    chosen = out[b].astype('string').where(valid_b, out[a].astype('string'))
    out['bc1'] = chosen.where((valid_a | valid_b) & ~conflict, pd.NA)
    out['plot_identity_conflict'] = conflict
    out['mean_purity'] = pd.concat([_numeric(out, name+'_colony_purity')
                                    for name in (array_a_name, array_b_name)], axis=1).mean(axis=1)
    return out


def unique_bc1_by_mean_purity(df: pd.DataFrame, array_a_name: str,
                               array_b_name: str, *, bc1_length=20) -> pd.DataFrame:
    """Rank unique callable identities; conflicting arrays remain unassessed.

    This is a unique-barcode utility, not a physical-position count. StageH does
    not use this deduplication before counting registered positions.
    """
    out = _observed_identity(df, array_a_name, array_b_name, bc1_length=bc1_length)
    return out.loc[out.bc1.notna()].sort_values('mean_purity', ascending=False, kind='stable').drop_duplicates('bc1')


def _plotting_implementation(matplotlib, seaborn) -> dict[str, object]:
    """Versions that determine the PNG bytes, recorded (not enforced) in H receipts.

    The H tables do not depend on plotting and the PNGs are hash-attested, so an
    upgrade does not invalidate; the receipt only says what drew the image. Backend
    and FreeType are read defensively because test stubs provide neither.
    """
    import sys
    backend = getattr(matplotlib, 'get_backend', None)
    ft2font = sys.modules.get('matplotlib.ft2font')
    return {'matplotlib': getattr(matplotlib, '__version__', None),
            'seaborn': getattr(seaborn, '__version__', None),
            'backend': backend() if callable(backend) else None,
            'freetype': getattr(ft2font, '__freetype_version__', None)}


def run_step_h(config: REDIConfig, annotated_cross_array_path: Path | None = None, *,
               min_counts: int | None = None, min_purity: float | None = None,
               min_replicates: int | None = None) -> dict[str, Path]:
    """Plot every configured pair, with counts and completion provenance.

    Standard configuration has no independent-replicate gate. Per-gene ORF and
    project-specific selection remain screen-driver responsibilities.
    """
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError as exc:
        raise ImportError('step_h requires matplotlib and seaborn') from exc
    if not config.array_pairs:
        raise ValueError('step_h needs at least one configured array pair')
    if annotated_cross_array_path is not None and len(config.array_pairs) != 1:
        raise ValueError('one annotated_cross_array_path requires exactly one array pair')
    pair_names = [(p.name_a,p.name_b) for p in config.array_pairs]
    if len(set(pair_names)) != len(pair_names):
        raise ValueError('duplicate array pairs')
    if len({a+'_'+b for a,b in pair_names}) != len(pair_names):
        raise ValueError('array-pair output filename collision')
    for pair in pair_names:
        if any(not n or '/' in n or '\\' in n for n in pair):
            raise ValueError('array names must be nonempty filename components')
    settings = config.plotting
    min_counts = settings.min_counts if min_counts is None else min_counts
    min_purity = settings.min_purity if min_purity is None else min_purity
    min_replicates = settings.min_replicates if min_replicates is None else min_replicates
    plots = config.paths.plots
    plots.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for pair in config.array_pairs:
        tag = f'{pair.name_a}_{pair.name_b}_' if len(config.array_pairs) > 1 else ''
        source = annotated_cross_array_path or config.paths.combined_data/f'{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv'
        table_path = plots/(tag+'bc1_registered_position_counts.tsv')
        prov = provenance('h_plot', table_path)
        prov.add_input(source, role='registered_cross_array_observations')
        df = keyed_frame(_read_key(source),
                         POSITION_COLUMNS, label='plot registered positions')
        df = _observed_identity(df, pair.name_a, pair.name_b, bc1_length=config.barcode.bc1_length)
        if config.bc1_reference_table is not None:
            prov.add_input(config.bc1_reference_table, role='design_annotation')
            ref = _read_key(config.bc1_reference_table)
            collisions = (set(df) & set(ref)) - {'bc1','counts','total_counts'}
            if collisions:
                raise ValueError(f'reference collides with measured plot columns: {sorted(collisions)}')
            df = annotate_with_bc1_reference(df, ref)
        filtered = filter_high_quality_clones(df, pair.name_a, pair.name_b,
                        min_counts=min_counts, min_purity=min_purity, min_replicates=min_replicates,
                        bc1_length=config.barcode.bc1_length)
        filtered = filtered.loc[filtered.bc1.notna()].copy()
        positions = filtered['bc1'].value_counts().rename_axis('bc1').reset_index(name='registered_positions')
        distribution = positions['registered_positions'].value_counts().sort_index()
        dist_path = plots/(tag+'bc1_clones_distribution.png')
        prov.set_parameters(array_pair=[pair.name_a,pair.name_b], min_counts=min_counts,
            min_observed_bc1_read_fraction=min_purity,min_independent_replicates=min_replicates,
            replicate_policy='no independent-replicate gate' if not min_replicates else 'explicit input independent-replicate count',
            counting_unit='registered array positions; repeated quadrants may be the same isolate',
            purity_interpretation='observed barcode read fraction, not cell purity',
            plotting_implementation=_plotting_implementation(matplotlib, sns))
        prov.set_counts(registered_positions=len(df), callable_identities=int(df.bc1.notna().sum()),
            conflicting_identities=int(df.plot_identity_conflict.sum()), passing_positions=len(filtered),
            distinct_passing_barcodes=len(positions))
        tables=[(table_path,positions),(plots/(tag+'filtered_registered_positions.tsv'),filtered)]
        with tempfile.TemporaryDirectory(prefix='redi_plot_',dir=plots) as staging:
            staging=Path(staging);image=staging/'distribution.png'
            plt.figure(figsize=(10,6))
            try:
                if len(distribution):sns.barplot(x=distribution.index,y=distribution.values)
                plt.xlabel('Registered positions per bc1')
                plt.ylabel('Number of barcodes')
                plt.title('Registered-position distribution' if len(distribution) else 'No registered positions pass the declared filters')
                plt.xticks(rotation=45);plt.tight_layout();plt.savefig(image)
            finally:
                plt.close()
            files=[(image,dist_path,None)]
            pair_outputs={tag+'bc1_clones_distribution':dist_path,tag+'registered_position_counts':table_path,
                          tag+'filtered_registered_positions':tables[1][0]}
            if 'sublibrary_ID' in filtered:
                staged_sublib=staging/'sublibrary.png'
                sublib_path=plots/(tag+'sublibrary_ID_counts_barchart.png')
                _plot_sublibrary_counts(df,filtered,staged_sublib)
                files.append((staged_sublib,sublib_path,None))
                pair_outputs[tag+'sublibrary_ID_counts_barchart']=sublib_path
            if set(outputs) & set(pair_outputs):
                raise ValueError('array-pair output filename collision')
            publish_tables(tables,prov,files=files)
            outputs.update(pair_outputs)
    return outputs


def _plot_sublibrary_counts(all_df: pd.DataFrame, filtered_df: pd.DataFrame, out_path: Path) -> None:
    import matplotlib.pyplot as plt
    def counts(frame):
        return frame['sublibrary_ID'].replace('',pd.NA).fillna('UNASSIGNED').value_counts().sort_index()
    data=pd.DataFrame({'registered_positions':counts(all_df),
                       'passing_registered_positions':counts(filtered_df)}).fillna(0).astype(int)
    try:
        if not data.empty:data.plot(kind='bar',color=['#1f77b4','#ff7f0e'])
        else:plt.figure(figsize=(10,6))
        plt.xlabel('sublibrary_ID (ambiguous/missing annotation is UNASSIGNED)')
        plt.ylabel('Registered positions');plt.title('Registered positions by sublibrary')
        plt.tight_layout();plt.legend(title='Filter',bbox_to_anchor=(1.05,1),loc='upper left')
        plt.savefig(out_path,bbox_inches='tight')
    finally:
        plt.close()
