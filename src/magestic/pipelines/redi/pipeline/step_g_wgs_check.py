"""Step G: compare WGS and REDI observations at registered physical positions."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from magestic.pipelines.wgs.readers.wgs_calls import load_wgs_bc1_calls
from ..core.config import REDIConfig
from ..core.identity import keyed_frame
from ..core.io import provenance, publish_tables
from ..core.wgs_check import compare_to_wgs
from ..readers.wgs_outcome import load_wgs_plate_key
from ..readers.keyfiles import _read_key


def run_step_g(config: REDIConfig, cross_array_path: Path | None = None,
               mat_alpha_array_barcode_filter: int | None = None) -> Path:
    """One row per registered sample and pair, over every configured array pair.

    No historical plate is selected by default. Registration must come from
    independent dated physical plate/transfer records. Both WGS emit and legacy
    tables pass through the shared loader, retaining native calls. An optional
    namespace table explicitly maps source_sample to sample_name.
    """
    if config.wgs_outcome_table is None:
        raise ValueError('config.wgs_outcome_table is required for step_g')
    if not config.array_pairs:
        raise ValueError('step_g requires an explicit array pair')
    if cross_array_path is not None and len(config.array_pairs) != 1:
        raise ValueError('one cross_array_path requires exactly one configured pair')
    names=[(p.name_a,p.name_b) for p in config.array_pairs]
    if len(names)!=len(set(names)):
        raise ValueError('duplicate array pairs')
    combined=config.paths.combined_data
    combined.mkdir(parents=True,exist_ok=True)
    out=combined/'WGS_bc1_call_outer_merge_with_REDI_df_and_plate_key.tsv'
    prov=provenance('g_wgs_check',out)
    key_path=config.paths.keyfiles/'WGS_plate_key.tsv'
    fmt_path=config.paths.keyfiles/'384_PIXL_rearray_format.tsv'
    prov.add_input(key_path,role='independent_physical_registration')
    if fmt_path.is_file():prov.add_input(fmt_path,role='explicit_sample_number_format')
    plate=load_wgs_plate_key(key_path,fmt_path if fmt_path.is_file() else None)
    plate=plate.rename(columns={'MAT_alpha_array_plate_384':'REDI_plate_color_384',
                               'row_384':'REDI_row_384','col_384':'REDI_col_384'})
    if 'REDI_plate_color_384' not in plate and 'plate_384' in plate:
        plate=plate.rename(columns={'plate_384':'REDI_plate_color_384'})
    namespace=None
    if config.wgs_sample_namespace is not None:
        prov.add_input(config.wgs_sample_namespace,role='sample_namespace')
        ns=keyed_frame(_read_key(config.wgs_sample_namespace),['source_sample'],label='sample namespace')
        ns=keyed_frame(ns,['sample_name'],label='sample namespace')
        namespace=dict(zip(ns.source_sample,ns.sample_name))
    prov.add_input(config.wgs_outcome_table,role='wgs_observations')
    calls=load_wgs_bc1_calls(config.wgs_outcome_table,sample_namespace=namespace)
    calls=calls.rename(columns={'bc1':'bc1_from_wgs'})
    unregistered=set(calls.sample_name)-set(plate.sample_name)
    if unregistered:
        raise ValueError(f'WGS calls lack physical registration for {sorted(unregistered)}')
    registered_before=len(plate)
    if mat_alpha_array_barcode_filter is not None:
        wanted=str(mat_alpha_array_barcode_filter)
        plate=plate[plate.MAT_alpha_array_1536_barcode.astype(str)==wanted]
        calls=calls[calls.sample_name.isin(plate.sample_name)]
    if plate.empty:
        raise ValueError('no physical WGS positions selected')
    tables=[]
    for pair in config.array_pairs:
        path=cross_array_path or combined/f'{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv'
        prov.add_input(path,role=f'redi_observations:{pair.name_a}:{pair.name_b}')
        cross=_read_key(path)
        result=compare_to_wgs(cross,calls,plate,pair.name_a,pair.name_b,config.barcode.bc1_length)
        result.insert(0,'array_pair',f'{pair.name_a}:{pair.name_b}')
        tables.append(result)
    result=pd.concat(tables,ignore_index=True)
    prov.set_parameters(array_pairs=names,bc1_length=config.barcode.bc1_length,
                        physical_plate_filter=mat_alpha_array_barcode_filter,
                        interpretation='paired observations; accuracy requires independent expected identity')
    counts={'registered_samples_before_filter':registered_before,'registered_samples':len(plate),
            'wgs_observation_rows':len(calls),'array_pairs':len(tables),'rows_out':len(result)}
    for pair,table in zip(config.array_pairs,tables):
        for name in (pair.name_a,pair.name_b):
            prefix=f'{pair.name_a}:{pair.name_b}:{name}'
            counts[prefix+':callable']=int(table[f'wgs_and_{name}_callable'].sum())
            counts[prefix+':agree']=int(table[f'wgs_and_{name}_bc1_agreement'].fillna(False).sum())
    prov.set_counts(**counts)
    publish_tables([(out,result)],prov)
    return out
