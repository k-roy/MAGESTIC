"""Compare explicitly selected assays at registered physical positions."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import re
import pandas as pd
from ..core.config import REDIConfig
from ..core.cross_array import compare_arrays
from ..core.identity import keyed_frame
from ..core.io import provenance,publish_tables
from ..readers.keyfiles import _read_key


def _bool_column(values,label):
    text=values.astype('string').str.strip().str.casefold()
    if not text.isin(['true','false','1','0','']).all():raise ValueError(f'invalid boolean metadata in {label}')
    return text.map({'true':True,'1':True,'false':False,'0':False}).astype('boolean')


def run_step_f(config: REDIConfig,annotated_path: Path | None=None) -> list[Path]:
    if not config.array_pairs:raise ValueError('Stage F requires an explicit array pair')
    combined=config.paths.combined_data
    annotated_path=annotated_path or combined/'annotated_top_REDI_bc_MAGESTIC_bc1_df.tsv'
    base=provenance('f_cross_array',combined/'cross_array_run.tsv')
    base.add_input(annotated_path,role='independent_assay_observations')
    df=_read_key(annotated_path)
    if 'purpose' not in df:raise ValueError('cross-array observations require assay purpose')
    names=[f'{p.name_a}_{p.name_b}' for p in config.array_pairs]
    if len(names)!=len(set(names)):raise ValueError('duplicate array-pair output names')
    geometry_path=config.paths.keyfiles/'1536_PIXL_rearray_format.tsv'
    geometry=None
    if geometry_path.is_file():
        base.add_input(geometry_path,role='physical_plate_geometry')
        geometry=keyed_frame(_read_key(geometry_path),['plate_384','row_384','col_384'],label='1536/384 geometry')
        geometry=keyed_frame(geometry,['row_1536','col_1536'],label='1536 geometry')
        geometry=geometry[['plate_384','row_384','col_384','row_1536','col_1536']]
    pending=[]
    for pair in config.array_pairs:
        if any(not re.fullmatch('[A-Za-z0-9_.-]+',name) for name in (pair.name_a,pair.name_b)):
            raise ValueError('array names must be safe nonempty output identifiers')
        selected=pair.purpose_b
        if selected is None:
            if pair.purpose_b_not:raise ValueError('purpose_b_not is ambiguous; declare purpose_b explicitly')
            alternatives=sorted(set(df.purpose)-{pair.purpose_a})
            if len(alternatives)!=1:raise ValueError('ambiguous second assay purpose; declare purpose_b explicitly')
            selected=alternatives
        elif isinstance(selected,str):selected=[selected]
        if not isinstance(selected,list) or not selected or any(not isinstance(x,str) or not x for x in selected):
            raise ValueError('purpose_b must name one or more explicit purposes')
        if len(selected)!=len(set(selected)) or pair.purpose_a in selected:
            raise ValueError('array selectors overlap or repeat the same purpose')
        unknown=set([pair.purpose_a,*selected])-set(df.purpose)
        if unknown:raise ValueError(f'array selector has missing purposes: {sorted(unknown)}')
        mask_a=df.purpose.eq(pair.purpose_a);mask_b=df.purpose.isin(selected)
        selected_a_before=int(mask_a.sum())
        if type(pair.require_validation_a) is not bool:raise ValueError('require_validation_a must be boolean')
        if pair.require_validation_a:
            if 'validated_by_yT100_mating' not in df:
                raise ValueError('required assay validation is missing validated_by_yT100_mating; '
                                 'set require_validation_a=False explicitly for unvalidated comparisons')
            mask_a=mask_a & _bool_column(df.validated_by_yT100_mating,'validated_by_yT100_mating').fillna(False)
            if 'confirmed_contaminated_by_yT100_mating' in df:
                mask_a=mask_a & _bool_column(df.confirmed_contaminated_by_yT100_mating,'confirmed_contaminated_by_yT100_mating').eq(False).fillna(False)
        out=combined/f'{pair.name_a}_{pair.name_b}_merged_REDI_df.tsv'
        prov=deepcopy(base);prov.output=str(out)
        merged=compare_arrays(df,mask_a,mask_b,pair.name_a,pair.name_b,
                              min_counts=config.purity.min_counts,bc1_length=config.barcode.bc1_length)
        if geometry is not None:
            merged=merged.assign(plate_384=merged.REDI_plate_color_384,row_384=merged.REDI_row_384,col_384=merged.REDI_col_384)
            merged=merged.merge(geometry,on=['plate_384','row_384','col_384'],how='left',validate='many_to_one',indicator=True)
            if merged['_merge'].ne('both').any():raise ValueError('compared physical positions lack geometry registration')
            merged=merged.drop(columns='_merge')
        prov.set_parameters(purpose_a=pair.purpose_a,purpose_b=selected,
            require_validation_a=pair.require_validation_a,min_counts=config.purity.min_counts,
            bc1_length=config.barcode.bc1_length,interpretation='this is not an accuracy result')
        prov.set_counts(input_rows=len(df),selected_a_before_validation=selected_a_before,
            selected_a=int(mask_a.sum()),selected_b=int(mask_b.sum()),
            rows_not_selected=int((~(mask_a|mask_b)).sum()),physical_positions=len(merged),
            jointly_callable=int(merged.bc1_agreement.notna().sum()),agree=int(merged.bc1_agreement.fillna(False).sum()))
        pending.append((out,merged,prov))
    for out,merged,prov in pending:publish_tables([(out,merged)],prov)
    return [out for out,_,_ in pending]
