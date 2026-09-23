"""Attach reference candidates without multiplying physical observations."""
from __future__ import annotations
import json
import numpy as np
import pandas as pd


def annotate_with_bc1_reference(top_clone_df: pd.DataFrame,bc1_reference_df: pd.DataFrame,
                                redi_prefix: str='REDI_bc_bc1_',
                                aa_position_col: str='amino_acid_position') -> pd.DataFrame:
    """One output row per input observation; ambiguous designs stay explicit.

    Multiple distinct reference rows for a barcode are retained as JSON candidate
    records. Their scalar annotations remain missing rather than choosing a row.
    Count/denominator names are stable regardless of the reference table schema.
    """
    top=top_clone_df.copy();ref=bc1_reference_df.copy()
    for label,frame in (('REDI',top),('reference',ref)):
        if 'bc1' not in frame or frame.columns.duplicated().any():
            raise ValueError(f'{label}: bc1 required and column names must be unique')
    reserved={'bc1_reference_match_count','bc1_reference_ambiguous','bc1_reference_candidates_json'}
    if reserved & (set(top)|set(ref)):
        raise ValueError('reserved reference annotation columns already present')
    if any(c.startswith(redi_prefix) for c in ref):
        raise ValueError('reference has reserved REDI observation column collision')
    observations={'sample_name','REDI_bc','REDI_library','REDI_plate_color_384',
        'REDI_row_384','REDI_col_384','MAT_alpha_array_1536_barcode',
        'MAT_a_array_1536_barcode','REDI_plate_barcode','gDNA_plate_name',
        'registration_status','purpose','outer_primers','inner_primers','sample_number',
        'validated_by_yT100_mating','confirmed_contaminated_by_yT100_mating',
        'all_bc1_counts_json','observed_bc1_read_fraction','identity_ambiguous',
        'colony_purity','purity_basis','top_bc1','edit_distance','top_barcode_count',
        'near_top_barcode_counts','legacy_filtered_top_n_total_counts','legacy_filtered_top_n_purity'}
    if observations & set(ref):
        raise ValueError(f'reference has reserved observation column collision: {sorted(observations & set(ref))}')
    if aa_position_col in ref:
        original=ref[aa_position_col].replace('',pd.NA)
        numeric=pd.to_numeric(original,errors='coerce')
        if ((original.notna() & numeric.isna()).any() or
            (numeric.dropna()%1!=0).any() or not np.isfinite(numeric.dropna()).all()):
            raise ValueError(f'{aa_position_col} must be an integer or missing')
        ref[aa_position_col]=numeric.astype('Int64')
    if ref['bc1'].isna().any() or ref['bc1'].astype('string').str.strip().eq('').any():
        raise ValueError('reference contains missing barcode key')
    # Identical rows represent redundant annotation, not independent observations.
    ref=ref.drop_duplicates()
    records=[]
    annotation=[c for c in ref if c!='bc1']
    for barcode,group in ref.groupby('bc1',sort=False):
        ambiguous=len(group)>1
        record={'bc1':barcode,**{c:pd.NA if ambiguous else group.iloc[0][c] for c in annotation}}
        candidates=json.loads(group.to_json(orient='records'))
        record.update(bc1_reference_match_count=len(group),bc1_reference_ambiguous=ambiguous,
                      bc1_reference_candidates_json=json.dumps(candidates,sort_keys=True,separators=(',',':')))
        records.append(record)
    collapsed=pd.DataFrame(records,columns=['bc1',*annotation,*sorted(reserved)])
    overlaps=(set(top)&set(ref))-{'bc1','counts','total_counts'}
    rename={c:redi_prefix+c for c in overlaps}
    for c in ('counts','total_counts'):
        if c in top:rename[c]=redi_prefix+c
    if (set(rename.values()) & (set(top)-set(rename))) or set(rename.values()) & set(collapsed):
        raise ValueError('reserved REDI column collision after annotation rename')
    merged=top.rename(columns=rename).merge(collapsed,on='bc1',how='left',validate='many_to_one')
    merged['bc1_reference_match_count']=merged['bc1_reference_match_count'].fillna(0).astype(int)
    merged['bc1_reference_ambiguous']=merged['bc1_reference_ambiguous'].fillna(False).astype(bool)
    if aa_position_col in merged:merged[aa_position_col]=merged[aa_position_col].astype('Int64')
    return merged
