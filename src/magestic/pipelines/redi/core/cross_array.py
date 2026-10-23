"""Compare independently selected arrays at explicit physical positions.

Pool complete barcode distributions before applying minimum support. Top-only
legacy observations can be compared singly, but cannot reconstruct a pooled
barcode distribution. Missing/tied/low-support calls remain unassessed.
"""
from __future__ import annotations
from collections import Counter
import json
import pandas as pd
from .identity import POSITION_COLUMNS, keyed_frame, valid_barcode, count_integer

PLATE_GROUP_COLS=['bc1','plate_and_position',*POSITION_COLUMNS]
PLATE_AGG={'REDI_bc_bc1_counts':'sum','REDI_bc_bc1_total_counts':'sum'}


def _attach_plate_and_position(df: pd.DataFrame) -> pd.DataFrame:
    out=keyed_frame(df,POSITION_COLUMNS,label='REDI physical positions',unique=False)
    out['plate_and_position']=out[list(POSITION_COLUMNS)].agg('_'.join,axis=1)
    return out


def _integer(value, label):
    return count_integer(value,label)


def _unique_object(pairs):
    values={}
    for key,value in pairs:
        if key in values:raise ValueError('duplicate barcode')
        values[key]=value
    return values


def _legacy_ambiguity(value) -> bool:
    """Respect an explicit tie flag, including all-string TSV roundtrips."""
    if pd.isna(value) or value == '':
        return False  # Legacy tables can omit the flag; their limit stays labelled.
    text = str(value).strip().casefold()
    if text not in ('true', 'false', '1', '0'):
        raise ValueError('invalid legacy identity_ambiguous flag')
    return text in ('true', '1')


def _collapse_array(df: pd.DataFrame,array_name: str,min_counts: int=10,
                    bc1_length: int=20) -> pd.DataFrame:
    records=[]
    required={'bc1','REDI_bc_bc1_counts','REDI_bc_bc1_total_counts'}
    if not required<=set(df):
        raise ValueError(f'cross-array calls missing {sorted(required-set(df))}')
    if {'sample_name','REDI_bc'}<=set(df):
        keyed_frame(df,['sample_name','REDI_bc'],label=f'{array_name} physical observations')
    for keys,group in df.groupby(list(POSITION_COLUMNS),sort=False):
        counts=Counter();total=0;legacy=None
        for _,row in group.iterrows():
            declared=_integer(row['REDI_bc_bc1_total_counts'],'declared denominator')
            maximum=_integer(row['REDI_bc_bc1_counts'],'top counts')
            if maximum>declared:
                raise ValueError('top counts exceed denominator')
            raw=row.get('all_bc1_counts_json')
            if pd.isna(raw) or raw=='':
                if len(group)!=1:
                    raise ValueError('pooling multiple observations requires all_bc1_counts_json distribution')
                legacy=(row['bc1'],maximum,declared,_legacy_ambiguity(row.get('identity_ambiguous')))
            else:
                try:
                    pairs=json.loads(raw,object_pairs_hook=_unique_object)
                    if not isinstance(pairs,dict) or not pairs:
                        raise ValueError('distribution must be a nonempty JSON object')
                    values={}
                    for barcode,n in pairs.items():
                        if barcode in values or not isinstance(barcode,str):
                            raise ValueError('duplicate or nonstring barcode')
                        if not all(base in 'ACGT' for base in barcode) or not barcode:
                            raise ValueError('invalid barcode')
                        values[barcode]=_integer(n,'distribution')
                    if sum(values.values())!=declared or max(values.values())!=maximum:
                        raise ValueError('distribution does not reconcile with declared denominator/top count')
                except (TypeError,ValueError) as exc:
                    raise ValueError(f'invalid all_bc1_counts_json distribution: {exc}') from exc
                counts.update(values)
            total+=declared
        if legacy is not None:
            winner,maximum,total,tied=legacy
        else:
            maximum=max(counts.values(),default=0)
            winners=[bc for bc,n in counts.items() if n==maximum and n>0]
            tied=len(winners)>1
            winner=winners[0] if len(winners)==1 else pd.NA
        valid=bool(valid_barcode(pd.Series([winner]),bc1_length).iloc[0])
        callable_call=valid and maximum>=min_counts and total>0 and not tied
        record=dict(zip(POSITION_COLUMNS,keys))
        record.update({f'{array_name}_bc1':winner if callable_call else pd.NA,
            f'{array_name}_REDI_bc_bc1_counts':maximum,
            f'{array_name}_REDI_bc_bc1_total_counts':total,
            f'{array_name}_colony_purity':maximum/total if total else float('nan'),
            f'{array_name}_identity_ambiguous':tied,
            f'{array_name}_callable':callable_call,
            f'{array_name}_observations':len(group),
            f'{array_name}_purity_basis':'observed_bc1_read_fraction' if legacy is None else 'legacy_top_only_read_fraction'})
        records.append(record)
    fields=['bc1','REDI_bc_bc1_counts','REDI_bc_bc1_total_counts','colony_purity',
            'identity_ambiguous','callable','observations','purity_basis']
    result=pd.DataFrame(records,columns=[*POSITION_COLUMNS,*[f'{array_name}_{f}' for f in fields]])
    # An outer physical-position join introduces missing values. Keep Python
    # integers so pandas cannot round exact counts via an int64→float64 cast.
    for field in ('REDI_bc_bc1_counts','REDI_bc_bc1_total_counts','observations'):
        column=f'{array_name}_{field}'
        result[column]=pd.Series([record[column] for record in records],dtype=object)
    return result


def compare_arrays(annotated_top_bc1_df: pd.DataFrame,array_a_mask: pd.Series,
                   array_b_mask: pd.Series,array_a_name: str,array_b_name: str,
                   min_counts: int=10,bc1_length: int=20) -> pd.DataFrame:
    """One row per union of physical positions; callable denominator is explicit."""
    if not array_a_name or not array_b_name or array_a_name==array_b_name:
        raise ValueError('two distinct array names are required')
    if type(min_counts) is not int or min_counts<1:
        raise ValueError('min_counts must be a positive integer')
    df=_attach_plate_and_position(annotated_top_bc1_df)
    for mask in (array_a_mask,array_b_mask):
        if (not isinstance(mask,pd.Series) or not mask.index.equals(df.index)
                or mask.isna().any() or not pd.api.types.is_bool_dtype(mask.dtype)):
            raise ValueError('array masks must be nonmissing booleans aligned to input rows')
    if (array_a_mask & array_b_mask).any():
        raise ValueError('array masks overlap; the same observation cannot validate itself')
    if {'sample_name','REDI_bc'}<=set(df):
        # Disjoint row masks are insufficient when one physical observation was
        # copied under two different purpose labels. Validate their union first.
        keyed_frame(df[array_a_mask|array_b_mask],['sample_name','REDI_bc'],
                    label='independent cross-array observations')
    left=_collapse_array(df[array_a_mask],array_a_name,min_counts,bc1_length)
    right=_collapse_array(df[array_b_mask],array_b_name,min_counts,bc1_length)
    merged=left.merge(right,on=list(POSITION_COLUMNS),how='outer',validate='one_to_one')
    merged['plate_and_position']=merged[list(POSITION_COLUMNS)].agg('_'.join,axis=1)
    both=(merged[f'{array_a_name}_callable'].astype('boolean').fillna(False)
          & merged[f'{array_b_name}_callable'].astype('boolean').fillna(False))
    merged[f'{array_a_name}_bc1_and_{array_b_name}_bc1_both_present']=both
    merged['bc1_agreement']=(merged[f'{array_a_name}_bc1'].astype('string')==merged[f'{array_b_name}_bc1'].astype('string')).astype('boolean').where(both,pd.NA)
    return merged
