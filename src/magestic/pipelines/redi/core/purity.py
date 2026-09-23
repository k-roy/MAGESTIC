"""Observed barcode read fractions, with explicit legacy heuristic diagnostics.

Every observed barcode contributes to the primary denominator. Near neighbors are
not declared sequencing errors from distance alone. The old distance-filtered,
top-N measure remains separately named for comparison, never as cell purity.
"""
from __future__ import annotations
from collections.abc import Iterable
import json
import numpy as np
import pandas as pd
from rapidfuzz.distance.Levenshtein import distance as lev_distance
from .identity import keyed_frame, count_integer


def _validated_counts(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    data=keyed_frame(df,group_cols,label='REDI count groups',unique=False)
    if 'bc1' not in data or 'counts' not in data:
        raise ValueError('REDI purity requires bc1 and counts columns')
    if not data['bc1'].astype('string').str.fullmatch('[ACGT]+',na=False).all():
        raise ValueError('REDI purity requires nonmissing ACGT barcode observations')
    counts=data['counts'].map(count_integer)
    if (counts>np.iinfo(np.int64).max).any():
        raise ValueError('counts must be finite nonnegative integers')
    data['counts']=counts.astype('int64')
    return data


def top_bc1_per_redi_bc(df: pd.DataFrame, group_cols: Iterable[str]) -> pd.DataFrame:
    """Aggregate identical barcodes first; a tied identity remains uncalled."""
    keys=list(group_cols)
    result=compute_purity(df,group_cols=keys)
    return result[keys+['top_bc1']]


def _edit_distance_vec(top: np.ndarray, observed: np.ndarray) -> list[float]:
    return [np.nan if pd.isna(t) or pd.isna(o) else lev_distance(str(t),str(o))
            for t,o in zip(top,observed)]


def filter_edit_distance(df: pd.DataFrame, drop_distances: Iterable[int]=(1,2)) -> pd.DataFrame:
    """Explicit legacy diagnostic filter; not the primary read denominator."""
    return df[~df['edit_distance'].isin(set(drop_distances))].copy()


def compute_purity(df: pd.DataFrame, group_cols: Iterable[str]=('sample_name','REDI_bc'),
                   top_n: int=4, edit_distance_drop: Iterable[int]=(1,2)) -> pd.DataFrame:
    """One row per registered sample/REDI barcode with all-read accounting.

    ``colony_purity`` is retained as a compatibility column and equals the
    largest observed bc1 read fraction; it is not an estimate of cell fraction.
    ``counts`` is the maximum aggregated barcode count, ``total_counts`` counts
    every observed barcode, and tied/no-support identities have null bc1.
    The complete count distribution is retained for correct downstream pooling.
    """
    keys=list(group_cols)
    if type(top_n) is not int or top_n<1:
        raise ValueError('top_n must be a positive integer')
    drops=tuple(edit_distance_drop)
    if any(type(n) is not int or n<1 for n in drops):
        raise ValueError('edit_distance_drop must contain positive integers')
    data=_validated_counts(df,keys)
    derived={'bc1','counts','top_bc1','edit_distance','total_counts','colony_purity',
             'observed_bc1_read_fraction','purity_basis','identity_ambiguous',
             'top_barcode_count','near_top_barcode_counts','all_bc1_counts_json',
             'legacy_filtered_top_n_total_counts','legacy_filtered_top_n_purity'}
    metadata=[c for c in data.columns if c not in set(keys)|derived]
    records=[]
    grouper=keys[0] if len(keys)==1 else keys
    # A group is one physical observation. Conflicting metadata cannot be
    # resolved by DataFrame.first(), which may synthesize a row from many inputs.
    for _,group in data.groupby(grouper,sort=False,dropna=False):
        conflicts=[c for c in metadata if group[c].nunique(dropna=False)>1]
        if conflicts:
            raise ValueError(f'conflicting physical/group metadata: {conflicts}')
        rec={c:group.iloc[0][c] for c in keys+metadata}
        # Python integer sums avoid an overflowing int64 aggregate denominator.
        counts={str(b):sum(int(v) for v in rows['counts'])
                for b,rows in group.groupby('bc1',sort=True)}
        ordered=sorted(counts.items(),key=lambda item:(-item[1],item[0]))
        total=sum(counts.values());maximum=ordered[0][1]
        winners=[b for b,n in ordered if n==maximum and n>0]
        unique=len(winners)==1
        winner=winners[0] if unique else pd.NA
        fraction=maximum/total if total else np.nan
        near=legacy_total=legacy_fraction=np.nan
        if unique:
            near=sum(n for b,n in ordered if lev_distance(winner,b) in drops)
            retained=[n for b,n in ordered if lev_distance(winner,b) not in drops]
            legacy_total=sum(retained[:top_n])
            legacy_fraction=maximum/legacy_total if legacy_total else np.nan
        rec.update(bc1=winner,top_bc1=winner,counts=maximum,total_counts=total,
                   edit_distance=0 if unique else np.nan,colony_purity=fraction,
                   observed_bc1_read_fraction=fraction,purity_basis='observed_bc1_read_fraction',
                   identity_ambiguous=len(winners)>1,top_barcode_count=len(winners),
                   near_top_barcode_counts=near,legacy_filtered_top_n_total_counts=legacy_total,
                   legacy_filtered_top_n_purity=legacy_fraction,
                   all_bc1_counts_json=json.dumps(counts,sort_keys=True,separators=(',',':')))
        records.append(rec)
    columns=keys+metadata+[c for c in ('bc1','top_bc1','counts','total_counts','edit_distance',
        'colony_purity','observed_bc1_read_fraction','purity_basis','identity_ambiguous',
        'top_barcode_count','near_top_barcode_counts','legacy_filtered_top_n_total_counts',
        'legacy_filtered_top_n_purity','all_bc1_counts_json') if c not in keys+metadata]
    return pd.DataFrame(records,columns=columns)
