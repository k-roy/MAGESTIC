"""A scalar annotation tuple must come from one actual GDB observation."""
import pandas as pd
from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources


def test_missing_guide_does_not_borrow_a_different_library_annotation():
    bdb = pd.DataFrame([dict(bc1='G'*20, donor='A'*129, bc0='A'*10,
                             counts=10, pcr_replicate_id='p1')])
    gdb = pd.DataFrame([
        dict(donor='A'*129, bc0='A'*10, guide=None, library_ID='L1', scaffold='S1', counts=1000),
        dict(donor='A'*129, bc0='A'*10, guide='T'*20, library_ID='L2', scaffold='S2', counts=10)])
    result = merge_data_sources(df_2x300=bdb, df_gdb=gdb).iloc[0]
    assert (result['library_id'], result['scaffold'], result['counts_gdb']) == ('L1', 'S1', 1000)
    assert pd.isna(result['guide'])
    assert result['gdb_ambiguous']
    assert result['total_counts'] == 10 and result['total_pcr_replicates'] == 1
