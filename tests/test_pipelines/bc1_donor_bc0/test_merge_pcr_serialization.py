"""The core merge must consume the numbered workflow's actual PCR-ID encoding."""
import numpy as np
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.link_collapse import decode_masks
from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources


@pytest.mark.parametrize('with_scalar', [False, True])
def test_numbered_producer_id_encoding_is_unioned(with_scalar):
    ids = decode_masks(np.array([[3], [6]], dtype=np.uint64), ['p1', 'p2', 'p3'])
    assert list(ids) == ['p1;p2', 'p2;p3']
    rows = pd.DataFrame(dict(bc1=['G'*20]*2, bc0=['A'*10]*2, donor=['A'*129]*2,
                             counts=[10, 20], pcr_replicate_ids=ids))
    if with_scalar:
        rows['num_PCR_replicates'] = 2
    result = merge_data_sources(df_2x300=rows).iloc[0]
    assert result['total_counts'] == 30
    assert result['total_pcr_replicates'] == 3
    assert result['pcr_support_exact']
    assert result['pcr_replicate_ids'] == 'p1;p2;p3'


def test_singular_identity_is_not_parsed_as_a_plural_list():
    rows = pd.DataFrame(dict(bc1=['G'*20], bc0=['A'*10], donor=['A'*129],
                             counts=[10], pcr_replicate_id=['one;opaque;identity']))
    # Existing semicolon output cannot represent this identity losslessly.
    # Refuse it rather than publishing something Step08 would count as3PCRs.
    with pytest.raises(ValueError, match='PCR identity'):
        merge_data_sources(df_2x300=rows)


@pytest.mark.parametrize('encoded', ['["p1", "p2"]', ' ["p1", "p2"]', '\n["p1", "p2"]\t'])
def test_json_plural_list_accepts_json_whitespace(encoded):
    rows = pd.DataFrame(dict(bc1=['G'*20], bc0=['A'*10], donor=['A'*129],
                             counts=[10], pcr_replicate_ids=[encoded]))
    result = merge_data_sources(df_2x300=rows).iloc[0]
    assert result['total_pcr_replicates'] == 2
    assert result['pcr_replicate_ids'] == 'p1;p2'


def test_opaque_json_looking_id_cannot_change_identity_on_roundtrip():
    rows = pd.DataFrame(dict(bc1=['G'*20], bc0=['A'*10], donor=['A'*129],
                             counts=[10], pcr_replicate_id=['["p1", "p2"]']))
    try:
        first = merge_data_sources(df_2x300=rows)
    except ValueError as error:
        assert 'PCR identity' in str(error)
        return  # Explicit refusal is valid; silently changing one identity to2 is not.
    second = merge_data_sources(df_2x300=first)
    assert first.iloc[0]['total_pcr_replicates'] == second.iloc[0]['total_pcr_replicates'] == 1
