"""Large integer counts and malformed distributions must not silently round."""
import json
import pandas as pd
import pytest
from magestic.pipelines.redi.core.purity import compute_purity
from magestic.pipelines.redi.core.cross_array import compare_arrays

BIG=9007199254740993
BC='A'*20

def test_purity_decimal_string_count_stays_exact():
    data=pd.DataFrame({'sample_name':['s'],'REDI_bc':['T'*20],'bc1':[BC],'counts':[str(BIG)+'.0']})
    result=compute_purity(data).iloc[0]
    assert result['counts']==BIG and result['total_counts']==BIG

def fixture(raw=None):
    return pd.DataFrame([{'MAT_alpha_array_1536_barcode':'6156','REDI_plate_color_384':'RED_384',
                         'REDI_row_384':'A','REDI_col_384':'1','sample_name':sample,'REDI_bc':'T'*20,
                         'bc1':BC,'REDI_bc_bc1_counts':BIG,'REDI_bc_bc1_total_counts':BIG,
                         'all_bc1_counts_json':raw or json.dumps({BC:BIG})} for sample in ('earlier','later')])

def test_cross_array_large_integer_read_count_stays_exact():
    data=fixture();result=compare_arrays(data,data.sample_name=='earlier',data.sample_name=='later','a','b')
    assert result.a_REDI_bc_bc1_counts.iloc[0]==BIG
    assert result.a_REDI_bc_bc1_total_counts.iloc[0]==BIG

def test_distribution_must_be_json_object_not_list_of_pairs():
    data=fixture(json.dumps([[BC,BIG]]))
    with pytest.raises(ValueError,match='distribution'):
        compare_arrays(data,data.sample_name=='earlier',data.sample_name=='later','a','b')

def test_boolean_read_count_is_not_one_molecule():
    with pytest.raises(ValueError,match='counts'):
        compute_purity(pd.DataFrame({'sample_name':['s'],'REDI_bc':['T'*20],'bc1':[BC],'counts':[True]}))
