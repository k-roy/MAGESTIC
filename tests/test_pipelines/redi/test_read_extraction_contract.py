"""Independent molecular-boundary fixtures; no observed calls used as truth."""
import pandas as pd
import pytest
from magestic.pipelines.redi.readers.fastq import extract_bc1_from_seq, extract_redi_bc_vectorized

BC='ACGTATGCCATTCAGGTAAC'
PREFIX='CAGGTCATGCTC'
SUFFIX='CCCTAGGATACG'
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def call(seq,**kwargs):
    return extract_bc1_from_seq(seq,fwd_bc1_prefix=PREFIX,fwd_bc1_suffix=SUFFIX,**kwargs)

@pytest.mark.parametrize('suffix',['CCCTAGGATACG','ATAACAGGGT'])
@pytest.mark.parametrize('reverse',[False,True])
def test_independently_declared_cassette_orientation(suffix,reverse):
    molecule=PREFIX+BC+suffix
    assert extract_bc1_from_seq(rc(molecule) if reverse else molecule,
        fwd_bc1_prefix=PREFIX,fwd_bc1_suffix=suffix)==BC

@pytest.mark.parametrize('sequence',['T'*100,PREFIX+BC+'G'*20,'C'*16+BC,
                                   PREFIX+BC[:11],PREFIX+BC[:8]+'N'+BC[9:]+SUFFIX])
def test_no_call_without_complete_declared_boundary(sequence):
    assert call(sequence) is None

def test_first_invalid_anchor_does_not_hide_later_valid_molecule():
    assert call(PREFIX+'T'*40+PREFIX+BC+SUFFIX)==BC

def test_two_disagreeing_molecules_refuse_identity():
    assert call(PREFIX+BC+SUFFIX+'GGG'+PREFIX+'T'*20+SUFFIX) is None

def test_opposite_orientations_must_not_disagree():
    assert call(PREFIX+BC+SUFFIX+'GGG'+rc(PREFIX+'T'*20+SUFFIX)) is None

def test_repeated_identical_evidence_preserves_identity():
    assert call((PREFIX+BC+SUFFIX)*2)==BC

@pytest.mark.parametrize('size',[19,20,21])
def test_length_tolerance_requires_suffix(size):
    barcode='A'*size
    assert call(PREFIX+barcode+SUFFIX)==barcode

def test_explicit_fixed_offset_requires_opt_in():
    assert call('C'*16+BC+'GG',fallback_prefix_length=16) is None
    assert call('C'*16+BC+'GG',fallback_prefix_length=16,allow_positional_fallback=True)==BC

@pytest.mark.parametrize('kwargs',[{'bc1_length':0},{'bc1_length':20.5},
                                 {'length_tolerance':-1},{'length_tolerance':True}])
def test_invalid_dimensions_fail(kwargs):
    with pytest.raises(ValueError):call(PREFIX+BC+SUFFIX,**kwargs)

def test_prefix_only_mode_remains_anchored_and_validated():
    assert extract_bc1_from_seq(PREFIX+BC,fwd_bc1_prefix=PREFIX)==BC
    assert extract_bc1_from_seq('T'*100,fwd_bc1_prefix=PREFIX) is None

def vector(sequences,prefixes,lengths):
    return extract_redi_bc_vectorized(pd.Series(sequences,dtype=object),
        pd.Series(prefixes,dtype=object),pd.Series(lengths,dtype=object))

def test_empty_vector_preserves_index():
    seq=pd.Series([],index=pd.Index([],name='observations'),dtype=object)
    result=extract_redi_bc_vectorized(seq,seq.copy(),seq.copy())
    assert result.empty and result.index.equals(seq.index)

def test_vector_requires_matching_indices():
    with pytest.raises(ValueError,match='align|index'):
        extract_redi_bc_vectorized(pd.Series(['ACGT'],index=[1]),
            pd.Series(['A'],index=[2]),pd.Series([3],index=[1]))

@pytest.mark.parametrize('length',[20.5,-1,float('inf'),True,'twenty'])
def test_vector_refuses_invalid_lengths(length):
    with pytest.raises(ValueError):vector([PREFIX+BC],[PREFIX],[length])

def test_disabled_extraction_requires_empty_anchor_and_zero_length():
    assert vector(['ACGT'],[None],[0]).isna().all()
    with pytest.raises(ValueError):vector([PREFIX+BC],[PREFIX],[0])

def test_vector_case_invalid_bases_short_reads_and_ambiguity():
    result=vector([PREFIX.lower()+BC.lower(),PREFIX+'N'*20,PREFIX+BC[:5],
                   PREFIX+BC+PREFIX+'T'*20],[PREFIX]*4,[20]*4)
    assert result.iloc[0]==BC
    assert result.iloc[1:].isna().all()

def test_vector_integral_spreadsheet_length_is_supported():
    assert vector([PREFIX+BC],[PREFIX],[20.0]).iloc[0]==BC
