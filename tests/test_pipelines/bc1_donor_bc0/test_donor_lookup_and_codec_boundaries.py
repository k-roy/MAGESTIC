"""One physical donor-reference file counts once; step 08 refuses mixed PCR-id encodings."""
import importlib.util
import os
import sys
from pathlib import Path

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.readers import data_loader
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import load_donor_bc0_lookup


@pytest.mark.parametrize('alias', ['symlink', 'hardlink'])
def test_aliased_reference_file_is_read_once(tmp_path, alias):
    """Two keyfile rows naming one physical file (symlink or hardlink) give 5 reads, not 10."""
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir()
    keys.mkdir()
    first = data / '20250101_plate_1A_1_complete_donor_bc0_counts.tsv'
    second = data / '20250101_plate_1B_1_complete_donor_bc0_counts.tsv'
    pd.DataFrame([dict(donor_bc0_fragment='C' * 30 + 'G' * 20 + 'A' * 10,
                       donor='A' * 99 + 'C' * 30, bc0='A' * 10, counts=5)]).to_csv(
        first, sep='\t', index=False)
    if alias == 'symlink':
        second.symlink_to(first)
    else:
        os.link(first, second)
    rows = [dict(description='donor-bc0', sequencing_date='20250101', outer_primers=outer,
                 sample_number=1, trafo_ID='trafo_3') for outer in ('plate_1A', 'plate_1B')]
    pd.DataFrame(rows).to_csv(keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    candidates = load_donor_bc0_lookup(data, keys)
    assert len(candidates) == 1
    assert int(candidates.iloc[0]['counts']) == 5


def _run_numbered08(tmp_path, monkeypatch, identities):
    inputs = []
    for index, identity in enumerate(identities):
        path = tmp_path / f'run{index}.tsv'
        pd.DataFrame([dict(bc1='G' * 20, donor='A' * 129, bc0='C' * 10, counts=7,
                           num_PCR_replicates=1, pcr_replicate_ids=identity)]).to_csv(
            path, sep='\t', index=False)
        inputs.append(str(path))
    source = Path(data_loader.__file__).parents[1] / 'snakemake/08_integrate_multiple_runs.py'
    spec = importlib.util.spec_from_file_location('test_codec_numbered08', source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    out = tmp_path / 'integrated.tsv'
    monkeypatch.setattr(sys, 'argv', [str(source), '--input-files', *inputs,
                                      '--output', str(out), '--n-jobs', '1'])
    module.main()
    return pd.read_csv(out, sep='\t')


METADATA = dict(gDNA_plate_name='prepA', inner_primers='innerA', outer_primers='plateA',
                plate_polymerase='Q5', condition='before', amplicon_level='trafo', sample_number=1)


def test_resequenced_pcr_with_one_encoding_is_one_replicate(tmp_path, monkeypatch):
    identity = data_loader.build_pcr_replicate_id(METADATA)
    result = _run_numbered08(tmp_path, monkeypatch, [identity, identity])
    assert result['counts'].tolist() == [14]
    assert int(result.iloc[0]['num_PCR_replicates']) == 1


def test_mixed_legacy_and_v1_encodings_are_refused(tmp_path, monkeypatch):
    legacy = '_'.join(str(METADATA[field]) for field in data_loader.PCR_REPLICATE_ID_FIELDS)
    with pytest.raises(ValueError, match='mix pcr:v1: and legacy'):
        _run_numbered08(tmp_path, monkeypatch,
                        [data_loader.build_pcr_replicate_id(METADATA), legacy])
