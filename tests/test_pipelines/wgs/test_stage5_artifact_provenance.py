"""Every published join table needs its own output and sample-scope receipt."""
import hashlib
import json
from pathlib import Path

import pytest

from magestic.pipelines.wgs.config import WGSConfig, WGSPaths
from magestic.pipelines.wgs.pipeline import step_2_donor_precall as s2
from magestic.pipelines.wgs.pipeline import step_4_reconcile as s4
from magestic.pipelines.wgs.pipeline import step_5_join_outcomes as s5


def declared_join(tmp_path):
    manifest = tmp_path / 'samples.tsv'
    manifest.write_text('sample\tbam\treference\ns2\ttwo.bam\tr.fa\ns1\tone.bam\tr.fa\n')
    cfg = WGSConfig('declared', '20260920', WGSPaths(tmp_path, manifest,
        tmp_path/'intermediate', tmp_path/'processed', tmp_path/'keys'))
    a = s2.table_path(cfg); a.parent.mkdir(parents=True)
    a.write_text('sample_name\tediting_outcome\ns1\tdesigned_variant_installed\ns2\tdesigned_variant_installed\n')
    b = s4.table_path(cfg)
    b.write_text('sample\tlabel\tconfidence\ns1\tSINGLE_EDIT\thigh\ns2\tNO_EDIT\thigh\n')
    return cfg, a, b


@pytest.mark.parametrize('artifact', ['outcome_join.tsv', 'outcome_disagreements.tsv', 'outcome_agreement_matrix.tsv'])
def test_stage5_each_artifact_receipt_binds_inputs_and_output(tmp_path, artifact):
    cfg, a, b = declared_join(tmp_path)
    outputs = s5.run_step_5(cfg)
    output = next(p for p in outputs if p.name == artifact)
    receipt_path = Path(str(output) + '.provenance.json')
    assert receipt_path.is_file(), artifact + ' lacks producer provenance'
    receipt = json.loads(receipt_path.read_text())
    actual_inputs = {Path(row['path']).resolve(): row['md5'] for row in receipt['inputs']}
    for path in (a, b, cfg.paths.manifest):
        assert actual_inputs[path.resolve()] == hashlib.md5(path.read_bytes()).hexdigest()
    actual_outputs = {Path(row['path']).resolve(): row['md5'] for row in receipt['outputs']}
    assert actual_outputs[output.resolve()] == hashlib.md5(output.read_bytes()).hexdigest()
