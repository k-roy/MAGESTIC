"""Cache completion requires the exact five Stage A artifacts, not five records."""
import gzip
import json
from pathlib import Path

import pytest

from magestic.pipelines.redi.core.config import REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline import step_a_trim_counts as stage


def test_repeated_counts_provenance_cannot_replace_missing_companion_outputs(tmp_path, monkeypatch):
    paths = REDIPaths(tmp_path, tmp_path/'raw', tmp_path/'intermediate',
                      tmp_path/'processed', tmp_path/'keys')
    paths.raw_fastq.mkdir()
    config = REDIConfig('cache_fixture', '20250627', paths)
    sample = '20250627_plate_4A_sample_1'
    for mate in (1, 2):
        with gzip.open(paths.raw_fastq/f'{sample}_R{mate}.fastq.gz', 'wt') as handle:
            handle.write(f'@same/{mate}\nACGTACGT\n+\nIIIIIIII\n')

    def trim(first, second, out_first, out_second, **kwargs):
        for source, target in ((first, out_first), (second, out_second)):
            with gzip.open(source, 'rt') as handle:
                Path(target).write_text(handle.read())

    def merge(first, second, output, **kwargs):
        Path(output).write_text(Path(first).read_text())

    monkeypatch.setattr(stage, '_bbduk', trim)
    monkeypatch.setattr(stage, '_bbmerge', merge)
    monkeypatch.setattr(stage, '_tool_identity', lambda: {'controlled_test_tool': 'unchanged'})
    stage.run_step_a_one_sample(sample, config)
    counts = paths.counts_tables/f'{sample}_R1_R2_counts.tsv'
    sidecar = Path(str(counts)+'.provenance.json')
    saved = json.loads(sidecar.read_text())
    main = next(item for item in saved['outputs'] if item['path'] == str(counts))
    companions = [Path(item['path']) for item in saved['outputs'] if item['path'] != str(counts)]
    assert len(companions) == 4
    for path in companions:
        path.unlink()
    saved['outputs'] = [main.copy() for _ in range(5)]
    sidecar.write_text(json.dumps(saved))
    with pytest.raises(FileExistsError, match='complete provenance'):
        stage.run_step_a_one_sample(sample, config, overwrite=False)
