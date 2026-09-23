"""Execute the array-task shell body; only interactive shell setup is omitted."""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

import magestic.pipelines.wgs as wgs


def invoke(tmp_path, manifest_text, *, index='1', stale=False, foreign=False):
    manifest = tmp_path / 'manifest with spaces.tsv'
    manifest.write_text(manifest_text)
    config_manifest = manifest
    if foreign:
        config_manifest = tmp_path / 'other.tsv'
        config_manifest.write_text('sample\tbam\treference\nother\tb.bam\tr.fa\n')
    config = tmp_path / 'config with spaces.yaml'
    config.write_text(yaml.safe_dump({'cohort_name': 'declared', 'sequencing_date': '20260920',
        'paths': {'project_root': str(tmp_path), 'manifest': str(config_manifest)}}))
    calls = tmp_path / 'called.json'
    bin_dir = tmp_path / 'bin'; bin_dir.mkdir()
    cli = bin_dir / 'magestic-wgs-bc1'
    cli.write_text('#!' + sys.executable + '\nimport json,os,sys\nfrom pathlib import Path\n'
        "Path(os.environ['CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n")
    cli.chmod(0o755)
    if stale:
        out = tmp_path / 'processed' / 'bc1_calls'; out.mkdir(parents=True)
        for suffix in ('bc1_calls.tsv', 'bc1_tracks.tsv'):
            (out / ('s.' + suffix)).write_text('invalid stale output\n')
    template = Path(wgs.__file__).parent / 'slurm/template.sbatch'
    text = template.read_text()
    assert text.count('source ~/.bashrc\n') == 1
    # The user's startup file is not part of the application under test.
    text = text.replace('source ~/.bashrc\n', '', 1)
    env = {**os.environ, 'CONFIG': str(config), 'MANIFEST': str(manifest),
        'SLURM_ARRAY_TASK_ID': index, 'PROCESSED_DIR': str(tmp_path / 'processed'),
        'PATH': str(bin_dir) + os.pathsep + os.environ['PATH'], 'PYTHON': sys.executable,
        'CAPTURE': str(calls)}
    result = subprocess.run(['bash', '-c', text], env=env, cwd=tmp_path,
                            text=True, capture_output=True, timeout=30)
    arguments = json.loads(calls.read_text()) if calls.exists() else None
    return result, arguments


@pytest.mark.parametrize('manifest_text,expected', [
    ('sample\tbam\treference\ns\tb.bam\tr.fa\n', 's'),
    ('sample\tbam\treference\nwell A1\tb.bam\tr.fa\n', 'well A1'),
    ('sample\tbam\treference\n# comment\n\ns\tb.bam\tr.fa\n', 's'),
    ('bam\treference\tsample\nb.bam\tr.fa\ts\n', 's'),
])
def test_array_selection_matches_manifest_parser(tmp_path, manifest_text, expected):
    result, args = invoke(tmp_path, manifest_text)
    assert result.returncode == 0, result.stderr
    assert args is not None and args[args.index('--sample')+1] == expected


def test_stale_output_names_do_not_bypass_python_cache_validation(tmp_path):
    result, args = invoke(tmp_path, 'sample\tbam\treference\ns\tb.bam\tr.fa\n', stale=True)
    assert result.returncode == 0, result.stderr
    assert args is not None, 'stage entrypoint must validate content/identity before reuse'


@pytest.mark.parametrize('index', ['0', '-1', '2', 'bad'])
def test_invalid_array_task_never_invokes_caller(tmp_path, index):
    result, args = invoke(tmp_path, 'sample\tbam\treference\ns\tb.bam\tr.fa\n', index=index)
    assert result.returncode != 0
    assert args is None


def test_config_and_shell_cannot_select_different_manifests(tmp_path):
    result, args = invoke(tmp_path, 'sample\tbam\treference\ns\tb.bam\tr.fa\n', foreign=True)
    assert result.returncode != 0
    assert args is None
