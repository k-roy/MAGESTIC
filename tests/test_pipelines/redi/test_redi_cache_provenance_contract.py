"""REDI cache/provenance hypotheses H1-H3 as before-tests, staged 2026-09-20.

Design and source citations: arms/redi_cache_provenance_before_tests.md.
Copy into tests/test_pipelines/redi/ and run from the workspace root.

Naming: ``test_hyp_*`` cases state the expected behavior that a hypothesis says is
missing; a FAILURE on current source reproduces that hypothesis. ``test_control_*``
cases change something current code already detects; they must PASS on current
source, which shows each instrument can see a change. An assertion message that
starts with ``PRECONDITION`` is an instrument failure, never a reproduction.

No product file is modified. H1 changes a dependency only inside a copy of the
imported ``magestic`` package under ``tmp_path``, and runs REDI in a child
interpreter whose ``magestic`` import is asserted to resolve to that copy. H2 never
executes BBTools or Java: ``_bbduk``/``_bbmerge`` are stubbed exactly as in
test_stage_a_contract.py, and the real ``_tool_identity()`` fingerprints a
test-only installation on PATH. H3 stubs the plotting modules as
test_ghi_followthrough_contract.py does; it asserts receipts, not rendering.
"""
import gzip
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import ModuleType

import pandas as pd
import pytest

import magestic
from magestic.pipelines.redi.core.config import REDIArrayPair, REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline import step_a_trim_counts as stage_a
from magestic.pipelines.redi.pipeline import step_h_plot as stage_h
from magestic.pipelines.redi.pipeline.run_pipeline import verify_completed_artifacts


# The copied package in H1 is taken from whatever magestic this process imported, so
# importing any other tree (for example the venv's editable install) would silently
# test the wrong source. Fail at collection instead. Valid only once this file sits in
# <checkout>/tests/test_pipelines/redi/.
_CHECKOUT_PACKAGE = Path(__file__).resolve().parents[3]/'src'/'magestic'
if Path(magestic.__file__).resolve().parent != _CHECKOUT_PACKAGE:
    raise RuntimeError(
        f'PRECONDITION: imported magestic is {Path(magestic.__file__).resolve().parent}, not the '
        f'checkout this test file belongs to ({_CHECKOUT_PACKAGE}); run under '
        'scripts/run_locked_v2.py, which sets PYTHONPATH to the isolated src')


# A real code change (a new function), appended only to a file in a copied tree.
MUTATION = ('\n\ndef _redi_cache_provenance_mutation_probe():\n'
            '    """Test-only change made in a copied package tree."""\n'
            '    return "changed"\n')


# --------------------------------------------------------------------------------------
# H1: a change to a non-REDI module that REDI executes must invalidate completed REDI
# artifacts. Stage G runs for real in a child interpreter importing a copied package.
# --------------------------------------------------------------------------------------

_STAGE_G_CHILD = r'''
import json, sys
from pathlib import Path
import pandas as pd
import magestic
from magestic.pipelines.redi.core.config import REDIArrayPair, REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline.step_g_wgs_check import run_step_g
from magestic.pipelines.redi.pipeline.run_pipeline import verify_completed_artifacts
root = Path(sys.argv[1])
paths = REDIPaths(root, root/'raw', root/'intermediate', root/'processed', root/'keys')
paths.keyfiles.mkdir(parents=True)
paths.combined_data.mkdir(parents=True)
config = REDIConfig('independent', '20250627', paths,
                    array_pairs=[REDIArrayPair('old', 'new', 'earlier', purpose_b='later')])
barcode = 'A' * 20
cross = pd.DataFrame([dict(MAT_alpha_array_1536_barcode='6156', REDI_plate_color_384='RED',
    REDI_row_384='A', REDI_col_384='1', num_independent_replicates=2,
    old_bc1=barcode, new_bc1=barcode, old_REDI_bc_bc1_counts=30, new_REDI_bc_bc1_counts=30,
    old_colony_purity=1.0, new_colony_purity=1.0)])
cross.to_csv(paths.combined_data/'old_new_merged_REDI_df.tsv', sep='\t', index=False)
positions = ['MAT_alpha_array_1536_barcode', 'REDI_plate_color_384', 'REDI_row_384', 'REDI_col_384']
cross[positions].assign(sample_name='registered').to_csv(
    paths.keyfiles/'WGS_plate_key.tsv', sep='\t', index=False)
config.wgs_outcome_table = root/'calls.tsv'
pd.DataFrame({'sample': ['registered'], 'barcode': [barcode]}).to_csv(
    config.wgs_outcome_table, sep='\t', index=False)
output = run_step_g(config)
verify_completed_artifacts([output])  # a freshly completed artifact must be accepted
print(json.dumps({'magestic': str(Path(magestic.__file__).resolve()), 'output': str(output)}))
'''

_GATE_CHILD = r'''
import json, sys
from pathlib import Path
import magestic
from magestic.pipelines.redi.pipeline.run_pipeline import verify_completed_artifacts
try:
    verify_completed_artifacts([sys.argv[1]])
except ValueError as exc:
    result = {'accepted': False, 'reason': str(exc)}
else:
    result = {'accepted': True}
result['magestic'] = str(Path(magestic.__file__).resolve())
print(json.dumps(result))
'''


def _copied_package(tmp_path):
    """Copy the imported magestic package (the checkout under test) into tmp_path."""
    source = Path(magestic.__file__).resolve().parent
    root = tmp_path/'copied_src'
    shutil.copytree(source, root/'magestic',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '.snakemake'))
    return root.resolve()


def _child(tmp_path, package_root, script, *args):
    env = os.environ.copy()
    # The copy comes first; the rest of this interpreter's path keeps pandas/yaml/TRACE
    # resolvable exactly as in the parent test process.
    env['PYTHONPATH'] = os.pathsep.join([str(package_root), *(str(p) for p in sys.path if p)])
    env['PYTHONDONTWRITEBYTECODE'] = '1'
    env['MPLCONFIGDIR'] = str(tmp_path/'mpl')
    result = subprocess.run([sys.executable, '-c', script, *(str(a) for a in args)],
                            cwd=tmp_path, env=env, text=True, capture_output=True, timeout=300)
    assert result.returncode == 0, 'PRECONDITION: child failed\n' + result.stdout + '\n' + result.stderr
    payload = json.loads(result.stdout.strip().splitlines()[-1])
    imported = Path(payload['magestic'])
    assert imported.is_relative_to(package_root), (
        f'PRECONDITION: child imported {imported}, not the copied tree under {package_root}; '
        'a comparison against the real checkout would be vacuous')
    return payload


def _stage_g_gate_after_change(tmp_path, relative):
    package = _copied_package(tmp_path)
    target = package/'magestic'/relative
    assert target.is_file(), f'PRECONDITION: {relative} is not part of the imported package'
    produced = _child(tmp_path, package, _STAGE_G_CHILD, tmp_path/'screen')
    target.write_text(target.read_text() + MUTATION)
    return _child(tmp_path, package, _GATE_CHILD, produced['output'])


def test_control_h1_redi_module_change_invalidates_completed_stage_g(tmp_path):
    after = _stage_g_gate_after_change(tmp_path, 'pipelines/redi/core/wgs_check.py')
    assert after['accepted'] is False, (
        'CONTROL FAILED: a change inside pipelines/redi no longer invalidates a completed Stage G')


@pytest.mark.parametrize('relative', [
    'utils/provenance.py',
    'pipelines/wgs/readers/wgs_calls.py',
], ids=['shared_provenance_module', 'stage_g_wgs_calls_reader'])
def test_hyp_h1_non_redi_dependency_change_invalidates_completed_stage_g(tmp_path, relative):
    after = _stage_g_gate_after_change(tmp_path, relative)
    assert after['accepted'] is False, (
        f'H1 REPRODUCED: magestic/{relative} changed, yet the completion gate that a warm '
        'workflow runs accepted the previously completed Stage G artifact')


# The Stage A skip-if-valid cache (_cache_matches) is the other consumer of the same
# digest. BBTools is stubbed exactly as in test_stage_a_cache_contract.py, including a
# constant _tool_identity, so only the source change can distinguish the two runs.
_STAGE_A_CHILD = r'''
import gzip, hashlib, json, sys
from pathlib import Path
import magestic
from magestic.pipelines.redi.core.config import REDIConfig, REDIPaths
from magestic.pipelines.redi.pipeline import step_a_trim_counts as stage
root, mode = Path(sys.argv[1]), sys.argv[2]
paths = REDIPaths(root, root/'raw', root/'intermediate', root/'processed', root/'keys')
config = REDIConfig('fixture', '20250627', paths)
prefix = '20250627_plate_4A_sample_1'
def trim(first, second, out_first, out_second, **kwargs):
    for source, target in ((first, out_first), (second, out_second)):
        with gzip.open(source, 'rt') as handle:
            Path(target).write_text(handle.read())
def merge(first, second, output, **kwargs):
    Path(output).write_text(Path(first).read_text())
stage._bbduk, stage._bbmerge = trim, merge
stage._tool_identity = lambda: {'controlled_test_tool': 'unchanged'}
sidecar = paths.counts_tables/f'{prefix}_R1_R2_counts.tsv.provenance.json'
if mode == 'produce':
    paths.raw_fastq.mkdir(parents=True)
    molecule = 'CAGGTCATGCTC' + 'ACGTATGCCATTCAGGTAAC' + 'CCCTAGGATACG'
    for mate in (1, 2):
        with gzip.open(paths.raw_fastq/f'{prefix}_R{mate}.fastq.gz', 'wt') as handle:
            for index in range(3):
                handle.write(f'@read{index}/{mate}\n{molecule}\n+\n' + 'I'*len(molecule) + '\n')
    stage.run_step_a_one_sample(prefix, config)
before = hashlib.sha256(sidecar.read_bytes()).hexdigest()
try:
    stage.run_step_a_one_sample(prefix, config, overwrite=False)
except FileExistsError:
    reused = False
else:
    reused = hashlib.sha256(sidecar.read_bytes()).hexdigest() == before
print(json.dumps({'magestic': str(Path(magestic.__file__).resolve()), 'reused': reused}))
'''


def _stage_a_cache_after_change(tmp_path, relative):
    package = _copied_package(tmp_path)
    target = package/'magestic'/relative
    assert target.is_file(), f'PRECONDITION: {relative} is not part of the imported package'
    screen = tmp_path/'screen'
    produced = _child(tmp_path, package, _STAGE_A_CHILD, screen, 'produce')
    assert produced['reused'] is True, 'PRECONDITION: an unchanged Stage A cache was not reusable'
    target.write_text(target.read_text() + MUTATION)
    return _child(tmp_path, package, _STAGE_A_CHILD, screen, 'reuse')


def test_control_h1_redi_module_change_invalidates_stage_a_cache(tmp_path):
    after = _stage_a_cache_after_change(tmp_path, 'pipelines/redi/pipeline/step_a_trim_counts.py')
    assert after['reused'] is False, (
        'CONTROL FAILED: Stage A silently reused its cache after its own module changed')


def test_hyp_h1_shared_provenance_change_invalidates_stage_a_cache(tmp_path):
    after = _stage_a_cache_after_change(tmp_path, 'utils/provenance.py')
    assert after['reused'] is False, (
        'H1 REPRODUCED: magestic/utils/provenance.py (which supplies md5_file and the receipt) '
        'changed, yet Stage A silently reused its overwrite=False cache')


# --------------------------------------------------------------------------------------
# H2: the warm-run completion gate must notice a changed BBTools/Java implementation
# even when the launcher script bytes are unchanged.
# --------------------------------------------------------------------------------------

MOLECULE = 'CAGGTCATGCTC' + 'ACGTATGCCATTCAGGTAAC' + 'CCCTAGGATACG'


def _fake_bbtools(tmp_path, monkeypatch):
    """A test-only installation laid out as _tool_identity() and _run() expect.

    Nothing here is ever executed: the launchers and java exit 97 if run.
    """
    install = tmp_path/'fake_bbtools'
    (install/'current'/'jgi').mkdir(parents=True)
    (install/'jni').mkdir()
    for name in ('bbduk.sh', 'bbmerge.sh'):
        launcher = install/name
        launcher.write_text('#!/bin/sh\n# test-only launcher; never executed\nexit 97\n')
        launcher.chmod(0o755)
    (install/'current'/'jgi'/'BBDuk.class').write_bytes(b'\xca\xfe\xba\xbe test-only BBDuk v1')
    (install/'current'/'jgi'/'BBMerge.class').write_bytes(b'\xca\xfe\xba\xbe test-only BBMerge v1')
    (install/'jni'/'libbbtoolsjni.test').write_bytes(b'test-only JNI library v1')
    jdk = tmp_path/'fake_jdk'/'bin'
    jdk.mkdir(parents=True)
    java = jdk/'java'
    java.write_text('#!/bin/sh\n# test-only java; never executed\nexit 97\n')
    java.chmod(0o755)
    monkeypatch.setenv('PATH', os.pathsep.join([str(install), str(jdk), os.environ.get('PATH', '')]))
    return {'compiled_class': install/'current'/'jgi'/'BBDuk.class',
            'jni_library': install/'jni'/'libbbtoolsjni.test',
            'java_runtime': java,
            'launcher': install/'bbduk.sh'}


def _completed_stage_a(tmp_path, monkeypatch):
    components = _fake_bbtools(tmp_path, monkeypatch)
    paths = REDIPaths(tmp_path, tmp_path/'raw', tmp_path/'intermediate',
                      tmp_path/'processed', tmp_path/'keys')
    paths.raw_fastq.mkdir()
    config = REDIConfig('fixture', '20250627', paths)
    prefix = '20250627_plate_4A_sample_1'
    for mate in (1, 2):
        with gzip.open(paths.raw_fastq/f'{prefix}_R{mate}.fastq.gz', 'wt') as handle:
            for index in range(3):
                handle.write(f'@read{index}/{mate}\n{MOLECULE}\n+\n' + 'I'*len(MOLECULE) + '\n')

    def trim(first, second, out_first, out_second, **kwargs):
        for source, target in ((first, out_first), (second, out_second)):
            with gzip.open(source, 'rt') as handle:
                Path(target).write_text(handle.read())

    def merge(first, second, output, **kwargs):
        Path(output).write_text(Path(first).read_text())

    monkeypatch.setattr(stage_a, '_bbduk', trim)
    monkeypatch.setattr(stage_a, '_bbmerge', merge)
    stage_a.run_step_a_one_sample(prefix, config)
    counts = paths.counts_tables/f'{prefix}_R1_R2_counts.tsv'
    recorded = json.loads(Path(str(counts)+'.provenance.json').read_text())['parameters']['tool_implementation']
    assert recorded == stage_a._tool_identity(), 'PRECONDITION: tool identity is not stable between calls'
    verify_completed_artifacts([counts])  # a freshly completed Stage A artifact is accepted
    return config, prefix, counts, components, recorded


def _change(path):
    path.write_bytes(path.read_bytes() + b'\n# changed by test\n')  # keeps the executable bit


@pytest.mark.parametrize('component', ['compiled_class', 'jni_library', 'java_runtime'])
def test_hyp_h2_warm_gate_rejects_stage_a_after_tool_implementation_change(tmp_path, monkeypatch, component):
    config, prefix, counts, components, recorded = _completed_stage_a(tmp_path, monkeypatch)
    _change(components[component])
    assert stage_a._tool_identity() != recorded, (
        f'PRECONDITION: Stage A _tool_identity() does not see the {component} change')
    try:
        verify_completed_artifacts([counts])
    except ValueError:
        return
    pytest.fail(f'H2 REPRODUCED: the {component} changed with launcher bytes unchanged; Stage A itself '
                'would refuse this cache, but the warm-run completion gate accepted it')


def test_control_h2_warm_gate_rejects_stage_a_after_launcher_change(tmp_path, monkeypatch):
    config, prefix, counts, components, recorded = _completed_stage_a(tmp_path, monkeypatch)
    _change(components['launcher'])
    assert stage_a._tool_identity() != recorded, 'PRECONDITION: launcher change not seen by _tool_identity()'
    try:
        verify_completed_artifacts([counts])
    except ValueError:
        return
    pytest.fail('CONTROL FAILED: a changed launcher (an attested Stage A input) was accepted by the gate')


def test_control_h2_stage_a_itself_refuses_cache_after_class_change(tmp_path, monkeypatch):
    config, prefix, counts, components, recorded = _completed_stage_a(tmp_path, monkeypatch)
    _change(components['compiled_class'])
    current = stage_a._tool_identity()
    assert current != recorded, 'PRECONDITION: class change not seen by _tool_identity()'
    try:
        stage_a.run_step_a_one_sample(prefix, config, overwrite=False)
    except FileExistsError:
        return  # refused the stale cache
    saved = json.loads(Path(str(counts)+'.provenance.json').read_text())['parameters']['tool_implementation']
    assert saved == current, 'CONTROL FAILED: Stage A reused a cache recorded under another tool implementation'


# --------------------------------------------------------------------------------------
# H3: the Stage H receipt must record the plotting implementation that produced it.
# --------------------------------------------------------------------------------------

PLOT_VERSIONS = {'matplotlib': 'redi-h3-matplotlib-version-sentinel',
                 'seaborn': 'redi-h3-seaborn-version-sentinel'}


def _override_distribution_versions(monkeypatch, overrides):
    """Make importlib.metadata.version agree with the patched module __version__."""
    real = importlib.metadata.version

    def version(name):
        return overrides[name] if name in overrides else real(name)
    monkeypatch.setattr(importlib.metadata, 'version', version)


def _plotting_stub(monkeypatch):
    captured = []
    matplotlib = ModuleType('matplotlib')
    pyplot = ModuleType('matplotlib.pyplot')
    sns = ModuleType('seaborn')
    matplotlib.__version__ = PLOT_VERSIONS['matplotlib']
    sns.__version__ = PLOT_VERSIONS['seaborn']
    matplotlib.use = lambda *a, **k: None
    matplotlib.pyplot = pyplot
    for name in ('figure', 'xlabel', 'ylabel', 'title', 'xticks', 'tight_layout', 'close', 'legend'):
        setattr(pyplot, name, lambda *a, **k: None)
    pyplot.savefig = lambda path, **kw: Path(path).write_bytes(b'plot-data-test-only')
    sns.barplot = lambda *, x, y, **kw: captured.append((list(x), list(y)))
    for name, module in [('matplotlib', matplotlib), ('matplotlib.pyplot', pyplot), ('seaborn', sns)]:
        monkeypatch.setitem(sys.modules, name, module)
    _override_distribution_versions(monkeypatch, dict(PLOT_VERSIONS))
    return captured


def _stage_h_receipt(tmp_path, captured):
    paths = REDIPaths(tmp_path, tmp_path/'raw', tmp_path/'intermediate',
                      tmp_path/'processed', tmp_path/'keys')
    paths.keyfiles.mkdir()
    paths.combined_data.mkdir(parents=True)
    config = REDIConfig('independent', '20250627', paths,
                        array_pairs=[REDIArrayPair('old', 'new', 'earlier', purpose_b='later')])
    barcode = 'A' * 20
    pd.DataFrame([dict(MAT_alpha_array_1536_barcode='6156', REDI_plate_color_384='RED',
                       REDI_row_384='A', REDI_col_384=str(column), num_independent_replicates=2,
                       old_bc1=barcode, new_bc1=barcode, old_REDI_bc_bc1_counts=30,
                       new_REDI_bc_bc1_counts=30, old_colony_purity=1.0, new_colony_purity=1.0)
                  for column in (1, 2)]).to_csv(
        paths.combined_data/'old_new_merged_REDI_df.tsv', sep='\t', index=False)
    stage_h.run_step_h(config)
    assert captured == [([2], [1])], 'PRECONDITION: Stage H did not plot through the stubbed libraries'
    return json.loads((paths.plots/'bc1_registered_position_counts.tsv.provenance.json').read_text())


def _recorded_strings(node):
    if isinstance(node, dict):
        for key, value in node.items():
            yield str(key)
            yield from _recorded_strings(value)
    elif isinstance(node, (list, tuple)):
        for value in node:
            yield from _recorded_strings(value)
    else:
        yield str(node)


def _records(receipt, value):
    return any(value in text for text in _recorded_strings(receipt))


def test_control_h3_stage_h_receipt_records_the_pandas_version_in_effect(tmp_path, monkeypatch):
    captured = _plotting_stub(monkeypatch)
    sentinel = 'redi-h3-pandas-version-sentinel'
    monkeypatch.setattr(pd, '__version__', sentinel)
    _override_distribution_versions(monkeypatch, {**PLOT_VERSIONS, 'pandas': sentinel})
    receipt = _stage_h_receipt(tmp_path, captured)
    assert _records(receipt, sentinel), 'CONTROL FAILED: the H receipt no longer records the pandas version'


@pytest.mark.parametrize('library', ['matplotlib', 'seaborn'])
def test_hyp_h3_stage_h_receipt_records_the_plotting_library_version_in_effect(tmp_path, monkeypatch, library):
    captured = _plotting_stub(monkeypatch)
    receipt = _stage_h_receipt(tmp_path, captured)
    assert _records(receipt, PLOT_VERSIONS[library]), (
        f'H3 REPRODUCED: Stage H plotted with {library} {PLOT_VERSIONS[library]!r}, '
        'and its completion receipt does not record that version')
