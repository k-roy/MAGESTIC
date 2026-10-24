"""Exercise the guide-link CLI with arguments from its shipped workflow rule."""
import ast
import json
from pathlib import Path
import re
import shlex
import subprocess
import sys
from types import SimpleNamespace

import magestic
import pandas as pd


def test_shipped_guide_link_arguments_produce_expected_link(tmp_path):
    scripts = Path(magestic.__file__).parent/'pipelines/bc1_donor_bc0/snakemake'
    rule = (scripts/'Snakefile.rules').read_text().split('rule link_03_guide_donor_bc0:', 1)[1].split('rule map_04_oligos:', 1)[0]
    shell_literal = re.search(r'shell\(("""[\s\S]*?""")\)', rule).group(1)
    shell_template = ast.literal_eval(shell_literal)
    source, output = tmp_path/'collapsed.tsv', tmp_path/'linked.tsv'
    guide_dir = tmp_path/'gdb'
    guide_dir.mkdir()
    pd.DataFrame([dict(bc1='A'*20, donor='C'*129, bc0='G'*20, total_counts=13)]).to_csv(source, sep='\t', index=False)
    pd.DataFrame([dict(library_ID='V_test', guide='T'*20,
                       donor_bc0_fragment='C'*30+'A'*20+'G'*10,
                       total_counts=21, SPS='A'*20)]).to_csv(
                           guide_dir/'V_test_guide_donor_bc0_counts.tsv', sep='\t', index=False)
    command = shell_template.format(PYTHON=sys.executable,
        params=SimpleNamespace(script=scripts/'03_link_guide_donor_bc0.py', guide_donor_bc0_dir=guide_dir),
        input=SimpleNamespace(collapsed=source), output=SimpleNamespace(linked=output), log=tmp_path/'link.log')
    argv = shlex.split(command.split('2>&1', 1)[0].replace('\\\n', ' '))
    result = subprocess.run(argv, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout+'\n'+result.stderr
    observed = pd.read_csv(output, sep='\t')
    assert len(observed) == 1
    assert observed.iloc[0].bc1 == 'A'*20
    assert observed.iloc[0].guide == 'T'*20
    assert observed.iloc[0].total_counts == 13
    assert observed.iloc[0].guide_link_counts == 21
    assert observed.iloc[0].guide_link_n_guides == 1
    sidecar = json.loads(Path(str(output)+'.provenance.json').read_text())
    assert sidecar['counts']['links_guide_linked'] == 1
    assert sidecar['counts']['links_in'] == sidecar['counts']['links_out'] == 1
    assert all(row.get('md5') for row in sidecar['inputs'])
