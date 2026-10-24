"""Importing a utility or requesting its help must not create output trees."""
import os
import subprocess
import sys

import pytest


@pytest.mark.parametrize('help_requested', [False, True])
def test_import_and_help_have_no_output_directory_side_effect(tmp_path, help_requested):
    oak = tmp_path/'unused_output_root'
    env = dict(os.environ, OAK=str(oak))
    command = 'import magestic.utils.harmonize_oligos as h; '
    if help_requested:
        command += "import sys; sys.argv=['magestic-harmonize-oligos','--help']; h._cli()"
    else:
        command += "print('IMPORT_OK')"
    result = subprocess.run([sys.executable, '-c', command], env=env,
                            capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    assert ('usage:' if help_requested else 'IMPORT_OK') in result.stdout
    assert not oak.exists(), 'import/help created the configured output tree'
