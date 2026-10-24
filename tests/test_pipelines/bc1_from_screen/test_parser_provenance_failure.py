"""An empty extraction is valid; missing required provenance is not completion."""
import gzip
import importlib.util
import json
from pathlib import Path
import sys

import pandas as pd
import pytest
import magestic
from magestic.utils.provenance import StepProvenance


@pytest.mark.parametrize("sidecar_fails", [False, True])
def test_empty_parser_completion_requires_sidecar(tmp_path, monkeypatch, sidecar_fails):
    source = Path(magestic.__file__).parent / "pipelines/bc1_from_screen/snakemake/02_parse_bc1_reads.py"
    spec = importlib.util.spec_from_file_location("screen_parser_provenance_test", source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    fastq = tmp_path / "empty.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        handle.write("")
    output = tmp_path / "counts.tsv"
    monkeypatch.setattr(sys, "argv", [str(source), "--input", str(fastq), "--output", str(output), "--sample-name", "empty"])
    if sidecar_fails:
        def fail_write(self):
            raise OSError("test mandatory sidecar failure")
        monkeypatch.setattr(StepProvenance, "write", fail_write)
    status = module.main()
    # Retain recoverable counts regardless of sidecar failure; never certify it.
    counts = pd.read_csv(output, sep="\t")
    assert counts.empty and len(counts.columns) == 12
    assert list(counts.columns[:3]) == ["bc1", "sample_name", "counts"]
    sidecar = Path(str(output) + ".provenance.json")
    if sidecar_fails:
        assert status != 0
        assert not sidecar.exists()
    else:
        assert status == 0
        receipt = json.loads(sidecar.read_text())
        assert receipt["counts"]["total_reads"] == 0
        assert receipt["counts"]["lineages_written"] == 0
