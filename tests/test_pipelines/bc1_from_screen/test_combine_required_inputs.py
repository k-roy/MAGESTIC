"""Combining counts must account for every discovered input, including empty wells."""
import hashlib
import importlib.util
import json
from pathlib import Path

import magestic
import pandas as pd
import pytest
from magestic.utils.bc1_extraction import EMITTING_TIERS


def load_combiner():
    path = Path(magestic.__file__).parent / "pipelines/bc1_from_screen/snakemake/04_combine_counts.py"
    spec = importlib.util.spec_from_file_location("screen_combine_required_inputs", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valid_table(path, empty=False):
    row = dict(bc1="A" * 20, sample_name="S", counts=12)
    row.update({f"n_{tier}": 12 if tier == "confirmed" else 0 for tier in EMITTING_TIERS})
    row.update(n_fwd=12, n_rev=0)
    pd.DataFrame([] if empty else [row], columns=list(row)).to_csv(path, sep="\t", index=False)


@pytest.mark.parametrize("failure", ["unreadable", "missing_column", "empty_missing_column"])
def test_bad_discovered_input_does_not_publish_partial_counts(tmp_path, failure):
    incoming = tmp_path / "input"
    incoming.mkdir()
    valid_table(incoming / "good_bc1_counts.tsv")
    bad = incoming / "bad_bc1_counts.tsv"
    if failure == "unreadable":
        bad.write_bytes(b"\xff\xfeinvalid")
    elif failure == "missing_column":
        bad.write_text("bc1\tsample_name\n" + "C" * 20 + "\tT\n")
    else:
        bad.write_text("bc1\tsample_name\n")
    output = tmp_path / "combined.tsv"
    sidecar = Path(str(output) + ".provenance.json")
    output.write_text("previous counts\n")
    sidecar.write_text("previous receipt\n")
    with pytest.raises(ValueError, match="bad_bc1_counts.tsv"):
        load_combiner().combine_count_files(incoming, output, n_workers=1)
    assert output.read_text() == "previous counts\n"
    assert sidecar.read_text() == "previous receipt\n"


def test_valid_and_empty_inputs_preserve_counts_and_record_both_files(tmp_path):
    incoming = tmp_path / "input"
    incoming.mkdir()
    valid_table(incoming / "good_bc1_counts.tsv")
    valid_table(incoming / "empty_bc1_counts.tsv", empty=True)
    output = tmp_path / "combined.tsv"
    load_combiner().combine_count_files(incoming, output, n_workers=1)
    result = pd.read_csv(output, sep="\t")
    assert result.to_dict("records") == [dict(bc1="A" * 20, sample_name="S", counts=12)]
    receipt = json.loads(Path(str(output) + ".provenance.json").read_text())
    observed = {Path(row["path"]).name: row["md5"] for row in receipt["inputs"]}
    for path in incoming.iterdir():
        assert observed[path.name] == hashlib.md5(path.read_bytes()).hexdigest()
    assert receipt["counts"]["input_files"] == 2
    assert receipt["counts"]["empty_input_files"] == 1
    assert receipt["counts"]["rows_in"] == receipt["counts"]["rows_out"] == 1
    assert receipt["counts"]["reads_in"] == receipt["counts"]["reads_out"] == 12


def test_entirely_empty_collection_is_explicitly_refused(tmp_path):
    incoming = tmp_path / "input"
    incoming.mkdir()
    valid_table(incoming / "empty_bc1_counts.tsv", empty=True)
    output = tmp_path / "combined.tsv"
    with pytest.raises(ValueError, match="All count files were empty"):
        load_combiner().combine_count_files(incoming, output, n_workers=1)
    assert not output.exists()
