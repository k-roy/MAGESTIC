"""Actual05->06 boundary must support an ordinary run with no detected errors."""
import importlib.util
import json
from pathlib import Path
import sys

import magestic
import pandas as pd
import pytest


def load_step(filename):
    path = Path(magestic.__file__).parent / "pipelines/bc1_from_screen/snakemake" / filename
    spec = importlib.util.spec_from_file_location(filename.replace(".", "_"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("with_error", [False, True])
def test_actual_error_mapping_and_stats_are_readable_and_provenanced(tmp_path, monkeypatch, with_error):
    parent, child = "A" * 20, "A" * 19 + "T"
    rows = [dict(bc1=parent, sample_name="S", counts=500)]
    if with_error:
        rows.append(dict(bc1=child, sample_name="S", counts=5))
    counts = tmp_path / "counts.tsv"
    pd.DataFrame(rows).to_csv(counts, sep="\t", index=False)
    reference = tmp_path / "reference.tsv"
    pd.DataFrame(dict(bc1=[parent])).to_csv(reference, sep="\t", index=False)
    output, stats = tmp_path / "errors.tsv", tmp_path / "stats.tsv"
    monkeypatch.setattr(sys, "argv", ["05", "--counts", str(counts), "--bc1-reference", str(reference),
                                      "--output", str(output), "--stats-output", str(stats)])
    load_step("05_error_detection.py").main()
    mapping = pd.read_csv(output, sep="\t")
    assert list(mapping.columns) == ["error_bc1", "parent_bc1", "hamming_distance", "error_counts",
                                      "parent_counts", "count_ratio", "error_type"]
    assert len(mapping) == int(with_error)
    if with_error:
        assert mapping.iloc[0]["error_bc1"] == child
        assert mapping.iloc[0]["parent_bc1"] == parent
        assert mapping.iloc[0]["count_ratio"] == 100
    cleaned = tmp_path / "cleaned.tsv"
    monkeypatch.setattr(sys, "argv", ["06", "--counts", str(counts), "--error-mappings", str(output),
                                      "--output", str(cleaned)])
    load_step("06_absorb_error_counts.py").main()
    observed = pd.read_csv(cleaned, sep="\t")
    assert observed.to_dict("records") == [dict(bc1=parent, sample_name="S", counts=500 + 5 * with_error)]
    for path in (output, stats):
        receipt = json.loads(Path(str(path) + ".provenance.json").read_text())
        assert len(receipt["inputs"]) == 2
        assert {Path(row["path"]).name for row in receipt["inputs"]} == {"counts.tsv", "reference.tsv"}
        assert receipt["parameters"]["min_count_ratio"] == 100
        assert receipt["parameters"]["max_error_count"] == 50
        assert receipt["parameters"]["max_error_samples"] == 3
        assert receipt["counts"]["input_reads"] == 500 + 5 * with_error
        assert receipt["counts"]["output_rows"] == len(pd.read_csv(path, sep="\t"))


def test_empty_analysis_is_refused_explicitly_before_publishing(tmp_path, monkeypatch):
    counts, reference, output = tmp_path / "counts.tsv", tmp_path / "reference.tsv", tmp_path / "errors.tsv"
    counts.write_text("bc1\tsample_name\tcounts\n")
    reference.write_text("bc1\n" + "A" * 20 + "\n")
    monkeypatch.setattr(sys, "argv", ["05", "--counts", str(counts), "--bc1-reference", str(reference), "--output", str(output)])
    with pytest.raises(ValueError, match="No count rows"):
        load_step("05_error_detection.py").main()
    assert not output.exists()
