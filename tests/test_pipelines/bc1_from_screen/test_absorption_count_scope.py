"""Read mass must not depend on how a sample/barcode was split across input files."""
import importlib.util
import json
from pathlib import Path
import sys

import magestic
import pandas as pd
import pytest


def load_absorber():
    path = Path(magestic.__file__).parent / "pipelines/bc1_from_screen/snakemake/06_absorb_error_counts.py"
    spec = importlib.util.spec_from_file_location("absorption_count_scope", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("parent_counts", [[30], [10, 20]])
def test_absorption_is_invariant_to_split_parent_rows(parent_counts):
    parent, child = "A" * 20, "A" * 19 + "C"
    counts = pd.DataFrame([dict(bc1=parent, sample_name="S", counts=n) for n in parent_counts]
                          + [dict(bc1=child, sample_name="S", counts=3)])
    mapping = pd.DataFrame([dict(error_bc1=child, parent_bc1=parent)])
    result, log = load_absorber().absorb_error_counts(counts, mapping)
    assert result.to_dict("records") == [dict(bc1=parent, sample_name="S", counts=33)]
    assert log.iloc[0]["total_absorbed_counts"] == 3
    assert log.iloc[0]["n_error_bc1s"] == 1


@pytest.mark.parametrize("with_errors", [False, True])
def test_actual_absorber_writes_clean_counts_and_even_empty_log(tmp_path, monkeypatch, with_errors):
    parent, child = "A" * 20, "A" * 19 + "C"
    rows = [dict(bc1=parent, sample_name="S", counts=10), dict(bc1=parent, sample_name="S", counts=20)]
    if with_errors:
        rows.append(dict(bc1=child, sample_name="S", counts=3))
    counts, mappings = tmp_path / "counts.tsv", tmp_path / "errors.tsv"
    pd.DataFrame(rows).to_csv(counts, sep="\t", index=False)
    pd.DataFrame([dict(error_bc1=child, parent_bc1=parent)] if with_errors else [],
                 columns=["error_bc1", "parent_bc1"]).to_csv(mappings, sep="\t", index=False)
    output, log = tmp_path / "cleaned.tsv", tmp_path / "absorption_log.tsv"
    log.write_text("stale previous absorption log\n")
    monkeypatch.setattr(sys, "argv", ["06", "--counts", str(counts), "--error-mappings", str(mappings),
                                      "--output", str(output), "--log-output", str(log)])
    load_absorber().main()
    expected = 30 + 3 * with_errors
    assert pd.read_csv(output, sep="\t").to_dict("records") == [dict(bc1=parent, sample_name="S", counts=expected)]
    record = pd.read_csv(log, sep="\t")
    assert list(record.columns) == ["parent_bc1", "n_error_bc1s", "total_absorbed_counts",
                                    "n_sample_contributions", "error_bc1s"]
    assert len(record) == int(with_errors)
    for path in (output, log):
        receipt = json.loads(Path(str(path) + ".provenance.json").read_text())
        assert receipt["counts"]["reads_in"] == receipt["counts"]["reads_out"] == expected
        assert receipt["counts"]["absorbed_reads"] == 3 * with_errors
        assert len(receipt["inputs"]) == 2
        assert receipt["counts"]["rows_out"] == len(pd.read_csv(path, sep="\t"))


def test_absorption_keeps_samples_separate():
    parent, child = "A" * 20, "A" * 19 + "C"
    counts = pd.DataFrame([
        (parent, "S1", 10), (parent, "S1", 20), (child, "S1", 3),
        (parent, "S2", 7), (child, "S2", 2)], columns=["bc1", "sample_name", "counts"])
    mapping = pd.DataFrame([dict(error_bc1=child, parent_bc1=parent)])
    result, log = load_absorber().absorb_error_counts(counts, mapping)
    assert result.set_index("sample_name")["counts"].to_dict() == {"S1": 33, "S2": 9}
    assert result["bc1"].tolist() == [parent, parent]
    assert log.iloc[0]["total_absorbed_counts"] == 5
    assert log.iloc[0]["n_sample_contributions"] == 2


def test_absent_error_child_returns_normalized_counts_and_empty_log():
    counts = pd.DataFrame(dict(bc1=["A" * 20] * 2, sample_name=["S"] * 2, counts=[10, 20]))
    mapping = pd.DataFrame(dict(error_bc1=["C" * 20], parent_bc1=["A" * 20]))
    result, log = load_absorber().absorb_error_counts(counts, mapping)
    assert result["counts"].tolist() == [30]
    assert log.empty and "parent_bc1" in log.columns
