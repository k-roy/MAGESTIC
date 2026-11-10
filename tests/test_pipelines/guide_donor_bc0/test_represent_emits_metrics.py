"""Step 05 ("represent") must EMIT the representation metrics on a normal run.

The wiring in ``run_library_representation`` is deliberately guarded: a failure inside the metrics
block logs ERROR and writes ``representation_metrics.FAILED.txt`` rather than failing the step,
because the metrics are additive and no existing output depends on them. 🔴 That guard is also the
way the wiring could rot silently -- a renamed argument or a bad import would leave the step
reporting success with no metrics on disk. So these tests run the real step function on a synthetic
project and assert the positive evidence: the metrics file exists, the FAILED marker does not, and
the pre-existing step-05 outputs are still written exactly as before.
"""

import json
import logging
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from magestic.pipelines.guide_donor_bc0.core.library_representation import (
    METRICS_FAILED_NAME,
    METRICS_JSON_NAME,
    METRICS_TSV_NAME,
)
from magestic.pipelines.guide_donor_bc0.pipeline.run_pipeline import (
    run_library_representation,
)

PURITY_SUFFIX = "_guide_donor_bc0_purity_filtered.tsv"


def _make_project(tmp_path, *, n_oligos=12, library="V629"):
    """A minimal but real step-05 input set: one purity file and one oligo design file."""
    output_dir = tmp_path / "processed_data" / "guide_donor_bc0"
    purity_dir = output_dir / "guide_donor_bc0_purity_filtered"
    purity_dir.mkdir(parents=True)

    # A skewed count distribution with a singleton tail, which is what makes the floors interesting
    counts = [1, 1, 1, 1, 2, 2, 3, 5, 9, 21, 60, 400][:n_oligos]
    names = [f"oligo_{i:03d}" for i in range(len(counts))]
    oligo_seqs = [
        ("A" * 51) + ("ACGT" * 32)[:129] + "CTTGACTCGACAGTGAGAGC" for _ in names
    ]

    pd.DataFrame({
        "library_ID": library,
        "bc0": [f"BC{i:04d}" for i in range(len(counts))],
        "guide": ["ACGTACGTACGTACGTACGT"] * len(counts),
        "donor": [s[51:180] for s in oligo_seqs],
        "oligo_name": names,
        "grouped_guide_bc0_counts": counts,
    }).to_csv(purity_dir / f"{library}{PURITY_SUFFIX}", sep="\t", index=False)

    design = tmp_path / "designs_Twist_200mer.tsv"
    pd.DataFrame({"oligo_name": names, "oligo_seq": oligo_seqs}).to_csv(
        design, sep="\t", index=False)

    config = SimpleNamespace(
        project_dir=tmp_path,
        output_dir=output_dir,
        oligo_design_file=design,
        library_sps_map={},
        use_harmonized_table=False,
        metrics_rarefy_to=None,
    )
    return config, purity_dir / "summary_plots"


def test_represent_step_emits_metrics_on_a_normal_run(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    config, summary_dir = _make_project(tmp_path)

    assert run_library_representation(config) is True

    # the pre-existing step-05 outputs, unchanged
    assert (summary_dir / "all_libraries_oligo_name_filtered.tsv").exists()
    assert (summary_dir / "V629_oligo_name_filtered.tsv").exists()

    # 🔴 the positive evidence, and the absence of the silent-failure marker
    assert not (summary_dir / METRICS_FAILED_NAME).exists(), \
        "the metrics block raised; see the ERROR in the step log"
    assert (summary_dir / METRICS_TSV_NAME).exists()
    assert (summary_dir / METRICS_JSON_NAME).exists()
    sidecar = summary_dir / f"{METRICS_TSV_NAME}.provenance.json"
    assert sidecar.exists()
    # the guide_bc0 rows are computed from the purity file, so it must be a recorded input with its
    # md5: two count distributions can share one oligo total and so one oligo-table checksum
    recorded = {Path(rec["path"]).name: rec for rec in json.loads(sidecar.read_text())["inputs"]}
    purity_record = recorded[f"V629{PURITY_SUFFIX}"]
    assert purity_record["role"] == "purity_filtered_counts" and purity_record["md5"]
    assert "Representation metrics FAILED" not in caplog.text


def test_emitted_metrics_table_carries_the_depth_with_every_value(tmp_path):
    config, summary_dir = _make_project(tmp_path)
    assert run_library_representation(config) is True

    table = pd.read_csv(summary_dir / METRICS_TSV_NAME, sep="\t")
    assert len(table) > 0
    for col in ("library", "grain", "metric", "value", "depth", "min_count",
                "n_observations", "universe", "weighting", "rarefied"):
        assert col in table.columns, col
    # no row may carry a value without the context needed to quote it
    assert table["depth"].notna().all()
    assert table["min_count"].notna().all()
    assert table["n_observations"].notna().all()
    # both grains are present, because a floor means a different thing at each
    assert {"oligo", "guide_bc0"} <= set(table["grain"])
    # several floors, not one
    assert len(set(table["min_count"])) >= 4
    # and the statistics asked for are all there
    assert {"gini", "shannon_entropy", "pielou_evenness", "simpson_index",
            "inverse_simpson", "hill_q1", "hill_q2"} <= set(table["metric"])
    # un-rarefied by default, hence sample-specific
    assert not table["rarefied"].any()


def test_the_step_logs_the_gini_together_with_its_depth(tmp_path, caplog):
    caplog.set_level(logging.INFO)
    config, summary_dir = _make_project(tmp_path)
    run_library_representation(config)
    assert "gini=" in caplog.text
    assert "at depth=" in caplog.text
    assert "NOT rarefied: sample-specific" in caplog.text


def test_rarefy_to_is_reachable_from_the_config(tmp_path):
    """The knob is a real config field, not a getattr that can never be set."""
    from magestic.pipelines.guide_donor_bc0.config import GuideDonorBC0Config

    assert "metrics_rarefy_to" in GuideDonorBC0Config.__dataclass_fields__
    assert GuideDonorBC0Config.__dataclass_fields__["metrics_rarefy_to"].default is None

    config, summary_dir = _make_project(tmp_path)
    config.metrics_rarefy_to = 100
    assert run_library_representation(config) is True
    table = pd.read_csv(summary_dir / METRICS_TSV_NAME, sep="\t")
    assert table["rarefied"].any()
    assert set(table.loc[table["rarefied"], "rarefied_to"]) == {100}
    assert not (summary_dir / METRICS_FAILED_NAME).exists()


def test_a_failing_metrics_block_marks_itself_and_does_not_fail_the_step(tmp_path, monkeypatch,
                                                                        caplog):
    """The guard, tested from both sides: the step still succeeds, and it leaves evidence."""
    caplog.set_level(logging.ERROR)
    config, summary_dir = _make_project(tmp_path)

    import magestic.pipelines.guide_donor_bc0.core.library_representation as lr

    def boom(*a, **kw):
        raise RuntimeError("synthetic metrics failure")

    monkeypatch.setattr(lr, "compute_representation_metrics", boom)
    assert run_library_representation(config) is True
    assert (summary_dir / METRICS_FAILED_NAME).exists()
    assert "synthetic metrics failure" in (summary_dir / METRICS_FAILED_NAME).read_text()
    assert "Representation metrics FAILED" in caplog.text
    # and the pre-existing outputs are still written
    assert (summary_dir / "all_libraries_oligo_name_filtered.tsv").exists()
