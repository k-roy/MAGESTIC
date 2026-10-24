"""The packaged runner's downstream route (steps 10-18) may not report success without outputs.

UP-C012 / docs/ARCHITECTURE_POLICY.md section 7: the runner used to log "PIPELINE COMPLETED SUCCESSFULLY" with
zero outputs. `pipeline/downstream.py` runs the numbered scripts and passes a stage only when it exits 0 AND
every output it declares holds data. These tests drive that contract with stand-in scripts, so they say
nothing about the science of any step; the real-data run is recorded in planning/272.
"""
import json
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from magestic.pipelines.bc1_from_screen.pipeline import downstream as DS


def _route(tmp_path, **extra):
    inputs = tmp_path / "inputs"
    inputs.mkdir(exist_ok=True)
    for name in ("count_matrix", "sample_key"):
        (inputs / f"{name}.tsv").write_text("a\tb\n1\t2\n")
    extra.setdefault("split", {"transformation_scope": "pooled"})
    return DS.RouteConfig(screen="S", library="yLX", output_dir=tmp_path / "out",
                          inputs={"count_matrix": inputs / "count_matrix.tsv",
                                  "sample_key": inputs / "sample_key.tsv"}, **extra)


@pytest.fixture
def fake_scripts(tmp_path, monkeypatch):
    """Point the route at stand-in scripts; each test writes the body of 10_split_comparisons.py."""
    scripts = tmp_path / "scripts"
    scripts.mkdir()
    monkeypatch.setattr(DS, "SCRIPTS", scripts)

    def write(body):
        (scripts / DS.STAGES["split"]).write_text(textwrap.dedent(body))
    return write


SPLIT_WRITES = """
    import sys, pathlib
    out = pathlib.Path(sys.argv[sys.argv.index('--output-dir') + 1]); out.mkdir(parents=True, exist_ok=True)
    (out / 'comparison_manifest.tsv').write_text({payload!r})
"""


def test_a_stage_that_exits_zero_without_writing_fails(tmp_path, fake_scripts):
    fake_scripts("import sys; sys.exit(0)\n")
    with pytest.raises(DS.StageFailed, match="wrote no data"):
        DS.run_downstream(_route(tmp_path), ["split"])
    manifest = json.loads((tmp_path / "out" / "run_manifest.json").read_text())
    assert manifest["status"] == "failed" and manifest["stages"][0]["rc"] == 0


def test_a_header_only_table_is_not_an_output(tmp_path, fake_scripts):
    fake_scripts(SPLIT_WRITES.format(payload="library\tcomparison_id\n"))
    with pytest.raises(DS.StageFailed, match="wrote no data"):
        DS.run_downstream(_route(tmp_path), ["split"])


def test_a_failing_stage_reports_its_exit_code_and_log(tmp_path, fake_scripts):
    fake_scripts("import sys; print('step 10 refuses: no transformation scope'); sys.exit(3)\n")
    with pytest.raises(DS.StageFailed, match="exited 3") as err:
        DS.run_downstream(_route(tmp_path), ["split"])
    assert "no transformation scope" in str(err.value)


def test_a_stage_with_data_passes_and_the_manifest_counts_its_rows(tmp_path, fake_scripts):
    fake_scripts(SPLIT_WRITES.format(payload="library\tcomparison_id\nyLX\tc1\nyLX\tc2\n"))
    manifest = json.loads(DS.run_downstream(_route(tmp_path), ["split"]).read_text())
    assert manifest["status"] == "completed"
    (out,) = manifest["stages"][0]["outputs"]
    assert out["rows"] == 2 and out["bytes"] > 0
    assert "--pool-transformations" in manifest["stages"][0]["argv"]


def test_step10_has_no_default_scope(tmp_path, fake_scripts):
    fake_scripts(SPLIT_WRITES.format(payload="x\ny\n"))
    with pytest.raises(DS.StageFailed, match="transformation_scope"):
        DS.run_downstream(_route(tmp_path, split={}), ["split"])


def test_a_named_comparison_absent_from_step10_is_refused(tmp_path, fake_scripts):
    fake_scripts(SPLIT_WRITES.format(payload="library\tcomparison_id\treference_file\ttreatment_file\n"
                                              "yLX\tc1\tc1_r.tsv\tc1_t.tsv\n"))
    route = _route(tmp_path, comparisons=["c1", "c9"])
    with pytest.raises(DS.StageFailed, match="c9"):
        DS.run_downstream(route, ["split", "deseq2"])


def test_unknown_steps_and_mixed_entry_points_are_usage_errors(tmp_path):
    config = tmp_path / "route.yaml"
    config.write_text("screen: S\nlibrary: yLX\noutput_dir: %s\ninputs: {}\n" % (tmp_path / "out"))
    cli = [sys.executable, "-m", "magestic.pipelines.bc1_from_screen.pipeline.run_pipeline"]
    bad_step = subprocess.run([*cli, "--route-config", config, "--steps", "split,no_such_step"],
                              capture_output=True, text=True)
    assert bad_step.returncode == 2 and "no_such_step" in bad_step.stderr
    mixed = subprocess.run([*cli, "--route-config", config, "--project-dir", tmp_path],
                           capture_output=True, text=True)
    assert mixed.returncode == 2 and "--route-config" in mixed.stderr


def test_the_route_config_refuses_unknown_keys(tmp_path):
    config = tmp_path / "route.yaml"
    config.write_text("screen: S\nlibrary: yLX\noutput_dir: out\ninputs: {}\nsplitt: {}\n")
    with pytest.raises(ValueError, match="splitt"):
        DS.RouteConfig.from_yaml(config)


def test_every_stage_names_a_shipped_script():
    for name, script in DS.STAGES.items():
        assert (Path(DS.__file__).resolve().parent.parent / "snakemake" / script).exists(), name


def test_keep_going_runs_later_stages_and_still_fails(tmp_path, fake_scripts):
    fake_scripts("import sys; sys.exit(1)\n")
    (DS.SCRIPTS / DS.STAGES["rarefaction"]).write_text(textwrap.dedent("""
        import sys, pathlib
        out = pathlib.Path(sys.argv[sys.argv.index('--output') + 1]); out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text('depth\\tn\\n10\\t3\\n'); out.with_suffix('.summary.tsv').write_text('depth\\tn\\n10\\t3\\n')
    """))
    table = tmp_path / "s1_bc1_counts.tsv"
    table.write_text("bc1\tcounts\nA\t1\n")
    route = _route(tmp_path, rarefaction={"samples": {"s1": str(table)}})
    with pytest.raises(DS.StageFailed, match="1 stage"):
        DS.run_downstream(route, ["split", "rarefaction"], keep_going=True)
    manifest = json.loads((tmp_path / "out" / "run_manifest.json").read_text())
    assert manifest["status"] == "failed" and list(manifest["failed_stages"]) == ["split"]
    assert [s["stage"] for s in manifest["stages"]] == ["split", "rarefaction"]
    assert manifest["stages"][1]["outputs"][0]["rows"] == 1
