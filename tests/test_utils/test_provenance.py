"""StepProvenance writes the policy's minimum sidecar fields at run time."""
import json

from magestic.utils import provenance as pv


def test_sidecar_fields(tmp_path):
    inp = tmp_path / "in.tsv"
    inp.write_text("a\tb\n1\t2\n")
    out = tmp_path / "out.tsv"
    out.write_text("x\n")
    p = pv.StepProvenance(pipeline="t", step="s", output=out, argv=["prog", "--x", "1"])
    p.add_input(inp, role="in").set_parameters(threshold=3).set_counts(rows_in=1, rows_out=1)
    p.add_output(out, rows=1).add_decision(date="2026-09-10", who="Kevin", source="chat", quote="do it right")
    path = p.write()
    d = json.load(open(path))
    assert path == str(out) + ".provenance.json"
    assert d["inputs"][0]["md5"] == pv.md5_file(inp) and d["inputs"][0]["md5_kind"] == "raw"
    assert d["outputs"][0]["md5"] == pv.md5_file(out) and d["outputs"][0]["rows"] == 1
    assert d["command"]["argv"] == ["prog", "--x", "1"]
    assert {"commit", "dirty"} <= set(d["code"])
    assert {"host", "python", "packages"} <= set(d["environment"])
    assert d["parameters"] == {"threshold": 3}
    assert d["counts"] == {"rows_in": 1, "rows_out": 1}
    assert d["decisions"][0]["date"] == "2026-09-10" and d["decisions"][-1]["quote"] == "do it right"
    assert d["started"] <= d["finished"]


def test_run_provenance_indexes_step_sidecars(tmp_path):
    out = tmp_path / "o.tsv"
    out.write_text("x\n")
    sc = pv.StepProvenance(pipeline="t", step="s", output=out).add_output(out, rows=1).write()
    run = pv.write_run_provenance(tmp_path, pipeline="t", run_name="r", step_sidecars=[sc])
    d = json.load(open(run))
    assert d["steps"][0]["sidecar"] == sc and d["steps"][0]["outputs"][0]["rows"] == 1
    assert d["steps"][0]["sidecar_md5"] == pv.md5_file(sc)
