"""Stage 4 end to end on windows that declare NO intended allele (planning/268).

Kevin, 2026-10-24: "we should not require an intended edit in order to be able to
assign confidence to variant calls from the WGS". So a window without
`intended_variants` is scored from its read evidence (`core.edit_profile`), gated at
the designed base that frame-verified `wt_seq` + `design_donor` give. With no designed
base the window is REFUSED (Kevin 2026-10-25, ruling 4 of planning/284 §8; WITHDRAWN:
"and at the window minimum otherwise"). These tests drive that route through
`run_step_4`, so the designed-position wiring, the provenance flag in the shipped shard
and both refusals are exercised on the path the pipeline takes.

Fixtures come from the 205 lane's `test_reconcile_functional` (synthetic molecules, no
external data). Reference chrII:301-320 is `A`; the design is A>C at 301.
"""
from dataclasses import replace

import pytest

from magestic.pipelines.wgs.core.reconcile_driver import (
    FLAG_WINDOW_EVIDENCE,
    DesignedPositionsRequired,
)
from magestic.pipelines.wgs.pipeline import step_1_bc1_call as s1
from magestic.pipelines.wgs.pipeline import step_4_reconcile as s4
from magestic.pipelines.wgs.pipeline._common import StepIncompleteError

from . import test_reconcile_functional as rf

fixture_factory = rf.fixture_factory  # the 205 lane's synthetic molecule builder
read, read_rows, write_bam = rf.read, rf.read_rows, rf.write_bam

WT_SEQ = "A" * 20
DESIGN = "C" + "A" * 19


def _case(f, reads, *, wt_seq=True, **reconcile):
    """One well, explicit assignment, no plasmid track and no bridge, as in the 205 tests."""
    if wt_seq:
        f["windows"].write_text(
            "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\twt_seq\tdesign_donor\n"
            f"oligo_1\tchrII\t301\t320\t{WT_SEQ}\t{DESIGN}\n")
    else:
        f["windows"].write_text(
            "oligo_name\tchrom\tdonor_start_coord\tdonor_end_coord\noligo_1\tchrII\t301\t320\n")
    write_bam(f["bam"], f["refs"], f["records"] + reads)
    assignments = f["root"] / "assigned.tsv"
    assignments.write_text("sample\toligo_name\nwell_A1\toligo_1\n")
    config = replace(f["config"], reconcile=replace(
        f["config"].reconcile, oligo_assignments=assignments, **reconcile))
    s1.run_step_1(config)
    return config


@pytest.mark.parametrize("seqs,label,status", [
    (["C" + "A" * 19] * 4, "SINGLE_EDIT", "full"),
    (["A" * 20] * 4, "NO_EDIT", "wt"),
])
def test_an_undeclared_window_gets_a_confident_call(fixture_factory, seqs, label, status):
    f = fixture_factory(bridge=False, plasmid=False)
    config = _case(f, [read(f"r{i}", s) for i, s in enumerate(seqs)])
    row = read_rows(s4.run_step_4(config))[0]
    assert (row["label"], row["confidence"]) == (label, "high")
    assert row["edits"].endswith(":" + status)
    assert FLAG_WINDOW_EVIDENCE in row["flags"].split(",")
    assert row["allele_counts"] == "[]"  # no declared allele, so no per-allele counts
    assert row["n_unassessed"] == "0"


def test_a_hole_at_the_designed_base_is_low_conf_not_no_edit(fixture_factory):
    """The rest of the window is deep; only the designed base (301) is unread."""
    f = fixture_factory(bridge=False, plasmid=False)
    config = _case(f, [read(f"r{i}", "A" * 19, start=301) for i in range(6)])
    row = read_rows(s4.run_step_4(config))[0]
    assert (row["label"], row["confidence"]) == ("LOW_CONF", "low")
    assert row["edits"].endswith(":nocov") and row["n_nocov"] == "1"


def test_no_wt_seq_is_refused_by_default_and_by_score_sites_when_accepted(fixture_factory):
    """A wrong frame inverts read-evidence calls, so `require_frame_check` refuses a table
    with no `wt_seq` while any window declares no allele. `False` passes that check and
    records the gap in the audit sidecar, but the scorer then refuses the window, because
    no designed base can be derived (ruling 4). WITHDRAWN: "`False` runs, gated at the
    window minimum"."""
    f = fixture_factory(bridge=False, plasmid=False)
    edits = [read(f"r{i}", "C" + "A" * 19) for i in range(4)]
    config = _case(f, edits, wt_seq=False)
    with pytest.raises(StepIncompleteError, match="no `wt_seq` column"):
        s4.run_step_4(config)
    g = fixture_factory(bridge=False, plasmid=False)
    accepted = _case(g, edits, wt_seq=False, require_frame_check=False)
    with pytest.raises(DesignedPositionsRequired, match="oligo_1"):
        s4.run_step_4(accepted)
    sidecar = s4.audit_sidecar_path(accepted, "well_A1").read_text()
    assert '"status": "no_wt_seq_column"' in sidecar
