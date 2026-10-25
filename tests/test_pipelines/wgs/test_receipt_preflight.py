"""No WGS table without a truthful receipt (planning/268 §5; arm F, planning/284).

X-S1-PARTIAL-PUBLISH: on the plate-5031 rerun, a git-archive tree with no code sha ran
Stage 1 to completion, published a 96-row `bc1_calls.tsv`, and only then raised
`ProvenanceError`, so a complete-looking table stood with no receipt. X-PROV-ENV-SHA:
every plate-5031 receipt named the sha the driver exported, `9d654b8`, while other
trees ran, and nothing checked it. Every `run_step_N` now calls `preflight_provenance`
before it computes, and a receipt says whether its commit was verified.
"""
import json

import pytest

from magestic.pipelines.wgs.pipeline import step_1_bc1_call as s1
from magestic.pipelines.wgs.qc import provenance as P

from . import test_reconcile_functional as rf

fixture_factory = rf.fixture_factory  # the 205 lane's synthetic molecule builder

HEAD = "a" * 40


def _checkout(monkeypatch, *, commit=HEAD, dirty=False):
    monkeypatch.setattr(P, "git_state", lambda: {"repo": "/x/magestic", "commit": commit,
                                                 "dirty": dirty, "dirty_files": []})


def _no_checkout(monkeypatch):
    monkeypatch.setattr(P, "git_state", lambda: {"repo": None, "commit": None, "dirty": None,
                                                 "note": "not a git checkout"})


def test_no_sha_refuses_before_stage_1_publishes_anything(fixture_factory, monkeypatch):
    """BEFORE: the cohort table was written, then the receipt raised."""
    config = fixture_factory(bridge=False, plasmid=False)["config"]
    monkeypatch.delenv(P.BUILD_SHA_ENV, raising=False)
    monkeypatch.setattr(P, "_sha_from_file", lambda: "")
    _no_checkout(monkeypatch)
    with pytest.raises(P.ProvenanceError, match="X-S1-PARTIAL-PUBLISH"):
        s1.run_step_1(config)
    assert not s1.table_path(config).exists()
    assert not s1.shard_path(config, "well_A1").exists()


def test_the_same_run_publishes_with_its_receipt_when_the_sha_resolves(fixture_factory,
                                                                       monkeypatch):
    """Positive control: the refusal above is the sha, not the fixture."""
    config = fixture_factory(bridge=False, plasmid=False)["config"]
    monkeypatch.setenv(P.BUILD_SHA_ENV, HEAD)
    _no_checkout(monkeypatch)
    out = s1.run_step_1(config)
    receipt = json.loads(open(str(out) + P.SIDECAR_SUFFIX).read())
    assert receipt["code"]["commit"] == HEAD
    assert receipt["code"]["commit_source"] == P.BUILD_SHA_ENV
    assert receipt["code"]["commit_verified"] is False  # no checkout to check it against


def test_an_exported_sha_that_is_not_the_checkout_head_is_refused(monkeypatch):
    monkeypatch.setenv(P.BUILD_SHA_ENV, "9d654b8")
    _checkout(monkeypatch)
    with pytest.raises(P.ProvenanceError, match="X-PROV-ENV-SHA"):
        P.preflight_provenance("step_1_bc1_call")


@pytest.mark.parametrize("exported,dirty,verified", [
    (HEAD, False, True),       # the full sha of a clean checkout
    (HEAD[:7], False, True),   # an abbreviated sha of the same commit
    (HEAD, True, False),       # right commit, uncommitted changes: not verified
])
def test_an_exported_sha_is_verified_against_the_checkout(monkeypatch, exported, dirty, verified):
    monkeypatch.setenv(P.BUILD_SHA_ENV, exported)
    _checkout(monkeypatch, dirty=dirty)
    assert P.preflight_provenance("step_4_reconcile") == exported
    state = P.code_state()
    assert (state["commit"], state["commit_verified"]) == (exported, verified)
    assert state["checkout_commit"] == HEAD


def test_a_checkout_sha_records_its_own_dirty_state(monkeypatch):
    monkeypatch.delenv(P.BUILD_SHA_ENV, raising=False)
    monkeypatch.setattr(P, "_sha_from_file", lambda: "")
    _checkout(monkeypatch, dirty=True)
    state = P.code_state()
    assert (state["commit_source"], state["commit_verified"]) == ("git_state", False)
