"""Regression contracts for the NNS-set WGS fixes (SGTC planning/328, 363, 385).

Every TRACE fix made on the NRD1/NAB3/SEN1 NNS cohort carries its own unit test in TRACE. Two things were not
pinned down on the magestic side:

1. that the TRACE magestic pins actually carries those fixes, so a pin bump backwards (or an environment that
   imports a stale TRACE, as Sherlock's `software/trace` did at 7c4c10f in planning/363b section 8) fails here
   rather than silently running the aligned-base reader; and
2. `50a48cb` (Stage 4 passes `reconcile.allele_evidence` through to TRACE), which shipped without a test.

Failing-first evidence (planning/385): at magestic `1b10956` (TRACE pin 12b35c5) the pin contract fails on the
missing `local_realign` / `target_calls` modules and the `evidence` parameter; at HEAD it passes.
"""
import dataclasses
import inspect

import pytest
import yaml

from . import test_reconcile_functional as rf

fixture_factory = rf.fixture_factory  # the 205 lane's synthetic molecule builder


def test_pinned_trace_carries_every_nns_fix():
    from trace_crispr.wgs import reconcile as tr
    from trace_crispr.wgs import sv_plasmid as sv
    # 9209016 / 6f5e439: local realignment to the site haplotypes is the DEFAULT reader
    from trace_crispr.wgs import local_realign  # noqa: F401
    assert set(tr.EVIDENCE_MODES) == {"realign", "aligned"}
    params = inspect.signature(tr.score_edit_site).parameters
    assert params["evidence"].default == "realign"
    # a0e8ec0 / 4279260: v7 haplotype-level target calls
    from trace_crispr.wgs import target_calls
    assert {"calling", "target", "contamination_rate"} <= set(params)
    assert hasattr(target_calls, "CE_UNCONFIRMED_TARGET") and hasattr(target_calls, "HELPERS_ONLY")
    # 620b7c1: Half B reads the mate CIGAR (MC:Z) for a dovetailed fragment
    assert "mate_cigar" in {f.name for f in dataclasses.fields(sv.SamRecord)}
    # 4ea1861 / 3e373c6: contained fragments kept when the plasmids are purged; plasmid load bound
    assert "keep_contained" in inspect.signature(sv.classify_donor_reads).parameters
    assert callable(sv.plasmid_load) and callable(sv.plasmids_purged)


class _CapturingTrace:
    """Only what `score_sites` touches on the intended-allele route."""

    SITE_PARTIAL = "partial"
    SITE_NOCOV = "nocov"
    SITE_UNASSESSED = "unassessed"

    def __init__(self, site):
        self._site = site
        self.calls = []

    def build_edit_sites(self, oligos, windows, fasta=None):
        return [self._site]

    def score_edit_site(self, bam, fasta, site, **kwargs):
        from trace_crispr.wgs.reconcile import SiteCall
        self.calls.append(kwargs)
        return SiteCall(site=site, vaf=0.0, cov=0, status="nocov")


@pytest.mark.parametrize("mode,expected", [
    (None, None),            # omitted: TRACE's own default (realign) decides, nothing is forced
    ("aligned", "aligned"),  # parity with pre-328 runs
    ("realign", "realign"),
])
def test_stage4_driver_passes_allele_evidence_through(fixture_factory, monkeypatch, mode, expected):
    from magestic.pipelines.wgs.core import edit_profile
    import magestic.pipelines.wgs.core.reconcile_driver as mod
    from trace_crispr.wgs.reconcile import EditSite, IntendedAllele
    import pysam

    f = fixture_factory(records=[])
    fake = _CapturingTrace(EditSite("chrII", 301, 310, "oligo_1", (IntendedAllele(301, "A", "C"),)))
    monkeypatch.setattr(mod._trace, "reconcile", lambda: fake)
    monkeypatch.setattr(edit_profile, "one_sided_coverage_loss", lambda *a, **k: None)
    with pysam.AlignmentFile(str(f["bam"])) as bam, pysam.FastaFile(str(f["fa"])) as fa:
        mod.score_sites(bam, fa, ["oligo_1"], {"oligo_1": None}, allele_evidence=mode)
    assert len(fake.calls) == 1
    if expected is None:
        assert "evidence" not in fake.calls[0]
    else:
        assert fake.calls[0]["evidence"] == expected


@pytest.mark.parametrize("value,ok", [("realign", True), ("aligned", True), ("pileup", False)])
def test_allele_evidence_config_values(fixture_factory, value, ok):
    from magestic.pipelines.wgs.config import load_config
    f = fixture_factory(records=[])
    data = yaml.safe_load(f["config_path"].read_text())
    data.setdefault("reconcile", {})["allele_evidence"] = value
    path = f["root"] / f"ae_{value}.yaml"
    path.write_text(yaml.safe_dump(data))
    if ok:
        assert load_config(path).reconcile.allele_evidence == value
    else:
        with pytest.raises(ValueError, match="realign or aligned"):
            load_config(path)
