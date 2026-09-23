"""Two review findings on patches 40 and 41 that do not depend on the link-scope ruling (218 patch 44).

F3 (correctness panel): ``resolve_donor_conflicts`` checked only ``donor_bc0_fragment``, which holds the
smallest fragment a link was observed with, so a fragment that disagrees with the donor passed whenever
an agreeing one sorted first. Every observed fragment is now checked.

F6 (correctness panel): aggregation's total PCR support was not floored at its per-source counts. An
exact minimum cover never shrinks when identities are added, but the search's fallback lower bound can,
so the total could fall below what one source alone proves. merge already floors it.

Reviews: planning/218_magestic_production_corrections/reviews/merge_regressions_repair/panel_correctness/.
"""
import json

import pandas as pd

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources
from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter
from magestic.pipelines.bc1_donor_bc0.readers import data_loader
from magestic.pipelines.bc1_donor_bc0.readers.data_loader import build_pcr_replicate_id

SPS = 'CGATAGACGACTGGACAGCA'
BC1, BC0 = 'G' * 20, 'A' * 10
DONOR = 'T' * 100 + 'ACGTACGTACGTACGTACGTACGTACGTAC'
GOOD = DONOR[-30:] + SPS + BC0          # agrees with the donor's tail
SHORTER_GOOD = DONOR[-25:] + SPS + BC0  # also agrees, spelled with a shorter tail
BAD_AFTER = 'G' * 30 + SPS + BC0        # disagrees; sorts after GOOD
BAD_FIRST = 'A' * 30 + SPS + BC0        # disagrees; sorts before GOOD


def _row(count, pcr, fragment):
    return dict(counts=count, bc1=BC1, bc0=BC0, donor=DONOR, SPS=SPS, donor_bc0_fragment=fragment,
                pcr_replicate_id=pcr, trafo_ID='trafo_3')


def test_a_disagreeing_fragment_is_flagged_when_an_agreeing_one_sorts_first():
    merged = merge_data_sources(df_2x300=pd.DataFrame([_row(10, 'p1', GOOD), _row(10, 'p2', BAD_AFTER)]))
    assert len(merged) == 1
    assert set(json.loads(merged.iloc[0]['observed_fragments'])) == {GOOD, BAD_AFTER}
    assert bool(merged.iloc[0]['has_donor_conflict'])


def test_a_disagreeing_fragment_that_sorts_first_is_flagged():
    """Control: this case was already caught, because the checked fragment was the bad one."""
    merged = merge_data_sources(df_2x300=pd.DataFrame([_row(10, 'p1', GOOD), _row(10, 'p2', BAD_FIRST)]))
    assert bool(merged.iloc[0]['has_donor_conflict'])


def test_agreeing_fragments_of_different_lengths_are_not_a_conflict():
    """Control: two spellings of the donor's own tail agree with it."""
    merged = merge_data_sources(df_2x300=pd.DataFrame([_row(10, 'p1', GOOD), _row(10, 'p2', SHORTER_GOOD)]))
    assert not bool(merged.iloc[0]['has_donor_conflict'])


def _id(**fields):
    return build_pcr_replicate_id(dict(sample_number=1, **fields))


# A five-cycle of compatible identities (cover 3, independence number 2) plus W, which declares only the
# well and so is compatible with all five; with the exact-search size cut to 5 the union falls back.
CYCLE = [_id(gDNA_plate_name='1', plate_polymerase='2'), _id(inner_primers='1', condition='2'),
         _id(gDNA_plate_name='2', outer_primers='1'), _id(inner_primers='2', plate_polymerase='1'),
         _id(outer_primers='2', condition='1')]
W = _id()


def test_aggregation_total_is_never_below_a_per_source_count(monkeypatch):
    monkeypatch.setattr(data_loader, '_PCR_COVER_EXACT_MAX', 5)
    variant = DONOR[:50] + 'G' + DONOR[51:]
    full = pd.DataFrame([dict(counts=5, bc1=BC1, bc0=BC0, donor=DONOR, SPS=SPS, pcr_replicate_id=pcr)
                         for pcr in CYCLE[:3]] +
                        [dict(counts=2, bc1=BC1, bc0=BC0, donor=variant, SPS=SPS, pcr_replicate_id=pcr)
                         for pcr in CYCLE[3:]])
    short = pd.DataFrame([dict(counts=1, bc1=BC1, bc0=BC0, donor=variant, SPS=SPS, pcr_replicate_id=W,
                               donor_bc0_fragment=variant[-30:] + SPS + BC0)])
    result = aggregate_and_filter(merge_data_sources(df_2x150=short, df_2x300=full), min_reads=1,
                                  min_pcr_replicates=3, collapse_sequences=True, show_progress=False)
    assert len(result) == 1
    row = result.iloc[0]
    assert row['pcr_replicates_2x300'] == 3
    assert row['num_PCR_replicates'] >= row['pcr_replicates_2x300']
    assert row['support_tier'] == 'confident'
