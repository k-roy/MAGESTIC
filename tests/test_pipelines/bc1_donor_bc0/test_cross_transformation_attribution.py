"""Each BDB link reports how its bc1's reads divide among real transformations (218 patch 46).

Kevin, 2026-09-23: "When comparing two separate trafos, we simply want to make sure that the bc1 is
observed in one trafo with much higher read counts/purity than in the other. Indeed, the birthday
statistics predicts very low, but not negligible, chance of collision."

A link lives in one transformation (patch 45). Aggregation adds, for every row, the number of real
transformations its bc1 is seen in, the reads per real transformation, and the row's own
transformation's share of them. A real transformation is a 'trafo_N' scope; the loaders map a plate
('trafo_N-M') and a trafo-key yL to it. Multi-transformation pools such as yL442, 'empty' wells and
undeclared rows are neither counted nor in the denominator. No threshold is applied: "much higher" is
Kevin's number, so nothing is dropped or relabelled here.

Design: planning/218_magestic_production_corrections/reviews/merge_regressions_repair/DESIGN.md §6 item 3.
"""
import json

import numpy as np
import pandas as pd

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources
from magestic.pipelines.bc1_donor_bc0.pipeline.run_pipeline import aggregate_and_filter

BC1, BC0, SPS = 'G' * 20, 'A' * 10, 'CGATAGACGACTGGACAGCA'
DONOR = 'T' * 100 + 'ACGTACGTACGTACGTACGTACGTACGTAC'
COLUMNS = ('bc1_real_transformations', 'bc1_transformation_reads', 'transformation_read_share')


def link(count, pcr, trafo, donor=DONOR, bc1=BC1):
    return dict(counts=count, bc1=bc1, bc0=BC0, donor=donor, SPS=SPS, pcr_replicate_id=pcr, trafo_ID=trafo)


def aggregate(rows, min_reads=1):
    merged = merge_data_sources(df_2x300=pd.DataFrame(rows))
    return aggregate_and_filter(merged, min_reads=min_reads, min_pcr_replicates=2,
                                collapse_sequences=False, show_progress=False).set_index('trafo_ID')


def test_a_bc1_in_two_transformations_reports_each_share():
    """20250912 shape: trafo_3 (SpG) and trafo_4 (SpCas9) both carry cV586."""
    result = aggregate([link(90, 'p1', 'trafo_3'), link(10, 'p2', 'trafo_4')])
    assert result['bc1_real_transformations'].tolist() == [2, 2]
    assert result.loc['trafo_3', 'transformation_read_share'] == 0.9
    assert result.loc['trafo_4', 'transformation_read_share'] == 0.1
    assert {json.loads(value)['trafo_3'] for value in result['bc1_transformation_reads']} == {90}
    assert json.loads(result.loc['trafo_4', 'bc1_transformation_reads']) == {'trafo_3': 90, 'trafo_4': 10}


def test_pools_empty_wells_and_undeclared_rows_are_not_transformations():
    result = aggregate([link(30, 'p1', 'trafo_3'), link(50, 'p2', 'yL442'), link(5, 'p3', 'empty'),
                        link(7, 'p4', None)])
    assert result['bc1_real_transformations'].tolist() == [1, 1, 1, 1]
    assert result.loc['trafo_3', 'transformation_read_share'] == 1.0
    assert np.isnan(result.loc['yL442', 'transformation_read_share'])
    assert np.isnan(result.loc['empty', 'transformation_read_share'])
    assert json.loads(result.loc['yL442', 'bc1_transformation_reads']) == {'trafo_3': 30}


def test_the_share_sums_every_link_of_the_bc1_in_that_transformation():
    """Two donors of one bc1 in trafo_3 (a competing link) against one in trafo_4."""
    result = aggregate([link(40, 'p1', 'trafo_3'), link(10, 'p2', 'trafo_3', donor='C' * 130),
                        link(50, 'p3', 'trafo_4')])
    assert result.loc['trafo_3', 'transformation_read_share'].tolist() == [0.5, 0.5]
    assert result.loc['trafo_4', 'transformation_read_share'] == 0.5


def test_reads_below_min_reads_stay_in_the_denominator():
    """Like the purity denominator: the 1-read trafo_3 row is dropped, its reads still count."""
    result = aggregate([link(1, 'p1', 'trafo_3'), link(9, 'p2', 'trafo_4')], min_reads=2)
    assert list(result.index) == ['trafo_4']
    assert result.loc['trafo_4', 'bc1_real_transformations'] == 2
    assert result.loc['trafo_4', 'transformation_read_share'] == 0.9


def test_other_bc1s_do_not_enter_the_share():
    result = aggregate([link(20, 'p1', 'trafo_1'), link(80, 'p2', 'trafo_2', bc1='C' * 20)])
    assert result['bc1_real_transformations'].tolist() == [1, 1]
    assert result['transformation_read_share'].tolist() == [1.0, 1.0]


def test_attribution_applies_no_threshold():
    """A 1% share keeps its row, support and tier; only the columns report it."""
    result = aggregate([link(99, 'p1', 'trafo_3'), link(99, 'p2', 'trafo_3'),
                        link(1, 'p3', 'trafo_4'), link(1, 'p4', 'trafo_4')])
    assert result.loc['trafo_4', 'transformation_read_share'] == 0.01
    assert result.loc['trafo_4', 'num_PCR_replicates'] == 2
    assert result.loc['trafo_4', 'support_tier'] == 'confident'


def test_empty_aggregation_has_the_attribution_columns():
    result = aggregate_and_filter(merge_data_sources(), min_reads=1, min_pcr_replicates=2,
                                  collapse_sequences=False, show_progress=False)
    assert result.empty and set(COLUMNS) <= set(result.columns)


def test_rows_without_a_transformation_column_report_none():
    rows = [link(10, 'p1', None), link(5, 'p2', None)]
    for row in rows:
        del row['trafo_ID']
    merged = merge_data_sources(df_2x300=pd.DataFrame(rows))
    result = aggregate_and_filter(merged, min_reads=1, min_pcr_replicates=2,
                                  collapse_sequences=False, show_progress=False)
    assert result['bc1_real_transformations'].tolist() == [0]
    assert np.isnan(result['transformation_read_share'].iloc[0])
    assert result['bc1_transformation_reads'].tolist() == ['{}']
