"""Rule R* for the donor-bc0 bridge: preserve the candidates AND tier them by support.

Kevin, 2026-10-24 (planning/259 section 7, decision 2): when one donor_bc0_fragment is seen with more
than one donor, use NN 267's rule R* (planning/222_triumph_embl_library_design/267_bc1_oligo_assignment_audit.md
sections 20.4 and 23): keep every candidate (no winner-takes-row collapse, no unique-only filter), assign
the support leader, and report perfect-donor reads plus a floor of 10 reads and a lead of 3x as a
CONFIDENCE TIER, not as a hard filter. R* measured 8 correct / 0 wrong against the shipped 6 / 14.

R* was defined for the bc1 -> oligo lookup; this file pins its translation to the fragment -> donor bridge
(planning/272 section 2). Neither merged lane's route does this: master collapsed each fragment to one
winner and then left-joined, and the 218 lane assigned only a unique candidate.
"""
import json

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core import donor_recovery as DR
from magestic.pipelines.bc1_donor_bc0.readers import data_loader as DL

FRAGMENT = 'C' * 30 + 'G' * 20 + 'A' * 10
D1, D2, D3 = 'A' * 99 + 'C' * 30, 'T' * 99 + 'C' * 30, 'G' * 99 + 'C' * 30


def observation(n=1, **extra):
    return pd.DataFrame([dict(bc1='G' * 20, bc0='A' * 10, donor_bc0_fragment=FRAGMENT,
                              counts=20, pcr_replicate_id=f'p{i}', **extra) for i in range(n)])


def lookup(support, **extra):
    """support: list of (donor, reads)."""
    return pd.DataFrame(dict(donor_bc0_fragment=[FRAGMENT] * len(support),
                             donor=[d for d, _ in support], bc0=['A' * 10] * len(support),
                             counts=[n for _, n in support], **extra))


def test_rstar_thresholds_are_the_measured_ones():
    assert (DR.RSTAR_FLOOR_READS, DR.RSTAR_LEAD_RATIO) == (10, 3)


def test_a_clear_support_leader_is_assigned_instead_of_refused():
    """The 218 route left this row unassigned ('ambiguous'); master's route took D1 but dropped D2."""
    row = DL.recover_full_donors(observation(), lookup([(D2, 5), (D1, 50)])).iloc[0]
    assert row['donor'] == D1
    assert row['donor_recovery_status'] == 'recovered_by_support'
    assert row['donor_recovery_tier'] == 'confident'
    assert (row['donor_recovery_lead_reads'], row['donor_recovery_runner_up_reads']) == (50, 5)
    # Every candidate is preserved, ranked by support.
    assert json.loads(row['donor_recovery_candidate_support']) == [[D1, 50], [D2, 5]]
    assert set(json.loads(row['donor_recovery_candidates'])) == {D1, D2}


def test_a_weak_lead_is_assigned_but_only_provisionally():
    row = DL.recover_full_donors(observation(), lookup([(D1, 12), (D2, 6)])).iloc[0]
    assert row['donor'] == D1 and row['donor_recovery_status'] == 'recovered_by_support'
    assert row['donor_recovery_tier'] == 'provisional'       # 12 < 3 x 6


def test_the_floor_is_a_tier_not_a_filter():
    row = DL.recover_full_donors(observation(), lookup([(D1, 5)])).iloc[0]
    assert row['donor'] == D1 and row['donor_recovery_status'] == 'recovered_unique'
    assert row['donor_recovery_tier'] == 'provisional'       # 5 < 10
    row = DL.recover_full_donors(observation(), lookup([(D1, 10)])).iloc[0]
    assert row['donor_recovery_tier'] == 'confident'


def test_an_exact_tie_at_the_top_is_not_guessed():
    row = DL.recover_full_donors(observation(), lookup([(D1, 7), (D2, 7), (D3, 1)])).iloc[0]
    assert pd.isna(row['donor']) and row['donor_recovery_status'] == 'ambiguous'
    assert row['donor_recovery_tier'] == 'unassigned'
    assert json.loads(row['donor_recovery_candidate_support']) == [[D1, 7], [D2, 7], [D3, 1]]


def test_unknown_donor_reads_compete_but_cannot_win():
    clear = DL.recover_full_donors(observation(), lookup([(D1, 50), (None, 5)])).iloc[0]
    assert clear['donor'] == D1 and clear['donor_recovery_status'] == 'recovered_by_support'
    assert clear['donor_recovery_runner_up_reads'] == 5
    unclear = DL.recover_full_donors(observation(), lookup([(D1, 5), (None, 5)])).iloc[0]
    assert pd.isna(unclear['donor']) and unclear['donor_recovery_status'] == 'unresolved_reference'


def test_perfect_donor_reads_set_the_confident_tier():
    """S3: rank by perfect-donor reads for the tier. Assignment still follows all-read support (S1+S2+S6)."""
    perfect = {D2}.__contains__
    # D1 leads on all reads, but the perfect-donor leader is D2: assigned D1, not confident.
    row = DR.recover_ranked_donors(observation(), lookup([(D1, 40), (D2, 12)]), is_perfect_donor=perfect).iloc[0]
    assert row['donor'] == D1 and row['donor_recovery_tier'] == 'provisional'
    assert row['donor_recovery_perfect_reads'] == 0
    # D2 leads on all reads AND on perfect reads; the all-read lead (40 vs 30) is under 3x, the perfect one is not.
    row = DR.recover_ranked_donors(observation(), lookup([(D2, 40), (D1, 30)]), is_perfect_donor=perfect).iloc[0]
    assert row['donor'] == D2 and row['donor_recovery_tier'] == 'confident'
    assert row['donor_recovery_perfect_reads'] == 40
    # Without the predicate the same row is tiered on all reads, and the perfect-read column is not emitted
    # (a column that is constant on every row names nothing).
    row = DR.recover_ranked_donors(observation(), lookup([(D2, 40), (D1, 30)])).iloc[0]
    assert row['donor_recovery_tier'] == 'provisional'
    assert 'donor_recovery_perfect_reads' not in row.index


def test_recovery_never_multiplies_rows_or_adds_reference_reads():
    obs = observation(n=3)
    result = DL.recover_full_donors(obs, lookup([(D1, 50), (D2, 5), (D3, 2)]))
    assert len(result) == 3 and result['counts'].tolist() == [20, 20, 20]
    assert result['donor'].tolist() == [D1, D1, D1]


@pytest.mark.parametrize('seed', range(5))
def test_recovery_is_independent_of_candidate_order(seed):
    ref = lookup([(D1, 9), (D2, 30), (D3, 30), (D1, 22)])          # D1 totals 31 across two rows
    a = DL.recover_full_donors(observation(), ref)
    b = DL.recover_full_donors(observation(), ref.sample(frac=1.0, random_state=seed).reset_index(drop=True))
    pd.testing.assert_frame_equal(a, b)
    assert a.iloc[0]['donor'] == D1 and a.iloc[0]['donor_recovery_lead_reads'] == 31


def test_candidate_rows_with_a_blank_library_id_are_kept(tmp_path):
    """R* S1: a blank library_ID is missing data, not a reason to discard the evidence."""
    data, keys = tmp_path / 'data', tmp_path / 'keys'
    data.mkdir(); keys.mkdir()
    pd.DataFrame([dict(description='donor-bc0', sequencing_date='20250101', outer_primers='plate_1A',
                       sample_number=1, library_ID=None),
                  dict(description='donor-bc0', sequencing_date='20250101', outer_primers='plate_1B',
                       sample_number=2, library_ID='V429')]).to_csv(
        keys / 'bc1_donor_bc0_combined_key.tsv', sep='\t', index=False)
    lookup([(D1, 30)]).to_csv(data / '20250101_plate_1A_1_complete_donor_bc0_counts.tsv', sep='\t', index=False)
    lookup([(D2, 4)]).to_csv(data / '20250101_plate_1B_2_complete_donor_bc0_counts.tsv', sep='\t', index=False)
    candidates = DL.load_donor_bc0_lookup(data, keys)
    assert sorted(candidates['counts']) == [4, 30]
    assert candidates['library_ID'].isna().sum() == 1
    blank = DL.recover_full_donors(observation(library_ID=None), candidates).iloc[0]
    assert blank['donor'] == D1
    scoped = DL.recover_full_donors(observation(library_ID='V429'), candidates).iloc[0]
    assert scoped['donor'] == D2 and scoped['donor_recovery_status'] == 'recovered_unique'


# --- NN 250 / C17, carried onto the R* route ---------------------------------------------------------
# Zero recovery from two NON-EMPTY tables is a failure, never a result. The R* route labels rows instead
# of raising, so without this guard a drifted fragment key reads as a legitimate all-ambiguous cohort.
# "Recovered nothing" is defined on the fragment strings themselves: an all-ambiguous cohort, or one whose
# candidates are excluded by scope, still shares fragments with its observations and must not raise.

def _stage01_fragment():
    return D1[-30:] + 'G' * 20 + 'A' * 10


def test_c17_a_drifted_key_raises_on_the_rstar_route():
    obs = observation(); obs['donor_bc0_fragment'] = _stage01_fragment()
    drifted = lookup([(D1, 50)]); drifted['donor_bc0_fragment'] = D1[-30:] + 'A' * 10     # the 40-nt key
    with pytest.raises(DL.BridgeRecoveredNothing):
        DL.recover_full_donors(obs, drifted)


def test_c17_a_healthy_key_still_bridges():
    obs = observation(); obs['donor_bc0_fragment'] = _stage01_fragment()
    good = lookup([(D1, 50)]); good['donor_bc0_fragment'] = _stage01_fragment()
    assert DL.recover_full_donors(obs, good)['donor'].tolist() == [D1]


def test_c17_an_all_ambiguous_cohort_is_a_result_not_a_failure():
    result = DL.recover_full_donors(observation(), lookup([(D1, 5), (D2, 5)]))
    assert result['donor_recovery_status'].tolist() == ['ambiguous']


def test_c17_scope_exclusion_is_not_key_drift():
    result = DL.recover_full_donors(observation(trafo_ID='trafo_3'), lookup([(D1, 50)], trafo_ID=['trafo_4']))
    assert result['donor_recovery_status'].tolist() == ['no_positive_candidate']
