"""A merged BDB row must remain an observed link with its own support."""
import json

import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.merge import merge_data_sources

BC1 = "G" * 20
SPS = "CGAT" + "G" * 16


def full(donor, count, pcr="p1", bc0="A" * 10, **extra):
    return dict(bc1=BC1, bc0=bc0, donor=donor, counts=count, pcr_replicate_id=pcr, **extra)


@pytest.mark.parametrize("reverse", [False, True])
def test_competing_observations_never_make_unobserved_tuple(reverse):
    rows = [full("A" * 129, 100), full("C" * 129, 1, "p2", "C" * 10)]
    merged = merge_data_sources(df_2x300=pd.DataFrame(rows[::-1] if reverse else rows))
    observed = {(r.donor, r.bc0): (r.total_counts, r.total_pcr_replicates)
                for r in merged.itertuples()}
    assert observed == {("A" * 129, "A" * 10): (100, 1), ("C" * 129, "C" * 10): (1, 1)}


def test_exact_link_counts_sum_but_pcr_identities_union():
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full("A" * 129, 8), full("A" * 129, 12, "p2"), full("A" * 129, 4)]))
    assert len(merged) == 1
    assert merged.iloc[0]["total_counts"] == 24
    assert merged.iloc[0]["total_pcr_replicates"] == 2
    assert merged.iloc[0]["pcr_replicate_ids"] == "p1;p2"


def test_scope_and_sps_are_not_reconstructed_from_barcode():
    merged = merge_data_sources(df_2x300=pd.DataFrame([
        full("A" * 129, 10, trafo_ID="trafo_3", SPS=SPS),
        full("A" * 129, 20, "p2", trafo_ID="trafo_4", SPS=SPS)]))
    assert sorted(merged["total_counts"]) == [10, 20]
    assert set(merged["trafo_ID"]) == {"trafo_3", "trafo_4"}
    assert merged["total_pcr_replicates"].tolist() == [1, 1]


@pytest.mark.parametrize("same_pcr", [False, True])
def test_unique_explicit_short_geometry_attaches_once(same_pcr):
    donor = "T" * 99 + "ACGAT" + "A" * 25
    short = dict(bc1=BC1, bc0="A" * 10, donor_bc0_fragment=donor[-30:] + SPS + "A" * 10,
                 SPS=SPS, counts=10, pcr_replicate_id="p1")
    merged = merge_data_sources(pd.DataFrame([short]),
                               pd.DataFrame([full(donor, 20, "p1" if same_pcr else "p2", SPS=SPS)]))
    assert len(merged) == 1
    row = merged.iloc[0]
    assert (row["donor"], row["counts_2x150"], row["counts_2x300"]) == (donor, 10, 20)
    assert row["total_counts"] == 30
    assert row["total_pcr_replicates"] == (1 if same_pcr else 2)
    assert row["donor_bc0_fragment"] == short["donor_bc0_fragment"]


@pytest.mark.parametrize("reverse", [False, True])
def test_ambiguous_short_fragment_keeps_its_own_counts(reverse):
    rows = [full("A" * 99 + "G" * 30, 100, SPS=SPS),
            full("C" * 99 + "G" * 30, 10, "p2", SPS=SPS)]
    short = dict(bc1=BC1, bc0="A" * 10, donor_bc0_fragment="G" * 30 + SPS + "A" * 10,
                 SPS=SPS, counts=5, pcr_replicate_id="p3")
    result = merge_data_sources(pd.DataFrame([short]), pd.DataFrame(rows[::-1] if reverse else rows))
    assert len(result) == 3 and result["total_counts"].sum() == 115
    unknown = result[result["donor"].isna()].iloc[0]
    assert unknown["total_counts"] == 5 and unknown["short_link_candidates"] == 2
    assert sorted(result["total_pcr_replicates"]) == [1, 1, 1]


def test_unknown_fragment_geometry_is_not_a_first_motif_guess():
    short = dict(bc1=BC1, bc0="A" * 10, donor_bc0_fragment="CGAT" + "A" * 56,
                 counts=5, pcr_replicate_id="p1")
    result = merge_data_sources(pd.DataFrame([short]), pd.DataFrame([full("A" * 129, 100, "p2")]))
    assert len(result) == 2
    assert result[result["donor"].isna()]["total_counts"].tolist() == [5]
    assert result["total_counts"].sum() == 105


def test_missing_sequences_remain_observable():
    result = merge_data_sources(df_2x300=pd.DataFrame([
        full("A" * 129, 10), full(None, 90, "p2"), full("C" * 129, 11, "p3", None)]))
    assert len(result) == 3 and result["total_counts"].sum() == 111
    assert result["donor"].isna().sum() == result["bc0"].isna().sum() == 1


@pytest.mark.parametrize("bad", [None, -1, 0.5, float("inf"), "bad"])
def test_invalid_read_counts_are_refused(bad):
    with pytest.raises(ValueError, match="count"):
        merge_data_sources(df_2x300=pd.DataFrame([full("A" * 129, bad)]))


def test_anonymous_support_does_not_become_a_sum_of_independent_pcrs():
    rows = [full("A" * 129, 10), full("A" * 129, 20)]
    for row in rows:
        del row["pcr_replicate_id"]
        row["num_PCR_replicates"] = 2
    result = merge_data_sources(df_2x300=pd.DataFrame(rows)).iloc[0]
    assert result["total_counts"] == 30
    assert result["total_pcr_replicates"] == 2
    assert not result["pcr_support_exact"]


def test_zero_reads_cannot_add_pcr_support():
    result = merge_data_sources(df_2x300=pd.DataFrame([
        full("A" * 129, 10), full("A" * 129, 0, "p2")])).iloc[0]
    assert result["total_pcr_replicates"] == 1


def test_gdb_suffix_collision_does_not_assign_a_conflicting_full_donor():
    donor = "A" * 99 + "G" * 30
    gdb = pd.DataFrame([dict(donor="C" * 99 + "G" * 30, bc0="A" * 10, guide="C" * 20, counts=1000)])
    result = merge_data_sources(df_2x300=pd.DataFrame([full(donor, 10)]), df_gdb=gdb).iloc[0]
    assert result["donor"] == donor and result["total_counts"] == 10
    assert not result["has_gdb"] and pd.isna(result["guide"])


def test_gdb_match_adds_annotation_without_bdb_support():
    donor = "A" * 129
    gdb = pd.DataFrame([dict(donor=donor, bc0="A" * 10, guide="C" * 20, counts=1000)])
    result = merge_data_sources(df_2x300=pd.DataFrame([full(donor, 10)]), df_gdb=gdb).iloc[0]
    assert result["guide"] == "C" * 20 and result["counts_gdb"] == 1000
    assert result["total_counts"] == 10 and result["total_pcr_replicates"] == 1


def test_empty_merge_has_a_usable_link_schema():
    result = merge_data_sources()
    assert result.empty
    assert {"bc1", "donor", "bc0", "total_counts", "pcr_replicate_ids"} <= set(result.columns)


@pytest.mark.parametrize("reverse", [False, True])
def test_missing_full_sps_does_not_let_first_short_choose_its_scope(reverse):
    rows = [dict(bc1=BC1, bc0="A" * 10, donor_bc0_fragment="A" * 30 + sps + "A" * 10,
                 SPS=sps, counts=count, pcr_replicate_id=pcr)
            for sps, count, pcr in [(SPS, 5, "p2"), ("T" * 20, 7, "p3")]]
    result = merge_data_sources(pd.DataFrame(rows[::-1] if reverse else rows),
                               pd.DataFrame([full("A" * 129, 100)]))
    assert len(result) == 3 and result["total_counts"].sum() == 112
    assert result[result["donor"].notna()]["total_counts"].tolist() == [100]


def test_partial_gdb_lookup_cannot_erase_observed_short_link_ambiguity():
    d1, d2 = "A" * 99 + "G" * 30, "C" * 99 + "G" * 30
    short = dict(bc1=BC1, bc0="A" * 10, donor_bc0_fragment="G" * 30 + SPS + "A" * 10,
                 SPS=SPS, counts=5, pcr_replicate_id="p3")
    gdb = pd.DataFrame([dict(donor=d1, bc0="A" * 10, guide="T" * 20, counts=1000, SPS=SPS)])
    result = merge_data_sources(pd.DataFrame([short]), pd.DataFrame([
        full(d1, 100, SPS=SPS), full(d2, 10, "p2", SPS=SPS)]), gdb)
    unresolved = result[result["donor"].isna()]
    assert unresolved["total_counts"].tolist() == [5]
    assert result["total_counts"].sum() == 115
