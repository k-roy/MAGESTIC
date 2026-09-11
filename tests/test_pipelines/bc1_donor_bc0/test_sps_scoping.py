"""SPS scoping: a row is matched only against the subpool its own molecule came from (planning 34_/36_).

Ordered by Kevin 2026-09-10: "SPS scoping in 04: implement and re-run." Measured on the 2026-09-10
rebuild, 313,570 of 4,458,480 named rows with a readable SPS (7.0%) carried a name from a different
subpool than the read came from, and the transformation independently agreed with the SPS on 99.3%.

🔴 SPS names the SUBPOOL ONLY (policy §4): V450/V454 share one SPS, V451/V455 the other, so nothing here
may be read as a scaffold claim.
"""
import pandas as pd
import pytest

from magestic.pipelines.bc1_donor_bc0.core.matching import (
    UNSCOPED, assign_sps_scope, build_lookup_tables, load_sps_mapping, match_dataframe,
    match_dataframe_scoped)

SPS_A = "AAAACCCCGGGGTTTTAAAA"
SPS_B = "CCCCGGGGTTTTAAAACCCC"
DONOR = "A" * 60 + "TTTT" + "A" * 65     # ONE donor, designed in BOTH subpools with different edits

# the two subpools install different variants from an identical donor, so only the scope can decide
OLIGOS_A = pd.DataFrame({"oligo_name": ["LbCas12a_x"], "guide": ["A" * 20], "donor": [DONOR],
                         "individual_variants": ["10_A_C_SNP"]})
OLIGOS_B = pd.DataFrame({"oligo_name": ["impLbCas12a_y"], "guide": ["C" * 20], "donor": [DONOR],
                         "individual_variants": ["20_G_T_SNP"]})
TABLES = {"V454": build_lookup_tables(OLIGOS_A), "V455": build_lookup_tables(OLIGOS_B),
          None: build_lookup_tables(pd.concat([OLIGOS_A, OLIGOS_B], ignore_index=True))}


@pytest.fixture
def mapping_file(tmp_path):
    p = tmp_path / "sps.tsv"
    p.write_text("sps\tlibrary_id\n" + f"{SPS_A}\tV454\n{SPS_B}\tV455\n" + "GGGG\tV429\n")
    return p


def test_load_sps_mapping_uppercases_and_filters(mapping_file):
    m = load_sps_mapping(str(mapping_file))
    assert m == {SPS_A: "V454", SPS_B: "V455", "GGGG": "V429"}
    assert load_sps_mapping(str(mapping_file), ["V454", "V455"]) == {SPS_A: "V454", SPS_B: "V455"}


def test_load_sps_mapping_refuses_to_run_unscoped_while_claiming_to_be_scoped(mapping_file, tmp_path):
    """An empty map would turn scoping OFF without saying so -- the run would look scoped and not be."""
    with pytest.raises(ValueError, match="empty after filtering"):
        load_sps_mapping(str(mapping_file), ["V999"])
    bad = tmp_path / "bad.tsv"
    bad.write_text("sps\twrong_col\nAAAA\tV1\n")
    with pytest.raises(ValueError, match="lacks"):
        load_sps_mapping(str(bad))


def test_assign_sps_scope():
    df = pd.DataFrame({"SPS": [SPS_A, SPS_B.lower(), "TTTT", None]})
    sc = assign_sps_scope(df, {SPS_A: "V454", SPS_B: "V455"})
    assert sc.tolist()[:2] == ["V454", "V455"]          # case-insensitive
    assert pd.isna(sc.iloc[2]) and pd.isna(sc.iloc[3])  # unknown and absent both fall back
    assert assign_sps_scope(pd.DataFrame({"x": [1]}), {}).isna().all()   # no SPS column at all


def test_scoping_changes_the_name_on_a_donor_designed_in_both_subpools():
    """The whole point: unscoped, one donor present in both subpools resolves by sorted-first across
    subpools. Scoped, each row takes the design its own molecule could actually have come from."""
    df = pd.DataFrame({"donor": [DONOR, DONOR], "SPS": [SPS_A, SPS_B]})
    scope = assign_sps_scope(df, {SPS_A: "V454", SPS_B: "V455"})

    unscoped = match_dataframe(df.copy(), TABLES[None], donor_col="donor")
    assert unscoped["oligo_name"].tolist() == ["LbCas12a_x", "LbCas12a_x"]   # both take the same name
    assert not unscoped["is_unambiguous"].any()                              # and both are ambiguous

    scoped = match_dataframe_scoped(df.copy(), scope, TABLES.get, donor_col="donor")
    assert scoped["oligo_name"].tolist() == ["LbCas12a_x", "impLbCas12a_y"]
    assert scoped["is_unambiguous"].all()


def test_row_order_and_index_are_preserved():
    df = pd.DataFrame({"donor": [DONOR] * 5, "SPS": [SPS_B, SPS_A, SPS_B, "TTTT", SPS_A]},
                      index=[10, 11, 12, 13, 14])
    scope = assign_sps_scope(df, {SPS_A: "V454", SPS_B: "V455"})
    out = match_dataframe_scoped(df.copy(), scope, TABLES.get, donor_col="donor")
    assert out.index.tolist() == [10, 11, 12, 13, 14]
    assert out["oligo_name"].tolist() == ["impLbCas12a_y", "LbCas12a_x", "impLbCas12a_y",
                                          "LbCas12a_x", "LbCas12a_x"]
    assert out.loc[13, "oligo_name"] == "LbCas12a_x"      # unrecognised SPS -> full scope, ambiguous
    assert not bool(out.loc[13, "is_unambiguous"])


def test_unscoped_rows_go_to_the_full_lookup():
    seen = []

    def build(v):
        seen.append(v)
        return TABLES[v]

    df = pd.DataFrame({"donor": [DONOR, DONOR], "SPS": ["TTTT", SPS_A]})
    scope = assign_sps_scope(df, {SPS_A: "V454"})
    match_dataframe_scoped(df, scope, build, donor_col="donor")
    assert None in seen and "V454" in seen and UNSCOPED not in seen


def test_empty_frame():
    df = pd.DataFrame({"donor": pd.Series(dtype=str), "SPS": pd.Series(dtype=str)})
    out = match_dataframe_scoped(df, assign_sps_scope(df, {SPS_A: "V454"}), TABLES.get, donor_col="donor")
    assert len(out) == 0
