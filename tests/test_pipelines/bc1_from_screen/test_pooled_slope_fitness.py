"""Hand-checkable tests for the NRD1 pooled-slope fitness method (core.pooled_slope_fitness).

Real-data parity (all six tables of the NRD1 R script on three inputs, and the 23 Supplementary Table 1 hits) is
planning/380 and `test_nrd1_parity_on_oak` below, which runs only where Oak is mounted."""
import os

import numpy as np
import pandas as pd
import pytest

from magestic.pipelines.bc1_from_screen.core import pooled_slope_fitness as psf

GENS = (0, 4, 8)


def _rows(bc1, pos, mtype, treatment, values, replicate=1):
    return [{"bc1": bc1, "amino_acid_position": pos, "mutation_type": mtype, "treatment": treatment,
             "replicate": replicate, "generation": g, "value": v} for g, v in zip(GENS, values)]


def _long():
    r = []
    # edit 1: depleted with the drug only. log2FC 0, -2, -4 -> slope -0.5; control flat
    for b in ("e1a", "e1b"):
        r += _rows(b, 1, "missense", "T", [100, 25, 6.25]) + _rows(b, 1, "missense", "C", [100, 100, 100])
    # edit 2: two barcodes of different abundance; the pooled series is the g0-weighted mean
    r += _rows("e2a", 2, "missense", "T", [100, 50, 25]) + _rows("e2a", 2, "missense", "C", [100, 100, 100])
    r += _rows("e2b", 2, "missense", "T", [300, 300, 300]) + _rows("e2b", 2, "missense", "C", [300, 300, 300])
    # edit 3: enriched. log2FC 0, 1, 2 -> slope +0.25
    for b in ("e3a", "e3b"):
        r += _rows(b, 3, "synonymous", "T", [100, 200, 400]) + _rows(b, 3, "synonymous", "C", [100, 100, 100])
    # a barcode below the generation-0 floor is dropped
    r += _rows("low", 3, "synonymous", "T", [5, 5, 5]) + _rows("low", 3, "synonymous", "C", [5, 5, 5])
    return pd.DataFrame(r)


CFG = psf.PooledSlopeConfig(treatment="T", control="C", tail_prob=0.5, depletion_fc=0.1)


def test_pooled_series_slope_centering_and_hits():
    res = psf.pooled_slope_fitness(_long(), CFG)
    m = res["mean_log2FC_by_group"].set_index(["amino_acid_position", "treatment", "generation"])["mean_log2fc"]
    assert m[(2, "T", 4)] == pytest.approx((100 * -1 + 300 * 0) / 400)        # weights max(g0, 1)
    st = res["slope_table_log2FC"]
    t = st[st.treatment == "T"].set_index("amino_acid_position")
    assert t.loc[1, "slope"] == pytest.approx(-0.5) and t.loc[3, "slope"] == pytest.approx(0.25)
    assert t.loc[2, "slope"] == pytest.approx(-0.0625)
    med = np.median([-0.5, -0.0625, 0.25])
    assert t.loc[1, "slope_adj"] == pytest.approx(-0.5 - med)
    h = res["hits_with_support"]
    assert list(h.amino_acid_position) == [1] and int(h.n_support.iloc[0]) == 2
    assert set(st[st.treatment == "HIT"].amino_acid_position) == {1}
    assert "low" not in set(res["barcode_support_flags"].bc1)


def test_per_barcode_slope_uses_the_first_replicate_at_each_generation():
    d = _long()
    extra = pd.DataFrame(_rows("e1a", 1, "missense", "T", [100, 100, 100], replicate=2))
    res = psf.pooled_slope_fitness(pd.concat([d, extra], ignore_index=True), CFG)
    f = res["barcode_support_flags"].set_index("bc1")
    med = np.median([-0.5, -0.5, -0.25, 0.0, 0.25, 0.25])  # per-barcode slopes; replicate 1 decides e1a
    assert f.loc["e1a", "T"] == pytest.approx(-0.5 - med)
    pooled = res["slope_table_log2FC"].query("treatment == 'T' and amino_acid_position == 1")["slope"].iloc[0]
    assert pooled > -0.5  # the pooled series does include replicate 2


def test_truncation_at_three_points_below_the_floor():
    g = np.array([0, 4, 8, 12, 16, 20], float)
    y = np.array([0, -1, -8, -9, -10, -11], float)
    assert psf._cut_generation(g, y, 0.01, 3) == 16.0          # third consecutive point below log2(0.01)
    assert psf._cut_generation(g[[1, 2, 3]], y[[1, 2, 3]], 0.01, 3) == 12.0  # no generation 0: no trim


def test_control_none_gives_fitness_without_hit_calls():
    res = psf.pooled_slope_fitness(_long(), psf.PooledSlopeConfig(treatment="T", control=None))
    assert "hits_summary" not in res and "hits_with_support" not in res
    assert "hit" not in res["slope_table_log2FC"].columns


def test_refuses_without_generation_zero():
    d = _long()
    d["generation"] = d["generation"] - 25        # the 20220607 key style (planning/288 item 4)
    with pytest.raises(ValueError, match="generation-0"):
        psf.pooled_slope_fitness(d, CFG)


def test_count_matrix_adapter_cpm_order_and_refusals():
    counts = pd.DataFrame({"s0": [10, 30], "s4": [5, 45]}, index=["A", "B"])
    samples = pd.DataFrame({"sample": ["s0", "s4"], "treatment": ["T", "T"], "replicate": [1, 1], "generation": [0, 4]})
    keys = pd.DataFrame({"bc1": ["A", "B"], "variant_id": ["v1", "v2"]})
    long = psf.long_from_count_matrix(counts, samples, keys, pseudocount=1)
    a = long[long.bc1 == "A"].set_index("generation")["value"]
    assert a[0] == pytest.approx(10 / 40 * 1e6 + 1) and a[4] == pytest.approx(5 / 50 * 1e6 + 1)
    with pytest.raises(ValueError, match="duplicate"):
        psf.long_from_count_matrix(counts, samples, pd.concat([keys, keys]))
    with pytest.raises(ValueError, match="integer generation"):
        psf.long_from_count_matrix(counts, samples.assign(generation=["0", "x"]), keys)


def test_mutation_type_and_wide_parsing():
    wide = pd.DataFrame({"bc1": ["A", "B"], "gene_name": ["NRD1"] * 2, "amino_acid_position": ["5", "6"],
                         "ref_amino_acid": ["L", "M"], "mutant_amino_acid": ["*", "M"],
                         "multiple_donor_bc0_fragments_passing_filter": ["True", "False"],
                         "YPD_set_R1_Normalized_Count_0": [1.0, 2.0], "YPD noIAA_set_R2_Normalized_Count_4": [3.0, 4.0]})
    long = psf.long_from_wide(wide)
    assert set(long.bc1) == {"A"} and set(long.mutation_type) == {"PTC"}
    assert sorted(long.treatment) == ["YPD", "YPD_noIAA"] and sorted(long.replicate.tolist()) == [1, 2]


W = "/oak/stanford/groups/larsms/Users/kevinroy/projects/NNS/20250813_NRD1_MAGESTIC_screen_reps_1-3/reanalysis_289"


@pytest.mark.integration
@pytest.mark.skipif(not os.path.isdir(W), reason="needs Oak")
def test_nrd1_parity_on_oak():
    res = psf.pooled_slope_fitness(psf.long_from_wide(pd.read_csv(f"{W}/inputs/pc1_fill-zero_mdf-keep.csv")),
                                   psf.PooledSlopeConfig(tail_prob=0.05, depletion_fc=0.1))
    r = pd.read_csv(f"{W}/runs/pc1_fill-zero_mdf-keep__tail0.05_dep0.1/slope_table_log2FC.csv")
    mine = res["slope_table_log2FC"]
    j = mine.merge(r, on=["amino_acid_position", "mutation_type", "treatment"], suffixes=("", "_r"))
    assert len(j) == len(r) == len(mine)
    assert np.allclose(j.slope_adj, j.slope_adj_r, atol=1e-12, equal_nan=True)
    assert len(res["hits_summary"]) == 26 and len(res["hits_with_support"]) == 23


def test_both_directions_call_the_enriched_mirror_without_changing_depleted_calls():
    cfg = psf.PooledSlopeConfig(treatment="T", control="C", tail_prob=0.5, depletion_fc=0.3, direction="both")
    res = psf.pooled_slope_fitness(_long(), cfg)
    assert list(res["hits_with_support"].amino_acid_position) == [1]
    e = res["hits_with_support_enriched"]
    assert list(e.amino_acid_position) == [3] and int(e.n_support.iloc[0]) == 2
    base = psf.pooled_slope_fitness(_long(), psf.PooledSlopeConfig(treatment="T", control="C", tail_prob=0.5,
                                                                   depletion_fc=0.3))
    assert "hits_with_support_enriched" not in base
    pd.testing.assert_frame_equal(base["hits_with_support"].reset_index(drop=True),
                                  res["hits_with_support"].reset_index(drop=True))


def test_median_depth_normalization_keeps_a_read_floor():
    counts = pd.DataFrame({"s0": [10, 30], "s4": [20, 80]}, index=["A", "B"])
    samples = pd.DataFrame({"sample": ["s0", "s4"], "treatment": ["T", "T"], "replicate": [1, 1], "generation": [0, 4]})
    keys = pd.DataFrame({"bc1": ["A", "B"], "variant_id": ["v1", "v2"]})
    long = psf.long_from_count_matrix(counts, samples, keys, normalize="median_depth")
    a = long[long.bc1 == "A"].set_index("generation")["value"]
    assert a[0] == pytest.approx(10 / 40 * 70) and a[4] == pytest.approx(20 / 100 * 70)


# ---------------------------------------------------------------- the endpoint route (planning/384 §3, Kevin 2026-11-03)

def _two_point_long():
    """Generation 0 and ONE later generation: the shape that yields no slope at all (Dec 2021; yL437/yL442)."""
    r = []
    for b, vals in (("a1", [100, 25]), ("a2", [100, 25]), ("b1", [100, 100]), ("b2", [100, 100])):
        pos = 1 if b.startswith("a") else 2
        r += [{"bc1": b, "amino_acid_position": pos, "mutation_type": "missense", "treatment": t,
               "replicate": 1, "generation": g, "value": v}
              for t, series in (("T", vals), ("C", [100, 100])) for g, v in zip((0, 6), series)]
    return pd.DataFrame(r)


ECFG = psf.PooledSlopeConfig(treatment="T", control="C", tail_prob=0.5, depletion_fc=0.1, index="endpoint")


def test_endpoint_index_is_the_pooled_log2fc_at_the_last_generation():
    res = psf.pooled_slope_fitness(_long(), ECFG)
    t = res["slope_table_log2FC"].query("treatment == 'T'").set_index("amino_acid_position")
    # edit 1 reaches -4 at generation 8; edit 3 reaches +2; edit 2 is the g0-weighted mean (100*-2 + 300*0)/400
    assert t.loc[1, "slope"] == pytest.approx(-4.0) and t.loc[3, "slope"] == pytest.approx(2.0)
    assert t.loc[2, "slope"] == pytest.approx(-0.5)
    assert set(t.endpoint_gen.dropna()) == {8.0}
    med = np.median([-4.0, -0.5, 2.0])
    assert t.loc[1, "slope_adj"] == pytest.approx(-4.0 - med)
    assert list(res["hits_with_support"].amino_acid_position) == [1]


def test_endpoint_generation_selects_that_point_and_never_substitutes_a_neighbour():
    res = psf.pooled_slope_fitness(_long(), psf.PooledSlopeConfig(
        treatment="T", control="C", tail_prob=0.5, depletion_fc=0.1, index="endpoint", endpoint_generation=4))
    t = res["slope_table_log2FC"].query("treatment == 'T'").set_index("amino_acid_position")
    assert t.loc[1, "slope"] == pytest.approx(-2.0) and set(t.endpoint_gen.dropna()) == {4.0}
    missing = psf.pooled_slope_fitness(_long(), psf.PooledSlopeConfig(
        treatment="T", control="C", tail_prob=0.5, depletion_fc=0.1, index="endpoint", endpoint_generation=99))
    assert missing["slope_table_log2FC"].empty


def test_endpoint_rescues_a_two_timepoint_screen_the_slope_route_cannot_fit():
    d = _two_point_long()
    slope = psf.pooled_slope_fitness(d, psf.PooledSlopeConfig(treatment="T", control="C", tail_prob=0.5,
                                                              depletion_fc=0.1))
    assert slope["slope_table_log2FC"].slope.notna().sum() == 0      # min_points=3 on a 2-point series
    end = psf.pooled_slope_fitness(d, ECFG)
    t = end["slope_table_log2FC"].query("treatment == 'T'").set_index("amino_acid_position")
    assert t.loc[1, "slope"] == pytest.approx(-2.0) and t.loc[2, "slope"] == pytest.approx(0.0)
    assert t.loc[1, "slope_adj"] == pytest.approx(-1.0) and t.loc[2, "slope_adj"] == pytest.approx(1.0)
    assert int(end["hit_barcode_support_summary"].set_index("amino_acid_position").loc[1, "n_barcodes"]) == 2


def test_endpoint_refuses_without_any_generation_above_zero():
    d = _long()
    d = d[d.generation == 0]
    with pytest.raises(ValueError, match="generation above 0"):
        psf.pooled_slope_fitness(d, ECFG)


def test_config_rejects_a_bad_route_and_generation_zero_as_the_endpoint():
    with pytest.raises(ValueError, match="index must be"):
        psf.PooledSlopeConfig(index="trend")
    with pytest.raises(ValueError, match="generation 0 itself"):
        psf.PooledSlopeConfig(index="endpoint", endpoint_generation=0)


def test_slope_route_is_byte_identical_after_the_endpoint_change():
    """The endpoint route must not perturb the published path: same tables, same values, no new columns."""
    res = psf.pooled_slope_fitness(_long(), CFG)
    assert "endpoint_gen" not in res["slope_table_log2FC"].columns
    st = res["slope_table_log2FC"]
    t = st[st.treatment == "T"].set_index("amino_acid_position")
    assert t.loc[1, "slope"] == pytest.approx(-0.5) and t.loc[2, "slope"] == pytest.approx(-0.0625)
    assert list(res["hits_with_support"].amino_acid_position) == [1]
