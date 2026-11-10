"""magestic.utils.metrics: known-answer property tests, and the depth trap as arithmetic.

These statistics are the kind that are easy to get subtly wrong and never notice, because every
one of them returns a plausible number from wrong code. So the tests pin the known answers
(uniform -> gini 0, pielou 1; one dominant member -> gini -> 1, pielou -> 0; Hill(0) == richness;
Hill(2) == inverse Simpson), the identities between the statistics, and the behaviour that matters
scientifically: **the values move with depth and with the count floor**. The depth tests are
documentation as much as verification — a reader who only reads the test names learns the trap.

The parity test against ``guide_donor_bc0.core.library_representation.gini_index`` is the proof
that this module is additive: the shared kernel reproduces the pipeline-local function exactly, so
nothing that already used it changes.
"""

import math

import numpy as np
import pytest

from magestic.utils import metrics as M


# ---------------------------------------------------------------------------
# No bare numbers
# ---------------------------------------------------------------------------

def test_every_public_statistic_returns_a_metric_not_a_float():
    counts = [10, 20, 30, 40]
    for fn, args in (
        (M.gini, ()), (M.shannon_entropy, ()), (M.pielou_evenness, ()),
        (M.simpson_index, ()), (M.inverse_simpson, ()), (M.richness, ()),
        (M.hill_number, (1.0,)), (M.top_n_share, (2,)),
    ):
        m = fn(counts, *args)
        assert isinstance(m, M.Metric), fn
        assert isinstance(m.context, M.MetricContext)


def test_a_metric_cannot_be_coerced_to_a_bare_float():
    """The structural half of the contract: quoting the value takes an explicit .value."""
    m = M.gini([1, 2, 3])
    with pytest.raises(TypeError):
        float(m)


def test_formatting_a_metric_always_carries_the_depth_and_floor():
    m = M.gini([1, 2, 3, 94], min_count=1)
    s = f"{m}"
    assert "depth=100" in s
    assert "floor=1" in s
    assert "NOT rarefied: sample-specific" in s


def test_context_records_what_the_floor_removed():
    m = M.gini([0, 1, 1, 5, 50], min_count=2)
    assert m.context.n_observations == 2            # 5 and 50
    assert m.context.n_excluded_by_floor == 3       # 0, 1, 1
    assert m.context.counts_excluded_by_floor == 2
    assert m.context.depth == 55


# ---------------------------------------------------------------------------
# Known answers
# ---------------------------------------------------------------------------

def test_uniform_vector_gives_gini_zero_and_pielou_one():
    counts = [25] * 8
    assert M.gini(counts).value == pytest.approx(0.0, abs=1e-12)
    assert M.pielou_evenness(counts).value == pytest.approx(1.0, abs=1e-12)
    assert M.hill_number(counts, 1.0).value == pytest.approx(8.0)
    assert M.inverse_simpson(counts).value == pytest.approx(8.0)
    assert M.top_n_share(counts, 1).value == pytest.approx(1 / 8)


@pytest.mark.parametrize("n", [10, 100, 1000])
def test_maximally_concentrated_vector_pushes_gini_to_its_ceiling_and_pielou_to_zero(n):
    """One member holds almost everything: gini -> (n-1)/n, pielou -> 0.

    🔴 The ceiling is (n-1)/n, **not 1**, and that is a real property of the Gini coefficient on a
    finite vector, not an artefact: with all mass on one of n members the value is exactly
    (n-1)/n. So "Gini approaching 1" is only true as n grows, and a Gini of 0.90 on ten members is
    maximal inequality, not moderate. A reader comparing two libraries with different numbers of
    members is comparing against two different ceilings -- one more reason the context travels
    with the value.
    """
    ceiling = (n - 1) / n
    counts = np.ones(n, dtype=float)
    counts[0] = 10 ** 6
    g = M.gini(counts).value
    j = M.pielou_evenness(counts).value
    assert 0.99 * ceiling < g < ceiling
    assert 0.0 < j < 0.02
    # the ceiling is approached from below and never exceeded, however extreme the skew
    harder = np.r_[10 ** 9, np.ones(n - 1)]
    assert g < M.gini(harder).value < ceiling
    assert M.pielou_evenness(harder).value < j
    # all mass on exactly one of n members hits the ceiling exactly
    only_one = np.zeros(n); only_one[0] = 7.0
    assert M.gini(only_one, min_count=0).value == pytest.approx(ceiling)
    # and the ceiling rises toward 1 only with n
    assert M.gini(np.r_[1.0, np.zeros(10 * n)], min_count=0).value > ceiling


def test_gini_of_a_single_occupied_member_is_zero_without_zeros_and_n_minus_1_over_n_with_them():
    """🔴 The floor trap as arithmetic — the same library, two correct answers.

    With the designed zeros kept (min_count=0) one occupied member out of n scores (n-1)/n, which
    is nearly 1. Drop the zeros (min_count=1) and the same library scores 0, because the only
    surviving member is trivially 'evenly' distributed among itself. This is exactly the
    difference between ``library_representation.gini_index``'s documented behaviour and what its
    caller ``plot_gini_indices`` actually computes with ``x[x > 0]``.
    """
    n = 10
    counts = [0] * (n - 1) + [500]
    with_zeros = M.gini(counts, min_count=0, universe=M.UNIVERSE_DESIGNED)
    without = M.gini(counts, min_count=1, universe=M.UNIVERSE_DESIGNED)
    assert with_zeros.value == pytest.approx((n - 1) / n)
    assert without.value == pytest.approx(0.0, abs=1e-12)
    # the context, not the value, is what tells them apart
    assert with_zeros.context.zeros_included is True
    assert without.context.zeros_included is False
    assert "designed zero-count" in (without.context.note or "")


def test_hill_at_order_zero_equals_richness():
    counts = [1, 1, 2, 3, 5, 8, 13, 0, 0]
    rich = M.richness(counts, min_count=1)
    hill0 = M.hill_number(counts, 0.0, min_count=1)
    assert rich.value == 7
    assert hill0.value == pytest.approx(rich.value)


def test_inverse_simpson_equals_hill_at_order_two():
    rng = np.random.default_rng(2174)
    for _ in range(5):
        counts = rng.integers(1, 5000, size=200)
        assert M.inverse_simpson(counts).value == pytest.approx(
            M.hill_number(counts, 2.0).value, rel=1e-12)


def test_hill_at_order_one_is_exp_shannon_and_the_limit_is_continuous():
    counts = [4, 9, 16, 25, 36]
    h1 = M.hill_number(counts, 1.0).value
    assert h1 == pytest.approx(math.exp(M.shannon_entropy(counts).value), rel=1e-12)
    # q=1 is a special case in the implementation; it must agree with the two-sided limit
    near_below = M.hill_number(counts, 1.0 - 1e-6).value
    near_above = M.hill_number(counts, 1.0 + 1e-6).value
    assert h1 == pytest.approx(near_below, rel=1e-5)
    assert h1 == pytest.approx(near_above, rel=1e-5)


def test_hill_at_infinity_is_the_reciprocal_of_the_largest_proportion():
    counts = [50, 25, 25]
    assert M.hill_number(counts, math.inf).value == pytest.approx(2.0)


def test_hill_profile_is_non_increasing_in_q():
    rng = np.random.default_rng(7)
    counts = rng.integers(1, 1000, size=300)
    values = [M.hill_number(counts, q).value for q in (0, 0.5, 1, 2, 4, math.inf)]
    for a, b in zip(values, values[1:]):
        assert a >= b - 1e-9


def test_simpson_is_the_plug_in_sum_of_squared_proportions():
    counts = np.array([1.0, 1.0, 2.0])
    p = counts / counts.sum()
    assert M.simpson_index(counts).value == pytest.approx(float((p ** 2).sum()))


def test_top_n_share_is_the_sorted_head_and_saturates_at_one():
    counts = [1, 2, 3, 94]
    assert M.top_n_share(counts, 1).value == pytest.approx(0.94)
    assert M.top_n_share(counts, 2).value == pytest.approx(0.97)
    assert M.top_n_share(counts, 99).value == pytest.approx(1.0)


def test_richness_is_reported_at_its_floor_and_names_it():
    counts = [0, 1, 1, 1, 4, 9, 100]
    assert M.richness(counts, min_count=1).value == 6
    assert M.richness(counts, min_count=2).value == 3
    assert M.richness(counts, min_count=2).name == "richness_min2"
    assert M.richness(counts, min_count=2).context.weighting == M.WEIGHTING_DISTINCT


# ---------------------------------------------------------------------------
# Degenerate input
# ---------------------------------------------------------------------------

def test_empty_and_all_zero_vectors_return_nan_not_zero():
    for counts in ([], [0, 0, 0]):
        assert math.isnan(M.gini(counts, min_count=0).value)
        assert math.isnan(M.shannon_entropy(counts, min_count=0).value)
        assert math.isnan(M.pielou_evenness(counts, min_count=0).value)


def test_pielou_is_undefined_for_a_single_member_rather_than_one_or_zero():
    assert math.isnan(M.pielou_evenness([42]).value)


def test_negative_and_nonfinite_counts_are_refused():
    with pytest.raises(ValueError):
        M.gini([1, -1, 3])
    with pytest.raises(ValueError):
        M.gini([1, np.nan, 3])


# ---------------------------------------------------------------------------
# 🔴 Depth dependence — the trap this module exists to encode
# ---------------------------------------------------------------------------

def _skewed_library(n=2000, seed=2174):
    """A log-normal count distribution with a heavy singleton tail, like a real plasmid pool."""
    rng = np.random.default_rng(seed)
    counts = np.rint(rng.lognormal(mean=1.0, sigma=2.0, size=n)).astype(float)
    counts[counts < 0] = 0
    return counts


def test_subsampling_a_known_distribution_moves_every_statistic():
    """The same library at two depths gives different numbers for every statistic here.

    This is the error the module is built to prevent: a value quoted without its depth is not a
    property of the library. Nothing is asserted about the direction of gini (subsampling both
    loses rare members and flattens the survivors); what is asserted is that the numbers MOVE,
    and that richness and the effective-number measures move the way depth forces them to.
    """
    counts = _skewed_library()
    full_depth = int(counts.sum())
    shallow = M.rarefy_counts(counts, full_depth // 20, seed=11)

    rich_full = M.richness(counts, min_count=1).value
    rich_shallow = M.richness(shallow, min_count=1).value
    assert rich_shallow < rich_full                      # depth only ever adds members

    hill_full = M.hill_number(counts, 1.0, min_count=1).value
    hill_shallow = M.hill_number(shallow, 1.0, min_count=1).value
    assert hill_shallow < hill_full                       # effective number falls with depth

    gini_full = M.gini(counts, min_count=1).value
    gini_shallow = M.gini(shallow, min_count=1).value
    assert gini_full != pytest.approx(gini_shallow, abs=1e-3)

    pielou_full = M.pielou_evenness(counts, min_count=1).value
    pielou_shallow = M.pielou_evenness(shallow, min_count=1).value
    assert pielou_full != pytest.approx(pielou_shallow, abs=1e-3)

    # and the contexts make the comparison's invalidity visible without reading the values
    assert M.richness(counts).context.depth != M.richness(shallow).context.depth


def test_raising_the_floor_moves_gini_a_long_way_on_a_singleton_heavy_library():
    """Why several floors are reported rather than one.

    The singleton tail dominates the distinct-key denominator, so the unfloored value is largely a
    statement about the error tail, not the library.
    """
    counts = _skewed_library()
    singleton_fraction = float(np.mean(counts[counts > 0] == 1))
    assert singleton_fraction > 0.1          # the regime this project is in
    g1 = M.gini(counts, min_count=1).value
    g10 = M.gini(counts, min_count=10).value
    assert g1 > g10 + 0.05                   # a large, scientifically material move
    assert M.gini(counts, min_count=10).context.n_excluded_by_floor > 0


def test_rarefaction_is_seeded_reproducible_and_labelled():
    counts = _skewed_library()
    target = int(counts.sum()) // 10
    a = M.rarefied_metric(M.gini, counts, target, replicates=5, seed=99)
    b = M.rarefied_metric(M.gini, counts, target, replicates=5, seed=99)
    c = M.rarefied_metric(M.gini, counts, target, replicates=5, seed=100)
    assert a.value == pytest.approx(b.value)
    assert a.value != pytest.approx(c.value, abs=1e-9)
    assert a.context.is_rarefied
    assert a.context.rarefied_to == target
    assert a.context.rarefaction_seed == 99
    assert a.context.rarefaction_replicates == 5
    assert a.sd is not None and a.sd >= 0
    assert "RAREFIED" in f"{a}"


def test_two_libraries_at_different_depths_compare_only_after_rarefaction():
    """The use case: the same underlying pool sequenced 10x deeper looks less even unrarefied."""
    pool = _skewed_library(seed=5)
    deep = M.rarefy_counts(pool, int(pool.sum()), seed=1)
    shallow = M.rarefy_counts(pool, int(pool.sum()) // 10, seed=2)

    raw_gap = abs(M.hill_number(deep, 1.0).value - M.hill_number(shallow, 1.0).value)
    target = int(shallow.sum())
    rare_deep = M.rarefied_metric(M.hill_number, deep, target, replicates=10, seed=3, q=1.0)
    rare_shallow = M.rarefied_metric(M.hill_number, shallow, target, replicates=10, seed=3, q=1.0)
    rare_gap = abs(rare_deep.value - rare_shallow.value)
    assert rare_gap < raw_gap


def test_a_library_shallower_than_the_target_returns_nan_with_a_flag_never_the_raw_value():
    counts = [5, 5, 5]
    m = M.rarefied_metric(M.gini, counts, 10 ** 6, replicates=3, seed=1)
    assert math.isnan(m.value)
    assert m.n_replicates == 0
    assert "below the rarefaction target" in m.context.note
    assert m.context.rarefied_to == 10 ** 6


def test_rarefy_counts_subsamples_without_replacement_and_preserves_length():
    counts = np.array([100.0, 50.0, 0.0, 1.0])
    sub = M.rarefy_counts(counts, 50, seed=4)
    assert sub.shape == counts.shape
    assert sub.sum() == 50
    assert np.all(sub <= counts)              # without replacement: never more than was there
    assert sub[2] == 0
    with pytest.raises(ValueError):
        M.rarefy_counts(counts, 10 ** 6)


def test_rarefaction_refuses_fractional_counts_instead_of_rounding_them_away():
    """[0.4, 0.4, 0.4] totals 1.2, so depth 1 passes the depth check, and every colour rounds to 0."""
    with pytest.raises(ValueError, match="whole-number"):
        M.rarefy_counts([0.4, 0.4, 0.4], 1, seed=1)
    with pytest.raises(ValueError, match="whole-number"):
        M.rarefied_metric(M.gini, [2.5, 3.5, 4.0], 5, replicates=2, seed=1)
    # whole numbers stored as floats are counts, and still subsample to exactly the target
    assert M.rarefy_counts([4.0, 3.0, 3.0], 5, seed=1).sum() == 5


def test_a_rarefied_metric_still_reports_what_the_floor_removed():
    """Rarefying [1, 9] to its full depth with floor 2 drops the singleton in every replicate."""
    m = M.rarefied_metric(M.richness, [1, 9], 10, replicates=3, seed=5, min_count=2)
    assert m.context.depth == 9 and m.context.n_observations == 1
    assert m.context.n_excluded_by_floor == 1
    assert m.context.counts_excluded_by_floor == 1


def test_hill_number_stays_finite_and_bounded_at_large_finite_orders():
    """The direct power sum underflows to 0.0 at q=1000, which the exponent then turned into inf."""
    uniform = [7] * 10
    for q in (50.0, 1000.0, 1e6):
        assert M.hill_number(uniform, q).value == pytest.approx(10.0, rel=1e-9)
    skewed = [900, 50, 30, 20]
    floor = M.hill_number(skewed, math.inf).value            # 1 / max(p), the profile's lower bound
    values = [M.hill_number(skewed, q).value for q in (2.0, 100.0, 1000.0, 1e5)]
    assert all(math.isfinite(v) for v in values)
    for a, b in zip(values, values[1:]):
        assert a >= b - 1e-9                                  # still non-increasing in q
    assert values[-1] == pytest.approx(floor, rel=1e-3) and values[-1] >= floor - 1e-9
    # a large negative order overflows instead; the limit there is 1 / min(p)
    assert M.hill_number(skewed, -1000.0).value == pytest.approx(
        M.hill_number(skewed, -math.inf).value, rel=1e-2)


# ---------------------------------------------------------------------------
# 🔴 The distinct-key vs read-weighted denominator
# ---------------------------------------------------------------------------

def test_coverage_fractions_returns_both_denominators_from_one_call():
    """Reproduces the shape of the 2168 result: ~1 % distinct, ~90 % read-weighted.

    One abundant member matches and a thousand singleton error variants do not. Both numbers are
    correct and they differ by two orders of magnitude; the function returns both so neither can
    be quoted alone.
    """
    counts = np.r_[9000.0, np.ones(1000)]
    selected = np.r_[True, np.zeros(1000, dtype=bool)]
    out = M.coverage_fractions(counts, selected, grain="guide_bc0", label="design-matched")
    assert out["distinct_key_fraction"].value == pytest.approx(1 / 1001, rel=1e-6)
    assert out["read_weighted_fraction"].value == pytest.approx(9000 / 10000)
    assert out["distinct_key_fraction"].context.weighting == M.WEIGHTING_DISTINCT
    assert out["read_weighted_fraction"].context.weighting == M.WEIGHTING_READ
    assert set(out) == {"distinct_key_fraction", "read_weighted_fraction"}


def test_coverage_fractions_rejects_a_mask_of_the_wrong_length():
    with pytest.raises(ValueError):
        M.coverage_fractions([1, 2, 3], [True, False])


# ---------------------------------------------------------------------------
# Batch report
# ---------------------------------------------------------------------------

def test_representation_metrics_reports_at_every_requested_floor():
    counts = _skewed_library(n=300)
    report = M.representation_metrics(
        counts, min_counts=(1, 2, 5), labels={"library": "V629"})
    floors = {m.context.min_count for m in report}
    assert floors == {1, 2, 5}
    names = {m.name for m in report}
    assert {"gini", "shannon_entropy", "pielou_evenness", "simpson_index",
            "inverse_simpson", "hill_q0", "hill_q1", "hill_q2", "top_1_share"} <= names
    df = report.to_dataframe()
    assert len(df) == len(report)
    assert set(df["library"]) == {"V629"}
    for col in ("depth", "min_count", "n_observations", "universe", "grain", "rarefied"):
        assert col in df.columns
    # no row may carry a value without a depth
    assert df["depth"].notna().all()


def test_representation_metrics_adds_a_rarefied_copy_when_asked():
    counts = _skewed_library(n=300)
    target = int(counts.sum()) // 4
    report = M.representation_metrics(counts, min_counts=(1,), rarefy_to=target, replicates=3)
    df = report.to_dataframe()
    assert set(df["rarefied"]) == {True, False}
    assert set(df.loc[df["rarefied"], "rarefied_to"]) == {target}
    # the un-rarefied rows are the ones labelled sample-specific
    assert df.loc[~df["rarefied"], "rarefied_to"].isna().all()


def test_report_get_selects_by_name_and_floor():
    report = M.representation_metrics([1, 2, 3, 4, 100], min_counts=(1, 5))
    assert report.get("gini", min_count=5).context.min_count == 5
    assert report.get("not_a_metric") is None


# ---------------------------------------------------------------------------
# Additive-only proof: the shared kernel reproduces the pipeline-local function
# ---------------------------------------------------------------------------

def test_shared_gini_reproduces_the_existing_guide_donor_bc0_gini_index_exactly():
    """``library_representation.gini_index`` is unchanged; the new kernel agrees with it.

    Both the zeros-included form (the function's documented behaviour) and the zeros-dropped form
    its caller ``plot_gini_indices`` actually computes are checked, so the two conventions already
    in the package both round-trip.
    """
    from magestic.pipelines.guide_donor_bc0.core.library_representation import gini_index

    rng = np.random.default_rng(314)
    for _ in range(10):
        counts = rng.integers(0, 10000, size=500).astype(float)
        assert M.gini(counts, min_count=0).value == pytest.approx(
            gini_index(counts), rel=1e-12, nan_ok=True)
        nonzero = counts[counts > 0]
        assert M.gini(counts, min_count=1).value == pytest.approx(
            gini_index(nonzero), rel=1e-12, nan_ok=True)


def test_metrics_is_importable_from_the_utils_namespace():
    from magestic.utils import gini as gini_from_namespace
    m = gini_from_namespace([1, 1, 1])
    assert isinstance(m, M.Metric)
