"""The declared order interface: a ``fill`` arm, ranked trimming, ``pct``, and ``enabled``.

NN 510/512. Kevin, 2026-11-11: *"Ideally the order can be assembled simply by specifying the
number/% of scrambled controls, what categories to include ... and how to fill any remaining
slots (default is 2nd and 3rd shots on goal for most likely recut targets)."*

🔴 The risk this whole file is about: a ``fill`` arm is allowed to be trimmed to fit, and "an arm
that may be trimmed to fit" is uncomfortably close to the NN 505 failure, where the total stayed at
exactly 24,000 while a whole category vanished into another arm. The tests below pin the
difference — **every other arm is still checked exactly**, an under-supplied fill arm still fails,
and the trim is by a declared rule rather than by arrival order.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import (
    PoolGateError,
    RecipeError,
    build_pool,
    gate_composition,
    load_recipe_string,
)

# natural 40 declared + extra = fill + ptc 10 + scrambled 10 -> total 100, so fill = 40
FILL_TOML = """
[pool]
name            = "FILL_POOL"
total           = 100
arm_column      = "subpool"

[arm.natural]
kind   = "natural_variant"
target = 40

[arm.ptc]
kind   = "control_ptc"
target = 10

[arm.scrambled]
kind     = "control_nontargeting"
pct      = 10.0
gc_match = "designed_pool"

[arm.natural_extra]
kind    = "natural_extra_shots"
target  = "fill"
select  = "rank"
rank_by = ["-recut_risk", "shot_index"]
"""


#: 🔴 A realistic ~50% GC oligo, NOT "A"*200. The non-targeting generator matches the GC
#: DISTRIBUTION of the designed rows, so an all-A reference asks it for 0%-GC 200-mers that also
#: carry the fixed GCTCTTC site — impossible, and it spins to its retry ceiling on every call.
#: A degenerate fixture that makes a generator hang looks exactly like a hung test.
SEQ = "ACGT" * 50


def rows(arm, n, start=0, **kw):
    out = []
    for i in range(n):
        r = {"subpool": arm, "oligo_id": f"{arm}_{start + i:04d}", "oligo_seq": SEQ}
        r.update({k: (v[i] if isinstance(v, list) else v) for k, v in kw.items()})
        out.append(r)
    return out


def test_fill_target_is_derived_from_the_total():
    r = load_recipe_string(FILL_TOML)
    t = r.resolved_targets()
    assert t == {"natural": 40, "ptc": 10, "scrambled": 10, "natural_extra": 40}
    assert sum(t.values()) == 100
    assert r.fill_arm is not None and r.fill_arm.name == "natural_extra"
    assert r.remainder_arm is None
    assert r.derived_arm is r.fill_arm


def test_pct_resolves_against_the_total():
    r = load_recipe_string(FILL_TOML)
    assert r.arm("scrambled").target == 10        # 10% of 100
    assert r.arm("scrambled").pct == 10.0


def test_pct_and_target_together_are_refused():
    bad = FILL_TOML.replace("pct      = 10.0", "pct      = 10.0\ntarget   = 10")
    with pytest.raises(RecipeError, match="either 'target' or 'pct'"):
        load_recipe_string(bad)


def test_fill_without_rank_is_refused():
    bad = FILL_TOML.replace('select  = "rank"\nrank_by = ["-recut_risk", "shot_index"]', "")
    with pytest.raises(RecipeError, match="requires select='rank'"):
        load_recipe_string(bad)


def test_rank_without_rank_by_is_refused():
    bad = FILL_TOML.replace('rank_by = ["-recut_risk", "shot_index"]', "")
    with pytest.raises(RecipeError, match="needs 'rank_by'"):
        load_recipe_string(bad)


def test_a_fill_and_a_remainder_arm_together_are_refused():
    bad = FILL_TOML.replace('[arm.natural]\nkind   = "natural_variant"\ntarget = 40',
                            '[arm.natural]\nkind   = "natural_variant"\ntarget = "remainder"')
    with pytest.raises(RecipeError, match="exactly one arm may be derived"):
        load_recipe_string(bad)


# ---------------------------------------------------------------------------
# the ranked trim
# ---------------------------------------------------------------------------
def _build(extra_n, **kw):
    r = load_recipe_string(FILL_TOML)
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20, internal_cloning_site="GCTCTTCAAAA"[:11],
                             reverse_tail="T" * 20)
    arms = {
        "natural": rows("natural", 40),
        "ptc": rows("ptc", 10),
        "natural_extra": rows(
            "natural_extra", extra_n,
            recut_risk=[round(i / extra_n, 4) for i in range(extra_n)],
            shot_index=[2 + (i % 2) for i in range(extra_n)],
        ),
    }
    return build_pool(r, arms, architecture=arch,
                      row_template={"recut_risk": "", "shot_index": ""},
                      on_gate_failure="report", **kw)


def test_an_oversupplied_fill_arm_is_trimmed_to_its_target():
    pool, rep = _build(60)
    assert len(pool) == 100
    fill = [a for a in rep.arms if a.name == "natural_extra"][0]
    assert (fill.supplied, fill.built, fill.target, fill.trimmed) == (60, 40, 40, 20)
    assert rep.passed, [g.detail for g in rep.gates if g.failed]


def test_the_trim_keeps_the_highest_ranked_rows():
    """🔴 `-recut_risk` descending: the 40 kept must be the 40 riskiest, not the first 40."""
    pool, rep = _build(60)
    kept = {r["oligo_id"] for r in pool if r["subpool"] == "natural_extra"}
    # risk was i/60 for i in 0..59, so the riskiest 40 are i = 20..59
    assert kept == {f"natural_extra_{i:04d}" for i in range(20, 60)}
    fill = [a for a in rep.arms if a.name == "natural_extra"][0]
    # the boundary is recorded, so the cut is readable afterwards
    assert fill.extra["rank_by"] == ["-recut_risk", "shot_index"]
    assert fill.extra["dropped"] == 20
    assert float(fill.extra["last_kept"]["recut_risk"]) >= \
        float(fill.extra["first_dropped"]["recut_risk"])


def test_the_trim_is_numeric_not_lexicographic():
    """"0.9" must outrank "0.12"; a string sort would invert them."""
    r = load_recipe_string(FILL_TOML.replace("total           = 100", "total           = 61"))
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    # total 61 => fill target = 61 - 40 - 10 - 6 (10% of 61 = 6) = 5
    arms = {
        "natural": rows("natural", 40),
        "ptc": rows("ptc", 10),
        "natural_extra": rows("natural_extra", 6,
                              recut_risk=["0.9", "0.12", "0.85", "0.2", "0.3", "0.05"],
                              shot_index=2),
    }
    pool, rep = build_pool(r, arms, architecture=arch,
                           row_template={"recut_risk": "", "shot_index": ""},
                           on_gate_failure="report")
    kept = {r["oligo_id"] for r in pool if r["subpool"] == "natural_extra"}
    # the five kept are everything except the smallest, 0.05
    assert "natural_extra_0005" not in kept
    assert "natural_extra_0000" in kept and "natural_extra_0002" in kept


def test_the_trim_is_deterministic_across_input_order():
    """The same designs in a different order must yield the SAME pool."""
    import random
    pool_a, _ = _build(60)
    r = load_recipe_string(FILL_TOML)
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    extra = rows("natural_extra", 60,
                 recut_risk=[round(i / 60, 4) for i in range(60)],
                 shot_index=[2 + (i % 2) for i in range(60)])
    random.Random(7).shuffle(extra)
    pool_b, _ = build_pool(r, {"natural": rows("natural", 40), "ptc": rows("ptc", 10),
                               "natural_extra": extra},
                           architecture=arch,
                           row_template={"recut_risk": "", "shot_index": ""},
                           on_gate_failure="report")
    ids = lambda p: sorted(x["oligo_id"] for x in p if x["subpool"] == "natural_extra")
    assert ids(pool_a) == ids(pool_b)


def test_ties_fall_back_to_the_id_so_a_rebuild_reproduces():
    r = load_recipe_string(FILL_TOML)
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    # every row identical on both rank keys -> only the id breaks the tie
    extra = rows("natural_extra", 50, recut_risk=0.5, shot_index=2)
    pool, _ = build_pool(r, {"natural": rows("natural", 40), "ptc": rows("ptc", 10),
                             "natural_extra": extra}, architecture=arch,
                         row_template={"recut_risk": "", "shot_index": ""},
                         on_gate_failure="report")
    kept = sorted(x["oligo_id"] for x in pool if x["subpool"] == "natural_extra")
    assert kept == [f"natural_extra_{i:04d}" for i in range(40)]


def test_a_missing_rank_value_sorts_last_and_does_not_crash():
    r = load_recipe_string(FILL_TOML)
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    extra = rows("natural_extra", 41, recut_risk=[""] + [0.5] * 40, shot_index=2)
    pool, _ = build_pool(r, {"natural": rows("natural", 40), "ptc": rows("ptc", 10),
                             "natural_extra": extra}, architecture=arch,
                         row_template={"recut_risk": "", "shot_index": ""},
                         on_gate_failure="report")
    kept = {x["oligo_id"] for x in pool if x["subpool"] == "natural_extra"}
    assert "natural_extra_0000" not in kept       # the blank one is the one dropped


# ---------------------------------------------------------------------------
# 🔴 what a fill arm must NOT be allowed to hide
# ---------------------------------------------------------------------------
def test_an_undersupplied_fill_arm_still_fails_the_gate():
    """A fill arm may be trimmed DOWN; it may never be padded or excused."""
    pool, rep = _build(30)
    assert len(pool) == 90
    g = [x for x in rep.gates if x.gate == "G-COMPOSITION"][0]
    assert g.failed
    assert "natural_extra" in g.detail and "90" in g.detail


def test_a_loss_in_another_arm_is_not_absorbed_by_the_fill_arm():
    """🔴 THE NN 505 TEST. ptc short by 4 must fail even though fill could cover the slots."""
    r = load_recipe_string(FILL_TOML)
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    pool, rep = build_pool(
        r,
        {"natural": rows("natural", 40), "ptc": rows("ptc", 6),          # 4 short
         "natural_extra": rows("natural_extra", 60, recut_risk=0.5, shot_index=2)},
        architecture=arch, row_template={"recut_risk": "", "shot_index": ""},
        on_gate_failure="report",
    )
    g = [x for x in rep.gates if x.gate == "G-COMPOSITION"][0]
    assert g.failed
    assert "'ptc'" in g.detail and "target 10, built 6" in g.detail
    # and the fill arm did NOT grow to cover it
    assert sum(1 for x in pool if x["subpool"] == "natural_extra") == 40


# ---------------------------------------------------------------------------
# enabled = false
# ---------------------------------------------------------------------------
def test_a_disabled_category_is_excluded_and_reported_not_silent():
    toml = FILL_TOML.replace('[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
                             '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10\nenabled = false')
    r = load_recipe_string(toml)
    assert "ptc" not in r.arm_names
    assert r.disabled_arms == ("ptc",)
    # the fill arm absorbs ptc's slots, because the recipe says ptc is not wanted
    assert r.resolved_targets() == {"natural": 40, "scrambled": 10, "natural_extra": 50}
    from magestic.design.pool import OligoArchitecture
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    pool, rep = build_pool(r, {"natural": rows("natural", 40),
                               "natural_extra": rows("natural_extra", 60, recut_risk=0.5,
                                                     shot_index=2)},
                           architecture=arch,
                           row_template={"recut_risk": "", "shot_index": ""},
                           on_gate_failure="report")
    assert len(pool) == 100 and rep.passed
    # 🔴 the exclusion is stated in the report; this is the whole difference from NN 505
    assert any("DISABLE" in n.upper() and "ptc" in n for n in rep.notes), rep.notes


def test_supplying_a_disabled_arm_is_still_an_error():
    """A category the recipe switched off must not sneak back in from the caller."""
    toml = FILL_TOML.replace('[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
                             '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10\nenabled = false')
    r = load_recipe_string(toml)
    from magestic.design.pool import OligoArchitecture, UndeclaredArmError
    arch = OligoArchitecture(forward_priming="A" * 20,
                             internal_cloning_site="GCTCTTCAAAA"[:11], reverse_tail="T" * 20)
    with pytest.raises(UndeclaredArmError, match="ptc"):
        build_pool(r, {"natural": rows("natural", 40), "ptc": rows("ptc", 10),
                       "natural_extra": rows("natural_extra", 50, recut_risk=0.5, shot_index=2)},
                   architecture=arch, on_gate_failure="report")
