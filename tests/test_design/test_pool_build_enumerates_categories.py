"""``build_pool`` enumerates the recipe's categories itself (NN 506 §3, NN 507 §1).

**The bug this pins.** Library A ran a fifth merge step that produced 252 ``syn_recoded``
recoding-helper controls. The polished-frame pipeline packs four arms, the fifth step was
never carried across, and the arm vanished from a pool that was uploaded, announced to
collaborators at EMBL and had its checksum quoted in email. A collaborator found it by eye.

Nothing could have caught it, because the composition existed only as a sequence of
``--arm`` flags typed into an sbatch script: the arm was never named anywhere a checker
could read.

**The contract.** 🔴 ``build_pool`` enumerates the categories from the **recipe**, not from
whatever the caller passed. A ``kind`` the recipe names but the caller does not supply is an
**error, not an omission**. The symmetric rule holds too — an arm the caller supplies that
the recipe does not declare is an error — so a category can neither leave nor enter a pool
without being written down.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import (
    GeneratedArmSuppliedError,
    MissingCategoryError,
    PoolGateError,
    UndeclaredArmError,
    build_pool,
)

from . import _pool_fixtures as fx


def test_the_happy_path_builds_and_every_gate_passes():
    pool, report = build_pool(
        fx.recipe(),
        fx.designed_arms(),
        architecture=fx.ARCHITECTURE,
        row_template={"gene": "", "primer_pair": "revP1",
                      "recut_cfd_max": "", "recut_cfd_sum": "", "recut_cfd_n": ""},
    )

    assert len(pool) == 100
    assert report.passed
    assert report.total_built == report.total_target == 100
    built = {a.name: a.built for a in report.arms}
    assert built == {"natural": 45, "ptc": 10, "syn": 10,
                     "missense_recoded": 20, "syn_recoded": 5, "scrambled": 10}


def test_omitting_a_declared_category_raises_rather_than_shipping_without_it():
    """🔴 The NN 505 regression, planted: the caller forgets ``syn_recoded``.

    The old pipeline would have packed the four arms it was given and produced a pool. This
    must refuse to build.
    """
    arms = fx.designed_arms()
    del arms["syn_recoded"]

    with pytest.raises(MissingCategoryError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)

    message = str(exc.value)
    assert "syn_recoded" in message
    assert "control_recoding_helper" in message
    assert "target 5" in message
    assert "ERROR, not an omission" in message


def test_omitting_several_categories_names_all_of_them():
    arms = fx.designed_arms()
    del arms["syn_recoded"]
    del arms["ptc"]

    with pytest.raises(MissingCategoryError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)

    assert "syn_recoded" in str(exc.value)
    assert "ptc" in str(exc.value)
    assert "declares 2 categories" in str(exc.value)


def test_an_empty_supplied_category_is_a_gate_failure_not_a_silent_absence():
    """Supplying an arm with zero rows is different from omitting it — and still fails."""
    arms = fx.designed_arms(syn_recoded=0, natural=50)

    with pytest.raises(PoolGateError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)

    assert "G-COMPOSITION" in str(exc.value)
    assert "syn_recoded" in str(exc.value)


def test_supplying_an_undeclared_arm_raises():
    arms = fx.designed_arms()
    arms["mystery_controls"] = fx.make_rows("mystery_controls", 5)

    with pytest.raises(UndeclaredArmError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)

    assert "mystery_controls" in str(exc.value)
    assert "not declared" in str(exc.value)


def test_supplying_rows_for_a_generated_arm_raises():
    """Two things would otherwise decide that arm's composition."""
    arms = fx.designed_arms()
    arms["scrambled"] = fx.make_rows("scrambled", 10)

    with pytest.raises(GeneratedArmSuppliedError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)

    assert "scrambled" in str(exc.value)


def test_a_generated_arm_without_an_architecture_raises_rather_than_skipping():
    with pytest.raises(Exception) as exc:
        build_pool(fx.recipe(), fx.designed_arms())
    assert "OligoArchitecture" in str(exc.value)


def test_drift_that_keeps_the_total_right_still_fails_the_build():
    """🔴 End to end: the 454-oligo drift shape, through ``build_pool``."""
    arms = fx.designed_arms(ptc=6, syn=7, missense_recoded=18, natural=54)

    with pytest.raises(PoolGateError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)

    message = str(exc.value)
    assert "G-COMPOSITION" in message
    for arm in ("ptc", "syn", "missense_recoded"):
        assert arm in message


def test_on_gate_failure_report_returns_the_verdicts_instead_of_raising():
    """A caller may collect the verdicts, but the failure is still recorded, never dropped."""
    arms = fx.designed_arms(ptc=6, natural=49)
    pool, report = build_pool(
        fx.recipe(), arms, architecture=fx.ARCHITECTURE, on_gate_failure="report"
    )

    assert not report.passed
    failed = [g.gate for g in report.gates if g.failed]
    assert "G-COMPOSITION" in failed
    assert report.to_dict()["passed"] is False


def test_the_report_carries_target_built_delta_and_drop_reasons():
    pool, report = build_pool(
        fx.recipe(),
        fx.designed_arms(),
        architecture=fx.ARCHITECTURE,
        drops={"natural": {"no_window_contains_the_edit": 622,
                           "inside_another_Verified_CDS": 18}},
        row_template={"gene": "", "primer_pair": "revP1",
                      "recut_cfd_max": "", "recut_cfd_sum": "", "recut_cfd_n": ""},
    )

    natural = next(a for a in report.arms if a.name == "natural")
    assert natural.target == 45 and natural.built == 45 and natural.delta == 0
    assert natural.drops["no_window_contains_the_edit"] == 622
    scrambled = next(a for a in report.arms if a.name == "scrambled")
    assert scrambled.generated is True
    assert scrambled.extra["ks_statistic_vs_reference"] is not None


def test_the_arm_column_is_stamped_and_is_not_called_subpool_idx():
    """🔴 'subpool' is overloaded; the module writes a CATEGORY, never a synthesis index."""
    pool, _ = build_pool(
        fx.recipe(),
        fx.designed_arms(),
        architecture=fx.ARCHITECTURE,
        row_template={"gene": "", "primer_pair": "revP1",
                      "recut_cfd_max": "", "recut_cfd_sum": "", "recut_cfd_n": ""},
    )

    assert {row["subpool"] for row in pool} == {
        "natural", "ptc", "syn", "missense_recoded", "syn_recoded", "scrambled",
    }
    assert all("subpool_idx" not in row for row in pool)


def test_a_declared_select_policy_trims_and_records_it():
    """Trimming is possible only because a human declared it, and it is counted."""
    from magestic.design.pool import load_recipe_string

    recipe = load_recipe_string(
        fx.RECIPE_TOML.replace(
            '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
            '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10\nselect = "head"',
        )
    )
    arms = fx.designed_arms(ptc=14)

    pool, report = build_pool(recipe, arms, architecture=fx.ARCHITECTURE,
                              row_template={"gene": "", "primer_pair": "revP1",
                                            "recut_cfd_max": "", "recut_cfd_sum": "",
                                            "recut_cfd_n": ""})

    ptc = next(a for a in report.arms if a.name == "ptc")
    assert ptc.supplied == 14 and ptc.built == 10 and ptc.trimmed == 4
    assert report.passed


def test_the_default_select_policy_does_not_trim():
    """🔴 Without a declared policy an over-supplied arm FAILS; it is never quietly cut."""
    from magestic.design.pool import PoolGateError

    arms = fx.designed_arms(ptc=14)
    with pytest.raises(PoolGateError) as exc:
        build_pool(fx.recipe(), arms, architecture=fx.ARCHITECTURE)
    assert "target 10, built 14" in str(exc.value)
