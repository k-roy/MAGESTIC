"""``pool_recipe.toml`` is the single source of truth, and it validates itself (NN 506 §2).

**The bug this pins.** The intended composition of the TRIUMPH pool existed only in prose —
in an email, in a design document, and in a sequence of ``--arm`` flags typed into an sbatch
script. Nothing in the build could read it, so nothing could check it, and three regressions
shipped.

**The contract.** One file declares ``total``, the reference and its md5, and a block per
arm with its target. It is read by the **build** and by the **gate**, so the two cannot
drift. A recipe that could not describe a real pool fails at load, before a build consumes
hours: at most one remainder arm, no unknown category names, targets that add up, and no
parent pointing at an arm that was never declared.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import RecipeError, load_recipe, load_recipe_string

from . import _pool_fixtures as fx


def test_the_fixture_recipe_parses_and_derives_the_remainder():
    recipe = fx.recipe()
    assert recipe.name == "SELFTEST_POOL"
    assert recipe.total == 100
    assert recipe.arm_names == (
        "natural", "ptc", "syn", "missense_recoded", "syn_recoded", "scrambled",
    )
    targets = recipe.resolved_targets()
    assert targets == {"natural": 45, "ptc": 10, "syn": 10,
                       "missense_recoded": 20, "syn_recoded": 5, "scrambled": 10}
    assert sum(targets.values()) == recipe.total


def test_the_remainder_target_moves_when_a_declared_target_changes():
    """🔴 Why a remainder arm cannot absorb a loss: its target is derived, not observed."""
    recipe = load_recipe_string(fx.RECIPE_TOML.replace("target = 20", "target = 25"))
    assert recipe.resolved_targets()["natural"] == 40


def test_a_generated_arm_is_recognised_without_being_declared_generated():
    recipe = fx.recipe()
    assert recipe.arm("scrambled").is_generated
    assert recipe.arm("scrambled").generate == "scrambled"
    assert not recipe.arm("ptc").is_generated


def test_the_category_column_is_declared_and_defaults_to_arm_not_subpool():
    """🔴 'subpool' is overloaded, so the module's own default never says it."""
    assert fx.recipe().arm_column == "subpool"  # this screen opted in explicitly
    minimal = load_recipe_string(
        """
        [pool]
        name  = "DEFAULTS"
        total = 1

        [arm.only]
        kind   = "natural_variant"
        target = 1
        """
    )
    assert minimal.arm_column == "arm"


def test_two_remainder_arms_are_refused():
    bad = fx.RECIPE_TOML.replace(
        '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
        '[arm.ptc]\nkind   = "control_ptc"\ntarget = "remainder"',
    )
    with pytest.raises(RecipeError, match="exactly one arm may be 'remainder'"):
        load_recipe_string(bad)


def test_an_unknown_category_name_is_refused():
    bad = fx.RECIPE_TOML.replace('kind   = "control_ptc"', 'kind   = "control_ptcc"')
    with pytest.raises(RecipeError, match="unknown kind"):
        load_recipe_string(bad)


def test_a_missing_target_is_refused():
    bad = fx.RECIPE_TOML.replace(
        '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
        '[arm.ptc]\nkind   = "control_ptc"',
    )
    with pytest.raises(RecipeError, match="'target' is required"):
        load_recipe_string(bad)


def test_targets_that_do_not_add_up_are_refused_when_there_is_no_remainder():
    with pytest.raises(RecipeError, match="never pads to fit"):
        load_recipe_string(
            """
            [pool]
            name  = "SHORT"
            total = 100

            [arm.a]
            kind   = "natural_variant"
            target = 40

            [arm.b]
            kind   = "control_ptc"
            target = 40
            """
        )


def test_declared_targets_exceeding_the_total_are_refused():
    bad = fx.RECIPE_TOML.replace("target = 20", "target = 200")
    with pytest.raises(RecipeError, match="would need a negative count"):
        load_recipe_string(bad)


def test_a_parent_naming_an_undeclared_arm_is_refused():
    bad = fx.RECIPE_TOML.replace('parent = "missense_recoded"', 'parent = "ghost_arm"')
    with pytest.raises(RecipeError, match="parent 'ghost_arm' is not a declared arm"):
        load_recipe_string(bad)


def test_an_arm_cannot_be_its_own_parent():
    bad = fx.RECIPE_TOML.replace('parent = "missense_recoded"', 'parent = "syn_recoded"')
    with pytest.raises(RecipeError, match="cannot be its own parent"):
        load_recipe_string(bad)


def test_a_parent_key_without_a_parent_is_refused():
    with pytest.raises(RecipeError, match="parent_key is set but no parent"):
        load_recipe_string(
            """
            [pool]
            name  = "ORPHAN_KEY"
            total = 10

            [arm.only]
            kind       = "natural_variant"
            target     = 10
            parent_key = "gene"
            """
        )


def test_a_designed_kind_cannot_be_declared_generated():
    with pytest.raises(RecipeError, match="cannot be generated by the module"):
        load_recipe_string(
            """
            [pool]
            name  = "FAKE_GENERATOR"
            total = 10

            [arm.only]
            kind     = "natural_variant"
            target   = 10
            generate = "scrambled"
            """
        )


def test_an_unknown_select_policy_is_refused():
    bad = fx.RECIPE_TOML.replace(
        '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
        '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10\nselect = "whatever"',
    )
    with pytest.raises(RecipeError, match="select must be one of"):
        load_recipe_string(bad)


def test_a_recipe_with_no_arms_is_refused():
    with pytest.raises(RecipeError, match="declares nothing"):
        load_recipe_string('[pool]\nname = "EMPTY"\ntotal = 10\n')


def test_a_recipe_with_no_total_is_refused():
    with pytest.raises(RecipeError, match="total is required"):
        load_recipe_string(
            '[pool]\nname = "NO_TOTAL"\n\n[arm.a]\nkind = "natural_variant"\ntarget = 1\n'
        )


def test_the_schema_block_round_trips_columns_and_groups():
    schema = fx.recipe().schema
    assert "oligo_id" in schema.columns
    assert schema.require_groups == ("recut_cfd",)
    group = schema.groups["recut_cfd"]
    assert group.prefix == "recut_cfd_" and group.min_count == 3


def test_segments_parse_including_a_negative_index():
    segments = fx.recipe().segments
    assert segments["forward_priming"] == (0, 20)
    assert segments["donor"] == (51, -20)
    assert segments["reverse_priming"] == (-20, 0)


def test_a_malformed_segment_is_refused():
    bad = fx.RECIPE_TOML.replace("forward_priming = [0, 20]", "forward_priming = [0]")
    with pytest.raises(RecipeError, match="two-element"):
        load_recipe_string(bad)


def test_load_recipe_reads_a_real_file_and_records_its_source(tmp_path):
    path = tmp_path / "pool_recipe.toml"
    path.write_text(fx.RECIPE_TOML)
    recipe = load_recipe(path)
    assert recipe.source == path
    assert recipe.total == 100
