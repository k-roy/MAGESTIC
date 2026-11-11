"""G-SCHEMA: the output must carry every column the recipe names (NN 505, NN 506 §4).

**The bug this pins.** The TRIUMPH rebuild wrote its own 27-column table where the previous
pool had 51 columns. **24 columns vanished, including all 15 re-cut scores** — 47 % of the
schema. Nothing compared the shipped schema against an intended one, so the loss was found
only by diffing two pools after the pool had been uploaded and announced.

**The contract.** Every column in ``[schema] columns`` and every group in
``require_groups`` must be present in the output. A column-group is how the 15 re-cut score
columns are pinned without naming each one: ``prefix = "recut_cfd_"`` with
``min_count = 15`` fails when the group thins out, not only when it disappears.

🔴 A recipe that declares no schema reports ``not_applicable`` **explicitly** rather than
passing. A silent skip is the failure mode this package exists to end.
"""

from __future__ import annotations

from magestic.design.pool import gate_schema, load_recipe_string
from magestic.design.pool.gates import FAIL, NOT_APPLICABLE, PASS

from . import _pool_fixtures as fx


def test_a_full_schema_passes():
    result = gate_schema(fx.built_pool(), fx.recipe())
    assert result.status == PASS, result.detail


def test_a_dropped_column_is_caught():
    """Defect planted: the writer stops emitting a declared column."""
    pool = fx.built_pool()
    for row in pool:
        del row["gene"]

    result = gate_schema(pool, fx.recipe())

    assert result.status == FAIL
    assert "gene" in result.detail
    assert "missing" in result.detail


def test_the_recut_score_group_thinning_out_is_caught():
    """🔴 The NN 505 defect proper: the re-cut score columns disappear.

    The recipe requires at least three ``recut_cfd_*`` columns (15 in the real pool). The
    rebuild here keeps one, which a plain column list would not have noticed if the names
    had changed.
    """
    pool = fx.built_pool()
    for row in pool:
        del row["recut_cfd_sum"]
        del row["recut_cfd_n"]

    result = gate_schema(pool, fx.recipe())

    assert result.status == FAIL
    assert "recut_cfd" in result.detail
    assert "1 present, 3 required" in result.detail


def test_all_the_recut_scores_vanishing_is_caught():
    pool = fx.built_pool()
    for row in pool:
        for col in ("recut_cfd_max", "recut_cfd_sum", "recut_cfd_n"):
            del row[col]

    result = gate_schema(pool, fx.recipe())

    assert result.status == FAIL
    assert "0 present, 3 required" in result.detail


def test_a_rebuild_that_loses_half_the_schema_names_every_loss():
    """The 24-of-51 shape: several columns and a whole group go at once."""
    pool = fx.built_pool()
    for row in pool:
        for col in ("gene", "primer_pair", "recut_cfd_max", "recut_cfd_sum", "recut_cfd_n"):
            del row[col]

    result = gate_schema(pool, fx.recipe())

    assert result.status == FAIL
    assert "gene" in result.detail and "primer_pair" in result.detail
    assert "recut_cfd" in result.detail
    # the failing rows are the artifact a reader works from
    assert len(result.rows) >= 3


def test_columns_can_be_passed_explicitly_for_a_sparse_writer():
    """A writer's intended header is authoritative even if some rows omit a key."""
    pool = fx.built_pool()
    for row in pool:
        del row["gene"]

    header = ["oligo_id", "subpool", "gene", "oligo_seq", "primer_pair",
              "recut_cfd_max", "recut_cfd_sum", "recut_cfd_n"]
    assert gate_schema(pool, fx.recipe(), columns=header).status == PASS


def test_no_declared_schema_is_not_applicable_rather_than_a_silent_pass():
    recipe = load_recipe_string(
        """
        [pool]
        name  = "NO_SCHEMA"
        total = 10

        [arm.only]
        kind   = "natural_variant"
        target = 10
        """
    )
    result = gate_schema([{"oligo_id": "x", "arm": "only"}], recipe)

    assert result.status == NOT_APPLICABLE
    assert "declares no [schema]" in result.detail


def test_a_required_group_with_no_definition_is_a_failure_not_a_pass():
    """A require_groups entry is validated at load; this pins that it cannot slip through."""
    import pytest

    from magestic.design.pool import RecipeError

    with pytest.raises(RecipeError, match="require_groups names 'ghost'"):
        load_recipe_string(
            """
            [pool]
            name  = "GHOST_GROUP"
            total = 10

            [arm.only]
            kind   = "natural_variant"
            target = 10

            [schema]
            require_groups = ["ghost"]
            """
        )
