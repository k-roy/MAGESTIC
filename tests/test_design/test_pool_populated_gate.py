"""G-POPULATED: presence is not population.

🔴 The defect this gate exists for ALREADY SHIPPED. The TRIUMPH v6 pool carried `aa_change` in its
header and empty on all 1,821 of its natural missense designs, where Library A populated all 1,894
(NN 512 finding 1, reproduced independently by the NN 516 cross-version diff). G-SCHEMA passed it,
because the column was right there. The first test below plants exactly that defect and requires
the gate to catch it; the second is its positive control, the same pool with the column populated.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import (
    RecipeError,
    gate_populated,
    gate_schema,
    load_recipe_string,
)

TOML = """
[pool]
name       = "POP_POOL"
total      = 6
arm_column = "subpool"

[arm.natural]
kind   = "natural_variant"
target = "remainder"

[arm.ptc]
kind   = "control_ptc"
target = 2

[schema]
columns = ["oligo_id", "subpool", "consequence", "aa_change", "recut_cfd"]

[schema.populated.aa_on_missense]
columns = ["aa_change"]
where   = { consequence = "missense" }

[schema.populated.recut_on_targeting]
columns = ["recut_cfd"]
where   = { subpool = "natural" }
"""


def pool(aa_on_natural_missense: str, cfd: str = "0.5"):
    """Four natural rows (2 missense, 2 non-coding) and two ptc rows."""
    rows = []
    for i in range(2):
        rows.append({"oligo_id": f"n{i}", "subpool": "natural", "consequence": "missense",
                     "aa_change": aa_on_natural_missense, "recut_cfd": cfd})
    for i in range(2):
        rows.append({"oligo_id": f"c{i}", "subpool": "natural", "consequence": "non_coding",
                     "aa_change": "", "recut_cfd": cfd})
    for i in range(2):
        rows.append({"oligo_id": f"p{i}", "subpool": "ptc", "consequence": "missense",
                     "aa_change": "W9*", "recut_cfd": ""})
    return rows


def test_the_v6_defect_is_caught():
    """🔴 `aa_change` present but empty on every natural missense design."""
    r = load_recipe_string(TOML)
    g = gate_populated(pool(""), r)
    assert g.failed
    assert "aa_change" in g.detail
    # the numbers are in the detail, so the report says how bad it is
    assert "2 of 4" in g.detail or "2 of 4 rows" in g.detail, g.detail


def test_positive_control_the_same_pool_populated_passes():
    r = load_recipe_string(TOML)
    g = gate_populated(pool("S73N"), r)
    assert not g.failed, g.detail


def test_g_schema_does_NOT_catch_it_which_is_why_this_gate_exists():
    """The companion assertion: the defect sails past the presence check."""
    r = load_recipe_string(TOML)
    g = gate_schema(pool(""), r)
    assert not g.failed, "G-SCHEMA should pass — the column IS present; that is the hole"


def test_the_scope_is_respected_so_a_meaningless_column_is_not_demanded():
    """Non-coding designs legitimately have no aa_change, and ptc rows no recut_cfd here."""
    r = load_recipe_string(TOML)
    g = gate_populated(pool("S73N"), r)
    assert not g.failed
    # ptc rows have an empty recut_cfd, and the rule is scoped to subpool=natural, so it passes
    assert any(row[0] == "recut_on_targeting" for row in g.rows)


def test_a_partial_loss_is_caught_not_only_a_total_one():
    r = load_recipe_string(TOML)
    rows = pool("S73N")
    rows[0]["aa_change"] = ""        # one of the two missense naturals goes blank
    g = gate_populated(rows, r)
    assert g.failed and "aa_change" in g.detail


def test_min_fraction_allows_a_declared_tolerance():
    toml = TOML.replace(
        'columns = ["aa_change"]\nwhere   = { consequence = "missense" }',
        'columns = ["aa_change"]\nwhere   = { consequence = "missense" }\nmin_fraction = 0.5')
    r = load_recipe_string(toml)
    rows = pool("S73N")
    rows[0]["aa_change"] = ""        # 3 of 4 populated = 75% >= 50%
    assert not gate_populated(rows, r).failed
    # the scope (consequence == missense) is 4 rows: the 2 natural missense AND the 2 ptc.
    # Leaving 2 of 4 populated is EXACTLY 50% and must pass; 1 of 4 is 25% and must fail.
    for row in rows:
        row["aa_change"] = ""
    rows[4]["aa_change"] = "W9*"
    rows[5]["aa_change"] = "W9*"
    assert not gate_populated(rows, r).failed, "2/4 is exactly min_fraction, so it passes"
    rows[5]["aa_change"] = ""        # now 1 of 4 = 25%, below the declared 50%
    assert gate_populated(rows, r).failed


def test_a_rule_that_matches_no_rows_is_itself_a_defect():
    """🔴 A population rule that selects nothing checks nothing, and must not read as a pass."""
    toml = TOML.replace('where   = { consequence = "missense" }',
                        'where   = { consequence = "nonsense_typo" }')
    r = load_recipe_string(toml)
    g = gate_populated(pool("S73N"), r)
    assert g.failed
    assert "matched NO rows" in g.detail


def test_whitespace_only_is_not_populated():
    r = load_recipe_string(TOML)
    assert gate_populated(pool("   "), r).failed


def test_no_rules_is_not_applicable_rather_than_a_pass():
    toml = TOML.split("[schema.populated.aa_on_missense]")[0]
    r = load_recipe_string(toml)
    g = gate_populated(pool(""), r)
    assert g.status.lower() in ("n/a", "not_applicable", "skipped"), g.status
    assert not g.failed


def test_a_rule_without_columns_is_refused():
    toml = TOML.replace('columns = ["aa_change"]\nwhere   = { consequence = "missense" }',
                        'where   = { consequence = "missense" }')
    with pytest.raises(RecipeError, match="needs 'columns'"):
        load_recipe_string(toml)


def test_an_out_of_range_min_fraction_is_refused():
    toml = TOML.replace('columns = ["aa_change"]',
                        'columns = ["aa_change"]\nmin_fraction = 0')
    with pytest.raises(RecipeError, match="min_fraction"):
        load_recipe_string(toml)


def test_the_list_shorthand_works():
    toml = TOML.replace(
        '[schema.populated.aa_on_missense]\ncolumns = ["aa_change"]\n'
        'where   = { consequence = "missense" }',
        '[schema.populated]\nalways = ["oligo_id"]')
    r = load_recipe_string(toml)
    assert "always" in r.schema.populated
    assert not gate_populated(pool("S73N"), r).failed
