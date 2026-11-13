"""G-PARENT and G-COAMPLIFY: a control must sit beside the thing it controls (NN 506 §4).

**The bugs these pin.**

*G-PARENT.* The ``syn_recoded`` arm is a **recoding helper**: 252 synonymous-recode controls
that exist only to sit beside the recoded missense designs and show what the recoding itself
costs. A helper whose parent designs are absent controls for nothing — it is 252 oligos of
synthesis paying for an answer nobody can read. G-COMPOSITION catches the helper going
missing; this gate catches the opposite, a helper that survives a rebuild its parent did
not.

*G-COAMPLIFY.* An earlier library shipped controls that did not share a reverse priming site
with the designs they controlled, so the two were measured in different reactions and the
comparison the controls existed for could not be made.

**The contract.** Every control arm with a declared ``parent`` must find that parent present
in the pool — row by row when the arm declares a ``parent_key``, arm by arm otherwise. And
every control must share a reverse priming site with its parent, unless the whole pool has
one site. 🔴 An arm is exempt from G-COAMPLIFY only by declaring ``coamplify = false``; an
exemption is never inferred.
"""

from __future__ import annotations

from magestic.design.pool import gate_coamplify, gate_parent, load_recipe_string
from magestic.design.pool.gates import FAIL, NOT_APPLICABLE, PASS

from . import _pool_fixtures as fx


# ---------------------------------------------------------------------------
# G-PARENT
# ---------------------------------------------------------------------------

def test_a_correct_pool_passes_parent():
    result = gate_parent(fx.built_pool(), fx.recipe())
    assert result.status == PASS, result.detail
    assert "1 parented arm" in result.detail


def test_a_helper_whose_parent_vanished_is_caught():
    """Defect planted: ``missense_recoded`` is gone, the helpers survive."""
    pool = fx.built_pool(missense_recoded=0, natural=65)
    result = gate_parent(pool, fx.recipe())

    assert result.status == FAIL
    assert "syn_recoded" in result.detail
    assert "missense_recoded" in result.detail
    assert "controls for nothing" in result.detail


def test_an_absent_child_is_left_to_the_composition_gate():
    """G-PARENT does not double-report a missing arm; that is G-COMPOSITION's finding."""
    pool = fx.built_pool(syn_recoded=0, natural=50)
    result = gate_parent(pool, fx.recipe())

    assert result.status == PASS
    assert any(row[4] == "child absent" for row in result.rows)


def test_row_level_orphans_are_caught_when_parent_key_is_declared():
    """With a ``parent_key`` the gate checks each helper against a real parent design."""
    recipe = load_recipe_string(
        fx.RECIPE_TOML.replace(
            'target = 5\nparent = "missense_recoded"',
            'target = 5\nparent = "missense_recoded"\nparent_key = "gene"',
        )
    )
    pool = fx.built_pool()
    # start from a pool where every helper DOES have a parent design...
    parent_genes = sorted({r["gene"] for r in pool if r["subpool"] == "missense_recoded"})
    helpers = [r for r in pool if r["subpool"] == "syn_recoded"]
    for i, r in enumerate(helpers):
        r["gene"] = parent_genes[i % len(parent_genes)]
    # ...then plant three pointing at a gene no recoded design covers
    for r in helpers[:3]:
        r["gene"] = "YGHOST999W"

    result = gate_parent(pool, recipe)

    assert result.status == FAIL
    assert "3 of 5 rows" in result.detail
    assert "YGHOST999W" not in result.detail  # it reports the oligo ids, not the key
    assert "syn_recoded_00000" in result.detail


def test_row_level_check_passes_when_every_helper_has_a_parent_design():
    recipe = load_recipe_string(
        fx.RECIPE_TOML.replace(
            'target = 5\nparent = "missense_recoded"',
            'target = 5\nparent = "missense_recoded"\nparent_key = "gene"',
        )
    )
    pool = fx.built_pool()
    parent_genes = sorted({r["gene"] for r in pool if r["subpool"] == "missense_recoded"})
    for i, r in enumerate([r for r in pool if r["subpool"] == "syn_recoded"]):
        r["gene"] = parent_genes[i % len(parent_genes)]

    result = gate_parent(pool, recipe)
    assert result.status == PASS, result.detail
    assert "row-by-row" in result.detail


def test_no_parented_arm_is_not_applicable():
    recipe = load_recipe_string(
        """
        [pool]
        name  = "NO_PARENTS"
        total = 10

        [arm.only]
        kind   = "natural_variant"
        target = 10
        """
    )
    result = gate_parent([{"oligo_id": "x", "arm": "only"}], recipe)
    assert result.status == NOT_APPLICABLE


# ---------------------------------------------------------------------------
# G-COAMPLIFY
# ---------------------------------------------------------------------------

def test_a_single_primer_pool_passes_coamplify():
    result = gate_coamplify(fx.built_pool(), fx.recipe())
    assert result.status == PASS
    assert "one reverse priming site" in result.detail


def test_a_control_on_a_different_primer_from_its_parent_is_caught():
    """Defect planted: the helpers get their own priming site."""
    pool = fx.built_pool()
    for r in pool:
        if r["subpool"] == "syn_recoded":
            r["primer_pair"] = "revP2"

    result = gate_coamplify(pool, fx.recipe())

    assert result.status == FAIL
    assert "syn_recoded" in result.detail
    assert "shares no reverse priming site" in result.detail
    assert "missense_recoded" in result.detail


def test_a_parentless_control_cut_off_from_the_designed_arms_is_caught():
    """A control with no parent must still co-amplify with the designs it is compared to."""
    pool = fx.built_pool()
    for r in pool:
        if r["subpool"] == "ptc":
            r["primer_pair"] = "revP2"

    result = gate_coamplify(pool, fx.recipe())

    assert result.status == FAIL
    assert "ptc" in result.detail
    assert "shares no reverse priming site with the designed arms" in result.detail
    assert "different reaction" in result.detail


def test_an_exemption_must_be_declared_and_then_is_honoured():
    recipe = load_recipe_string(
        fx.RECIPE_TOML.replace(
            '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10',
            '[arm.ptc]\nkind   = "control_ptc"\ntarget = 10\ncoamplify = false',
        )
    )
    pool = fx.built_pool()
    for r in pool:
        if r["subpool"] == "ptc":
            r["primer_pair"] = "revP2"

    assert gate_coamplify(pool, recipe).status == PASS


def test_a_control_sharing_one_of_several_primers_with_its_parent_passes():
    """A two-primer pool is fine as long as every control reaches what it controls."""
    pool = fx.built_pool()
    for i, r in enumerate(pool):
        if r["subpool"] == "missense_recoded":
            r["primer_pair"] = "revP1" if i % 2 else "revP2"
        else:
            r["primer_pair"] = "revP1"
    # helpers all on revP2, which their parent also uses
    for r in pool:
        if r["subpool"] == "syn_recoded":
            r["primer_pair"] = "revP2"

    result = gate_coamplify(pool, fx.recipe())

    assert result.status == PASS, result.detail
    # the parentless controls reach the designed arms on revP1; the helper reaches its
    # parent on revP2
    by_arm = {row[0]: row for row in result.rows}
    assert by_arm["syn_recoded"][1] == "missense_recoded"
    assert by_arm["ptc"][1] == "<designed arms>"


def test_no_primer_column_is_not_applicable_rather_than_a_silent_pass():
    pool = fx.built_pool()
    for r in pool:
        del r["primer_pair"]

    result = gate_coamplify(pool, fx.recipe())
    assert result.status == NOT_APPLICABLE
    assert "primer_column" in result.detail


def test_a_helper_arm_may_declare_SEVERAL_parent_arms():
    """🔴 One helper arm can control more than one category of parent.

    v7 adds three recoded-alternative arms beside the original missense one, and the silent
    change for every one of them lives in the SAME helper arm. With a single parent name the
    gate reported every helper of the other categories as an orphan: 624 of 1,033 rows on the
    real build, while all four parent arms were present in the pool.
    """
    recipe = load_recipe_string(
        fx.RECIPE_TOML.replace(
            'target = 5\nparent = "missense_recoded"',
            'target = 5\nparent = ["missense_recoded", "ptc"]\nparent_key = "gene"',
        )
    )
    pool = fx.built_pool()
    # half the helpers point at a missense parent, half at a ptc parent
    mis = sorted({r["gene"] for r in pool if r["subpool"] == "missense_recoded"})
    ptc = sorted({r["gene"] for r in pool if r["subpool"] == "ptc"})
    helpers = [r for r in pool if r["subpool"] == "syn_recoded"]
    for i, r in enumerate(helpers):
        r["gene"] = (mis[i % len(mis)] if i % 2 == 0 else ptc[i % len(ptc)])

    result = gate_parent(pool, recipe)
    assert result.status == PASS, result.detail


def test_several_parent_arms_still_catch_a_genuine_orphan():
    """Widening the parent set must not weaken the check it exists for."""
    recipe = load_recipe_string(
        fx.RECIPE_TOML.replace(
            'target = 5\nparent = "missense_recoded"',
            'target = 5\nparent = ["missense_recoded", "ptc"]\nparent_key = "gene"',
        )
    )
    pool = fx.built_pool()
    mis = sorted({r["gene"] for r in pool if r["subpool"] == "missense_recoded"})
    helpers = [r for r in pool if r["subpool"] == "syn_recoded"]
    for i, r in enumerate(helpers):
        r["gene"] = mis[i % len(mis)]
    helpers[0]["gene"] = "YGHOST999W"          # in NEITHER parent arm

    result = gate_parent(pool, recipe)
    assert result.status == FAIL
    assert "1 of 5 rows" in result.detail
