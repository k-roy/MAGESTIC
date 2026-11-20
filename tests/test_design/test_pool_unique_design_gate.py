"""G-UNIQUE, tested against the incident that produced it.

NN 600. A collaborator at EMBL found 20 variants present twice in a finished 24,000-oligo pool, 40
oligos. Each pair carried the SAME guide and the SAME edit over byte-identical donor content, and
differed only in where the 129-nt window sat, by 1 to 3 nt.

Two near misses are encoded as tests here, because both are the kind of check that feels sufficient:

  * a SEQUENCE-identity gate passes, legitimately. Shift the window and the 200-mers really are
    different molecules. Design identity is a different question.
  * the key `(gene, guide_seq, ref_allele, alt_allele)` sweeps in every legitimate second and third
    shot at the same edit, which on the real pool meant 308 groups and 313 excess rows hiding 20.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import gate_unique_design, load_recipe_string

RECIPE = """
[pool]
name = "t"
total = 4
oligo_nt = 8
id_column = "oligo_id"
arm_column = "subpool"
sequence_column = "oligo_seq"

[invariants]
design_identity_key = ["chrom", "edit_pos", "ref_allele", "alt_allele", "guide_seq"]
design_identity_exempt_arms = ["scrambled"]

[arm.a]
kind = "natural_variant"
target = 4
"""

NO_KEY = RECIPE.replace(
    'design_identity_key = ["chrom", "edit_pos", "ref_allele", "alt_allele", "guide_seq"]\n', ""
).replace('design_identity_exempt_arms = ["scrambled"]\n', "")


def row(oid, chrom="chrI", pos="100", ref="A", alt="T", guide="G" * 20,
        arm="a", seq="ACGTACGT"):
    return {"oligo_id": oid, "chrom": chrom, "edit_pos": pos, "ref_allele": ref,
            "alt_allele": alt, "guide_seq": guide, "subpool": arm, "oligo_seq": seq}


def test_a_unique_pool_passes():
    r = load_recipe_string(RECIPE)
    res = gate_unique_design([row("o1"), row("o2", pos="200")], r)
    assert res.status == "pass", res.detail


def test_the_real_defect_is_caught_even_though_the_sequences_differ():
    """🔴 The incident: same edit, same guide, DIFFERENT 200-mer because the window shifted."""
    r = load_recipe_string(RECIPE)
    a = row("helper_x", seq="ACGTACGT")
    b = row("syn_x", seq="CGTACGTA")        # a genuinely different molecule
    assert a["oligo_seq"] != b["oligo_seq"], "the fixture must not be two identical sequences"
    res = gate_unique_design([a, b], r)
    assert res.status == "fail"
    assert "same guide" in res.detail or "same experiment twice" in res.detail


def test_an_extra_shot_at_the_same_edit_is_NOT_a_duplicate():
    """The 1,584 legitimate second shots differ in guide_seq and must fall out by construction."""
    r = load_recipe_string(RECIPE)
    res = gate_unique_design([row("first", guide="A" * 20), row("second", guide="C" * 20)], r)
    assert res.status == "pass", res.detail


def test_an_exempt_arm_is_skipped():
    r = load_recipe_string(RECIPE)
    res = gate_unique_design([row("s1", arm="scrambled"), row("s2", arm="scrambled")], r)
    assert res.status == "pass"
    assert "exempt" in res.detail


def test_an_undeclared_key_is_not_applicable_rather_than_invented():
    """🔴 The wrong key is how the defect survived, so this gate never guesses one."""
    r = load_recipe_string(NO_KEY)
    res = gate_unique_design([row("o1"), row("o1")], r)
    assert res.status == "not_applicable"


def test_a_key_naming_a_missing_column_fails_rather_than_passing_vacuously():
    r = load_recipe_string(RECIPE.replace('"guide_seq"]', '"no_such_column"]'))
    res = gate_unique_design([row("o1")], r)
    assert res.status == "fail"
    assert "cannot be read" in res.detail


def test_rows_with_an_incomplete_key_are_skipped_not_collided():
    """Two rows both missing edit_pos must not be read as duplicates of each other."""
    r = load_recipe_string(RECIPE)
    res = gate_unique_design([row("o1", pos=""), row("o2", pos="")], r)
    assert res.status == "pass"


def test_three_rows_sharing_one_key_are_reported_as_one_group():
    r = load_recipe_string(RECIPE)
    res = gate_unique_design([row("o1"), row("o2"), row("o3")], r)
    assert res.status == "fail"
    assert "1 group" in res.detail and "3 oligo" in res.detail


def test_exempt_arms_without_a_key_is_rejected_by_the_loader():
    from magestic.design.pool import RecipeError
    bad = NO_KEY.replace("[invariants]\n", '[invariants]\ndesign_identity_exempt_arms = ["x"]\n')
    with pytest.raises(RecipeError, match="applies to no check"):
        load_recipe_string(bad)
