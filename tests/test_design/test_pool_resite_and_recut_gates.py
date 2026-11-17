"""G-RESITE and G-RECUT, plus a list-valued parent key.

🔴 These tests are written against the REAL 2026-11-17 incident, not against the implementation.
A 24,000-oligo order was assembled carrying 63 donors with an internal SapI site (the cloning
enzyme cuts the insert) and 127 whose repair product is still a perfect substrate for its own guide
(the edit is installed and then cut back out). The design code had guards for both
(``remove_bspqi_sites_from_donor``, ``guide_disruption``) and both are correct — but they run
per-oligo at DESIGN time, and the codon unit rule rebuilds donors AFTERWARDS. Nothing re-asked the
question of the assembled pool. The pool shipped before that rebuild had zero of either.

Each test below fails on the pre-fix code and passes on the fixed code.
"""

from __future__ import annotations

import pytest

from magestic.design.pool.gates import (
    FAIL,
    NOT_APPLICABLE,
    PASS,
    gate_parent,
    gate_recut,
    gate_restriction_sites,
    run_gates,
)
from magestic.design.pool.recipe import load_recipe_string

SAPI_FWD = "GAAGAGC"
# the architecture the TRIUMPH pool declares: donor occupies [51, 180) of a 200 nt oligo
RECIPE = """
[pool]
name = "t"
total = 2
oligo_nt = 200
arm_column = "subpool"
[segments]
forward_priming = [0, 20]
guide = [20, 40]
cloning_site = [40, 51]
donor = [51, 180]
reverse_priming = [180, 200]
[arm.natural]
kind = "natural_variant"
target = 2
[cloning]
restriction_sites = ["GAAGAGC", "GCTCTTC"]
restriction_sites_per_oligo = 1
recut_pam = "NGG"
"""

NO_CLONING = RECIPE.split("[cloning]")[0]


def _oligo(donor: str, guide: str = "A" * 20) -> str:
    """Assemble a 200 nt oligo with ONE SapI site in the cloning linker, as the real pool has."""
    fwd = "C" * 20
    linker = "GTTTGAAGAGC"            # 11 nt, carries the single intended site
    assert len(linker) == 11
    donor = (donor + "T" * 129)[:129]
    rev = "G" * 20
    oligo = fwd + guide + linker + donor + rev
    assert len(oligo) == 200, len(oligo)
    return oligo


def _row(oligo: str, guide: str = "A" * 20, **kw):
    d = {"oligo_id": kw.pop("oligo_id", "o1"), "subpool": "natural",
         "oligo_seq": oligo, "guide_seq": guide}
    d.update(kw)
    return d


# --------------------------------------------------------------------------- G-RESITE

def test_resite_passes_a_clean_pool():
    recipe = load_recipe_string(RECIPE)
    pool = [_row(_oligo("ACGT"), oligo_id="a"), _row(_oligo("TTGC"), oligo_id="b")]
    r = gate_restriction_sites(pool, recipe)
    assert r.status == PASS, r.detail


def test_resite_FAILS_on_a_site_inside_the_donor():
    """🔴 The 63-oligo incident. The extra site sits in the donor BODY, not at a junction, so the
    design-time junction fixer could never have removed it."""
    recipe = load_recipe_string(RECIPE)
    donor = "ACGT" + SAPI_FWD + "ACGT"
    pool = [_row(_oligo(donor), oligo_id="bad1")]
    r = gate_restriction_sites(pool, recipe)
    assert r.status == FAIL
    assert "INSIDE the donor" in r.detail
    assert r.rows and r.rows[0][0] == "bad1"


def test_resite_FAILS_when_the_per_oligo_count_is_wrong():
    recipe = load_recipe_string(RECIPE)
    # two sites: the linker's plus one in the donor
    pool = [_row(_oligo("A" + SAPI_FWD), oligo_id="two")]
    r = gate_restriction_sites(pool, recipe)
    assert r.status == FAIL
    assert "expected 1" in str(r.rows) or "INSIDE the donor" in r.detail


def test_resite_is_not_applicable_when_undeclared():
    recipe = load_recipe_string(NO_CLONING)
    r = gate_restriction_sites([_row(_oligo("ACGT"))], recipe)
    assert r.status == NOT_APPLICABLE
    assert "declared" in r.detail


# --------------------------------------------------------------------------- G-RECUT

def test_recut_FAILS_when_the_donor_still_matches_guide_plus_pam():
    """🔴 The 127-oligo incident: the repair product is a perfect substrate, so the edit is cut
    back out. 90 of those rows declared recut_hazard = 0."""
    recipe = load_recipe_string(RECIPE)
    guide = "ACGTACGTACGTACGTACGT"
    pool = [_row(_oligo(guide + "TGG" + "AAAA"), guide=guide, oligo_id="recut1")]
    r = gate_recut(pool, recipe)
    assert r.status == FAIL
    assert "re-cut" in r.detail
    assert r.rows and r.rows[0][0] == "recut1"


def test_recut_catches_the_reverse_complement_strand():
    recipe = load_recipe_string(RECIPE)
    guide = "ACGTACGTACGTACGTACGT"
    site = guide + "TGG"
    rc = site.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
    pool = [_row(_oligo("TT" + rc), guide=guide, oligo_id="rc1")]
    r = gate_recut(pool, recipe)
    assert r.status == FAIL, r.detail


def test_recut_passes_a_disrupted_donor():
    recipe = load_recipe_string(RECIPE)
    guide = "ACGTACGTACGTACGTACGT"
    # the PAM is destroyed: TGG -> TAA
    pool = [_row(_oligo(guide + "TAA" + "AAAA"), guide=guide, oligo_id="ok1")]
    r = gate_recut(pool, recipe)
    assert r.status == PASS, r.detail


def test_recut_is_not_applicable_without_a_declared_pam():
    recipe = load_recipe_string(NO_CLONING)
    r = gate_recut([_row(_oligo("ACGT"))], recipe)
    assert r.status == NOT_APPLICABLE


# --------------------------------------------------------------------------- G-PARENT, list form

PARENTED = """
[pool]
name = "t"
total = 4
arm_column = "subpool"
[arm.parentarm]
kind = "recoded_synonymous"
target = 2
[arm.helper]
kind = "control_recoding_helper"
target = 2
parent = ["parentarm"]
parent_key = "recode_pair_key"
"""


def test_parent_accepts_a_list_valued_key_when_every_element_resolves():
    """🔴 The regression that blocked collapsing duplicate molecules: one helper legitimately
    serves TWO parents, so its key carries both. Resolving the value atomically read it as an
    orphan and failed a correct pool."""
    recipe = load_recipe_string(PARENTED)
    pool = [
        {"oligo_id": "p1", "subpool": "parentarm", "recode_pair_key": "K1"},
        {"oligo_id": "p2", "subpool": "parentarm", "recode_pair_key": "K2"},
        {"oligo_id": "h1", "subpool": "helper", "recode_pair_key": "K1;K2"},
        {"oligo_id": "h2", "subpool": "helper", "recode_pair_key": "K2"},
    ]
    r = gate_parent(pool, recipe)
    assert r.status == PASS, r.detail


def test_parent_still_fails_when_one_element_of_the_list_is_an_orphan():
    """The list form must be STRICTLY STRONGER, not a way through."""
    recipe = load_recipe_string(PARENTED)
    pool = [
        {"oligo_id": "p1", "subpool": "parentarm", "recode_pair_key": "K1"},
        {"oligo_id": "p2", "subpool": "parentarm", "recode_pair_key": "K2"},
        {"oligo_id": "h1", "subpool": "helper", "recode_pair_key": "K1;NOPE"},
        {"oligo_id": "h2", "subpool": "helper", "recode_pair_key": "K2"},
    ]
    r = gate_parent(pool, recipe)
    assert r.status == FAIL
    # G-PARENT reports arm-level rows, so the orphan COUNT is the evidence, not the id
    assert "1 of 2" in r.detail and "orphan" in str(r.rows)


# --------------------------------------------------------------------------- wiring + invariants

def test_run_gates_includes_the_two_new_gates():
    recipe = load_recipe_string(RECIPE)
    names = [g.gate for g in run_gates([_row(_oligo("ACGT"))], recipe)]
    assert "G-RESITE" in names and "G-RECUT" in names


def test_invariants_block_is_parsed_not_dropped():
    """It used to be read by NOTHING while its own comment claimed a dedupe would "fail loudly"."""
    recipe = load_recipe_string(RECIPE + "\n[invariants]\nduplicate_biology_rows = 311\n")
    assert dict(recipe.invariants)["duplicate_biology_rows"] == 311


def test_cloning_block_rejects_a_count_without_sites():
    with pytest.raises(Exception):
        load_recipe_string(
            NO_CLONING + "\n[cloning]\nrestriction_sites_per_oligo = 1\n"
        )
