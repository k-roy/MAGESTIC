"""G-GC, G-IDLEN and G-NONTARGET, each tested against the defect that produced it.

NN 580 §8 measured three residual classes on a 24,000-oligo pool whose other eight gates all
passed. Each one is encoded here with the NEAR MISS that let it ship, because in all three cases a
weaker check existed and returned green:

  * **G-GC** — 200 rows carried a ``gc_frac`` true of the molecule BEFORE its own edit, because the
    codon unit rule rebuilt the oligo and nothing recomputed the column. A schema gate sees the
    column present and populated; presence is not correctness.
  * **G-IDLEN** — 110 identifiers ran over 80 characters, to a maximum of 151, every one carrying a
    long indel allele inline. Identifier UNIQUENESS was gated and passed; length was never asked.
  * **G-NONTARGET** — five "non-targeting" controls had a genomic NGG site two mismatches from
    their guide. 🔴 An EXACT-match screen passed and proved nothing: the closest approach of those
    guides is two mismatches, never zero or one.
"""

from __future__ import annotations

import pytest

from magestic.design.pool import (
    ProtospacerIndex,
    gate_gc_fraction,
    gate_identifier_length,
    gate_nontargeting,
    load_recipe_string,
    make_nontargeting_controls,
    OligoArchitecture,
)
from magestic.design.pool.recipe import RecipeError

BASE = """
[pool]
name = "t"
total = 2
oligo_nt = 8
id_column = "oligo_id"
arm_column = "subpool"
sequence_column = "oligo_seq"

[invariants]
{inv}

[arm.a]
kind = "natural_variant"
target = 1

[arm.scrambled]
kind = "control_nontargeting"
target = 1
"""


def rec(inv):
    return load_recipe_string(BASE.format(inv=inv))


def row(oid="o1", seq="ACGTACGT", gc=0.5, arm="a", guide=None):
    r = {"oligo_id": oid, "subpool": arm, "oligo_seq": seq, "gc_frac": gc}
    if guide is not None:
        r["guide_seq"] = guide
    return r


# --------------------------------------------------------------------- G-GC

def test_gc_undeclared_is_not_applicable():
    assert gate_gc_fraction([row()], rec("")).status == "not_applicable"


def test_gc_agreeing_passes():
    g = gate_gc_fraction([row(seq="GGCCAATT", gc=0.5)], rec('gc_fraction_column = "gc_frac"'))
    assert g.status == "pass", g.detail


def test_gc_stale_value_from_a_rebuilt_oligo_fails():
    """The real shape: the sequence changed by one base and the column did not."""
    # 🔴 THE FIXTURE IS ITSELF A TRAP: "GGCCAATA" still has 4 G/C of 8, so it agrees with 0.5
    # and this test passed-as-pass on the first run. The sequence has to actually move the GC.
    pool = [row(oid="rebuilt", seq="GGCAAATT", gc=0.5)]      # true GC is 3/8 = 0.375
    g = gate_gc_fraction(pool, rec('gc_fraction_column = "gc_frac"'))
    assert g.status == "fail"
    assert "rebuilt" in str(g.rows)
    assert "0.37500" in str(g.rows)


def test_gc_naming_a_missing_column_fails_rather_than_passing_vacuously():
    g = gate_gc_fraction([{"oligo_id": "o", "subpool": "a", "oligo_seq": "ACGT"}],
                         rec('gc_fraction_column = "not_a_column"'))
    assert g.status == "fail"


def test_gc_blank_cell_is_skipped_not_failed():
    g = gate_gc_fraction([row(gc="")], rec('gc_fraction_column = "gc_frac"'))
    assert g.status == "pass"


# ------------------------------------------------------------------ G-IDLEN

def test_idlen_undeclared_is_not_applicable():
    assert gate_identifier_length([row()], rec("")).status == "not_applicable"


def test_idlen_under_the_cap_passes():
    assert gate_identifier_length([row(oid="x" * 80)],
                                  rec("max_oligo_id_chars = 80")).status == "pass"


def test_idlen_over_the_cap_fails_and_names_the_longest_first():
    pool = [row(oid="x" * 81), row(oid="y" * 151)]
    g = gate_identifier_length(pool, rec("max_oligo_id_chars = 80"))
    assert g.status == "fail"
    assert g.rows[0][1] == 151


def test_idlen_rejects_a_nonsense_cap_at_load():
    with pytest.raises(RecipeError):
        rec("max_oligo_id_chars = 0")


# -------------------------------------------------------------- G-NONTARGET

GEN = "ACGT" * 15 + "AAAACCCCGGGGTTTTACGT" + "AGG" + "TTAA" * 15
SCREEN = 'nontargeting_arms = ["scrambled"]\nnontargeting_max_genomic_mismatch = 2'


def idx(max_mm=2):
    return ProtospacerIndex.from_genome({"chrT": GEN}, 20, max_mm, "NGG")


def test_nontargeting_undeclared_is_not_applicable():
    assert gate_nontargeting([row()], rec(""), idx()).status == "not_applicable"


def test_nontargeting_declared_without_an_index_FAILS():
    """🔴 The whole point. A claim the build could not check is not a passing claim."""
    g = gate_nontargeting([row(arm="scrambled", guide="A" * 20)], rec(SCREEN), None)
    assert g.status == "fail"
    assert "never checked" in g.detail


def test_nontargeting_clean_guide_passes():
    g = gate_nontargeting([row(arm="scrambled", guide="TTTTTTTTTTTTTTTTTTTT")],
                          rec(SCREEN), idx())
    assert g.status == "pass", g.detail


def test_nontargeting_two_mismatch_site_fails():
    """The incident: 2 mismatches, PAM-distal, which an exact-match screen cannot see."""
    guide = "TTAACCCCGGGGTTTTACGT"          # 2 mismatches at offsets 0 and 1 of the real site
    g = gate_nontargeting([row(oid="scrambled_00004", arm="scrambled", guide=guide)],
                          rec(SCREEN), idx())
    assert g.status == "fail"
    assert "scrambled_00004" in str(g.rows)
    assert g.rows[0][3] == 2


def test_the_exact_match_near_miss_would_have_passed():
    """🔴 Proof the OLD check was insufficient, not merely unlucky."""
    guide = "TTAACCCCGGGGTTTTACGT"
    assert ProtospacerIndex.from_genome({"chrT": GEN}, 20, 0, "NGG").hits(guide) == []
    assert idx().hits(guide)          # the same guide, seen only once mismatches are allowed


def test_a_narrower_index_cannot_answer_the_declared_question():
    g = gate_nontargeting([row(arm="scrambled", guide="A" * 20)], rec(SCREEN), idx(max_mm=1))
    assert g.status == "fail"
    assert "narrower index" in g.detail


def test_only_the_declared_arms_are_screened():
    guide = "TTAACCCCGGGGTTTTACGT"
    g = gate_nontargeting([row(oid="designed", arm="a", guide=guide)], rec(SCREEN), idx())
    assert g.status == "pass"


# ------------------------------------- the screen runs at DRAW time, not only after the fact

def test_generator_rejects_a_candidate_the_screen_refuses():
    arch = OligoArchitecture(forward_priming="A" * 20, internal_cloning_site="GTTTGAAGAGC",
                             reverse_tail="C" * 20, total_length=200, guide_length=20)
    calls = {"n": 0}

    def reject_everything_once(guide, oligo):
        calls["n"] += 1
        return "planted" if calls["n"] <= 50 else None

    rows, rep = make_nontargeting_controls(5, [0.45] * 100, arch,
                                           reject=reject_everything_once)
    assert len(rows) == 5
    assert rep["screened"] is True
    assert rep["rejected_by_screen"]["planted"] == 50


def test_an_unscreened_generation_says_so_in_its_report():
    arch = OligoArchitecture(forward_priming="A" * 20, internal_cloning_site="GTTTGAAGAGC",
                             reverse_tail="C" * 20, total_length=200, guide_length=20)
    _, rep = make_nontargeting_controls(3, [0.45] * 100, arch)
    assert rep["screened"] is False


# ------------------------------------------------------------- recipe hygiene

def test_half_a_screen_is_refused_at_load():
    with pytest.raises(RecipeError):
        rec('nontargeting_arms = ["scrambled"]')
    with pytest.raises(RecipeError):
        rec("nontargeting_max_genomic_mismatch = 2")


# ------------------- G-IDLEN's scope: the name that ships, not the name in flight

def test_idlen_on_provisional_names_is_not_applicable_and_says_why():
    """🔴 Measured: run at build time this gate failed on 23,000 internal names, longest 209,
    none of which reaches the order file. A gate pointed at the wrong artifact is noise."""
    g = gate_identifier_length([row(oid="x" * 209)], rec("max_oligo_id_chars = 80"),
                               names_are_final=False)
    assert g.status == "not_applicable"
    assert "PROVISIONAL" in g.detail and "renamed pool that ships" in g.detail


def test_idlen_still_fails_when_the_names_are_final():
    g = gate_identifier_length([row(oid="x" * 209)], rec("max_oligo_id_chars = 80"),
                               names_are_final=True)
    assert g.status == "fail"


def test_run_gates_threads_the_scope_through():
    from magestic.design.pool import run_gates
    pool = [row(oid="x" * 209, arm="a", guide="T" * 20),
            row(oid="scrambled_00001", arm="scrambled", guide="T" * 20)]
    r = rec('gc_fraction_column = "gc_frac"\nmax_oligo_id_chars = 80')
    by = {g.gate: g for g in run_gates(pool, r, identifier_names_are_final=False)}
    assert by["G-IDLEN"].status == "not_applicable"
    by = {g.gate: g for g in run_gates(pool, r, identifier_names_are_final=True)}
    assert by["G-IDLEN"].status == "fail"


# ------------------------------------------------------------- G-RECODED

RF = 'recoded_flag_arms = ["a"]\nrecoded_flag_column = "recoded"'


def rrow(oid, arm, recoded):
    return {"oligo_id": oid, "subpool": arm, "oligo_seq": "ACGT", "recoded": recoded}


def test_recoded_flag_undeclared_is_not_applicable():
    from magestic.design.pool import gate_recoded_flag
    assert gate_recoded_flag([rrow("o", "a", "true")], rec("")).status == "not_applicable"


def test_recoded_flag_agreeing_both_ways_passes():
    from magestic.design.pool import gate_recoded_flag
    pool = [rrow("i", "a", "true"), rrow("o", "scrambled", "false")]
    assert gate_recoded_flag(pool, rec(RF)).status == "pass"


def test_recoded_flag_false_inside_a_recoded_arm_fails():
    """🔴 The shipped defect: an arm that is recoded by construction reading FALSE."""
    from magestic.design.pool import gate_recoded_flag
    g = gate_recoded_flag([rrow("synalt_1", "a", "false")], rec(RF))
    assert g.status == "fail" and "synalt_1" in str(g.rows)


def test_recoded_flag_true_outside_is_also_a_failure():
    """A one-way gate would pass a column that is true everywhere, which says nothing."""
    from magestic.design.pool import gate_recoded_flag
    g = gate_recoded_flag([rrow("i", "a", "true"), rrow("x", "scrambled", "true")], rec(RF))
    assert g.status == "fail" and "expected false" in str(g.rows)


def test_recoded_flag_blank_inside_is_a_failure_not_a_skip():
    from magestic.design.pool import gate_recoded_flag
    g = gate_recoded_flag([rrow("i", "a", "")], rec(RF))
    assert g.status == "fail" and "(blank)" in str(g.rows)


def test_recoded_flag_naming_an_absent_arm_fails():
    from magestic.design.pool import gate_recoded_flag
    g = gate_recoded_flag([rrow("o", "scrambled", "false")], rec(RF))
    assert g.status == "fail" and "absent from this pool" in g.detail


def test_half_a_recoded_flag_declaration_is_refused_at_load():
    with pytest.raises(RecipeError):
        rec('recoded_flag_arms = ["a"]')
    with pytest.raises(RecipeError):
        rec('recoded_flag_column = "recoded"')
