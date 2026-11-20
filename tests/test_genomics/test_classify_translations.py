"""classify_translations, and the 7 real designs that proved it was needed.

A project copy of this logic shipped a 24,000-oligo pool in which 7 designs that destroy a gene's
start methionine were labelled `missense`: its classifier took the first differing residue and had
no start-codon case. Measured against this module on the same 14,748 coding rows the two agreed on
14,741, and every one of the 7 disagreements was codon 1.

The real oligo identifiers are kept in the parametrisation so the regression is traceable to the
designs that found it.
"""

from __future__ import annotations

import pytest

from magestic.genomics.consequence import classify_translations


def test_a_synonymous_change_names_no_residue():
    v = classify_translations("MKV*", "MKV*", "A", "G")
    assert v.consequence == "synonymous" and v.aa_change == "" and v.aa_changes == ()


def test_an_ordinary_substitution_is_missense():
    v = classify_translations("MKV*", "MKA*", "T", "C")
    assert v.consequence == "missense" and v.aa_change == "V3A"


@pytest.mark.parametrize(
    "oligo_id,prot_ref,prot_alt,expected_aa",
    [
        ("natural_YER017C_chrV_191381_T_C", "MKV*", "VKV*", "M1V"),
        ("natural_YER050C_chrV_253978_C_A", "MKV*", "IKV*", "M1I"),
        ("natural_YNL184C_chrXIV_293605_C_T", "MKV*", "IKV*", "M1I"),
        ("recoded_YER017C_chrV_191379_CAT_TAC", "MKV*", "VKV*", "M1V"),
    ],
)
def test_destroying_the_start_codon_is_start_lost_not_missense(
    oligo_id, prot_ref, prot_alt, expected_aa
):
    """🔴 The shipped defect. Each of these read `missense` in the project's own classifier."""
    v = classify_translations(prot_ref, prot_alt, "A", "G")
    assert v.consequence == "start_lost", oligo_id
    assert v.aa_change == expected_aa


def test_a_gained_stop_outranks_a_destroyed_start():
    """Documented precedence, pinned: both happen in one MNP and nonsense wins."""
    v = classify_translations("MKV*", "V*V*", "AAA", "GTG")
    assert v.consequence == "nonsense" and v.aa_change.endswith("*")


def test_a_destroyed_start_outranks_a_destroyed_stop():
    v = classify_translations("MKV*", "VKVQ", "AAA", "GTG")
    assert v.consequence == "start_lost" and v.aa_change == "M1V"


def test_a_destroyed_stop_alone_is_stop_lost():
    v = classify_translations("MKV*", "MKVQ", "T", "C")
    assert v.consequence == "stop_lost" and v.aa_change.startswith("*")


def test_a_multi_residue_change_reports_every_one():
    v = classify_translations("MKVA*", "MRLA*", "AAA", "GGC")
    assert v.consequence == "missense"
    assert v.aa_changes == ("K2R", "V3L")


@pytest.mark.parametrize(
    "prot_ref,prot_alt,ref,alt,expected",
    [
        ("MKVA*", "MKA*", "AAA", "", "inframe_deletion"),
        ("MKVA*", "MKVGA*", "", "GGG", "inframe_insertion"),
        ("MKVA*", "MKVA*", "A", "", "frameshift"),
        ("MKVA*", "MKVA*", "", "GG", "frameshift"),
    ],
)
def test_length_changing_alleles(prot_ref, prot_alt, ref, alt, expected):
    """🔴 THE FIXTURE HAD TO BE FIXED, TWICE OVER THE SAME LESSON. The first version passed
    IDENTICAL translations for an in-frame insertion, which cannot happen: an insertion lengthens
    the protein. With both proteins equal the premature-stop check fires correctly and the answer is
    `nonsense`, so the test failed against correct code. A fixture that could not occur in the data
    tests nothing, and the earlier G-GC fixture in this suite failed the same way."""
    assert classify_translations(prot_ref, prot_alt, ref, alt).consequence == expected


def test_an_inframe_deletion_that_truncates_is_nonsense_not_a_deletion():
    """The premature-stop check has to allow for the length the indel itself removes."""
    v = classify_translations("MKVAQRST*", "MK*", "AAAAAA", "")
    assert v.consequence == "nonsense"


def test_codon_is_optional_and_changes_no_consequence():
    with_c = classify_translations("MKVA*", "MKA*", "AAA", "", codon=2)
    without = classify_translations("MKVA*", "MKA*", "AAA", "")
    assert with_c.consequence == without.consequence == "inframe_deletion"
    # only the RESIDUE NUMBER is codon-dependent; "1 residue(s)" is in both, so assert on the
    # phrase that actually carries the number rather than on the word
    assert " at residue 3" in with_c.detail
    assert " at residue " not in without.detail


def test_it_agrees_with_annotate_variant_on_the_same_translations():
    """One implementation means the two entry points cannot drift."""
    from magestic.genomics.consequence import annotate_variant, build_cds, translate

    seq = "ATGAAAGTTTAA" + "G" * 30
    model = {"G1": ("chrT", "+", [(0, 12)])}
    genome = {"chrT": seq}
    for pos0, ref, alt in ((0, "A", "G"), (4, "A", "G"), (3, "AAA", "TAA"), (10, "A", "G")):
        av = annotate_variant(genome, model, "G1", "chrT", pos0, ref, alt)
        edited = seq[:pos0] + alt + seq[pos0 + len(ref):]
        ct = classify_translations(
            translate(build_cds(seq, model["G1"][2], "+")),
            translate(build_cds(edited, model["G1"][2], "+")),
            ref, alt,
        )
        assert av.consequence == ct.consequence, (pos0, ref, alt)
        assert av.aa_change == ct.aa_change, (pos0, ref, alt)


# ------------- the unshifted-window bug, which this module ALSO had until 2026-11-20

def test_shift_exons_for_indel_moves_only_what_the_edit_moves():
    from magestic.genomics.consequence import shift_exons_for_indel as sh
    assert sh([(0, 30)], 10, +3) == [(0, 33)]
    assert sh([(0, 30), (50, 80)], 10, +3) == [(0, 33), (53, 83)]
    assert sh([(0, 30), (50, 80)], 60, +3) == [(0, 30), (50, 83)]
    assert sh([(0, 30)], 10, -3) == [(0, 27)]


def test_shift_declines_rather_than_guesses_at_a_boundary():
    from magestic.genomics.consequence import shift_exons_for_indel as sh
    assert sh([(0, 30), (50, 80)], 40, +3) is None      # between the exons
    assert sh([(0, 30)], 10, -40) is None               # removes more than the exon holds


def test_annotate_variant_now_lengthens_the_protein_for_an_insertion():
    """🔴 The bug: splicing the edited genome with the ORIGINAL exons kept the CDS the same width,
    so an insertion displaced downstream bases instead of extending the protein."""
    from magestic.genomics.consequence import annotate_variant, build_cds, translate
    seq = "ATG" + "AAA" * 10 + "TAA" + "G" * 30
    model = {"G1": ("chrT", "+", [(0, 36)])}
    v = annotate_variant({"chrT": seq}, model, "G1", "chrT", 6, "", "GGG")
    assert v.consequence == "inframe_insertion"
    # and the edited protein really is one residue longer
    from magestic.genomics.consequence import shift_exons_for_indel
    ex2 = shift_exons_for_indel(model["G1"][2], 6, +3)
    edited = seq[:6] + "GGG" + seq[6:]
    assert len(translate(build_cds(edited, ex2, "+"))) == len(
        translate(build_cds(seq, model["G1"][2], "+"))) + 1


def test_classify_translations_REFUSES_an_unshifted_window():
    """The caller did the splicing, so the caller can have done it wrong. Equal-length
    translations with a non-zero in-frame net is the signature, and it is refused, not guessed."""
    v = classify_translations("MKVA*", "MKVA*", "", "GGG")
    assert v.consequence == "non_coding"
    assert "UNSHIFTED" in v.detail and "annotate_variant" in v.detail


def test_a_frameshift_is_still_classified_without_the_guard_firing():
    """net %% 3 != 0 cannot be checked this way and must still answer."""
    assert classify_translations("MKVA*", "MKVA*", "A", "").consequence == "frameshift"
