"""find_re_site: the GDB-01 root fix (planning/pipeline_provenance_rebuild 30_; policy 2026-09-10).

The defect: gdb step 02 searched the SpCas9 cloning site first with an indel-tolerant fuzzy fallback, so on a
Cas12a read whose guide ends in G the SpCas9 pattern tied (90.9) one base early and the guide came out as
'T' + guide[:19]. These tests fail against the old find_anchor()-based search and pass against find_re_site().
"""
import random

import pytest

from magestic.utils.sequence_anchors import (
    RE_SITE_PATTERNS, ReSiteHit, find_anchor, find_re_site, hamming_mismatches)

SP_SITE = dict(RE_SITE_PATTERNS)["SpCas9"]
LB_SITE = dict(RE_SITE_PATTERNS)["LbCas12a"]
rng = random.Random(20260911)


def _seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _read(guide, site, leader="ACGTAT", donor=None):
    # the leader ends in T, as the real MAGESTIC 1.0/2.5 leader does; that T is the base the old
    # search prepended to every G-ending Cas12a guide
    donor = donor or _seq(129)
    return leader + guide + site + donor + _seq(20) + _seq(10) + "AGGAAAACAGAC"


def _guide_from(read, hit):
    return read[hit.pos - 20:hit.pos]


@pytest.mark.parametrize("last", list("ACGT"))
def test_cas12a_guide_any_last_base_is_exact(last):
    guide = _seq(19) + last
    read = _read(guide, LB_SITE)
    hit = find_re_site(read, end=150)
    assert hit == ReSiteHit(26, "LbCas12a", 0, ("LbCas12a",))
    assert _guide_from(read, hit) == guide


def test_old_search_reproduces_the_defect_on_g_ending_cas12a_guides():
    """Documents the mechanism: SpCas9-first fuzzy search lands one base early ONLY for G-ending guides."""
    for last, shifted in (("G", True), ("A", False), ("C", False), ("T", False)):
        guide = _seq(19) + last
        read = _read(guide, LB_SITE)
        old_pos = find_anchor(read, SP_SITE, end=150)          # what step 02 used to do first
        assert old_pos is not None                              # the fuzzy fallback "succeeds" on a Cas12a read
        assert (read[old_pos - 20:old_pos] == "T" + guide[:19]) is shifted
        assert (read[old_pos - 20:old_pos] == guide) is not shifted


def test_g_ending_cas12a_guide_with_a_substitution_in_the_site_still_lands_on_the_true_site():
    """The tie case a 'prefer the higher fuzzy score' fix would still get wrong: a 1-nt error in the last site
    base leaves the SpCas9 ratio at the shifted window equal to the Cas12a ratio at the true window.
    Substitution-only scoring separates them (1 vs 6 mismatches)."""
    guide = _seq(19) + "G"
    read = _read(guide, LB_SITE[:-1] + "A")
    hit = find_re_site(read, end=150)
    assert hit.pos == 26 and hit.kind == "LbCas12a" and hit.mismatches == 1
    assert _guide_from(read, hit) == guide


def test_spcas9_reads_unchanged():
    guide = _seq(20)
    read = _read(guide, SP_SITE)
    hit = find_re_site(read, end=150)
    assert hit == ReSiteHit(26, "SpCas9", 0, ("SpCas9",))
    assert _guide_from(read, hit) == guide
    read1 = _read(guide, "GTTTGAAGAGA")
    hit1 = find_re_site(read1, end=150)
    assert (hit1.pos, hit1.kind, hit1.mismatches) == (26, "SpCas9", 1)


def test_exact_beats_approximate_regardless_of_pattern_order():
    guide = _seq(19) + "G"
    read = _read(guide, LB_SITE)
    reversed_patterns = tuple(reversed(RE_SITE_PATTERNS))
    assert find_re_site(read, end=150, patterns=reversed_patterns).pos == 26
    assert find_re_site(read, end=150).pos == 26


def test_ambiguous_kind_when_both_patterns_tie_at_the_same_position():
    read = _read(_seq(20), "GTTCGAAGAGC")   # 1 mismatch from BOTH sites, same window
    hit = find_re_site(read, end=150)
    assert hit.pos == 26 and hit.kind == "ambiguous" and hit.mismatches == 1
    assert set(hit.kinds) == {"SpCas9", "LbCas12a"}


def test_two_substitutions_in_the_site_is_not_a_hit():
    guide = _seq(20)
    read = _read(guide, "GTTTGAAGTTC")      # 2 substitutions from SpCas9, 4 from LbCas12a
    assert find_re_site(read, end=150) is None


def test_no_site_at_all_and_degenerate_inputs():
    assert find_re_site("ACGT" * 40, end=150) is None
    assert find_re_site("", end=150) is None
    assert find_re_site(_seq(5)) is None


def test_an_insertion_in_the_site_is_scored_as_substitutions_not_as_a_shift():
    """Indels are deliberately NOT modelled: the site is a fixed-length, fixed-position anchor, and
    indel tolerance is exactly what let a 1-base shift tie with a 1-base substitution (GDB-01).

    A 1-nt insertion therefore never reads as an exact site. Whether the read is dropped or
    anchored one base over depends on the inserted base: here 'GTTTTGAAGAG' leaves the next window
    ('TTTTGAAGAGC') one substitution from the LbCas12a site, so the read is anchored there and the
    guide loses its last base. That is rare (it needs a synthesis/sequencing insertion inside the
    11-nt site) and it is labelled: mismatches=1 in the hit, so a caller can require mismatches==0.
    The pre-fix behaviour was the opposite of rare -- every G-ending Cas12a guide, always."""
    guide = _seq(20)
    read = _read(guide, "GTTTTGAAGAG")
    hit = find_re_site(read, end=150)
    assert hit is not None and hit.mismatches == 1
    assert hit.pos == 27 and read[hit.pos - 20:hit.pos] != guide
    # with approximate matching off, an insertion is simply not a site
    assert find_re_site(read, end=150, max_mismatches=0) is None


def test_earliest_exact_position_wins_across_patterns():
    read = "AC" + LB_SITE + _seq(10) + SP_SITE + _seq(30)
    assert find_re_site(read).pos == 2
    read = "AC" + SP_SITE + _seq(10) + LB_SITE + _seq(30)
    assert find_re_site(read).pos == 2


def test_search_window_is_respected():
    guide = _seq(20)
    read = _seq(200) + guide + LB_SITE + _seq(40)
    assert find_re_site(read, end=150) is None
    assert find_re_site(read).pos == 220


def test_hamming_mismatches():
    assert hamming_mismatches("ACGT", "ACGT") == 0
    assert hamming_mismatches("ACGT", "ACGA") == 1
    assert hamming_mismatches("ACGT", "TGCA", stop_after=1) == 2
    with pytest.raises(ValueError):
        hamming_mismatches("A", "AC")


def test_pattern_validation():
    with pytest.raises(ValueError):
        find_re_site("ACGT", patterns=())
    with pytest.raises(ValueError):
        find_re_site("ACGT", patterns=(("a", "ACG"), ("b", "ACGT")))
