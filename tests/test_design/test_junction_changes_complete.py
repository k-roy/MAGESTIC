"""
NN 234 bug 4 -- ``junction_changes`` must account for every donor base it changed.

Both junction-repair implementations accept an edit that lowers the restriction-site
count without reaching 1, keep it, and move on to the other end::

    for end_idx in (0, -1):
        for new_base in "ATCG":
            if n_new < n_sites:
                n_sites = n_new
                donor = candidate              # edit accepted and kept
            if n_sites == 1:
                junction_changes.append(ONE record)   # only the LAST end recorded
                return donor, junction_changes

So when an index-0 edit takes the count 3 -> 2 and an index -1 edit then takes it
2 -> 1, the returned donor differs from the input at TWO positions while
``junction_changes`` names ONE.  That list is written verbatim into the design info
table (``oligo_generation.py``, ``str(junction_changes)``), so the design record
understates what was changed at the bench.

The contract these tests pin down: the positions named in ``junction_changes`` are
exactly the positions where the returned donor differs from the input.
"""

import pytest

from magestic.design.guide_donor_functions import (
    remove_BspQI_sites_from_donor_in_final_oligo,
)
from magestic.design.oligo_generation import remove_bspqi_sites_from_donor

SITE = "GAAGAGC"      # BspQI / SapI
RC_SITE = "GCTCTTC"

FWD_PRIMER = "ACACTCTTTCCCTACACG"
GUIDE = "ACGTACGTACGTACGTACGT"

#: Holds the ONE intended site (RC_SITE) and ends in GAAGAG, so donor[0] == "C"
#: completes a second site across the cloning-site/donor junction.
INTERNAL_CLONING_SITE = RC_SITE + "ATC" + "GAAGAG"

#: Starts with AAGAGC, so donor[-1] == "G" completes a third site across the
#: donor/rev-primer junction.
REV_PRIMER = "AAGAGC" + "TTACCTACTGCA"

#: donor[0] == "C" and donor[-1] == "G"; the interior carries no site of its own.
DONOR = "C" + "ATCATCATCATCATCATCATCATCATCA" + "G"


def _n_sites(seq: str) -> int:
    upper = seq.upper()
    return upper.count(SITE) + upper.count(RC_SITE)


def _changed_positions(before: str, after: str) -> list:
    """Indices (negative for the last base) where the donor was changed."""
    assert len(before) == len(after), "junction repair must not change donor length"
    out = []
    for i, (b, a) in enumerate(zip(before, after)):
        if b != a:
            out.append(i - len(before) if i == len(before) - 1 else i)
    return out


def _positions_named(junction_changes) -> set:
    """The donor indices named by the change records."""
    named = set()
    for rec in junction_changes or []:
        for item in rec:
            if isinstance(item, int):
                named.add(item)
                break
    return named


# ---------------------------------------------------------------------------
# The fixture itself must set up a two-ended repair, or the tests prove nothing
# ---------------------------------------------------------------------------

def test_fixture_has_three_sites_needing_repair_at_both_ends():
    oligo = FWD_PRIMER + GUIDE + INTERNAL_CLONING_SITE + DONOR + REV_PRIMER
    assert _n_sites(oligo) == 3, (
        f"fixture must start with 3 sites (1 intended + 2 junction); got "
        f"{_n_sites(oligo)}"
    )
    # Each junction site really is fixable by a single terminal-base change.
    left_fixed = FWD_PRIMER + GUIDE + INTERNAL_CLONING_SITE + ("A" + DONOR[1:]) + REV_PRIMER
    right_fixed = FWD_PRIMER + GUIDE + INTERNAL_CLONING_SITE + (DONOR[:-1] + "T") + REV_PRIMER
    assert _n_sites(left_fixed) == 2, "changing donor[0] must remove exactly one site"
    assert _n_sites(right_fixed) == 2, "changing donor[-1] must remove exactly one site"
    assert _n_sites(DONOR) == 0, "the donor interior must carry no site of its own"


# ---------------------------------------------------------------------------
# oligo_generation.remove_bspqi_sites_from_donor
# ---------------------------------------------------------------------------

def test_oligo_generation_records_every_changed_position():
    new_donor, junction_changes = remove_bspqi_sites_from_donor(
        FWD_PRIMER, GUIDE, INTERNAL_CLONING_SITE, DONOR, REV_PRIMER,
    )

    assert new_donor is not None, "the fixture is repairable; repair must succeed"
    oligo = FWD_PRIMER + GUIDE + INTERNAL_CLONING_SITE + new_donor + REV_PRIMER
    assert _n_sites(oligo) == 1, f"repair must leave exactly 1 site; got {_n_sites(oligo)}"

    changed = _changed_positions(DONOR, new_donor)
    named = _positions_named(junction_changes)

    assert len(junction_changes) == len(changed), (
        f"the donor was changed at {len(changed)} position(s) {changed} but "
        f"junction_changes records {len(junction_changes)}: {junction_changes}. "
        f"The design info table under-reports the donor edits."
    )
    assert named == set(changed), (
        f"junction_changes names positions {sorted(named)} but the donor actually "
        f"changed at {sorted(changed)}"
    )


def test_oligo_generation_records_the_true_original_base():
    """Each record's 'from' base must be the base that was really there."""
    new_donor, junction_changes = remove_bspqi_sites_from_donor(
        FWD_PRIMER, GUIDE, INTERNAL_CLONING_SITE, DONOR, REV_PRIMER,
    )
    assert new_donor is not None
    for rec in junction_changes:
        rec = list(rec)
        idx = rec[rec.index("from") - 1]
        from_base = rec[rec.index("from") + 1]
        to_base = rec[rec.index("to") + 1]
        assert DONOR[idx] == from_base, (
            f"record says index {idx} went from {from_base!r}, but the input donor "
            f"has {DONOR[idx]!r} there"
        )
        assert new_donor[idx] == to_base, (
            f"record says index {idx} went to {to_base!r}, but the returned donor "
            f"has {new_donor[idx]!r} there"
        )


# ---------------------------------------------------------------------------
# guide_donor_functions.remove_BspQI_sites_from_donor_in_final_oligo
# ---------------------------------------------------------------------------

def test_guide_donor_functions_records_every_changed_position():
    new_donor, junction_changes = remove_BspQI_sites_from_donor_in_final_oligo(
        FWD_PRIMER, GUIDE, INTERNAL_CLONING_SITE, DONOR, REV_PRIMER, (),
    )

    assert new_donor is not None, "the fixture is repairable; repair must succeed"
    oligo = FWD_PRIMER + GUIDE + INTERNAL_CLONING_SITE + new_donor + REV_PRIMER
    assert _n_sites(oligo) == 1, f"repair must leave exactly 1 site; got {_n_sites(oligo)}"

    changed = _changed_positions(DONOR, new_donor)
    named = _positions_named(junction_changes)

    assert len(junction_changes) == len(changed), (
        f"the donor was changed at {len(changed)} position(s) {changed} but "
        f"junction_changes records {len(junction_changes)}: {junction_changes}"
    )
    assert named == set(changed), (
        f"junction_changes names positions {sorted(named)} but the donor actually "
        f"changed at {sorted(changed)}"
    )


# ---------------------------------------------------------------------------
# Positive controls
# ---------------------------------------------------------------------------

def test_clean_oligo_is_returned_unchanged_with_no_records():
    """An oligo already carrying exactly one site needs no repair and no records."""
    clean_cloning_site = RC_SITE + "ATCATC"
    clean_rev = "TTACCTACTGCA"
    clean_donor = "ATCATCATCATCATCATCATCATCATCATC"
    assert _n_sites(FWD_PRIMER + GUIDE + clean_cloning_site + clean_donor + clean_rev) == 1

    new_donor, junction_changes = remove_bspqi_sites_from_donor(
        FWD_PRIMER, GUIDE, clean_cloning_site, clean_donor, clean_rev,
    )
    assert new_donor == clean_donor
    assert junction_changes == []


def test_single_ended_repair_records_exactly_one_change():
    """One junction site -> one changed position -> one record."""
    cloning_site = RC_SITE + "ATC" + "GAAGAG"   # donor[0] == "C" completes site 2
    rev = "TTACCTACTGCA"                        # no site at the right junction
    donor = "C" + "ATCATCATCATCATCATCATCATCATCA" + "T"
    assert _n_sites(FWD_PRIMER + GUIDE + cloning_site + donor + rev) == 2

    new_donor, junction_changes = remove_bspqi_sites_from_donor(
        FWD_PRIMER, GUIDE, cloning_site, donor, rev,
    )
    assert new_donor is not None
    changed = _changed_positions(donor, new_donor)
    assert len(changed) == 1, f"expected one changed position; got {changed}"
    assert len(junction_changes) == 1, f"got {junction_changes}"
    assert _positions_named(junction_changes) == set(changed)
