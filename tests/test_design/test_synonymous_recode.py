"""Tests for :mod:`magestic.design.synonymous_recode`.

Every test either asserts a property the recode must have, or plants a case it must REFUSE.
A recode that silently changed a protein would be worse than the missing design it replaces, so
the refusals matter as much as the successes.
"""

from __future__ import annotations

import pytest

from magestic.design.synonymous_recode import (
    AA_TO_CODONS,
    CODON_TO_AA,
    recode_donor_to_remove_sites,
)

SITE = "GCTCTTC"          # BspQI / SapI
COMP = str.maketrans("ACGT", "TGCA")


def rc(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def translate(seq: str) -> str:
    return "".join(CODON_TO_AA.get(seq[i: i + 3], "?") for i in range(0, len(seq) - 2, 3))


def _plus_strand_gene(donor_body: str, pad: int = 30):
    """A single-block plus-strand CDS that starts at genomic position 1."""
    cds = "ATG" + donor_body + "TAA"
    # pad so the donor is a window inside the CDS
    cds = cds + "A" * ((3 - len(cds) % 3) % 3)
    blocks = [(1, len(cds))]
    return cds, blocks


def test_no_site_returns_donor_unchanged():
    cds, blocks = _plus_strand_gene("GAAGAAGAAGAAGAA")
    donor = cds
    out, changes = recode_donor_to_remove_sites(donor, 1, blocks, "+", cds, [SITE])
    assert out == donor
    assert changes == []


def test_removes_a_site_and_preserves_the_protein():
    # put a BspQI site in frame, spanning codons, inside a CDS
    body = "GAA" + SITE + "GAAGAAGAA"
    cds, blocks = _plus_strand_gene(body)
    donor = cds
    assert SITE in donor
    before = translate(cds)

    out, changes = recode_donor_to_remove_sites(donor, 1, blocks, "+", cds, [SITE])
    assert out is not None, "a synonymous fix should exist here"
    assert SITE not in out.upper() and rc(SITE) not in out.upper()
    assert changes, "a successful recode must report what it changed"
    # the protein is untouched
    assert translate(out) == before
    # and only whole codons moved
    for pos, frm, to in changes:
        assert frm != to
        assert 1 <= pos <= len(cds)


def test_does_not_introduce_a_new_site():
    body = "GAA" + SITE + "GAAGCTCTTGAA"
    cds, blocks = _plus_strand_gene(body)
    out, _ = recode_donor_to_remove_sites(cds, 1, blocks, "+", cds, [SITE])
    if out is not None:
        assert out.upper().count(SITE) == 0
        assert out.upper().count(rc(SITE)) == 0


def test_refuses_to_touch_protected_edit_bases():
    body = "GAA" + SITE + "GAAGAAGAA"
    cds, blocks = _plus_strand_gene(body)
    # protect every base of the site and its codons: nothing may be changed there
    protected = range(1, len(cds) + 1)
    out, changes = recode_donor_to_remove_sites(
        cds, 1, blocks, "+", cds, [SITE], protected=protected
    )
    assert out is None
    assert changes == []


def test_returns_none_when_no_synonymous_fix_exists():
    # ATG (Met) and TGG (Trp) have no synonym; build the site only out of those codons
    # so every candidate codon is unchangeable.
    cds = "ATG" + "ATGTGG" + "TAA"
    blocks = [(1, len(cds))]
    # force a site that overlaps only Met/Trp codons by asking for a pattern present there
    pseudo_site = "ATGTGG"
    out, changes = recode_donor_to_remove_sites(cds, 1, blocks, "+", cds, [pseudo_site])
    assert out is None
    assert changes == []


def test_minus_strand_gene_protein_preserved():
    body = "GAA" + SITE + "GAAGAAGAA"
    coding, _ = _plus_strand_gene(body)         # translation-orientation CDS
    genomic = rc(coding)                        # what sits on the plus strand
    blocks = [(1, len(genomic))]
    before = translate(coding)

    out, changes = recode_donor_to_remove_sites(
        genomic, 1, blocks, "-", coding, [SITE]
    )
    assert out is not None
    assert SITE not in out.upper() and rc(SITE) not in out.upper()
    # re-derive the coding sequence from the recoded plus-strand donor
    assert translate(rc(out)) == before
    assert changes


def test_every_amino_acid_with_a_synonym_is_recodable_in_principle():
    """Sanity on the table itself: the codon table and its inverse must agree."""
    for codon, amino in CODON_TO_AA.items():
        assert codon in AA_TO_CODONS[amino]
    assert len(CODON_TO_AA) == 64
    # Met and Trp are the only amino acids with exactly one codon
    singles = {aa for aa, cods in AA_TO_CODONS.items() if len(cods) == 1}
    assert singles == {"M", "W"}


@pytest.mark.parametrize("strand", ["+", "-"])
def test_reports_changes_as_genomic_positions(strand):
    body = "GAA" + SITE + "GAAGAAGAA"
    coding, _ = _plus_strand_gene(body)
    genomic = coding if strand == "+" else rc(coding)
    blocks = [(101, 100 + len(genomic))]
    out, changes = recode_donor_to_remove_sites(
        genomic, 101, blocks, strand, coding, [SITE]
    )
    if out is not None:
        for pos, _frm, _to in changes:
            assert 101 <= pos <= 100 + len(genomic), "positions must be genomic, not donor offsets"
