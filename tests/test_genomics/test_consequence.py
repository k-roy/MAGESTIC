"""
Coding-consequence annotation must never return a blank, and must handle every
variant class -- not only single-base substitutions.

The annotators used before `magestic.genomics.consequence` emitted a call only for
single-base substitutions. Everything else came out EMPTY. Measured on a shipped
24,000-oligo library (15,641 natural designs): 100 % of deletions (914/914) and 100 %
of insertions (859/859) had a blank consequence, plus 882 of 13,868 substitutions,
essentially all of them multi-base. Of the 2,655 blanks, 365 sat INSIDE a CDS, where
the consequence is fully determined by the sequence -- 257 substitutions (120 of them
4 nt spanning two codons), 75 frameshift indels and 33 in-frame indels.

Guide-donor design keys on the consequence: the synonymous-recode rescue, the PTC arm,
the missense arm, and which variants deserve extra guides all depend on it. A blank is
therefore a design defect and not a cosmetic gap.

The contract these tests pin down:
  - every call returns a member of CONSEQUENCES; `""` is impossible, by construction
  - a multi-codon MNP is classified AND reports every residue it changes
  - an indel is classified as in-frame or frameshift, and a 3n deletion is NOT nonsense
  - `non_coding` is a positive statement about an edit outside the CDS
  - a CDS table that parses to nothing raises, rather than making everything non_coding
"""

import pytest

from magestic.genomics.consequence import (
    CONSEQUENCES,
    VariantConsequence,
    annotate_overlapping,
    annotate_variant,
    build_cds,
    build_gene_index,
    genes_at,
    load_cds_table,
    translate,
)

# ATG AAA CCC GGG TTT TGG TAA  ->  M K P G F W *
PLUS = "ATGAAACCCGGGTTTTGGTAA"
START = 4
GENOME = {"chrT": "GGGG" + PLUS + "GGGG"}
MODEL = {"PLUS": ("chrT", "+", [(START, START + len(PLUS))])}


def _call(pos_offset, ref, alt, **kw):
    return annotate_variant(
        GENOME, MODEL, "PLUS", "chrT", START + pos_offset, ref, alt, **kw
    )


def test_reference_translation_is_what_the_fixture_claims():
    assert translate(build_cds(GENOME["chrT"], MODEL["PLUS"][2], "+")) == "MKPGFW*"


def test_a_blank_consequence_cannot_be_constructed():
    with pytest.raises(ValueError):
        VariantConsequence("")


@pytest.mark.parametrize(
    "offset,ref,alt,expected",
    [
        (5, "A", "G", "synonymous"),           # AAA -> AAG
        (3, "A", "G", "missense"),             # AAA -> GAA
        (3, "AAA", "TAA", "nonsense"),         # a stop gained by an MNP
        (0, "A", "G", "start_lost"),
        (18, "T", "C", "stop_lost"),
        (4, "AAC", "TTG", "missense"),         # MNP across two codons
        (3, "AAA", "", "inframe_deletion"),    # 3n deletion is NOT nonsense
        (3, "A", "", "frameshift"),
        (3, "", "GGG", "inframe_insertion"),
        (3, "", "GG", "frameshift"),
    ],
)
def test_every_variant_class_is_classified(offset, ref, alt, expected):
    got = _call(offset, ref, alt)
    assert got.consequence == expected
    assert got.consequence in CONSEQUENCES
    assert got.consequence != ""


def test_multi_codon_mnp_reports_every_changed_residue():
    # AAC->TTG at PLUS[4:7] straddles codon 2 (AAA->ATT) and codon 3 (CCC->GCC).
    got = _call(4, "AAC", "TTG")
    assert got.aa_changes == ("K2I", "P3A")
    assert got.aa_change in got.aa_changes


def test_inframe_deletion_is_not_reported_as_nonsense():
    # A clean 3 nt deletion moves the terminal stop one residue earlier. Comparing stop
    # positions without allowing for that reports every in-frame deletion as nonsense.
    assert _call(3, "AAA", "").consequence == "inframe_deletion"


def test_non_coding_is_a_positive_answer_not_a_blank():
    got = annotate_variant(GENOME, MODEL, "PLUS", "chrT", 0, "G", "A")
    assert got.consequence == "non_coding"
    assert got.aa_change == ""
    absent = annotate_variant(GENOME, MODEL, "ABSENT", "chrT", START + 3, "A", "G")
    assert absent.consequence == "non_coding"


def test_minus_strand_is_handled_not_merely_accepted():
    rc = PLUS.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    genome = {"chrM": "GGGG" + rc + "GGGG"}
    model = {"MINUS": ("chrM", "-", [(4, 4 + len(rc))])}
    assert translate(build_cds(genome["chrM"], model["MINUS"][2], "-")) == "MKPGFW*"
    got = annotate_variant(
        genome, model, "MINUS", "chrM", 4 + len(rc) - 1 - 3, "T", "C"
    )
    assert got.consequence == "missense"
    assert got.aa_change == "K2E"


def test_a_wrong_reference_allele_raises_rather_than_annotating():
    with pytest.raises(ValueError):
        _call(3, "C", "G")


def test_overlapping_genes_each_get_their_own_answer():
    model = dict(MODEL)
    model["OTHER"] = ("chrT", "-", [(START, START + len(PLUS))])
    index = build_gene_index(model)
    assert genes_at(index, "chrT", START + 3) == ["OTHER", "PLUS"]
    results = annotate_overlapping(
        GENOME, model, index, "chrT", START + 3, "A", "G"
    )
    assert len(results) == 2
    assert {r.gene for r in results} == {"OTHER", "PLUS"}


def test_an_empty_cds_table_raises_instead_of_silently_non_coding(tmp_path):
    # The design-frame CDS tables are HEADERLESS. A parser that eats the first row as a
    # header yields zero intervals, which makes every position read as outside a CDS --
    # invisible unless the interval count is checked.
    path = tmp_path / "cds.tsv"
    path.write_text("chrom\tstart\tend\tstrand\tphase\tattrs\n")
    with pytest.raises(ValueError):
        load_cds_table(path)
