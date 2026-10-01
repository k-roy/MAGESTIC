"""Tests for magestic.design.strain_background.

Each test targets one of the four traps the module exists to close, plus the
coordinate arithmetic, which is the part an off-by-one would silently ruin: a design
near an indel would be built on the wrong sequence and fail at the bench, not in CI.
"""

from __future__ import annotations

import textwrap

import pytest

from magestic.design.strain_background import (
    CoordinateMap,
    StrainVariant,
    annotate_variants_with_features,
    build_strain_reference,
    load_strain_variants,
    project_variants,
    read_fasta,
    remap_gff,
    verify_designs_against_assembly,
    write_fasta,
)


# ---------------------------------------------------------------------------
# Trap 1 — dropping contigs
# ---------------------------------------------------------------------------


def test_read_fasta_keeps_mito_and_two_micron(tmp_path):
    """A reader that stops at the mitochondrial defline loses chrM and the 2-micron."""
    p = tmp_path / "ref.fsa"
    p.write_text(
        ">chrI nuclear\nACGTACGTAC\n"
        ">chrM mitochondrion\nTTTTTTTT\n"
        ">two_micron A364A\nGGGGCCCC\n"
    )
    seqs = read_fasta(p)
    assert set(seqs) == {"chrI", "chrM", "two_micron"}
    assert seqs["two_micron"] == "GGGGCCCC"


def test_build_reference_returns_every_contig_even_without_variants():
    ref = {"chrI": "ACGTACGTAC", "chrM": "TTTTTTTT", "two_micron": "GGGGCCCC"}
    out, _, _ = build_strain_reference(ref, [StrainVariant("chrI", 1, "A", "G")])
    assert set(out) == set(ref)
    assert out["chrM"] == ref["chrM"]
    assert out["two_micron"] == ref["two_micron"]


# ---------------------------------------------------------------------------
# Coordinate arithmetic
# ---------------------------------------------------------------------------


def test_snv_does_not_move_coordinates():
    ref = {"chrI": "AAAAACAAAAA"}
    out, cmap, rej = build_strain_reference(ref, [StrainVariant("chrI", 6, "C", "T")])
    assert out["chrI"] == "AAAAATAAAAA"
    assert rej == []
    for pos in (1, 5, 6, 7, 11):
        assert cmap.ref_to_strain("chrI", pos) == pos


def test_deletion_shifts_downstream_and_deleted_base_maps_to_none():
    #            123456789
    ref = {"chrI": "AAACCCGGG"}
    # Delete the three Cs: REF spans 4..6, ALT keeps only position 4's base.
    v = StrainVariant("chrI", 4, "CCC", "C")
    out, cmap, rej = build_strain_reference(ref, [v])
    assert out["chrI"] == "AAACGGG"
    assert rej == []
    assert cmap.ref_to_strain("chrI", 3) == 3          # upstream unchanged
    assert cmap.ref_to_strain("chrI", 4) == 4          # retained anchor base
    assert cmap.ref_to_strain("chrI", 5) is None       # deleted
    assert cmap.ref_to_strain("chrI", 6) is None       # deleted
    assert cmap.ref_to_strain("chrI", 7) == 5          # shifted by -2
    assert cmap.total_length_change("chrI") == -2


def test_insertion_shifts_downstream():
    ref = {"chrI": "AAACGGG"}
    v = StrainVariant("chrI", 4, "C", "CTT")
    out, cmap, _ = build_strain_reference(ref, [v])
    assert out["chrI"] == "AAACTTGGG"
    assert cmap.ref_to_strain("chrI", 3) == 3
    assert cmap.ref_to_strain("chrI", 5) == 7          # shifted by +2
    assert cmap.total_length_change("chrI") == 2


def test_coordinate_round_trip_across_mixed_indels():
    ref = {"chrI": "ACGT" * 25}  # 100 bp; 1-based pos p: p%4 -> 1=A 2=C 3=G 0=T
    variants = [
        StrainVariant("chrI", 10, "C", "A"),            # SNV       (pos 10 % 4 == 2 -> C)
        StrainVariant("chrI", 20, "TACG", "T"),         # 3 bp deletion
        StrainVariant("chrI", 50, "C", "CGGGG"),        # 4 bp insertion
    ]
    _, cmap, rej = build_strain_reference(ref, variants)
    assert rej == []
    for pos in range(1, 101):
        sp = cmap.ref_to_strain("chrI", pos)
        if sp is None:
            continue
        assert cmap.strain_to_ref("chrI", sp) == pos, f"round trip failed at {pos}"


# ---------------------------------------------------------------------------
# Trap 2 — two records at one position; Trap 3 — application order
# ---------------------------------------------------------------------------


def test_two_indels_on_one_contig_both_applied_correctly():
    """Left-to-right application would place the second indel at the wrong offset."""
    ref = {"chrI": "AAAACCCCGGGGTTTT"}
    variants = [
        StrainVariant("chrI", 5, "CCCC", "C"),      # -3
        StrainVariant("chrI", 13, "T", "TAAA"),     # +3
    ]
    out, cmap, rej = build_strain_reference(ref, variants)
    assert rej == []
    assert out["chrI"] == "AAAACGGGGTAAATTT"
    assert cmap.total_length_change("chrI") == 0


def test_overlapping_variants_are_rejected_not_applied():
    ref = {"chrI": "AAACCCGGG"}
    variants = [
        StrainVariant("chrI", 4, "CCC", "C"),
        StrainVariant("chrI", 5, "C", "T"),   # inside the deletion above
    ]
    out, _, rej = build_strain_reference(ref, variants)
    assert len(rej) == 1
    assert "overlaps" in rej[0].reason
    assert out["chrI"] == "AAACGGG"


def test_ref_allele_mismatch_is_rejected():
    """The cheapest detector for a VCF called against a different assembly."""
    ref = {"chrI": "AAACCCGGG"}
    out, _, rej = build_strain_reference(ref, [StrainVariant("chrI", 4, "G", "T")])
    assert len(rej) == 1
    assert "REF mismatch" in rej[0].reason
    assert out["chrI"] == ref["chrI"]


# ---------------------------------------------------------------------------
# Trap 4 — genotype is consulted
# ---------------------------------------------------------------------------


def _write_vcf(tmp_path, body: str):
    p = tmp_path / "v.vcf"
    p.write_text(textwrap.dedent(body).lstrip())
    return p


def test_genotype_filters_hom_ref_and_het(tmp_path):
    vcf = _write_vcf(tmp_path, """
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM22
        chrI\t10\t.\tA\tG\t60\tPASS\t.\tGT\t1/1
        chrI\t20\t.\tA\tG\t60\tPASS\t.\tGT\t0/0
        chrI\t30\t.\tA\tG\t60\tPASS\t.\tGT\t0/1
        chrI\t40\t.\tA\tG\t60\tPASS\t.\tGT\t./.
        """)
    variants, rej = load_strain_variants(vcf, "M22")
    assert [v.pos for v in variants] == [10]
    assert len(rej) == 3


def test_multiallelic_second_allele_is_selected(tmp_path):
    vcf = _write_vcf(tmp_path, """
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPW5
        chrI\t10\t.\tA\tG,T\t60\tPASS\t.\tGT\t2/2
        """)
    variants, _ = load_strain_variants(vcf, "PW5")
    assert len(variants) == 1
    assert variants[0].alt == "T"


def test_filter_and_qual_are_honored(tmp_path):
    vcf = _write_vcf(tmp_path, """
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
        chrI\t10\t.\tA\tG\t60\tLowQual\t.\tGT\t1/1
        chrI\t20\t.\tA\tG\t5\tPASS\t.\tGT\t1/1
        chrI\t30\t.\tA\tG\t60\tPASS\t.\tGT\t1/1
        """)
    variants, rej = load_strain_variants(vcf, "S", min_qual=30)
    assert [v.pos for v in variants] == [30]
    reasons = " ".join(r.reason for r in rej)
    assert "FILTER=LowQual" in reasons and "QUAL=5" in reasons


def test_two_records_at_same_position_both_survive_loading(tmp_path):
    """A position-keyed dict would keep only one of these."""
    vcf = _write_vcf(tmp_path, """
        ##fileformat=VCFv4.2
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
        chrI\t10\t.\tA\tG\t60\tPASS\t.\tGT\t1/1
        chrI\t10\t.\tAT\tA\t60\tPASS\t.\tGT\t1/1
        """)
    variants, _ = load_strain_variants(vcf, "S")
    assert len(variants) == 2


# ---------------------------------------------------------------------------
# GFF remapping
# ---------------------------------------------------------------------------


def test_remap_gff_shifts_feature_past_an_indel(tmp_path):
    ref = {"chrI": "A" * 100}
    _, cmap, _ = build_strain_reference(
        ref, [StrainVariant("chrI", 10, "AAAAA", "A")]  # -4
    )
    gff = tmp_path / "in.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chrI\t.\tgene\t1\t5\t.\t+\t.\tID=UPSTREAM\n"
        "chrI\t.\tgene\t50\t60\t.\t+\t.\tID=DOWNSTREAM\n"
    )
    out = tmp_path / "out.gff"
    counts = remap_gff(gff, cmap, out)
    lines = [l for l in out.read_text().splitlines() if not l.startswith("#")]
    assert lines[0].split("\t")[3:5] == ["1", "5"]       # upstream unmoved
    assert lines[1].split("\t")[3:5] == ["46", "56"]     # shifted by -4
    assert counts["written"] == 2


# ---------------------------------------------------------------------------
# Cross-strain projection
# ---------------------------------------------------------------------------


def test_project_variants_drops_alleles_the_recipient_already_has():
    """A variant the recipient already carries is not an edit; an oligo for it edits nothing."""
    ref = {"chrI": "AAACAAAGAAA"}
    recipient, cmap, _ = build_strain_reference(ref, [StrainVariant("chrI", 4, "C", "T")])
    assert recipient["chrI"] == "AAATAAAGAAA"

    donor = [
        StrainVariant("chrI", 4, "C", "T", sample="B"),   # recipient already T
        StrainVariant("chrI", 8, "G", "A", sample="B"),   # real edit
    ]
    projected, rej = project_variants(donor, cmap, recipient)
    assert [v.pos for v in projected] == [8]
    assert "already carries" in rej[0].reason
    assert projected[0].ref == "G"


def test_project_variants_shifts_coordinates_for_recipient_indels():
    ref = {"chrI": "AAACCCGGGTTT"}
    recipient, cmap, _ = build_strain_reference(
        ref, [StrainVariant("chrI", 4, "CCC", "C")]       # -2
    )
    projected, rej = project_variants(
        [StrainVariant("chrI", 10, "T", "A", sample="B")], cmap, recipient
    )
    assert rej == []
    assert projected[0].pos == 8          # 10 - 2
    assert projected[0].ref == recipient["chrI"][7]


# ---------------------------------------------------------------------------
# Feature annotation — never leave ORF empty
# ---------------------------------------------------------------------------


def test_noncoding_variants_get_a_real_label_not_an_empty_orf(tmp_path):
    """An empty ORF is dropped silently downstream by filter_variants."""
    gff = tmp_path / "a.gff"
    gff.write_text(
        "##gff-version 3\n"
        "chrI\t.\tgene\t1000\t2000\t.\t+\t.\tID=YAL001C;gene=TFC3\n"
    )
    variants = [
        StrainVariant("chrI", 1500, "A", "G"),   # inside the gene
        StrainVariant("chrI", 900, "A", "G"),    # 100 bp upstream (+ strand)
        StrainVariant("chrI", 50_000, "A", "G"), # distal
    ]
    out = annotate_variants_with_features(variants, gff, flank=1000)
    assert out[0].orf == "YAL001C"
    assert out[1].orf == "upstream:YAL001C"
    assert out[2].orf == "intergenic"
    assert all(v.orf for v in out), "no variant may be left with an empty ORF"


# ---------------------------------------------------------------------------
# Assembly verification
# ---------------------------------------------------------------------------


def test_verify_designs_accepts_revcomp_and_rejects_absent_sequence():
    assembly = {"ctg1": "AAAACCCCGGGGTTTTACGTACGTAAGGCCTT"}
    designs = [
        {"id": "ok_fwd", "guide": "CCCCGGGG",
         "left_homology_arm": "AAAACCCC", "right_homology_arm": "GGGGTTTT"},
        {"id": "ok_rev", "guide": _rc("ACGTACGT"),
         "left_homology_arm": _rc("TTTTACGT"), "right_homology_arm": _rc("ACGTAAGG")},
        {"id": "bad", "guide": "TTTTTTTTTTTT",
         "left_homology_arm": "AAAACCCC", "right_homology_arm": "GGGGTTTT"},
    ]
    passed, failed = verify_designs_against_assembly(designs, assembly)
    assert {d["id"] for d in passed} == {"ok_fwd", "ok_rev"}
    assert len(failed) == 1
    assert failed[0]["id"] == "bad"
    assert "guide" in failed[0]["failure_reason"]


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# ---------------------------------------------------------------------------
# Round trip through disk
# ---------------------------------------------------------------------------


def test_write_then_read_fasta_round_trips(tmp_path):
    seqs = {"chrI": "ACGT" * 40, "two_micron": "GGCC" * 10}
    p = tmp_path / "o.fa"
    write_fasta(seqs, p, width=60)
    assert read_fasta(p) == seqs
