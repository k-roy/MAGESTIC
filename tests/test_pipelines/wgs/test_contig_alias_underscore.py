"""GATE 3: contig names containing '_' must survive TRACE's positional variant-id parse.

The bug being guarded, verbatim from `trace_crispr.wgs.donor_outcome.precall_sample`:

    v_id = f"{v_chr}_{v_pos}_{v_ref}_to_{v_alt[0]}"
    ...
    sp = v_id.split("_")
    start_pos = int(sp[1]) - offset - 1   # <- sp[1] is only the position if chrom has no '_'
    ref_len = len(sp[2])
    combine[start_pos] = sp[4]

So the test reproduces that parse directly on real composite-reference contig names and
asserts it breaks before the alias and works after it.
"""

from __future__ import annotations

import pytest

from magestic.pipelines.wgs.core.contig_alias import ContigAlias, alias_rows

# Real contigs of yKR1033_1036_post_T2T_full.fasta, the reference the WGS
# reference-mapping policy mandates for this cohort.
REAL_CONTIGS = [
    "chrI",
    "chrII",
    "chrMito",
    "unplaced_contig_34",
    "blank_MAGESTIC3_pF835_guide_donor_plasmid",
    "native_2micron",
    # the forms the spec calls out for other cohorts
    "plasmid_o1",
    "integrated_plasmid_o1",
]


def trace_positional_parse(v_id: str) -> tuple[int, int, str]:
    """Exactly what precall_sample does. Raises ValueError on an underscored contig."""
    sp = v_id.split("_")
    return int(sp[1]), len(sp[2]), sp[4]


def make_id(chrom: str, pos: int = 1234, ref: str = "AC", alt: str = "G") -> str:
    return f"{chrom}_{pos}_{ref}_to_{alt}"


# --------------------------------------------------------------------------- the bug


@pytest.mark.parametrize(
    "chrom",
    [
        "unplaced_contig_34",
        "blank_MAGESTIC3_pF835_guide_donor_plasmid",
        "native_2micron",
        "integrated_plasmid_o1",
    ],
)
def test_trace_parse_really_does_break_on_underscored_contigs(chrom: str) -> None:
    """Without the alias, the parse raises — this is the gate's empirical basis."""
    with pytest.raises(ValueError):
        trace_positional_parse(make_id(chrom))


@pytest.mark.parametrize("chrom", ["chrI", "chrII", "chrMito"])
def test_trace_parse_is_fine_on_safe_contigs(chrom: str) -> None:
    pos, ref_len, alt = trace_positional_parse(make_id(chrom))
    assert (pos, ref_len, alt) == (1234, 2, "G")


# ------------------------------------------------------------------------- the fix


def test_aliased_ids_parse_correctly_under_traces_own_parse() -> None:
    alias = ContigAlias.from_names(REAL_CONTIGS)
    for chrom in REAL_CONTIGS:
        v_id = make_id(alias.encode(chrom), pos=1234, ref="AC", alt="G")
        pos, ref_len, alt = trace_positional_parse(v_id)
        assert (pos, ref_len, alt) == (1234, 2, "G"), chrom


def test_alias_is_bijective_and_underscore_free() -> None:
    alias = ContigAlias.from_names(REAL_CONTIGS)
    for chrom in REAL_CONTIGS:
        assert alias.decode(alias.encode(chrom)) == chrom
        assert "_" not in alias.encode(chrom)
    # one alias per name, no merging
    assert len({alias.encode(c) for c in REAL_CONTIGS}) == len(set(REAL_CONTIGS))


def test_safe_names_are_left_untouched_so_old_outputs_are_unchanged() -> None:
    alias = ContigAlias.from_names(REAL_CONTIGS)
    for safe in ["chrI", "chrII", "chrMito"]:
        assert alias.encode(safe) == safe
    assert set(alias.aliased) == {
        "unplaced_contig_34",
        "blank_MAGESTIC3_pF835_guide_donor_plasmid",
        "native_2micron",
        "plasmid_o1",
        "integrated_plasmid_o1",
    }


def test_noop_on_an_already_safe_reference() -> None:
    alias = ContigAlias.from_names(["chrI", "chrII", "chrMito"])
    assert alias.is_noop
    assert alias.aliased == ()


def test_alias_cannot_shadow_a_real_contig() -> None:
    """A reference that already contains 'ctg0001' must not have it stolen."""
    alias = ContigAlias.from_names(["ctg0001", "plasmid_o1", "chrI"])
    assert alias.encode("ctg0001") == "ctg0001"
    assert alias.encode("plasmid_o1") != "ctg0001"
    assert alias.decode(alias.encode("plasmid_o1")) == "plasmid_o1"


def test_deterministic_across_input_order() -> None:
    a = ContigAlias.from_names(REAL_CONTIGS)
    b = ContigAlias.from_names(list(reversed(REAL_CONTIGS)))
    assert a.encode_map == b.encode_map


def test_decode_variant_id_round_trips_including_the_LC_prefix() -> None:
    alias = ContigAlias.from_names(REAL_CONTIGS)
    for chrom in REAL_CONTIGS:
        enc = make_id(alias.encode(chrom), pos=77, ref="A", alt="T")
        assert alias.decode_variant_id(enc) == make_id(chrom, pos=77, ref="A", alt="T")
        # precall_sample prefixes low-confidence calls with "(LC)"
        assert alias.decode_variant_id(f"(LC){enc}") == f"(LC){make_id(chrom, 77, 'A', 'T')}"


def test_decode_variant_id_leaves_an_unexpected_shape_alone() -> None:
    alias = ContigAlias.from_names(REAL_CONTIGS)
    for junk in ["", "chrI", "chrI_100", "chrI_100_A_G", "no-underscores-at-all"]:
        assert alias.decode_variant_id(junk) == junk


def test_alias_rows_copies_and_does_not_mutate_the_caller() -> None:
    alias = ContigAlias.from_names(REAL_CONTIGS)
    rows = [{"chrom": "native_2micron", "pos": 5}, {"chrom": "chrI", "pos": 6}]
    out = alias_rows(rows, alias)
    assert rows[0]["chrom"] == "native_2micron"  # original untouched
    assert "_" not in str(out[0]["chrom"])
    assert out[1]["chrom"] == "chrI"
    assert out[0]["pos"] == 5


def test_multi_base_ref_and_alt_survive() -> None:
    """ref/alt are nucleotide strings of any length; the right-anchored decode must cope."""
    alias = ContigAlias.from_names(REAL_CONTIGS)
    enc = make_id(alias.encode("native_2micron"), pos=9, ref="ACTACTAC", alt="AA")
    pos, ref_len, alt = trace_positional_parse(enc)
    assert (pos, ref_len, alt) == (9, 8, "AA")
    assert alias.decode_variant_id(enc) == make_id("native_2micron", 9, "ACTACTAC", "AA")
