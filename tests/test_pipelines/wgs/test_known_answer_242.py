"""
NN 242 — known-answer reproductions for the three WGS trace/edit-calling defects,
plus the planted-truth confusion assertions that back the per-class rates.

Scope: `magestic` @ 9d654b8 over the PINNED TRACE sha 7c4c10f1 (the local clone
`~/repos/trace`, whose HEAD *is* that sha; the sha exists in no GitHub ref).

Each test either REPRODUCES a defect (and is marked `xfail(strict=True)` so the
shared tree is not left red while the defect stands) or PINS the measured
behaviour so a future fix cannot silently change a reported rate.

Full method and numbers:
  planning/222_triumph_embl_library_design/242_wgs_trace_known_answer.md
  planning/222_triumph_embl_library_design/242_known_answer_sim.py

Author: NN 242 (wgs-knownanswer), 2026-10-01
"""

from __future__ import annotations

from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")
trace_wgs = pytest.importorskip("trace_crispr.wgs")

from trace_crispr.wgs.donor_outcome import DonorSpec, precall_sample  # noqa: E402
from trace_crispr.wgs.reconcile import (  # noqa: E402
    CONTAMINATION,
    EditSite,
    LOW_CONF,
    NO_EDIT,
    SINGLE_EDIT,
    SITE_NOCOV,
    SITE_WT,
    SiteCall,
    classify_site,
    reconcile_cell,
)
from trace_crispr.wgs.sv_plasmid import load_blacklist  # noqa: E402


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #
def _site(vaf: float, cov: int, *, chrom="chrXV", start=1000, end=1010) -> SiteCall:
    return SiteCall(
        site=EditSite(chrom=chrom, start=start, end=end, oligo="o1"),
        vaf=vaf,
        cov=cov,
        status=classify_site(vaf, cov),
    )


def _cell(sites, **kw):
    kw.setdefault("n_genome_bc", 1)
    kw.setdefault("n_plasmid_bc", 1)
    return reconcile_cell("s", list(sites), **kw)


def _write_vcf(path: Path, records, *, gvcf: bool = False, chrom="chrXV", window=(1000, 1010),
               ref_seq: str | None = None):
    """Minimal VCF/gVCF in the FORMAT layout TRACE's positional parser assumes
    (`GT:AD:...`, i.e. AD in slot 2 — what `core/vcf_adapter.py` remaps 4 -> 2 to produce).

    `records` are (pos, ref, alt, ad_ref, ad_alt). When `gvcf` is set, a <NON_REF>
    reference block is emitted for every position in `window`, which is what gives
    `precall_sample` the `flag_covered` evidence that the `unedited` branch requires.

    `ref_seq` (1-based: `ref_seq[pos - 1]`) supplies each block's REF. planning/268: the
    pinned TRACE validates every gVCF REF against the reference, so the former constant
    `A` made a fixture that no real caller could write.
    """
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={chrom}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
    ]
    if gvcf:
        for p in range(window[0], window[1] + 1):
            ref = ref_seq[p - 1] if ref_seq else "A"
            lines.append(f"{chrom}\t{p}\t.\t{ref}\t<NON_REF>\t.\t.\t.\tGT:AD\t0/0:30,0")
    for pos, ref, alt, adr, ada in records:
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:AD\t0/1:{adr},{ada}")
    path.write_text("\n".join(lines) + "\n")
    return path


# =========================================================================== #
# DEFECT 1 — Stage 3's blacklist is inert and fails silently                   #
# =========================================================================== #
def test_blacklist_parser_accepts_contigs_that_match_nothing(tmp_path):
    """REPRODUCTION (documents current behaviour).

    `load_blacklist` is permissive by design (it ports Shengdi's `split("\\t")`).
    A blacklist written in `chr15` style therefore loads cleanly against a
    `chrXV`-style reference and filters NOTHING, with no error and no warning.
    """
    bl_path = tmp_path / "blacklist.txt"
    bl_path.write_text("chr15\t100\t200\tREDI_cassette\nchr2\t300\t400\n")
    bl = load_blacklist(str(bl_path))

    # It parsed happily...
    assert set(bl) == {"chr15", "chr2"}
    # ...and nothing it names exists on a roman-numeral reference, so every lookup
    # is a no-op. This is the silent class: a 0% filtered rate is indistinguishable
    # from a correctly-configured blacklist that found nothing to filter.
    assert bl.get("chrXV") is None
    assert bl.get("chrII") is None


def test_mismatched_blacklist_silently_filters_nothing_through_the_real_filter(tmp_path):
    """REPRODUCTION through the ACTUAL filter, not just the parser.

    Two blacklists describe the SAME cassette region: one in `chrXV` naming (matching
    the reference) and one in `chr15` naming (the committed copy). A read inside the
    region is filtered by the first and NOT by the second -- and the second raises
    nothing, so a run configured with it reports a clean zero.
    """
    from trace_crispr.wgs.sv_plasmid import _is_mapped_to_cassette

    right = tmp_path / "right.txt"
    right.write_text("chrXV\t100\t200\tREDI_cassette\n")
    wrong = tmp_path / "wrong.txt"
    wrong.write_text("chr15\t100\t200\tREDI_cassette\n")

    bl_right = load_blacklist(str(right))
    bl_wrong = load_blacklist(str(wrong))

    read_coord = 150  # squarely inside the cassette
    assert _is_mapped_to_cassette("chrXV", read_coord, bl_right) is True
    # The SAME read, the SAME region, a differently-named blacklist -> not filtered,
    # and nothing anywhere raised. This is the silent class.
    assert _is_mapped_to_cassette("chrXV", read_coord, bl_wrong) is False


def test_magestic_step3_validates_blacklist_existence_but_not_its_contigs(tmp_path):
    """Pins WHERE the missing check belongs. magestic's Stage 3 does guard the
    blacklist -- but only `is None or not Path(...).exists()`
    (`step_3_sv_plasmid.py:258-262`). It never compares the blacklist's contig names
    against the BAM header or the reference, which is the check that would have
    caught the `chr15` copy."""
    src = Path(
        __import__(
            "magestic.pipelines.wgs.pipeline.step_3_sv_plasmid", fromlist=["x"]
        ).__file__
    ).read_text()
    assert "sv_plasmid.blacklist is" in src  # the existence guard is present
    # ...and no contig reconciliation is present anywhere in the step
    assert "references" not in src.split("blacklist")[1][:400]
    for token in ("blacklist_contigs", "unknown_contig", "contig_mismatch"):
        assert token not in src, f"unexpected: {token} suggests the check now exists"


def test_blacklist_coordinate_bounds_are_strict_both_sides(tmp_path):
    """A second silent edge: `start < coord < end` is STRICT, so a read at exactly
    the blacklist boundary is NOT filtered. Pins the off-by-one so a later fix is
    a deliberate change rather than a drift in a reported rate."""
    from trace_crispr.wgs.sv_plasmid import _is_mapped_to_cassette

    bl_path = tmp_path / "bl.txt"
    bl_path.write_text("chrXV\t100\t200\n")
    bl = load_blacklist(str(bl_path))
    assert _is_mapped_to_cassette("chrXV", 150, bl) is True
    assert _is_mapped_to_cassette("chrXV", 100, bl) is False  # boundary, NOT filtered
    assert _is_mapped_to_cassette("chrXV", 200, bl) is False  # boundary, NOT filtered


# =========================================================================== #
# DEFECT 2 — GATE 2: `unedited` is unobservable without a gVCF                 #
# =========================================================================== #
def _spec() -> DonorSpec:
    # A 21-base donor window whose design differs from WT at one position.
    wt = "ACGTACGTACGTACGTACGTA"
    design = wt[:10] + "T" + wt[11:]
    return DonorSpec(
        sample="S1",
        chrom="chrXV",
        donor_start=1000,
        donor_end=1020,
        donor_wt=wt,
        donor_design=design,
        donor_step2=design,
    )


def test_unedited_requires_a_gvcf_and_is_otherwise_unknown(tmp_path):
    """CHARACTERIZATION of GATE 2, the single most consequential finding for the
    22.6% "no edit" figure.

    A genuinely unedited clone, scored with a variant-only VCF (which is all the
    cohort has), is called `unknown` -- never `designed_variant_unedited`. The
    `unedited` branch is reached only when `flag_covered` clears
    `quantification_window_size * quantification_window_cov_perc`, and
    `flag_covered` is incremented ONLY inside the gVCF pass. So `unedited` is
    UNOBSERVABLE without a gVCF, not rare.
    """
    spec = _spec()
    chrom_seq = "N" * 999 + spec.donor_wt + "N" * 100

    plain = _write_vcf(tmp_path / "plain.vcf", [])          # truly unedited, no gVCF
    empty_g = _write_vcf(tmp_path / "empty.g.vcf", [])      # a gVCF with NO ref blocks
    out_no_gvcf = precall_sample(spec, plain, empty_g, chrom_seq)

    full_g = _write_vcf(tmp_path / "full.g.vcf", [], gvcf=True, window=(1000, 1020),
                        ref_seq=chrom_seq)
    out_with_gvcf = precall_sample(spec, plain, full_g, chrom_seq)

    assert out_no_gvcf.editing_outcome == "unknown", (
        "without gVCF ref blocks an unedited clone must fall through to `unknown`"
    )
    assert out_with_gvcf.editing_outcome == "designed_variant_unedited", (
        "with gVCF ref blocks over the window the same clone is callable as unedited"
    )
    # The two differ ONLY in whether a gVCF was supplied. Same clone, same VCF.
    assert out_no_gvcf.editing_outcome != out_with_gvcf.editing_outcome


@pytest.mark.xfail(
    strict=True,
    reason="DEFECT 242-2: `unknown` and `designed_variant_unedited` are distinct "
    "outcomes, but nothing in the emitted row records WHICH gate produced an "
    "`unknown` -- no-coverage, no-gVCF, or genuinely-ambiguous. A rate computed "
    "over `unknown` therefore cannot be apportioned.",
)
def test_unknown_should_record_why_it_is_unknown(tmp_path):
    """DESIRED: the outcome string should distinguish `unknown(no_gvcf)` from
    `unknown(low_coverage)`. It does not."""
    spec = _spec()
    chrom_seq = "N" * 999 + spec.donor_wt + "N" * 100
    plain = _write_vcf(tmp_path / "p.vcf", [])
    empty_g = _write_vcf(tmp_path / "e.g.vcf", [])
    out = precall_sample(spec, plain, empty_g, chrom_seq)
    assert out.editing_outcome != "unknown", f"bare `unknown` with no reason: {out.editing_outcome}"


# =========================================================================== #
# DEFECT 3 — split("_") crashes on an underscore-bearing contig                #
# =========================================================================== #
def test_precall_sample_handles_an_underscore_bearing_contig(tmp_path):
    """FORMER REPRODUCTION, now the closed state (planning/268).

    TRACE `7c4c10f` re-parsed its own `f"{chrom}_{pos}_{ref}_to_{alt}"` id with a
    positional `split("_")`, so every mandated composite contig (`plasmid_<oligo>`,
    `integrated_plasmid_<oligo>`, ...) raised `ValueError` mid-loop. The pinned TRACE
    builds ids from structured fields. WITHDRAWN: "precall_sample crashes on an
    underscore-bearing contig". Asserted instead: the SAME evidence gives the SAME
    outcome under either name, and the shipped ids carry the real contig name.
    """
    wt = "ACGTACGTACGTACGTACGTA"
    design = wt[:10] + "T" + wt[11:]
    chrom_seq = "N" * 999 + wt + "N" * 100

    outs = {}
    for name in ("chrXV", "integrated_plasmid_o1"):
        spec = DonorSpec("S1", name, 1000, 1020, wt, design, design)
        vcf = _write_vcf(tmp_path / f"{name}.vcf", [(1010, "G", "T", 0, 30)], chrom=name)
        g = _write_vcf(tmp_path / f"{name}.g.vcf", [], chrom=name)
        outs[name] = precall_sample(spec, vcf, g, chrom_seq)
    assert outs["chrXV"].editing_outcome == outs["integrated_plasmid_o1"].editing_outcome
    for name, out in outs.items():
        assert all(v.startswith(name + "_") for v in out.designed_variants), out


def test_single_underscore_contig_is_also_handled(tmp_path):
    """One underscore was enough to crash TRACE `7c4c10f` (`plasmid_o1_1010_C_to_T` ->
    sp[1] == 'o1'). WITHDRAWN with the defect, which the pinned TRACE closes."""
    wt = "ACGTACGTACGTACGTACGTA"
    design = wt[:10] + "T" + wt[11:]
    chrom_seq = "N" * 999 + wt + "N" * 100
    spec = DonorSpec("S3", "plasmid_o1", 1000, 1020, wt, design, design)
    vcf = _write_vcf(tmp_path / "p.vcf", [(1010, "G", "T", 0, 30)], chrom="plasmid_o1")
    g = _write_vcf(tmp_path / "p.g.vcf", [], chrom="plasmid_o1")
    out = precall_sample(spec, vcf, g, chrom_seq)
    assert all(v.startswith("plasmid_o1_") for v in out.designed_variants), out


def test_pipeline_qc_guard_for_the_underscore_contig_exists_and_fires():
    """magestic's `qc/pipeline_qc.py` claims to guard GATE 3 loudly. Pin that the
    guard is reachable and names the contig, so the guard cannot rot unnoticed."""
    from magestic.pipelines.wgs.qc import pipeline_qc

    src = Path(pipeline_qc.__file__).read_text()
    assert "split('_')" in src or 'split("_")' in src
    assert "integrated_plasmid" in src


# =========================================================================== #
# PLANTED-TRUTH CONFUSION — the numbers behind Kevin's requested rates         #
# =========================================================================== #
@pytest.mark.parametrize(
    "label,vaf,cov,expect_site,expect_cell",
    [
        ("clean_edit", 1.00, 30, "full", SINGLE_EDIT),
        ("mixed_50", 0.50, 30, "partial", CONTAMINATION),
        ("mixed_30", 0.30, 30, "partial", CONTAMINATION),
        # 🔴 the artifact: 0.10 <= tau_wt (0.15), so a 10%-mixed well is WT -> NO_EDIT
        ("mixed_10", 0.10, 30, SITE_WT, NO_EDIT),
        ("no_edit", 0.00, 30, SITE_WT, NO_EDIT),
        ("cov_1x", 1.00, 1, SITE_NOCOV, LOW_CONF),
        ("cov_3x", 1.00, 3, SITE_NOCOV, LOW_CONF),
        ("cov_5x", 1.00, 5, "full", SINGLE_EDIT),
        ("zero_cov", 0.00, 0, SITE_NOCOV, LOW_CONF),
    ],
)
def test_planted_truth_maps_to_the_measured_call(label, vaf, cov, expect_site, expect_cell):
    s = _site(vaf, cov)
    assert s.status == expect_site, label
    assert _cell([s]).label == expect_cell, label


def test_a_ten_percent_mixed_well_is_reported_as_unedited():
    """🔴 THE FINDING. tau_wt = 0.15, so any well whose true edit fraction is at or
    below 15% is classified SITE_WT and reconciled as NO_EDIT -- byte-identical to a
    genuinely unedited clone. A mixed well is not an unedited well, but this
    pipeline cannot tell them apart below 15%.

    Consequence: the real cohort's 22.6% "no edit" is an UPPER bound on genuinely
    unedited clones; it also absorbs every low-fraction mixture.
    """
    mixed = _cell([_site(0.10, 30)])
    truly_unedited = _cell([_site(0.00, 30)])
    assert mixed.label == NO_EDIT
    assert truly_unedited.label == NO_EDIT
    assert mixed.label == truly_unedited.label
    # and nothing in the record flags the mixture
    assert "PARTIAL_VAF" not in mixed.flags


def test_a_double_edit_is_recovered_as_multi_edit():
    """Two FULL designed sites -> MULTI_EDIT. The double-edit rate IS supportable."""
    assert _cell([_site(1.0, 30), _site(1.0, 30, start=1500, end=1510)]).label == "MULTI_EDIT"


def test_a_double_edit_with_one_mixed_site_is_hidden_by_contamination():
    """Decision order puts PARTIAL above MULTI_EDIT, so a double edit in which one
    site is mixed is reported as CONTAMINATION, not as a double edit. Pins the
    precedence, which is a real confusion channel for the double-edit rate."""
    cell = _cell([_site(1.0, 30), _site(0.5, 30, start=1500, end=1510)])
    assert cell.label == CONTAMINATION
    assert "PARTIAL_VAF" in cell.flags


def test_reconcile_cannot_distinguish_an_indel_from_a_designed_snv_on_real_bams(tmp_path):
    """REPRODUCTION of defect 242-4 through the real pileup, not by comparing two
    identical objects.

    Two BAMs are built over the SAME donor window: one carrying a designed SNV at
    100%, one carrying a 2-base cut-site insertion at 100%. `window_edit_vaf`
    measures only the MAX NON-REFERENCE fraction, so both return VAF 1.0 and both
    reconcile to SINGLE_EDIT. The lesion TYPE is not recoverable at this stage, so
    an indel is reported as a correct edit.
    """
    import importlib.util as _ilu

    sim = Path(
        "/Users/kevinroy/work/Stanford/SGTC_genome_editing_group/planning/"
        "222_triumph_embl_library_design/242_known_answer_sim.py"
    )
    if not sim.exists():  # pragma: no cover - the planning tree is not always present
        pytest.skip("242_known_answer_sim.py (the BAM builder) not available")
    spec_ = _ilu.spec_from_file_location("nn242sim", sim)
    m = _ilu.module_from_spec(spec_)
    spec_.loader.exec_module(m)

    fa = m.build_reference(tmp_path)
    labels = {}
    for label, mode in (("snv", "snv"), ("indel", "indel_ins")):
        bam_p = tmp_path / f"{label}.bam"
        m.write_well_bam(fa, bam_p, depth=30, vaf=1.0, mode=mode)
        with pysam.FastaFile(str(fa)) as fasta, pysam.AlignmentFile(str(bam_p), "rb") as bam:
            vaf, cov = m.window_edit_vaf(bam, fasta, "chrXV", m.WIN_START, m.WIN_END)
        labels[label] = (round(vaf, 3), cov, _cell([_site(vaf, cov)]).label)

    assert labels["snv"][0] == 1.0 and labels["indel"][0] == 1.0
    # 🔴 Identical call for two different lesions. This is the defect.
    assert labels["snv"][2] == SINGLE_EDIT
    assert labels["indel"][2] == SINGLE_EDIT
    assert labels["snv"][2] == labels["indel"][2]


# =========================================================================== #
# STAGE 2 — the stage that ACTUALLY produced the real cohort's 22.6%           #
# =========================================================================== #
# In `235_per_well_table.tsv`, editing_outcome_final == "NE" is a 115/115 EXACT
# row-wise match to editing_outcome == "designed_variant_unedited", emitted only at
# donor_outcome.py:354. reconcile_cell never produced those labels.
_WT21 = "ACGTACGTACGTACGTACGTA"
_DESIGN21 = _WT21[:10] + "T" + _WT21[11:]
_DS, _DE = 1000, 1020
_EDIT = 1010
_CHROMSEQ = "N" * (_DS - 1) + _WT21 + "N" * 200


def _s2(chrom="chrXV"):
    return DonorSpec("S", chrom, _DS, _DE, _WT21, _DESIGN21, _DESIGN21)


def _s2_vcf(path: Path, records, *, chrom="chrXV", ref_blocks=None, end_block=None, gt="0/1"):
    # Block REF from `_CHROMSEQ`: the pinned TRACE validates every gVCF REF (planning/268).
    # The callers' variant records were written C>T at 1010, where `_WT21` carries G (the
    # design is G>T); TRACE `7c4c10f` never checked. They now state the real REF base.
    out = ["##fileformat=VCFv4.2", f"##contig=<ID={chrom}>",
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1"]
    if end_block:
        out.append(f"{chrom}\t{end_block[0]}\t.\t{_CHROMSEQ[end_block[0] - 1]}\t<NON_REF>\t.\t.\t"
                   f"END={end_block[1]}\tGT:AD\t0/0:30,0")
    for p, adr in ref_blocks or []:
        out.append(f"{chrom}\t{p}\t.\t{_CHROMSEQ[p - 1]}\t<NON_REF>\t.\t.\t.\tGT:AD\t0/0:{adr},0")
    for pos, ref, alt, adr, ada in records:
        out.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\tGT:AD\t{gt}:{adr},{ada}")
    path.write_text("\n".join(out) + "\n")
    return path


_FULLCOV = [(p, 30) for p in range(_DS, _DE + 1)]


@pytest.mark.parametrize("ad_ref,ad_alt", [(0, 30), (15, 15), (21, 9), (27, 3), (29, 1)])
def test_stage2_has_no_allele_fraction_threshold_at_all(tmp_path, ad_ref, ad_alt):
    """🔴 THE STAGE-2 FINDING, as it stands on the pinned TRACE (planning/268).

    TRACE `7c4c10f` bucket-sorted purely on the REF/ALT strings, so a well with ONE alt
    read in thirty (3%) was `designed_variant_installed`, byte-identical to a clean
    30/30 clone. WITHDRAWN: "a 3% heterozygous call ships as installed". The pinned
    TRACE selects alleles by GT: a heterozygous genotype in a clonal donor is refused as
    `unknown;ambiguous_genotype` at every fraction, and a homozygous-alt genotype is
    installed at every fraction. What survives, and is pinned: Stage 2 still has NO
    allele-fraction threshold of its own. AD never changes the call, so the upstream
    caller's genotype is the controlling threshold (`268_rescue_probe_results.tsv`).
    """
    g = _s2_vcf(tmp_path / "g.vcf", [], ref_blocks=_FULLCOV)
    het = _s2_vcf(tmp_path / "het.vcf", [(_EDIT, "G", "T", ad_ref, ad_alt)])
    hom = _s2_vcf(tmp_path / "hom.vcf", [(_EDIT, "G", "T", ad_ref, ad_alt)], gt="1")
    assert precall_sample(_s2(), het, g, _CHROMSEQ).editing_outcome == "unknown;ambiguous_genotype"
    assert precall_sample(_s2(), hom, g, _CHROMSEQ).editing_outcome == "designed_variant_installed"


def test_stage2_calls_unedited_at_gvcf_reference_depth_one(tmp_path):
    """🔴 `ne_threshold` defaults to 1, so a single reference read per position is
    enough to call a well confidently `unedited`. Stage 2 applies NO depth floor of
    its own; the 5x floor in the published figure was applied afterwards by
    `235_plates_2348_2349_early_pam.py` on `var_cov`. 17 of the 115 raw NE rows had
    var_cov < 5 and were removed by that floor (98 + 17 = 115)."""
    v = _s2_vcf(tmp_path / "v.vcf", [])
    for depth, expect in ((0, "unknown"), (1, "designed_variant_unedited"),
                          (3, "designed_variant_unedited")):
        g = _s2_vcf(tmp_path / f"g{depth}.vcf", [], ref_blocks=[(p, depth) for p in range(_DS, _DE + 1)])
        assert precall_sample(_s2(), v, g, _CHROMSEQ).editing_outcome == expect, depth


def test_stage2_separates_an_indel_and_an_imprecise_repair_from_a_correct_edit(tmp_path):
    """Stage 2 DOES discriminate lesion type, unlike reconcile: an imprecise repair
    and a cut-site indel both come back `designed_variant_unedited;spontaneous_mutations`,
    carrying a suffix that distinguishes them from a bare `designed_variant_unedited`.

    This is why imprecise repairs and indels are NOT hidden inside the real 22.6%:
    all 115 raw NE rows are bare `designed_variant_unedited` with no suffix."""
    g = _s2_vcf(tmp_path / "g.vcf", [], ref_blocks=_FULLCOV)
    # A clonal haploid lesion is GT `1`; the pinned TRACE refuses a heterozygous
    # genotype in a clonal donor as ambiguous (planning/268), so the fixture states it.
    for label, rec in (
        ("imprecise", (_EDIT, "G", "A", 0, 30)),
        ("insertion", (_EDIT, "G", "GAA", 0, 30)),
        ("deletion", (_EDIT - 1, "CG", "C", 0, 30)),
    ):
        v = _s2_vcf(tmp_path / f"{label}.vcf", [rec], gt="1")
        out = precall_sample(_s2(), v, g, _CHROMSEQ)
        assert out.editing_outcome == "designed_variant_unedited;spontaneous_mutations", label
        assert len(out.spontaneous_variants) == 1, label
    # the bare form, for contrast -- what the 115 real NE rows carry
    v0 = _s2_vcf(tmp_path / "none.vcf", [])
    assert precall_sample(_s2(), v0, g, _CHROMSEQ).editing_outcome == "designed_variant_unedited"


def test_stage2_misses_an_END_compressed_gvcf_block(tmp_path):
    """🔴 NEW DEFECT 242-5, not on the brief's list. `flag_covered` counts gVCF
    RECORDS by POS (donor_outcome.py:297) and ignores the `END` INFO key. Both
    `bcftools --gvcf` and GATK emit COMPRESSED reference blocks -- one record spanning
    many positions -- so a gVCF produced the normal way yields flag_covered == 1 and
    every genuinely unedited well comes back `unknown`.

    So GATE 2 is worse than documented: supplying a gVCF is not sufficient, it must
    be an UNCOMPRESSED per-position gVCF."""
    import math

    # The quantification window is CENTRED on the donor midpoint, not anchored at
    # donor_start: qwin_start = ceil(dend - (dend-dstart)/2 - qws/2) = 1005, qwin_end = 1014.
    qs = math.ceil(_DE - (_DE - _DS) / 2 - 10 / 2)
    qe = qs + 9
    assert (qs, qe) == (1005, 1014)

    v = _s2_vcf(tmp_path / "v.vcf", [])
    # Per-position records over EXACTLY the quantification window -> callable.
    g_pp = _s2_vcf(tmp_path / "g_pp.vcf", [], ref_blocks=[(p, 30) for p in range(qs, qe + 1)])
    assert precall_sample(_s2(), v, g_pp, _CHROMSEQ).editing_outcome == "designed_variant_unedited"

    # The SAME coverage expressed as ONE END-compressed block, starting INSIDE the
    # quantification window so there is no question of it being out of range -> `unknown`.
    g_end_in = _s2_vcf(tmp_path / "g_end_in.vcf", [], end_block=(qs, qe))
    assert precall_sample(_s2(), v, g_end_in, _CHROMSEQ).editing_outcome == "unknown"

    # And a block spanning the whole donor window, likewise.
    g_end_wide = _s2_vcf(tmp_path / "g_end_wide.vcf", [], end_block=(_DS, _DE))
    assert precall_sample(_s2(), v, g_end_wide, _CHROMSEQ).editing_outcome == "unknown"


def test_the_upstream_caller_is_recorded_and_is_variants_only_and_now_haploid():
    """Pins the recovered upstream command, which is the ACTUAL controlling threshold
    for where a mixed well lands (Stage 2 has none of its own -- see
    `test_stage2_has_no_allele_fraction_threshold_at_all`).

    `core/vcf_adapter.py` records it as `MPILEUP_CALL_COMMAND`:
      `bcftools mpileup -q 20 -Q 20 -d 500 -a FORMAT/AD,FORMAT/DP | bcftools call -mv`

    `-v` is VARIANTS ONLY, so no reference blocks are emitted anywhere -> GATE 2 is
    real, and `unedited` is unobservable from these VCFs. **That half is unchanged.**

    🔴 **UPDATED by NN 251 (defect 1), which is what this test was waiting for.** The
    original asserted `"--ploidy" not in src`: the recovered command used bcftools'
    DIPLOID default on haploid yeast, so a 10-30% mixture was emitted with GT 0/0 and
    `carries_alt()` dropped it. The recorded command now carries `--ploidy 1`.

    Two things NN 251 established that this test now pins, because both are easy to
    over-read:
      * `--ploidy 1` makes the GENOTYPE correct for a haploid; it does NOT recover
        minor alleles (a haploid model chooses between GT 0 and GT 1, so a 33% alt is
        still expected to be GT 0). The recovery is the separate AD-derived
        `min_allele_fraction` path, which ignores GT entirely.
      * `carries_alt`'s semantics are DELIBERATELY unchanged, so every existing
        output is reproducible; the rescue is additive and opt-in.
    """
    from magestic.pipelines.wgs.core import vcf_adapter

    src = Path(vcf_adapter.__file__).read_text()
    assert "call -mv" in src
    assert "mpileup -q 20 -Q 20 -d 500" in src
    assert "--ploidy 1" in vcf_adapter.MPILEUP_CALL_COMMAND
    assert vcf_adapter.CALLER_PLOIDY == 1
    # The limitation is recorded in the module, not only in a planning doc.
    #
    # 🔴 UPDATED by NN 255. This line asserted "UNTESTED AGAINST THE REAL BINARY",
    # which was the caveat's own opening claim. NN 252 F13 measured that premise to be
    # FALSE — bcftools 1.21 / htslib 1.21 is installed and `MPILEUP_CALL_COMMAND` runs
    # clean — so the caveat is now a MEASURED statement and the test pins the
    # measurement instead of the apology. The limits it records are unchanged in
    # substance and stronger in detail; NN 255's floor sweep added a third that
    # inverts the naive advice, so the test pins that too.
    assert "MEASURED against bcftools 1.21" in vcf_adapter.PLOIDY_CAVEAT
    assert "UNTESTED" not in vcf_adapter.PLOIDY_CAVEAT, (
        "the caveat's premise was measured false; do not reinstate it"
    )
    assert "does not recover minor alleles" in vcf_adapter.PLOIDY_CAVEAT
    # The third, measured limit: `-v` is what sets the private-mixture floor, and
    # `--ploidy 1` WITH `-v` makes that floor WORSE than the diploid default.
    assert "34/40" in vcf_adapter.PLOIDY_CAVEAT and "6/40" in vcf_adapter.PLOIDY_CAVEAT
    assert "--ploidy 1 -A" in vcf_adapter.MPILEUP_CALL_ALL_SITES_COMMAND
    assert "-v" not in vcf_adapter.MPILEUP_CALL_ALL_SITES_COMMAND.split("|")[1]
    # GT filtering is unchanged -- it is still what drops a 0/0 mixture by default.
    assert vcf_adapter.carries_alt("0/1:...") is True
    assert vcf_adapter.carries_alt("0/0:...") is False
    assert vcf_adapter.carries_alt("./.:...") is False
    # ...and the AD-derived fraction, which is the actual evidence, ignores GT.
    assert vcf_adapter.allele_fraction(["GT", "AD"], ["0/0", "20,10"]) == pytest.approx(
        1 / 3
    )
    assert vcf_adapter.allele_fraction(["GT", "AD"], ["0/0", "."]) is None
    assert vcf_adapter.allele_fraction(["GT", "DP"], ["0/0", "30"]) is None


@pytest.mark.parametrize(
    "chrom",
    ["plasmid_o1", "integrated_plasmid_o1", "native_2micron",
     "blank_MAGESTIC3_pF835_guide_donor_plasmid", "unplaced_contig_34"],
)
def test_gate3_is_closed_on_every_mandated_composite_contig(tmp_path, chrom):
    """GATE 3 across the full set of underscore-bearing contig names that
    `core/contig_alias.py` itself lists as mandated by the reference standard.

    planning/268: WITHDRAWN "every mandated composite contig crashes Stage 2". Under the
    pinned TRACE each one gives exactly the outcome the underscore-free `chrXV` gives
    for the same evidence, with its real name in every shipped id."""
    def run(name):
        v = _s2_vcf(tmp_path / f"{name}.v.vcf", [(_EDIT, "G", "T", 0, 30)], chrom=name)
        g = _s2_vcf(tmp_path / f"{name}.g.vcf", [], chrom=name,
                    ref_blocks=[(p, 30) for p in range(_DS, _DE + 1)])
        return precall_sample(_s2(name), v, g, _CHROMSEQ)

    control, out = run("chrXV"), run(chrom)
    assert out.editing_outcome == control.editing_outcome
    ids = list(out.designed_variants) + list(out.spontaneous_variants)
    assert all(v.startswith(chrom + "_") for v in ids), ids


def test_stage_2_needs_no_contig_alias_now_that_gate_3_is_closed():
    """🔴 **UPDATED by planning/268 -- this test's tripwire fired as designed, again.**

    NN 251 wired `core/contig_alias.py` into step_2 (encode every name TRACE reads,
    decode the ids before publishing) because TRACE `7c4c10f` crashed on underscores.
    The pinned TRACE builds ids from structured fields
    (`test_gate3_is_closed_on_every_mandated_composite_contig`), so the 205 lane dropped
    GATE 3 and the merge retired the alias wiring. Re-aliasing would rewrite per-sample
    VCFs after their producer receipts were written. WITHDRAWN: "the mitigation is wired
    into Stage 2". Pinned instead: step_2 neither aliases nor gates, and the module
    survives unchanged for anyone who needs the bijection.
    """
    from magestic.pipelines.wgs.core import contig_alias
    from magestic.pipelines.wgs.pipeline import step_2_donor_precall as s2

    src = Path(s2.__file__).read_text()
    assert "core.contig_alias import" not in src and "alias_vcf_contigs(" not in src
    assert not hasattr(s2, "GATE_3") and not hasattr(s2, "_gate_3_alias")
    assert hasattr(contig_alias, "ContigAlias")
