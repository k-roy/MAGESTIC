"""
NN 234 bug 2 -- ``--minimum-homology`` must be enforced, or the design rejected.

``oligo_generation.construct_donor`` sized the homology arms as::

    left_arm  = (donor_length - len(alt)) // 2
    right_arm = donor_length - len(alt) - left_arm

without consulting ``minimum_homology``.  For an insertion as long as the donor both
arms are **0** and the emitted "donor" is the bare insertion, which cannot recombine.
For an insertion longer than the donor the arms go negative and the genome slices come
back empty, giving the same result.  The sibling EMBL design run emitted 4 oligos at
207 nt with zero homology arms this way.

The contract these tests pin down: a design that cannot satisfy the homology floor is
REJECTED with a recorded reason, never emitted.
"""

import pytest

from magestic.design.oligo_generation import construct_donor

CHROM = "chrI"
CHROM_LENGTH = 6000
GUIDE_LENGTH = 20
DONOR_LENGTH = 207
MINIMUM_HOMOLOGY = 30
VARIANT_COORD = 3000  # comfortably mid-chromosome: no boundary clipping

# A genome with no BspQI/SapI site, so restriction-site handling never fires and the
# only thing under test is the homology-arm arithmetic.
GENOME_SEQ = {CHROM: "ACGT" * (CHROM_LENGTH // 4)}
PLUS_BG = {CHROM: {}}
MINUS_BG = {CHROM: {}}


def _construct(alt: str, *, donor_length: int = DONOR_LENGTH, ref: str = "A"):
    return construct_donor(
        CHROM,
        VARIANT_COORD,
        ref,
        alt,
        CHROM_LENGTH,
        [[str(VARIANT_COORD), ref, alt, "INDEL"]],
        GUIDE_LENGTH,
        donor_length,
        MINIMUM_HOMOLOGY,
        GENOME_SEQ,
        PLUS_BG,
        MINUS_BG,
    )


def _is_rejection(result) -> bool:
    return result is not None and result[0] is None


def _rejection_reason(result) -> str:
    """The reason a rejection carries. Empty string if it carries none."""
    if len(result) >= 3 and isinstance(result[2], str):
        return result[2]
    return ""


# ---------------------------------------------------------------------------
# The library-killing case: an insertion exactly as long as the donor
# ---------------------------------------------------------------------------

def test_insertion_as_long_as_donor_is_rejected_with_a_reason():
    """An insertion that leaves no room for arms must be rejected, not emitted."""
    alt = "AT" * 103 + "A"  # 207 nt, same as DONOR_LENGTH, no RE site
    assert len(alt) == DONOR_LENGTH

    result = _construct(alt)

    assert _is_rejection(result), (
        "construct_donor emitted a design for a 207-nt insertion into a 207-nt "
        f"donor. Returned left arm {result[4]!r} ({len(result[4])} nt) and right arm "
        f"{result[5]!r} ({len(result[5])} nt); donor is {len(result[8])} nt of bare "
        "insertion with ZERO homology and cannot recombine."
    )
    reason = _rejection_reason(result)
    assert "homology" in reason.lower(), (
        f"rejection must record a homology reason; got {reason!r}"
    )


def test_insertion_longer_than_donor_is_rejected_with_a_reason():
    """len(alt) > donor_length drives both arms negative -- also a rejection."""
    alt = "AT" * 130  # 260 nt > DONOR_LENGTH
    assert len(alt) > DONOR_LENGTH

    result = _construct(alt)

    assert _is_rejection(result), (
        f"construct_donor emitted a design for a {len(alt)}-nt insertion into a "
        f"{DONOR_LENGTH}-nt donor (arms computed negative)."
    )
    assert "homology" in _rejection_reason(result).lower()


def test_insertion_one_base_over_the_floor_is_rejected():
    """The boundary: an insertion leaving arms of minimum_homology - 1 is rejected."""
    # donor_length - len(alt) == 2 * MINIMUM_HOMOLOGY - 1  ->  one arm is 29 nt
    alt = "AT" * ((DONOR_LENGTH - (2 * MINIMUM_HOMOLOGY - 1)) // 2)
    result = _construct(alt)
    assert _is_rejection(result), (
        f"a {len(alt)}-nt insertion leaves arms below the {MINIMUM_HOMOLOGY}-nt floor "
        "and must be rejected"
    )


# ---------------------------------------------------------------------------
# Positive controls: the guard must not reject designs that DO satisfy the floor
# ---------------------------------------------------------------------------

def test_snp_still_produces_a_donor_with_both_arms_above_the_floor():
    """A plain SNP must still be accepted, with both arms at or above the floor."""
    result = _construct("G", ref="A")

    assert not _is_rejection(result), f"a plain SNP must not be rejected; got {result!r}"
    left_arm, right_arm = result[4], result[5]
    assert len(left_arm) >= MINIMUM_HOMOLOGY, f"left arm {len(left_arm)} nt"
    assert len(right_arm) >= MINIMUM_HOMOLOGY, f"right arm {len(right_arm)} nt"
    assert len(result[8]) == DONOR_LENGTH


def test_insertion_exactly_at_the_floor_is_accepted():
    """An insertion leaving exactly minimum_homology on each side is accepted."""
    alt = "AT" * ((DONOR_LENGTH - 2 * MINIMUM_HOMOLOGY) // 2)
    result = _construct(alt)

    assert not _is_rejection(result), (
        f"a {len(alt)}-nt insertion leaves arms of exactly {MINIMUM_HOMOLOGY} nt and "
        f"must be accepted; got {result!r}"
    )
    assert len(result[4]) >= MINIMUM_HOMOLOGY
    assert len(result[5]) >= MINIMUM_HOMOLOGY


# ---------------------------------------------------------------------------
# No emitted member may ever carry an under-length arm
# ---------------------------------------------------------------------------

def test_no_accepted_donor_ever_has_an_arm_below_the_floor():
    """Sweep insertion sizes: every ACCEPTED design satisfies the floor."""
    offenders = []
    for n in range(0, DONOR_LENGTH + 40, 7):
        alt = "AT" * (n // 2) if n else "G"
        result = _construct(alt)
        if _is_rejection(result):
            continue
        left_arm, right_arm = result[4], result[5]
        if len(left_arm) < MINIMUM_HOMOLOGY or len(right_arm) < MINIMUM_HOMOLOGY:
            offenders.append((len(alt), len(left_arm), len(right_arm), len(result[8])))

    assert not offenders, (
        "construct_donor emitted designs whose arms are below the "
        f"{MINIMUM_HOMOLOGY}-nt floor -- (len_alt, left, right, donor): {offenders}"
    )


def test_variant_near_chromosome_end_does_not_emit_a_short_arm():
    """Boundary clipping must not silently produce an under-length arm either."""
    near_end = CHROM_LENGTH - 20
    result = construct_donor(
        CHROM, near_end, "A", "G", CHROM_LENGTH,
        [[str(near_end), "A", "G", "SNP"]],
        GUIDE_LENGTH, DONOR_LENGTH, MINIMUM_HOMOLOGY,
        GENOME_SEQ, PLUS_BG, MINUS_BG,
    )
    if _is_rejection(result):
        assert _rejection_reason(result), "a rejection must carry a reason"
        return
    assert len(result[4]) >= MINIMUM_HOMOLOGY, f"left arm {len(result[4])} nt"
    assert len(result[5]) >= MINIMUM_HOMOLOGY, f"right arm {len(result[5])} nt"


# ---------------------------------------------------------------------------
# End-to-end smoke test on the SAME code path: the rejection must reach the
# untargetable output through design_oligos_from_vcf, carrying its own reason.
# ---------------------------------------------------------------------------

def test_design_oligos_from_vcf_records_the_homology_rejection(tmp_path):
    """The reason must survive the caller and land in the untargetable TSV."""
    from pathlib import Path

    from magestic.design.nuclease import get_nuclease
    from magestic.design.oligo_generation import design_oligos_from_vcf

    tmp_path = Path(tmp_path)
    nuclease = get_nuclease("SpG")
    rev_primers = ["A" * 20]
    donor_length = nuclease.donor_length_for_oligo(200, len(rev_primers[0]))

    # An insertion long enough that no 30-nt arm can fit in a 129-nt donor.
    long_alt = "AT" * ((donor_length - 2 * MINIMUM_HOMOLOGY + 2) // 2)
    assert len(long_alt) > donor_length - 2 * MINIMUM_HOMOLOGY

    vcf = tmp_path / "processed.tsv"
    vcf.write_text(
        "CHROM\tPOS\tREF\tALT\tvariant_description\tSAMPLES\tORF\tNAME\n"
        f"{CHROM}\t{VARIANT_COORD}\tA\t{long_alt}\t"
        f"{VARIANT_COORD}_A_{long_alt}_INDEL\tS1\tYBR001C\tGENE1\n"
    )

    guides_header = "chrom\tPAM_coord\tstrand\tfull_guide\tfull_PAM\n"
    allowed = tmp_path / "allowed.tsv"
    allowed.write_text(guides_header)
    disallowed = tmp_path / "disallowed.tsv"
    disallowed.write_text(guides_header)

    out_dir = tmp_path / "out"
    design_oligos_from_vcf(
        processed_vcf_file=vcf,
        genome_seq=GENOME_SEQ,
        nuclease=nuclease,
        allowed_guides_file=allowed,
        disallowed_guides_file=disallowed,
        plus_strand_bedgraph=PLUS_BG,
        minus_strand_bedgraph=MINUS_BG,
        rev_primer_sequences=rev_primers,
        output_oligo_pool_file=out_dir / "pool.tsv",
        output_info_file=out_dir / "info.tsv",
        output_untargetable_file=out_dir / "untargetable.tsv",
        oligo_length=200,
        minimum_homology=MINIMUM_HOMOLOGY,
        nuclease_label="SpG",
    )

    pool = (out_dir / "pool.tsv").read_text().splitlines()
    assert len(pool) == 1, f"no oligo may be emitted for this variant; got {pool}"

    untargetable = (out_dir / "untargetable.tsv").read_text().splitlines()
    assert len(untargetable) == 2, f"expected one rejection row; got {untargetable}"
    reason = untargetable[1].split("\t")[-1]
    assert "below_minimum_homology" in reason, (
        f"the untargetable row must carry the homology reason, not a mislabelled "
        f"restriction-site reason; got {reason!r}"
    )
