"""
NN 234 bug 4, residual defect -- junction_changes must name EVERY base where the
emitted donor differs from the DESIGNED donor, whichever code path changed it.

Measured on the three final EMBL ORDER pools (60,662 oligos): **203 terminal donor
base changes exist, 54 are recorded in junction_changes, 149 are not.**  Recording
every accepted edit inside ``remove_bspqi_sites_from_donor`` was necessary but not
sufficient, because that function is not the only route by which an emitted donor
comes to differ from the designed one.

The measured mechanism (54/54 variant groups, all three pools): ``donor`` is
reassigned in the enclosing scope of ``design_oligos_from_vcf``::

    new_donor, junction_changes = remove_bspqi_sites_from_donor(..., donor, ...)
    donor = new_donor          # persists into the NEXT guide of the SAME variant

so the second and later guides of one variant are handed an already-repaired donor.
That donor already carries exactly one restriction site, so the repair function
returns at its early exit with an empty list -- correctly, since *it* changed
nothing -- and the inherited change goes unrecorded.  In every one of the 54 groups
the first-emitted row is recorded and the rest are not, and all rows of a group
carry the identical change.

The fixture below is a REAL variant group from the shipped M22 ORDER pool,
``chrII:48676 G>T`` (ORF YBL090W), whose four oligos are rows 236-239 of
``final/M22/M22_ORDER_resp_ccm_info.tsv``.  Its designed donor begins ``TCTTCAAC``,
which together with the cloning site's trailing ``GAGC`` forms a spurious
``GCTCTTC``; the repair mutates donor index 0 from ``T`` to ``A``.  All four emitted
oligos carry that change.  Only the first records it.
"""

from pathlib import Path

import pytest

from magestic.design.nuclease import get_nuclease
from magestic.design.oligo_generation import design_oligos_from_vcf

# ---------------------------------------------------------------------------
# Real values, copied out of final/M22/M22_ORDER_resp_ccm_info.tsv
# ---------------------------------------------------------------------------

CHROM = "chrII"
VARIANT_COORD = 48676
REF = "G"
ALT = "T"
ORF = "YBL090W"
NAME = "YBL090W_id003"
SAMPLES = "M22"

DONOR_START_COORD = 48611
DONOR_END_COORD = 48740

#: ``genome[donor_start_coord : donor_end_coord + len(ref)]`` -- 130 nt.
GENOMIC_TARGET = (
    "CTTTTTTAAGTTCGCTTCTTTCCCAAGATTGAAAATAAAAAATAAGTGAGACCAGGTCATTAAAGG"
    "TCCCTTGAAGAGTAAAAAATCTGTCAAATATCACCCATAAGAGCATAACCATAATGTTGAAGAG"
)

#: The donor as STORED (this variant's donor_strand is '-'), designed vs emitted.
DESIGNED_DONOR = (
    "TCTTCAACATTATGGTTATGCTCTTATGGGTGATATTTGACAGATTTTTTACTCTTCAAGGGACA"
    "TTTAATGACCTGGTCTCACTTATTTTTTATTTTCAATCTTGGGAAAGAAGCGAACTTAAAAAAG"
)
EMITTED_DONOR = (
    "ACTTCAACATTATGGTTATGCTCTTATGGGTGATATTTGACAGATTTTTTACTCTTCAAGGGACA"
    "TTTAATGACCTGGTCTCACTTATTTTTTATTTTCAATCTTGGGAAAGAAGCGAACTTAAAAAAG"
)

#: (PAM_coord, PAM_strand, guide, PAM) for the four real oligos, in emitted order.
#: The first is the one whose junction_changes is non-empty in the shipped pool.
REAL_GUIDES = [
    (48685, "+", "GGTCATTAAAGGTCCCTTGA", "AGAG"),
    (48674, "+", "AAGTGAGACCAGGTCATTAA", "AGGT"),
    (48670, "-", "ACTCTTCAAGGGACCTTTAA", "TGAC"),
    (48682, "+", "CCAGGTCATTAAAGGTCCCT", "TGAA"),
]

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def test_fixture_matches_the_shipped_pool():
    """The fixture must be the real row, or the tests below prove nothing."""
    assert len(GENOMIC_TARGET) == (DONOR_END_COORD + len(REF)) - DONOR_START_COORD
    assert len(DESIGNED_DONOR) == len(EMITTED_DONOR) == 129

    # The designed donor really is the reference window with REF swapped for ALT.
    left = (VARIANT_COORD - 1) - DONOR_START_COORD
    assert left == 64
    plus = _rc(DESIGNED_DONOR)
    assert plus[:left] == GENOMIC_TARGET[:left]
    assert plus[left] == ALT
    assert plus[left + 1:] == GENOMIC_TARGET[left + len(REF):left + len(REF) + 64]

    # Exactly one difference, at stored index 0, T -> A.
    diffs = [i for i, (d, e) in enumerate(zip(DESIGNED_DONOR, EMITTED_DONOR)) if d != e]
    assert diffs == [0]
    assert (DESIGNED_DONOR[0], EMITTED_DONOR[0]) == ("T", "A")

    # And the change is real: the designed donor forms a spurious GCTCTTC against
    # the cloning site's trailing GAGC, while the emitted one does not.
    icl = get_nuclease("SpG").internal_cloning_site.upper()
    assert "GCTCTTC" in (icl + DESIGNED_DONOR).upper()
    assert "GCTCTTC" not in (icl + EMITTED_DONOR).upper()


# ---------------------------------------------------------------------------
# The invariant, end to end through the real design entry point
# ---------------------------------------------------------------------------

def _build_inputs(tmp_path: Path):
    """A chromosome carrying the real reference window at its real coordinates."""
    chrom_len = DONOR_END_COORD + len(REF) + 200
    seq = (
        "A" * DONOR_START_COORD
        + GENOMIC_TARGET
        + "A" * (chrom_len - DONOR_START_COORD - len(GENOMIC_TARGET))
    )
    genome_seq = {CHROM: seq}
    assert genome_seq[CHROM][DONOR_START_COORD:DONOR_END_COORD + len(REF)] == GENOMIC_TARGET

    # donor_strand is '-' for this variant, so minus-strand coverage must dominate.
    window = range(DONOR_START_COORD, DONOR_END_COORD)
    plus_bg = {CHROM: {}}
    minus_bg = {CHROM: {c: 1 for c in window}}

    vcf = tmp_path / "processed.tsv"
    vcf.write_text(
        "CHROM\tPOS\tREF\tALT\tvariant_description\tSAMPLES\tORF\tNAME\n"
        f"{CHROM}\t{VARIANT_COORD}\t{REF}\t{ALT}\t"
        f"{VARIANT_COORD}_{REF}_{ALT}_SNP\t{SAMPLES}\t{ORF}\t{NAME}\n"
    )

    allowed = tmp_path / "allowed.tsv"
    rows = ["chrom\tPAM_coord\tstrand\tfull_guide\tfull_PAM"]
    for pam_coord, pam_strand, guide, pam in REAL_GUIDES:
        rows.append(f"{CHROM}\t{pam_coord}\t{pam_strand}\t{guide}\t{pam}")
    allowed.write_text("\n".join(rows) + "\n")

    disallowed = tmp_path / "disallowed.tsv"
    disallowed.write_text("chrom\tPAM_coord\tstrand\tfull_guide\tfull_PAM\n")

    return genome_seq, plus_bg, minus_bg, vcf, allowed, disallowed


def _run_design(tmp_path: Path):
    genome_seq, plus_bg, minus_bg, vcf, allowed, disallowed = _build_inputs(tmp_path)
    out_dir = tmp_path / "out"
    design_oligos_from_vcf(
        processed_vcf_file=vcf,
        genome_seq=genome_seq,
        nuclease=get_nuclease("SpG"),
        allowed_guides_file=allowed,
        disallowed_guides_file=disallowed,
        plus_strand_bedgraph=plus_bg,
        minus_strand_bedgraph=minus_bg,
        rev_primer_sequences=["A" * 20],
        output_oligo_pool_file=out_dir / "pool.tsv",
        output_info_file=out_dir / "info.tsv",
        output_untargetable_file=out_dir / "untargetable.tsv",
        oligo_length=200,
        minimum_homology=30,
        max_guides_per_donor=4,
        nuclease_label="SpG",
    )
    import csv

    with open(out_dir / "info.tsv", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _designed_donor_for(row) -> str:
    """Rebuild the designed donor from the row's own reference window."""
    gt = row["genomic_target"]
    ref, alt = row["ref"], row["alt"]
    left = (int(row["variant_coord"]) - 1) - int(row["donor_start_coord"])
    stored = row["donor"]
    plus_len = len(stored)
    right = plus_len - left - len(alt)
    plus = gt[:left] + alt + gt[left + len(ref):left + len(ref) + right]
    return plus if row["donor_strand"] == "+" else _rc(plus)


def _positions_named(junction_changes_str: str, donor_len: int) -> set:
    import ast

    s = junction_changes_str.strip()
    recs = ast.literal_eval(s) if s else []
    named = set()
    for rec in recs:
        for item in rec:
            if isinstance(item, int):
                named.add(item % donor_len)
                break
    return named


def test_real_variant_group_reproduces_in_the_pipeline(tmp_path):
    """Sanity: the fixture really does yield several oligos sharing one donor."""
    rows = _run_design(tmp_path)
    assert len(rows) >= 2, (
        f"the fixture must emit several guides for one variant, or the leak across "
        f"guides cannot be exercised; got {len(rows)} rows"
    )
    assert len({r["donor"] for r in rows}) == 1, (
        "all guides of this variant share one donor -- that sharing is the defect"
    )
    assert rows[0]["donor"] == EMITTED_DONOR, (
        f"the pipeline must reproduce the shipped donor.\n"
        f"  got      {rows[0]['donor']}\n  expected {EMITTED_DONOR}"
    )


def test_every_changed_donor_base_is_named_in_junction_changes(tmp_path):
    """THE INVARIANT: no emitted oligo may carry an unrecorded donor change."""
    rows = _run_design(tmp_path)

    offenders = []
    for i, row in enumerate(rows):
        designed = _designed_donor_for(row)
        emitted = row["donor"]
        assert len(designed) == len(emitted)
        changed = {
            j for j, (d, e) in enumerate(zip(designed, emitted)) if d.upper() != e.upper()
        }
        named = _positions_named(row["junction_changes"], len(emitted))
        missing = changed - named
        if missing:
            offenders.append(
                (i, row["PAM_coord"], sorted(missing), row["junction_changes"])
            )

    assert not offenders, (
        f"{len(offenders)} of {len(rows)} emitted oligos carry a donor base that "
        f"differs from the designed donor but is NOT named in junction_changes. "
        f"(row_idx, PAM_coord, unrecorded positions, junction_changes): {offenders}"
    )


def test_recorded_change_names_the_true_designed_base(tmp_path):
    """A record must state the base that was designed, not an intermediate one."""
    import ast

    rows = _run_design(tmp_path)
    for row in rows:
        designed = _designed_donor_for(row)
        emitted = row["donor"]
        s = row["junction_changes"].strip()
        for rec in (ast.literal_eval(s) if s else []):
            rec = list(rec)
            idx = rec[rec.index("from") - 1]
            from_base = rec[rec.index("from") + 1]
            to_base = rec[rec.index("to") + 1]
            assert designed[idx].upper() == from_base.upper(), (
                f"record says index {idx} came from {from_base!r}, but the DESIGNED "
                f"donor has {designed[idx]!r} there"
            )
            assert emitted[idx].upper() == to_base.upper(), (
                f"record says index {idx} became {to_base!r}, but the emitted donor "
                f"has {emitted[idx]!r} there"
            )


def test_unchanged_donor_records_nothing(tmp_path):
    """Positive control: an oligo needing no repair must record no change."""
    rows = _run_design(tmp_path)
    for row in rows:
        designed = _designed_donor_for(row)
        if designed == row["donor"]:
            assert _positions_named(row["junction_changes"], len(designed)) == set(), (
                f"nothing changed, yet junction_changes says "
                f"{row['junction_changes']!r}"
            )
