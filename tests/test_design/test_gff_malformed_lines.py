"""
NN 234 bug 1 — a malformed GFF line must skip the record, not re-emit the previous one.

`GTF_GFF_manipulation.load_ORF_coordinates` used ``except: print(line)`` with no
``continue``, so after a parse failure the loop body fell through and re-appended the
*previous* feature's coordinates.  Every line of an embedded ``##FASTA`` block is such a
failure (those lines do not start with ``#``), which is how a collaborator's notebook
emitted ~152,000 duplicate GFF lines where 17 were expected.

These tests assert the invariant for all three implementations in the repo, plus that
the malformed-line count is *reported* rather than printed into stdout.
"""

import logging
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


# ---------------------------------------------------------------------------
# Fixture: a GFF with (a) two real genes, one of them two-exon, (b) a truncated
# line, and (c) an embedded FASTA block.
# ---------------------------------------------------------------------------

_N_FASTA_SEQ_LINES = 5

GFF_WITH_MALFORMED_LINE_AND_FASTA = (
    "##gff-version 3\n"
    "chrI\tSGD\tCDS\t335\t649\t.\t+\t0\t"
    "Parent=YAL069W_mRNA;Name=YAL069W_CDS;orf_classification=Dubious\n"
    "chrI\tSGD\tCDS\t1000\t1200\t.\t+\t0\t"
    "Parent=YAL068C_mRNA;Name=YAL068C_CDS;orf_classification=Verified\n"
    "chrI\tSGD\tCDS\t1300\t1500\t.\t+\t1\t"
    "Parent=YAL068C_mRNA;Name=YAL068C_CDS;orf_classification=Verified\n"
    # (b) truncated record: only 4 tab fields, so the 9-way unpack raises
    "chrI\tSGD\tCDS\t2000\n"
    # (c) embedded FASTA block -- the directive starts with '#' and is skipped,
    #     but the defline and every sequence line are parse failures
    "##FASTA\n"
    ">chrI\n"
    + ("ACGTACGTACGTACGTACGT\n" * _N_FASTA_SEQ_LINES)
)

#: The only two ORFs in the fixture, and how many CDS exons each really has.
EXPECTED_EXON_COUNTS = {"YAL069W": 1, "YAL068C": 2}


def _write_gff(tmp_path) -> Path:
    p = Path(tmp_path) / "malformed_with_fasta.gff"
    p.write_text(GFF_WITH_MALFORMED_LINE_AND_FASTA)
    return p


class _LogCapture(logging.Handler):
    """Collect records emitted on a named logger (works without pytest's caplog)."""

    def __init__(self, logger_name: str):
        super().__init__(level=logging.DEBUG)
        self._logger = logging.getLogger(logger_name)
        self.records: list = []

    def emit(self, record):  # noqa: D102
        self.records.append(record)

    def __enter__(self):
        self._prev_level = self._logger.level
        self._logger.setLevel(logging.DEBUG)
        self._logger.addHandler(self)
        return self

    def __exit__(self, *exc):
        self._logger.removeHandler(self)
        self._logger.setLevel(self._prev_level)
        return False

    @property
    def text(self) -> str:
        return "\n".join(r.getMessage() for r in self.records)


def _assert_each_feature_exactly_once(orf_dict, label: str) -> None:
    assert set(orf_dict) == set(EXPECTED_EXON_COUNTS), (
        f"{label}: ORF set is {sorted(orf_dict)}, expected "
        f"{sorted(EXPECTED_EXON_COUNTS)}"
    )
    for orf, expected in EXPECTED_EXON_COUNTS.items():
        got = orf_dict[orf]
        assert len(got) == expected, (
            f"{label}: {orf} has {len(got)} exon tuples, expected {expected}. "
            f"A malformed line re-emitted the previous feature. Got: {got}"
        )
        # Every tuple must be a distinct, correctly typed record.
        assert len(set(got)) == len(got), f"{label}: {orf} has duplicate tuples: {got}"
        for chrom, start, end, strand in got:
            assert isinstance(start, int) and isinstance(end, int), (
                f"{label}: {orf} has non-integer coordinates {(start, end)} -- a "
                f"failed int() conversion fell through into the record"
            )
            assert chrom == "chrI" and strand == "+"


# ---------------------------------------------------------------------------
# The loose root-level script -- the copy that was never fixed
# ---------------------------------------------------------------------------

def test_root_script_skips_malformed_gff_lines(tmp_path):
    """GTF_GFF_manipulation.load_ORF_coordinates must emit each feature once."""
    import GTF_GFF_manipulation

    orf_dict = GTF_GFF_manipulation.load_ORF_coordinates(str(_write_gff(tmp_path)))
    _assert_each_feature_exactly_once(orf_dict, "GTF_GFF_manipulation")


def test_root_script_reports_malformed_line_count(tmp_path):
    """The failure must be visible as a counted report, not printed into stdout."""
    import GTF_GFF_manipulation

    with _LogCapture("GTF_GFF_manipulation") as cap:
        GTF_GFF_manipulation.load_ORF_coordinates(str(_write_gff(tmp_path)))

    # 1 truncated record + 1 defline + N sequence lines
    n_malformed = 1 + 1 + _N_FASTA_SEQ_LINES
    assert cap.records, (
        "load_ORF_coordinates reported nothing for the malformed lines; a parse "
        "failure must be counted and reported, not silent (or printed to stdout)"
    )
    assert str(n_malformed) in cap.text, (
        f"expected the report to carry the malformed-line count {n_malformed}; "
        f"got: {cap.text!r}"
    )


# ---------------------------------------------------------------------------
# The two packaged copies
# ---------------------------------------------------------------------------

def test_design_package_skips_malformed_gff_lines(tmp_path):
    from magestic.design.guide_donor_functions import load_ORF_coordinates

    orf_dict = load_ORF_coordinates(str(_write_gff(tmp_path)))
    _assert_each_feature_exactly_once(orf_dict, "design.guide_donor_functions")


def test_design_package_reports_malformed_line_count(tmp_path):
    from magestic.design import guide_donor_functions

    with _LogCapture(guide_donor_functions.__name__) as cap:
        guide_donor_functions.load_ORF_coordinates(str(_write_gff(tmp_path)))

    n_malformed = 1 + 1 + _N_FASTA_SEQ_LINES
    assert str(n_malformed) in cap.text, (
        f"expected a counted report of {n_malformed} malformed lines; got: {cap.text!r}"
    )


def test_genomics_package_skips_malformed_gff_lines(tmp_path):
    from magestic.genomics.annotation import load_ORF_coordinates

    orf_dict = load_ORF_coordinates(str(_write_gff(tmp_path)))
    _assert_each_feature_exactly_once(orf_dict, "genomics.annotation")


def test_genomics_package_reports_malformed_line_count(tmp_path):
    from magestic.genomics import annotation

    with _LogCapture(annotation.__name__) as cap:
        annotation.load_ORF_coordinates(str(_write_gff(tmp_path)))

    n_malformed = 1 + 1 + _N_FASTA_SEQ_LINES
    assert str(n_malformed) in cap.text, (
        f"expected a counted report of {n_malformed} malformed lines; got: {cap.text!r}"
    )


# ---------------------------------------------------------------------------
# Positive control: a clean GFF must still load, and must report nothing
# ---------------------------------------------------------------------------

CLEAN_GFF = (
    "##gff-version 3\n"
    "chrI\tSGD\tCDS\t335\t649\t.\t+\t0\t"
    "Parent=YAL069W_mRNA;Name=YAL069W_CDS;orf_classification=Dubious\n"
    "chrI\tSGD\tCDS\t1000\t1200\t.\t+\t0\t"
    "Parent=YAL068C_mRNA;Name=YAL068C_CDS;orf_classification=Verified\n"
    "chrI\tSGD\tCDS\t1300\t1500\t.\t+\t1\t"
    "Parent=YAL068C_mRNA;Name=YAL068C_CDS;orf_classification=Verified\n"
)


def test_clean_gff_loads_and_reports_no_malformed_lines(tmp_path):
    """Positive control -- the guard must not reject well-formed records."""
    import GTF_GFF_manipulation
    from magestic.design.guide_donor_functions import (
        load_ORF_coordinates as design_loader,
    )
    from magestic.genomics.annotation import load_ORF_coordinates as genomics_loader

    p = Path(tmp_path) / "clean.gff"
    p.write_text(CLEAN_GFF)

    for label, loader in (
        ("GTF_GFF_manipulation", GTF_GFF_manipulation.load_ORF_coordinates),
        ("design.guide_donor_functions", design_loader),
        ("genomics.annotation", genomics_loader),
    ):
        _assert_each_feature_exactly_once(loader(str(p)), label)

    with _LogCapture("GTF_GFF_manipulation") as cap:
        GTF_GFF_manipulation.load_ORF_coordinates(str(p))
    assert "malformed" not in cap.text.lower(), (
        f"clean GFF must not report malformed lines; got: {cap.text!r}"
    )
