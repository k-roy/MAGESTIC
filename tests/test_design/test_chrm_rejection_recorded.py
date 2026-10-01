"""
NN 234 bug 3 -- a dropped contig must be counted and reported, never silent.

``vcf_processing.load_haplotype_vcf`` dropped non-nuclear contigs with a bare
``continue``::

    if "M" in chrom or chrom not in YEAST_CHROMS:
        continue

No counter, no reason, nothing in any output -- mitochondrial variants simply
disappeared.  That matters for the EMBL panel specifically, which is respiration and
central carbon metabolism, so the mitochondrial genome is squarely in scope.

The contract these tests pin down: every dropped row is counted under an explicit
reason, and a caller can opt chrM in rather than edit the module.
"""

import logging
from pathlib import Path

import pytest

from magestic.design import vcf_processing
from magestic.design.vcf_processing import load_haplotype_vcf

HEADER = "\t".join([
    "CHROM", "POS", "REF", "ALT", "SAMPLE", "HAPLOTYPE", "PID", "genotypes",
    "all_ALT_alleles", "ORF", "NAME",
]) + "\n"


def _row(chrom, pos, ref, alt, sample="S1", hap="1", orf="YBR001C", name="GENE1"):
    return "\t".join([
        chrom, str(pos), ref, alt, sample, hap, "PID1", "0|1", alt, orf, name,
    ]) + "\n"


#: chrI is nuclear and kept. The other three are the silently-dropped classes:
#: the mitochondrial genome under both spellings, and the 2-micron plasmid.
TSV_BODY = (
    _row("chrI", 1000, "A", "G")
    + _row("chrM", 2000, "C", "T", orf="Q0045", name="COX1")
    + _row("chrmt", 2500, "G", "A", orf="Q0105", name="COB")
    + _row("2-micron", 500, "T", "C", orf="REP1", name="REP1")
)


def _write_tsv(tmp_path, body=TSV_BODY) -> Path:
    p = Path(tmp_path) / "haplotype.tsv"
    p.write_text(HEADER + body)
    return p


class _LogCapture(logging.Handler):
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


def _n_variants(variants_by_sample) -> int:
    return sum(
        len(pos_lst)
        for s in variants_by_sample.values()
        for h in s.values()
        for c in h.values()
        for pos_lst in c.values()
    )


# ---------------------------------------------------------------------------
# The drop must be counted
# ---------------------------------------------------------------------------

def test_dropped_contigs_are_counted_by_reason(tmp_path):
    """load_haplotype_vcf must hand back a count of what it dropped and why."""
    skip_counts: dict = {}
    result = load_haplotype_vcf(_write_tsv(tmp_path), skip_counts=skip_counts)

    # The nuclear variant still loads.
    assert _n_variants(result) == 1

    assert skip_counts, (
        "three rows (chrM, chrmt, 2-micron) were dropped and nothing recorded it; "
        "silent loss is the failure class strain_background.py exists to prevent"
    )
    total_skipped = sum(v for v in skip_counts.values() if isinstance(v, int))
    assert total_skipped == 3, (
        f"expected 3 dropped rows accounted for, got {total_skipped}: {skip_counts}"
    )
    # The reason must name the cause, not just be a bare total.
    assert any(
        "contig" in str(k).lower() or "chrom" in str(k).lower() for k in skip_counts
    ), f"the reason key must name the cause; got keys {list(skip_counts)}"


def test_dropped_contigs_are_reported_in_the_log(tmp_path):
    """The drop must also be visible without the caller opting in to counting."""
    with _LogCapture(vcf_processing.__name__) as cap:
        load_haplotype_vcf(_write_tsv(tmp_path))

    assert cap.records, "dropping 3 of 4 rows must not be silent"
    # Look only at messages that are actually about skipped rows -- the "Loaded N
    # variants from <path>" line can carry a stray digit from the temp path.
    skip_msgs = [
        r.getMessage() for r in cap.records
        if any(w in r.getMessage().lower() for w in ("skip", "drop", "reject"))
    ]
    assert skip_msgs, (
        f"no message reported the dropped rows; got {[r.getMessage() for r in cap.records]!r}"
    )
    assert any("3" in m for m in skip_msgs), (
        f"expected the dropped count (3) in the report; got {skip_msgs!r}"
    )


# ---------------------------------------------------------------------------
# chrM must be *processable*, not only rejectable
# ---------------------------------------------------------------------------

def test_chrm_can_be_explicitly_included(tmp_path):
    """A caller must be able to opt chrM in without editing the module."""
    skip_counts: dict = {}
    result = load_haplotype_vcf(
        _write_tsv(tmp_path),
        allowed_chroms=list(vcf_processing.YEAST_CHROMS) + ["chrM", "chrmt"],
        skip_counts=skip_counts,
    )

    assert _n_variants(result) == 3, (
        f"chrI + chrM + chrmt should load when chrM is allowed; got "
        f"{_n_variants(result)} ({result})"
    )
    chroms_loaded = {c for s in result.values() for h in s.values() for c in h}
    assert "chrM" in chroms_loaded, f"chrM missing from {chroms_loaded}"
    assert "chrmt" in chroms_loaded, f"chrmt missing from {chroms_loaded}"
    # The 2-micron is still outside the allowed set, and is still accounted for.
    assert sum(v for v in skip_counts.values() if isinstance(v, int)) == 1


# ---------------------------------------------------------------------------
# Positive controls
# ---------------------------------------------------------------------------

def test_nuclear_only_file_loads_with_no_drops(tmp_path):
    """A nuclear-only file must load fully and record no drops."""
    body = _row("chrI", 1000, "A", "G") + _row("chrII", 4000, "C", "T")
    skip_counts: dict = {}
    result = load_haplotype_vcf(_write_tsv(tmp_path, body), skip_counts=skip_counts)

    assert _n_variants(result) == 2
    assert sum(v for v in skip_counts.values() if isinstance(v, int)) == 0, (
        f"a nuclear-only file must record no drops; got {skip_counts}"
    )


def test_short_and_sv_rows_are_also_counted(tmp_path):
    """The other two silent drops in the same loop must be accounted for too."""
    body = (
        _row("chrI", 1000, "A", "G")
        + _row("chrI", 1500, "A", "*")      # structural variant
        + "chrI\t1600\tA\tG\n"               # short row: fewer than 11 columns
    )
    skip_counts: dict = {}
    result = load_haplotype_vcf(_write_tsv(tmp_path, body), skip_counts=skip_counts)

    assert _n_variants(result) == 1
    assert sum(v for v in skip_counts.values() if isinstance(v, int)) == 2, (
        f"expected the SV row and the short row to be counted; got {skip_counts}"
    )
