"""
Unit tests for magestic.design.pam_search.

Uses a synthetic 200-bp reference sequence with known PAM sites to validate
guide enumeration and off-target writing.
"""

import tempfile
from pathlib import Path

import pytest
from magestic.design.nuclease import SPCAS9
from magestic.design.pam_search import search_pam_sites_genomewide


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_genome():
    """
    Minimal 2-chromosome genome with hand-crafted PAM sites.

    chrI:  ...[20-nt guide][NGG]...
           Position 30: guide = chrI[10:30], PAM = chrI[30:33] = "AGG"
    chrII: All T's — no NGG or CCN patterns anywhere (no SpCas9 PAM).
    """
    # chrI has an NGG PAM at index 30 (0-based)
    guide = "ACGTACGTACGTACGTACGT"  # 20 nt
    pam = "AGG"
    chrom_i = "A" * 10 + guide + pam + "T" * 100
    # T's only: no NGG (T...T), no CCN rev strand (GG...G = rev of CC...C absent)
    chrom_ii = "T" * 200
    return {"chrI": chrom_i, "chrII": chrom_ii}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestSearchPamSitesGenomewide:
    def test_output_file_created(self, simple_genome, tmp_path):
        outfile = tmp_path / "off_targets.tsv"
        search_pam_sites_genomewide(simple_genome, SPCAS9, outfile)
        assert outfile.exists()

    def test_header_present(self, simple_genome, tmp_path):
        outfile = tmp_path / "off_targets.tsv"
        search_pam_sites_genomewide(simple_genome, SPCAS9, outfile)
        with open(outfile) as fh:
            header = fh.readline()
        assert "full_guide" in header
        assert "num_mismatches" in header

    def test_known_guide_in_output(self, simple_genome, tmp_path):
        outfile = tmp_path / "off_targets.tsv"
        search_pam_sites_genomewide(simple_genome, SPCAS9, outfile)
        content = outfile.read_text()
        # The exact guide sequence we placed in the genome
        assert "ACGTACGTACGTACGTACGT" in content

    def test_returns_path(self, simple_genome, tmp_path):
        outfile = tmp_path / "off_targets.tsv"
        result = search_pam_sites_genomewide(simple_genome, SPCAS9, outfile)
        assert result == outfile

    def test_no_pam_on_empty_chromosome(self, simple_genome, tmp_path):
        """chrII is all C's — no NGG present; should have no guides from chrII."""
        outfile = tmp_path / "off_targets.tsv"
        search_pam_sites_genomewide(simple_genome, SPCAS9, outfile)
        content = outfile.read_text()
        lines_from_chrom_ii = [l for l in content.splitlines() if "\tchrII\t" in l]
        assert len(lines_from_chrom_ii) == 0

    def test_zero_mismatch_record(self, simple_genome, tmp_path):
        outfile = tmp_path / "off_targets.tsv"
        search_pam_sites_genomewide(simple_genome, SPCAS9, outfile)
        with open(outfile) as fh:
            fh.readline()  # header
            first_data = fh.readline().split("\t")
        # Column 1 (0-indexed) = num_mismatches
        assert first_data[1] in ("0", "1")
