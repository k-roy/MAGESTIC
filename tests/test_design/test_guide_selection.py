"""
Unit tests for magestic.design.guide_selection.
"""

import pytest
from magestic.design.guide_selection import assign_allowable_guides, load_allowed_guides


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_off_target_file(tmp_path, records):
    """
    Write a minimal off-target analysis file.

    Each record is a dict with keys:
      full_query_guide, num_mismatches, mismatch_position, pam_type,
      truncated_target, chrom, pam_coord, strand, full_target_guide, full_target_pam
    """
    header = "\t".join([
        "full_guide", "num_mismatches", "mismatch_position", "PAM_type",
        "truncated_target", "chrom", "idx", "strand", "full_guide", "full_PAM",
    ])
    lines = [header]
    for r in records:
        lines.append("\t".join(str(r[k]) for k in [
            "full_query_guide", "num_mismatches", "mismatch_position", "pam_type",
            "truncated_target", "chrom", "pam_coord", "strand",
            "full_target_guide", "full_target_pam",
        ]))
    outfile = tmp_path / "off_targets.tsv"
    outfile.write_text("\n".join(lines) + "\n")
    return outfile


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestAssignAllowableGuides:
    def test_single_target_allowed(self, tmp_path):
        """A guide with exactly one perfect-match target → allowed."""
        records = [
            dict(full_query_guide="AAAAAAAAAAAAAAAAAAAA",
                 num_mismatches=0, mismatch_position="NA", pam_type="AGG",
                 truncated_target="AAAAAAAAAAAAAAAAAAA",
                 chrom="chrI", pam_coord=100, strand="+",
                 full_target_guide="AAAAAAAAAAAAAAAAAAAA", full_target_pam="AGG"),
        ]
        off_target = _make_off_target_file(tmp_path, records)
        allowed = tmp_path / "allowed.tsv"
        disallowed = tmp_path / "disallowed.tsv"
        n_allowed, n_disallowed = assign_allowable_guides(off_target, allowed, disallowed)
        assert n_allowed == 1
        assert n_disallowed == 0

    def test_two_targets_disallowed(self, tmp_path):
        """A guide with two perfect-match targets → disallowed."""
        records = [
            dict(full_query_guide="GGGGGGGGGGGGGGGGGGGG",
                 num_mismatches=0, mismatch_position="NA", pam_type="CGG",
                 truncated_target="GGGGGGGGGGGGGGGGGGG",
                 chrom="chrI", pam_coord=200, strand="+",
                 full_target_guide="GGGGGGGGGGGGGGGGGGGG", full_target_pam="CGG"),
            dict(full_query_guide="GGGGGGGGGGGGGGGGGGGG",
                 num_mismatches=0, mismatch_position="NA", pam_type="CGG",
                 truncated_target="GGGGGGGGGGGGGGGGGGG",
                 chrom="chrII", pam_coord=500, strand="-",
                 full_target_guide="GGGGGGGGGGGGGGGGGGGG", full_target_pam="CGG"),
        ]
        off_target = _make_off_target_file(tmp_path, records)
        allowed = tmp_path / "allowed.tsv"
        disallowed = tmp_path / "disallowed.tsv"
        n_allowed, n_disallowed = assign_allowable_guides(off_target, allowed, disallowed)
        assert n_allowed == 0
        assert n_disallowed == 2

    def test_seed_mismatch_not_counted(self, tmp_path):
        """Mismatch at position ≤ 10 (seed) is NOT counted as effective off-target.

        Guide A has perfect match at site-A (count=1 → allowed).
        Guide A also has a 1-mismatch to site-X at position 5 (≤10, NOT counted).
        So site-X has count=0 → disallowed (neither site is a UNIQUE effective
        target of guide A due to seed proximity).
        """
        records = [
            # Perfect match: guide=CCCC... → site-A. Count for site-A: 1.
            dict(full_query_guide="CCCCCCCCCCCCCCCCCCCC",
                 num_mismatches=0, mismatch_position="NA", pam_type="TGG",
                 truncated_target="CCCCCCCCCCCCCCCCCCC",
                 chrom="chrI", pam_coord=300, strand="+",
                 full_target_guide="CCCCCCCCCCCCCCCCCCCC", full_target_pam="TGG"),
            # 1-mismatch at seed position 5 → NOT counted for site-X
            dict(full_query_guide="CCCCCCCCCCCCCCCCCCCC",
                 num_mismatches=1, mismatch_position=5, pam_type="TGG",
                 truncated_target="CCCCaCCCCCCCCCCCCCCC",
                 chrom="chrIII", pam_coord=700, strand="+",
                 full_target_guide="CCCCaCCCCCCCCCCCCCCC", full_target_pam="TGG"),
        ]
        off_target = _make_off_target_file(tmp_path, records)
        allowed = tmp_path / "allowed.tsv"
        disallowed = tmp_path / "disallowed.tsv"
        n_allowed, n_disallowed = assign_allowable_guides(off_target, allowed, disallowed)
        # site-A (count=1) → allowed; site-X (count=0) → disallowed
        assert n_allowed == 1
        assert n_disallowed == 1

    def test_outside_seed_mismatch_counted(self, tmp_path):
        """Mismatch at position > 10 IS counted → a site targeted by 2 guides → disallowed.

        Guide A has perfect match at site-A (count for A += 1).
        Guide B also has a 1-mismatch record to site-A at mismatch_position=15 (counted).
        Site-A now has count=2 → disallowed.
        """
        records = [
            # Guide A perfect match at site-A
            dict(full_query_guide="AAAAAAAAAAAAAAAAAAAT",
                 num_mismatches=0, mismatch_position="NA", pam_type="GGG",
                 truncated_target="AAAAAAAAAAAAAAAAAAT",
                 chrom="chrI", pam_coord=400, strand="+",
                 full_target_guide="AAAAAAAAAAAAAAAAAAAT", full_target_pam="GGG"),
            # Guide B has off-target at site-A, mismatch outside seed (pos 15)
            dict(full_query_guide="GAAAAAAAAAAAAAAAAAAAT",
                 num_mismatches=1, mismatch_position=15, pam_type="GGG",
                 truncated_target="AAAAAAAAAAAAAAAAAaAT",
                 chrom="chrI", pam_coord=400, strand="+",
                 full_target_guide="AAAAAAAAAAAAAAAAAAAT", full_target_pam="GGG"),
        ]
        off_target = _make_off_target_file(tmp_path, records)
        allowed = tmp_path / "allowed.tsv"
        disallowed = tmp_path / "disallowed.tsv"
        n_allowed, n_disallowed = assign_allowable_guides(off_target, allowed, disallowed)
        # site-A has count=2 → both records for site-A go to disallowed
        assert n_allowed == 0
        assert n_disallowed == 2

    def test_allowed_guides_header(self, tmp_path):
        records = [
            dict(full_query_guide="AAAAAAAAAAAAAAAAAAAT",
                 num_mismatches=0, mismatch_position="NA", pam_type="AGG",
                 truncated_target="AAAAAAAAAAAAAAAAAAT",
                 chrom="chrV", pam_coord=50, strand="+",
                 full_target_guide="AAAAAAAAAAAAAAAAAAAT", full_target_pam="AGG"),
        ]
        off_target = _make_off_target_file(tmp_path, records)
        allowed = tmp_path / "allowed.tsv"
        disallowed = tmp_path / "disallowed.tsv"
        assign_allowable_guides(off_target, allowed, disallowed)
        with open(allowed) as fh:
            header = fh.readline()
        assert "chrom" in header
        assert "full_guide" in header


@pytest.mark.unit
class TestLoadAllowedGuides:
    def test_returns_set(self, tmp_path):
        allowed_file = tmp_path / "allowed.tsv"
        allowed_file.write_text(
            "chrom\tPAM_coord\tstrand\tfull_guide\tfull_PAM\n"
            "chrI\t100\t+\tACGTACGTACGTACGTACGT\tAGG\n"
        )
        result = load_allowed_guides(allowed_file)
        assert isinstance(result, set)
        assert "ACGTACGTACGTACGTACGT" in result

    def test_empty_file(self, tmp_path):
        allowed_file = tmp_path / "allowed.tsv"
        allowed_file.write_text("chrom\tPAM_coord\tstrand\tfull_guide\tfull_PAM\n")
        result = load_allowed_guides(allowed_file)
        assert len(result) == 0
