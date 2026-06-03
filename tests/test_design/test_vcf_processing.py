"""
Unit tests for magestic.design.vcf_processing.
"""

import io
import tempfile
from pathlib import Path

import pytest

from magestic.design.vcf_processing import (
    load_haplotype_vcf,
    combine_linked_variants,
    filter_variants,
    load_gene_list,
    write_processed_vcf,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINIMAL_HAPLOTYPE_TSV = """\
CHROM\tPOS\tREF\tALT\tSAMPLE\tHAPLOTYPE\tPID\tgenotypes\tall_ALT_alleles\tORF\tNAME
chrI\t100\tA\tT\tS1\tH1\t.\tA|T\tT\tYAL001C\tTFC3
chrI\t200\tGG\tAA\tS1\tH1\t.\tGG|AA\tAA\tYAL002W\tVPS8
chrII\t500\tC\tG\tS2\tH1\t.\tC|G\tG\tYBL001C\tECM15
"""

GENOME_SEQ = {
    "chrI": "A" * 300,
    "chrII": "C" * 600,
}


@pytest.fixture
def haplotype_tsv(tmp_path):
    f = tmp_path / "haplotype.tsv"
    f.write_text(MINIMAL_HAPLOTYPE_TSV)
    return f


@pytest.fixture
def simple_genome():
    return GENOME_SEQ


# ---------------------------------------------------------------------------
# load_haplotype_vcf
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLoadHaplotypeVcf:
    def test_returns_dict(self, haplotype_tsv):
        result = load_haplotype_vcf(haplotype_tsv)
        assert isinstance(result, dict)

    def test_sample_key(self, haplotype_tsv):
        result = load_haplotype_vcf(haplotype_tsv)
        assert "S1" in result
        assert "S2" in result

    def test_variant_loaded(self, haplotype_tsv):
        result = load_haplotype_vcf(haplotype_tsv)
        assert "H1" in result["S1"]
        assert "chrI" in result["S1"]["H1"]
        # POS=100, REF=A, ALT=T
        assert 100 in result["S1"]["H1"]["chrI"]

    def test_snp_classification(self, haplotype_tsv):
        result = load_haplotype_vcf(haplotype_tsv)
        var_tuple = result["S1"]["H1"]["chrI"][100][0]
        pos, ref, alt, desc, orf, name = var_tuple
        assert "SNP" in desc

    def test_mnp_classification(self, haplotype_tsv):
        result = load_haplotype_vcf(haplotype_tsv)
        # POS=200, REF=GG, ALT=AA → MNP
        var_tuple = result["S1"]["H1"]["chrI"][200][0]
        pos, ref, alt, desc, orf, name = var_tuple
        assert "MNP" in desc


# ---------------------------------------------------------------------------
# combine_linked_variants
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestCombineLinkedVariants:
    def test_single_variants_present(self, haplotype_tsv, simple_genome):
        variants_by_sample = load_haplotype_vcf(haplotype_tsv)
        combined = combine_linked_variants(variants_by_sample, simple_genome)
        # Each individual variant should appear
        assert len(combined) >= 3

    def test_combined_key_format(self, haplotype_tsv, simple_genome):
        variants_by_sample = load_haplotype_vcf(haplotype_tsv)
        combined = combine_linked_variants(variants_by_sample, simple_genome)
        for key in combined:
            assert len(key) == 7  # (chrom, pos, ref, alt, desc, orf, name)


# ---------------------------------------------------------------------------
# filter_variants
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestFilterVariants:
    def _make_combined(self):
        return {
            ("chrI", 100, "A", "T", "100_A_T_SNP", "YAL001C", "TFC3"): ["S1"],
            ("chrI", 200, "GG", "AA", "200_GG_AA_MNP", "background_variant", "background_variant"): ["S1"],
            ("chrII", 500, "C", "G", "500_C_G_SNP", "YBL001C", "ECM15"): ["S2"],
        }

    def test_exclude_background(self):
        combined = self._make_combined()
        filtered = filter_variants(combined, exclude_background=True)
        assert all(k[5] != "background_variant" for k in filtered)

    def test_keep_background_when_disabled(self):
        combined = self._make_combined()
        filtered = filter_variants(combined, exclude_background=False)
        assert any(k[5] == "background_variant" for k in filtered)

    def test_gene_list_filter(self):
        combined = self._make_combined()
        filtered = filter_variants(combined, gene_list={"YAL001C"}, exclude_background=False)
        assert len(filtered) == 1
        assert list(filtered.keys())[0][5] == "YAL001C"

    def test_exclude_gene_names(self):
        combined = self._make_combined()
        filtered = filter_variants(combined, exclude_gene_names={"TFC3"}, exclude_background=False)
        assert all(k[6] != "TFC3" for k in filtered)


# ---------------------------------------------------------------------------
# load_gene_list
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestLoadGeneList:
    def test_single_column(self, tmp_path):
        f = tmp_path / "genes.txt"
        f.write_text("YAL001C\nYBL001C\n")
        genes = load_gene_list(f)
        assert "YAL001C" in genes
        assert "YBL001C" in genes

    def test_two_column(self, tmp_path):
        f = tmp_path / "genes.tsv"
        f.write_text("YAL001C\tTFC3\nYBL001C\tECM15\n")
        genes = load_gene_list(f)
        assert "TFC3" in genes
        assert "ECM15" in genes
        assert "YAL001C" in genes


# ---------------------------------------------------------------------------
# write_processed_vcf
# ---------------------------------------------------------------------------

@pytest.mark.unit
class TestWriteProcessedVcf:
    def test_creates_file(self, haplotype_tsv, simple_genome, tmp_path):
        variants_by_sample = load_haplotype_vcf(haplotype_tsv)
        combined = combine_linked_variants(variants_by_sample, simple_genome)
        outfile = tmp_path / "processed.tsv"
        write_processed_vcf(combined, outfile)
        assert outfile.exists()

    def test_header_columns(self, haplotype_tsv, simple_genome, tmp_path):
        variants_by_sample = load_haplotype_vcf(haplotype_tsv)
        combined = combine_linked_variants(variants_by_sample, simple_genome)
        outfile = tmp_path / "processed.tsv"
        write_processed_vcf(combined, outfile)
        with open(outfile) as fh:
            header = fh.readline()
        assert "CHROM" in header
        assert "REF" in header
        assert "ALT" in header
