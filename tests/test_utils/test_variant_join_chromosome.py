"""Variant annotations join designs by chromosome and allele string (register U09).

Allele strings give a position but no chromosome. The Bloom codon-change table uses
559546_G_A_SNP for TOR1 R45K (YJR066W, chrX) and for CAT5 G62G (YOR125C, chrXV); the
shipped common registry (md5 2dce503d...) carries four CAT5 designs labelled R45K.
"""
import pandas as pd
import pytest

from magestic.utils import harmonize_oligos as h

STRING = '559546_G_A_SNP'


def test_explicit_chromosome_selects_the_matching_annotation():
    oligos = pd.DataFrame({'oligo_name': ['chrII_design'], 'chrom': ['chrII'],
                           'individual_variants': ['100_A_G']})
    variants = pd.DataFrame({'chrom': ['chrI', 'chrII'], 'individual_variants': ['100_A_G', '100_A_G'],
                             'annotation': ['WRONG_chrI', 'RIGHT_chrII']})
    merged = h.merge_variants_info(oligos, variants)
    assert merged['annotation'].tolist() == ['RIGHT_chrII']


def test_chromosome_is_derived_from_the_annotation_orf():
    oligos = pd.DataFrame({'oligo_name': ['tor1_design', 'cat5_design'], 'chrom': ['chrX', 'chrXV'],
                           'individual_variants': [STRING, STRING]})
    variants = pd.DataFrame({'systematic_ORF_name': ['YJR066W', 'YOR125C'], 'gene_name': ['TOR1', 'CAT5'],
                             'aa_nums': [45, 62], 'ref_aas': ['R', 'G'], 'alt_aas': ['K', 'G'],
                             'individual_variants': [STRING, STRING]})
    merged = h.merge_variants_info(oligos, variants)
    assert merged[['oligo_name', 'gene_name', 'alt_aas']].values.tolist() == [
        ['tor1_design', 'TOR1', 'K'], ['cat5_design', 'CAT5', 'G']]


@pytest.mark.parametrize('name,chrom', [('YAL001C', 'chrI'), ('YJR066W', 'chrX'), ('YOR125C', 'chrXV'),
                                        ('YPR200C', 'chrXVI'), ('YOL083W', 'chrXV'),
                                        ('Q0010', None), (float('nan'), None), ('TPK2', None)])
def test_chromosome_from_systematic_name(name, chrom):
    assert h.chromosome_from_systematic_name(name) == chrom


def test_duplicate_annotation_keys_are_refused_not_multiplied():
    oligos = pd.DataFrame({'oligo_name': ['d1'], 'chrom': ['chrXV'], 'individual_variants': ['166088_G_A_SNP']})
    variants = pd.DataFrame({'systematic_ORF_name': ['YOL083W', 'YOL084W'], 'gene_name': ['ATG34', 'PHM7'],
                             'individual_variants': ['166088_G_A_SNP', '166088_G_A_SNP']})
    with pytest.raises(ValueError, match='not unique'):
        h.merge_variants_info(oligos, variants)


def test_string_only_join_refuses_a_string_on_two_chromosomes():
    oligos = pd.DataFrame({'oligo_name': ['a', 'b'], 'chrom': ['chrX', 'chrXV'],
                           'individual_variants': [STRING, STRING]})
    variants = pd.DataFrame({'individual_variants': [STRING], 'annotation': ['one']})
    with pytest.raises(ValueError, match='more than one chromosome'):
        h.merge_variants_info(oligos, variants)


def test_string_only_join_still_works_when_unambiguous():
    oligos = pd.DataFrame({'oligo_name': ['a'], 'chrom': ['chrX'], 'individual_variants': [STRING]})
    variants = pd.DataFrame({'individual_variants': [STRING], 'annotation': ['one']})
    assert h.merge_variants_info(oligos, variants)['annotation'].tolist() == ['one']


def test_unplaceable_annotation_is_refused():
    variants = pd.DataFrame({'systematic_ORF_name': ['not_an_orf'], 'individual_variants': [STRING]})
    with pytest.raises(ValueError, match='no chromosome'):
        h.add_variant_chromosome(variants)
