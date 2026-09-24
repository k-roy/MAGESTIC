"""Step 06 attributes a donor to a design only when it carries every designed allele (218 patch 53).

Registers BDB-H13 and the BDB-H14 family.
- H13: variants were compared by position only, so a designed A>C with an observed A>G gave no delta
  and was classed as the designed edit.
- H14: with no designed variant, every observed variant was a delta at infinite distance, classed
  CAN_ATTRIBUTE.
- H14 family, beyond the audit's narrow H14 claim: a donor with no variants at all (NO_VARIANTS),
  or missing some designed alleles (NO_DELTA), was classed attributable although it lacks the edit.

The designed variants come from the harmonized key without a chromosome ("25340_A_C_SNP"), and indels
are VCF-anchored like the observed ones. On the real references, 91.7-91.8% of close matches carry
every designed allele (probes/step06_ground_truth.py). A designed entry with ref equal to alt is a
no-op and is not required.
"""
import importlib.util
import sys
from pathlib import Path

import pytest

from magestic.pipelines.bc1_donor_bc0.readers import data_loader

ATTRIBUTABLE = {'tier_1_low_ed_attributable', 'tier_3_mod_ed_attributable'}


def _step06():
    path = Path(data_loader.__file__).parents[1] / 'snakemake/06_annotate_unmatched_donors.py'
    spec = importlib.util.spec_from_file_location('test_step06_rescue_identity', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope='module')
def step06():
    return _step06()


def classify(step06, observed, designed):
    info = step06.OligoInfo('oligo_1', 'G' * 20, 'A' * 129, 'A' * 60, '')
    return step06.classify_rescue_status(observed, designed, info)


def test_a_different_allele_at_the_designed_site_is_not_the_designed_edit(step06):
    observed, designed = [('chrI', 61, 'A', 'G')], [('', 61, 'A', 'C')]
    assert step06.compute_delta_variants(observed, designed) == observed
    result = classify(step06, observed, designed)
    assert result['rescue_status'] == 'DESIGNED_VARIANT_ABSENT' and result['designed_variants_missing'] == 1


def test_no_designed_variant_is_an_unknown_design(step06):
    assert classify(step06, [('chrI', 61, 'A', 'G')], [])['rescue_status'] == 'DESIGN_UNKNOWN'


def test_a_donor_without_variants_lacks_the_designed_allele(step06):
    result = classify(step06, [], [('', 61, 'A', 'C')])
    assert result['rescue_status'] == 'DESIGNED_VARIANT_ABSENT' and result['designed_variants_missing'] == 1


def test_a_donor_carrying_part_of_the_design_is_not_attributed(step06):
    result = classify(step06, [('chrI', 61, 'A', 'C')], [('', 61, 'A', 'C'), ('', 64, 'G', 'T')])
    assert result['rescue_status'] == 'DESIGNED_VARIANT_ABSENT' and result['designed_variants_missing'] == 1


def test_the_exact_designed_alleles_are_attributable(step06):
    """Controls: the designed SNV alone, and with an extra variant 25 bp away (beyond 20 bp)."""
    assert classify(step06, [('chrI', 61, 'A', 'C')], [('', 61, 'A', 'C')])['rescue_status'] == 'NO_DELTA'
    extra = classify(step06, [('chrI', 61, 'A', 'C'), ('chrI', 86, 'G', 'T')], [('', 61, 'A', 'C')])
    assert extra['rescue_status'] == 'CAN_ATTRIBUTE' and extra['closest_delta_to_designed'] == 25


def test_an_anchored_indel_matches_its_design(step06):
    observed, designed = [('chrI', 100, 'AT', 'A')], [('', 100, 'AT', 'A')]
    assert classify(step06, observed, designed)['rescue_status'] == 'NO_DELTA'


def test_a_no_op_designed_entry_is_not_required(step06):
    observed = [('chrI', 61, 'A', 'C')]
    assert classify(step06, observed, [('', 61, 'A', 'C'), ('', 70, 'G', 'G')])['rescue_status'] == 'NO_DELTA'
    assert classify(step06, observed, [('', 70, 'G', 'G')])['rescue_status'] == 'DESIGN_UNKNOWN'


@pytest.mark.parametrize('status', ['DESIGN_UNKNOWN', 'DESIGNED_VARIANT_ABSENT'])
@pytest.mark.parametrize('edit_distance', [1, 7])
def test_the_new_statuses_are_never_attributable(step06, status, edit_distance):
    annotation = step06.AnnotationResult(bc1='A' * 20, closest_oligo_name='oligo_1',
                                         closest_donor_edit_distance=edit_distance, rescue_status=status)
    assert step06.assign_confidence_tier(annotation) not in ATTRIBUTABLE
