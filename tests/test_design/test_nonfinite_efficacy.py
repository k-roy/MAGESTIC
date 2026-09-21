"""Missing model predictions must not become the highest calibrated efficacy."""
import math
import pytest

from magestic.design.domain_insertion.guides import load_score_predictions, guide_quality
from magestic.design.domain_insertion.score import design_score, p_eff


@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), float('-inf')])
def test_nonfinite_efficacy_is_unscored(invalid):
    assert design_score(invalid, 0, 0, 0, cut_distance=0, log=False) == 0.0
    assert design_score(invalid, 0, 0, 0, cut_distance=0) == float('-inf')


@pytest.mark.parametrize('invalid', [float('nan'), float('inf'), float('-inf')])
def test_calibration_refuses_nonfinite_input(invalid):
    with pytest.raises(ValueError, match='finite'):
        p_eff(invalid)


@pytest.mark.parametrize('invalid', ['NaN', 'inf', '-inf'])
def test_score_loader_does_not_create_prediction_from_nonfinite_unedited(tmp_path, invalid):
    path = tmp_path/'score.tsv'
    path.write_text('guide\tSCORE.NE.minmax\tSCORE.SV.minmax\n'
                    +'A'*20+'\t'+invalid+'\t0\n'
                    +'C'*20+'\t0.2\t0.1\n')
    scores = load_score_predictions(str(path))
    assert 'A'*20 not in scores
    assert guide_quality('A'*20, scores) is None
    assert scores['C'*20]['efficacy'] == pytest.approx(0.8)
    assert guide_quality('C'*20, scores) == pytest.approx(0.7)


@pytest.mark.parametrize('value', [0., 0.59, 0.6, 0.9, 0.95, 1.])
def test_finite_calibration_and_score_stay_rankable(value):
    assert 0 < p_eff(value) <= 1
    raw = design_score(value, 0.1, 0, 0, cut_distance=0, log=False)
    assert 0 < raw < 1
    assert math.exp(design_score(value, 0.1, 0, 0, cut_distance=0)) == pytest.approx(raw)
