"""SV score imputation must preserve the row receiving each chromosome's score."""
import importlib
import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def harmonizer(tmp_path, monkeypatch):
    # Existing module-level directory setup is a separately tracked CLI issue.
    monkeypatch.setenv('OAK', str(tmp_path/'oak'))
    return importlib.import_module('magestic.utils.harmonize_oligos')


@pytest.mark.parametrize('labels', [[0, 1, 2], [8, 3, 8]])
def test_interleaved_chromosomes_keep_row_identity(harmonizer, labels):
    variants = pd.DataFrame({'chrom': ['chrI', 'chrII', 'chrI'], 'PAM_coord': [10, 50, 100],
                             'SV_prediction_score': [np.nan]*3, 'identity': ['first', 'second', 'third']}, index=labels)
    reference = pd.DataFrame({'chr': ['I', 'I', 'I', 'I', 'II'], 'pos': [9, 11, 99, 101, 50],
                              'SCORE.SV.minmax': [0.1, 0.1, 0.9, 0.9, 0.5]})
    result = harmonizer.impute_missing_sv_scores(variants, reference)
    assert result.SV_prediction_score.tolist() == pytest.approx([0.1, 0.5, 0.9])
    assert result.index.tolist() == labels
    assert result.identity.tolist() == ['first', 'second', 'third']


def test_unknown_coordinates_stay_missing_and_observed_scores_stay_fixed(harmonizer):
    variants = pd.DataFrame({'chrom': ['chrI', None, 'chrII', 'chrI', 'chrI'],
                             'PAM_coord': [10, 10, 10, np.nan, 10],
                             'SV_prediction_score': [np.nan, np.nan, np.nan, np.nan, 0.7]})
    reference = pd.DataFrame({'chr': ['I'], 'pos': [10], 'SCORE.SV.minmax': [0.2]})
    result = harmonizer.impute_missing_sv_scores(variants, reference)
    assert result.iloc[0].SV_prediction_score == pytest.approx(0.2)
    assert result.iloc[1:4].SV_prediction_score.isna().all()
    assert result.iloc[4].SV_prediction_score == pytest.approx(0.7)
