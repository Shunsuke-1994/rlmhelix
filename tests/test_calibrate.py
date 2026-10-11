"""Calibrated P values of the paper's scan."""
import numpy as np

from rlmhelix import Calibrator


def test_pvalues_match_scan(windows):
    p = Calibrator.load().pvalues(windows.score, list(windows.sequence))
    np.testing.assert_allclose(p, windows.p_value, rtol=1e-9, atol=1e-300)
