"""End-to-end checks with the language model (GPU, bf16), against values of the paper.

Single sequences were scored in batches of 16 for the paper; other batch sizes change the bf16 rounding of the
backbone, so scores here agree to about 0.1% rather than exactly.
"""
import json

import numpy as np
import pytest

from conftest import DATA
from rlmhelix import Calibrator, rlmhelix_score
from rlmhelix.msa import read_stockholm, score_alignment

pytestmark = pytest.mark.model


def test_window_scores(windows, student):
    maps = student.predict(list(windows.sequence))
    scores = [rlmhelix_score(m, s).score for m, s in zip(maps, windows.sequence)]
    np.testing.assert_allclose(scores, windows.score, rtol=5e-3)
    p = Calibrator.load().pvalues(scores, list(windows.sequence))
    assert ((p <= 0.05) == (windows.p_value <= 0.05)).all()


def test_alignment_score(student):
    expected = json.loads((DATA / "msa_expected.json").read_text())["RF00031"]
    result = score_alignment(read_stockholm(DATA / "RF00031.sto"), student)
    for key, value in expected.items():
        assert result[key] == pytest.approx(value, abs=1e-4), key
