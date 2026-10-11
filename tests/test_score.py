"""Structure scoring without the language model."""
import numpy as np

from conftest import DATA
from rlmhelix.score import decode_pairs, nussinov_pairs, rlmhelix_score, stems


def test_stored_map_gives_stored_score():
    data = np.load(DATA / "window_apc.npz")
    result = rlmhelix_score(data["apc"], str(data["sequence"]))
    assert abs(result.score - float(data["expected_score"])) < 1e-6
    assert result.structure == str(data["expected_structure"])


def test_hairpin_is_decoded():
    sequence = "GGGGAAAACCCC"
    weights = np.zeros((12, 12), dtype=np.float32)
    for k in range(4):
        weights[k, 11 - k] = weights[11 - k, k] = 1.0
    assert nussinov_pairs(weights, sequence) == [(0, 11), (1, 10), (2, 9), (3, 8)]
    assert [len(s) for s in stems(decode_pairs(weights, sequence))] == [4]


def test_isolated_pairs_and_short_stems_are_not_scored():
    sequence = "GAAAACAAAAGGGAAAACCC"
    weights = np.zeros((20, 20), dtype=np.float32)
    weights[0, 5] = weights[5, 0] = 1.0                      # isolated pair, removed
    for i, j in [(10, 19), (11, 18), (12, 17)]:              # three-pair stem, below minL4
        weights[i, j] = weights[j, i] = 1.0
    bpp = np.ones((20, 20), dtype=np.float32)
    result = rlmhelix_score(weights, sequence, bpp=bpp)
    assert result.score == 0.0 and result.n_pairs == 3 and result.n_stems_scored == 0
