"""Score calibration.

Conditional generalized extreme value (GEV) model of the background distribution of log(1 + score), with location
and log-scale linear in 25 standardised features of length and mono- and dinucleotide composition and one shared
shape, fitted on 570,000 scores of dinucleotide-preserving shuffles (uniform, Dirichlet and genomic sequences).
The upper-tail probability is the calibrated P value.

The calibration was built for the 200-nt windows of the genome scan; it is anti-conservative for short sequences
(30-99 nt), so compare short sequences with their own shuffles instead.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

import numpy as np

FEATURES = (["L", "pct_A", "pct_C", "pct_G", "pct_U", "GC", "A_over_AU", "C_over_CG"]
            + [f"di_{a}{b}" for a in "ACGU" for b in "ACGU"] + ["log_L"])


def sequence_features(sequence: str) -> dict[str, float]:
    """Length and composition features used by the calibration."""
    sequence = sequence.upper().replace("T", "U")
    length = max(len(sequence), 1)
    counts = {base: sequence.count(base) for base in "ACGU"}
    au, cg = counts["A"] + counts["U"], counts["C"] + counts["G"]
    out = {"L": float(len(sequence)), **{f"pct_{base}": counts[base] / length for base in "ACGU"},
           "GC": cg / length, "A_over_AU": counts["A"] / au if au else 0.0,
           "C_over_CG": counts["C"] / cg if cg else 0.0}
    denominator = max(len(sequence) - 1, 1)
    pairs = [sequence[k:k + 2] for k in range(len(sequence) - 1)]
    for left in "ACGU":
        for right in "ACGU":
            out[f"di_{left}{right}"] = pairs.count(left + right) / denominator
    out["log_L"] = float(np.log(len(sequence) + 1))
    return out


@dataclass
class Calibrator:
    feature_mu: np.ndarray
    feature_sd: np.ndarray
    loc_beta: np.ndarray
    scale_beta: np.ndarray
    xi: float

    @classmethod
    def load(cls, path: str | Path | None = None) -> "Calibrator":
        path = path or Path(str(resources.files("rlmhelix").joinpath("models/calibrator_gev.json")))
        model = json.loads(Path(path).read_text())
        if model["feature_names"] != FEATURES or model["transform"] != "log1p":
            raise ValueError(f"{path}: unexpected calibrator")
        return cls(*(np.asarray(model[k], dtype=np.float64)
                     for k in ("feature_mu", "feature_sd", "loc_beta", "scale_beta")), float(model["xi"]))

    def pvalues(self, scores, sequences: list[str]) -> np.ndarray:
        """Calibrated upper-tail P values of RLMHelix scores for the given sequences."""
        from scipy.stats import genextreme

        x_raw = np.asarray([[f[k] for k in FEATURES] for f in map(sequence_features, sequences)], dtype=np.float64)
        x = np.column_stack([(x_raw - self.feature_mu) / self.feature_sd, np.ones(len(x_raw))])
        loc = x @ self.loc_beta
        scale = np.exp(np.clip(x @ self.scale_beta, -5, 8))
        p = genextreme.sf(np.log1p(np.asarray(scores, dtype=np.float64)), -self.xi, loc=loc, scale=scale)
        return np.clip(np.asarray(p, dtype=np.float64), 0, 1)
