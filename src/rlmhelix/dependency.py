"""Teacher dependency maps from an RNA language model.

Substitution-based nucleotide dependencies (3L + 1 forward passes for a sequence of length L): each position is
replaced by each of the three other nucleotides, and the dependency of position j on position i is the largest
absolute log2 odds ratio of the predicted nucleotide probabilities at j, over the substitutions at i and the four
nucleotides at j. Directional dependencies are symmetrised by their maximum and corrected by APC.

The teacher is slow and is not needed for scoring (the student predicts the APC map directly); it is included to
reproduce the training targets and the teacher–student comparisons.
"""
from __future__ import annotations

import contextlib

import numpy as np

BACKBONE = "multimolecule/rinalmo-giga"
# Hugging Face revision of the backbone used for the paper (also used by rlmhelix.student). Defined here, not in
# rlmhelix.student, so that this module does not import torch.
BACKBONE_REVISION = "2a71f6f98fb41dd2e6542a5e131d3778111d1468"
NUCLEOTIDES = "ACGU"


def apc(dependency: np.ndarray, min_loop_size: int = 3) -> np.ndarray:
    """Average product correction with means over pairs |i - j| >= min_loop_size + 1; negatives set to 0."""
    scores = np.asarray(dependency, dtype=np.float32)
    length = scores.shape[0]
    corrected = np.zeros_like(scores, dtype=np.float32)
    index = np.arange(length)
    valid = np.abs(np.subtract.outer(index, index)) >= min_loop_size + 1
    if not valid.any():
        return corrected
    masked = np.where(valid, scores, 0.0)
    counts = valid.sum(axis=1).astype(np.float32)
    row_mean = np.divide(masked.sum(axis=1), counts, out=np.zeros(length, dtype=np.float32), where=counts > 0)
    global_mean = float(masked[valid].mean())
    if global_mean == 0.0:
        corrected[valid] = masked[valid]
    else:
        corrected[valid] = masked[valid] - (np.outer(row_mean, row_mean) / global_mean)[valid]
    corrected = np.maximum(corrected, 0.0)
    np.fill_diagonal(corrected, 0.0)
    return corrected


class TeacherModel:
    """RiNALMo-giga masked-language-model probabilities over A, C, G, U (no masking).

    `revision` is the Hugging Face revision of `model_name` (by default the one used for the paper); pass None with
    another model.
    """

    def __init__(self, device: str = "cuda", dtype: str = "bf16", model_name: str = BACKBONE,
                 revision: str | None = BACKBONE_REVISION):
        import torch
        from multimolecule import RiNALMoForMaskedLM, RnaTokenizer

        self.device = device
        self.dtype = {"bf16": torch.bfloat16, "fp16": torch.float16, "fp32": None}[dtype]
        self.tokenizer = RnaTokenizer.from_pretrained(model_name, revision=revision)
        self.model = RiNALMoForMaskedLM.from_pretrained(model_name, revision=revision,
                                                        attn_implementation="sdpa").to(device).eval()
        vocab = self.tokenizer.get_vocab()
        self.ids = [vocab[base] for base in NUCLEOTIDES]

    def _autocast(self):
        import torch

        if self.dtype is None or self.device == "cpu":
            return contextlib.nullcontext()
        return torch.autocast(device_type=self.device, dtype=self.dtype)

    def probabilities(self, sequences: list[str], batch_size: int = 128) -> np.ndarray:
        """(N, L, 4) probabilities for sequences of equal length L."""
        import torch
        import torch.nn.functional as F

        chunks = []
        for start in range(0, len(sequences), batch_size):
            batch = [s.upper().replace("T", "U") for s in sequences[start:start + batch_size]]
            tokens = self.tokenizer(batch, return_tensors="pt", padding=True).to(self.device)
            with torch.no_grad(), self._autocast():
                logits = self.model(**tokens).logits
            chunks.append(F.softmax(logits.float()[:, 1:-1, :][:, :, self.ids], dim=-1).cpu().numpy())
        return np.concatenate(chunks, axis=0)


def dependency_map(sequence: str, model: TeacherModel, batch_size: int = 128, eps: float = 1e-7) -> np.ndarray:
    """Symmetrised (maximum over directions) dependency map, before APC."""
    sequence = sequence.upper().replace("T", "U")
    length = len(sequence)
    variants, meta = [sequence], []
    for i, ref in enumerate(sequence):
        if ref not in NUCLEOTIDES:
            continue
        for alt in NUCLEOTIDES:
            if alt != ref:
                variants.append(sequence[:i] + alt + sequence[i + 1:])
                meta.append(i)
    probs = np.clip(model.probabilities(variants, batch_size), eps, 1.0 - eps)
    log_odds_ref = np.log2(probs[0] / (1.0 - probs[0]))
    directed = np.zeros((length, length, 4), dtype=np.float32)
    for k, i in enumerate(meta):
        log_odds = np.log2(probs[k + 1] / (1.0 - probs[k + 1]))
        np.maximum(directed[i], np.abs(log_odds - log_odds_ref), out=directed[i])
    dep = directed.max(axis=2)
    np.fill_diagonal(dep, 0.0)
    return np.maximum(dep, dep.T)


def teacher_apc(sequence: str, model: TeacherModel, batch_size: int = 128) -> np.ndarray:
    """APC-corrected teacher dependency map (the student's training target)."""
    return apc(dependency_map(sequence, model, batch_size))
