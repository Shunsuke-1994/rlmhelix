"""Distilled student head.

A bilinear head on the frozen RNA language model (RiNALMo-giga) that outputs the APC-corrected dependency map in
one forward pass. Each position's 1,280-dimensional representation is layer-normalised and passed to a source and
a target tower (Linear(1280, 768) -> SiLU -> Dropout -> Linear(768, 128) -> SiLU); with outputs S and T,
A_ij = softplus((Z_ij + Z_ji) / 2), Z = S T^T / sqrt(128), with the diagonal set to 0.
"""
from __future__ import annotations

import contextlib
import math
import warnings
from importlib import resources
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from .dependency import BACKBONE, BACKBONE_REVISION  # multimolecule/rinalmo-giga at the revision of the paper

TRAINED_LENGTHS = (30, 300)


class StudentHead(nn.Module):
    """Bilinear head (parameter names match the trained checkpoint)."""

    def __init__(self, embed_dim: int = 1280, position_hidden_dim: int = 768, pair_dim: int = 128,
                 dropout: float = 0.1):
        super().__init__()
        self.pair_dim = pair_dim
        self.norm = nn.LayerNorm(embed_dim)

        def tower():
            return nn.Sequential(nn.Linear(embed_dim, position_hidden_dim), nn.SiLU(), nn.Dropout(dropout),
                                 nn.Linear(position_hidden_dim, pair_dim), nn.SiLU())

        self.source_tower = tower()
        self.target_tower = tower()

    def forward(self, hidden: torch.Tensor) -> torch.Tensor:
        """(B, L, 1280) hidden states -> (B, L, L) non-negative APC-corrected dependencies."""
        normalized = self.norm(hidden)
        directed = torch.einsum("blr,bmr->blm", self.source_tower(normalized), self.target_tower(normalized))
        directed = directed / math.sqrt(self.pair_dim)
        prediction = F.softplus(0.5 * (directed + directed.transpose(-1, -2)))
        diagonal = torch.eye(prediction.shape[-1], dtype=torch.bool, device=prediction.device)
        return prediction.masked_fill(diagonal[None], 0.0)


def default_weights() -> Path:
    return Path(str(resources.files("rlmhelix").joinpath("models/student_bilinear_apc.pt")))


def load_head(path: str | Path | None = None, device: str = "cpu") -> StudentHead:
    checkpoint = torch.load(path or default_weights(), map_location=device, weights_only=True)
    config = checkpoint["config"]
    head = StudentHead(config["embed_dim"], config["position_hidden_dim"], config["pair_dim"], config["dropout"])
    head.load_state_dict(checkpoint["state_dict"])
    return head.to(device).eval()


class Student:
    """RiNALMo-giga backbone plus the student head: APC-corrected dependency maps for many sequences.

    On a GPU the backbone and head run under bf16 autocast, as for the paper; on the CPU they run in float32,
    which gives slightly different values. `revision` is the Hugging Face revision of `backbone` (by default the one
    used for the paper); pass None with another backbone.
    """

    def __init__(self, device: str | None = None, dtype: str = "bf16", weights: str | Path | None = None,
                 backbone: str = BACKBONE, revision: str | None = BACKBONE_REVISION):
        from multimolecule import RiNALMoModel, RnaTokenizer

        self.device = device or ("cuda" if torch.cuda.is_available() else "cpu")
        self.dtype = None if self.device == "cpu" else {"bf16": torch.bfloat16, "fp16": torch.float16,
                                                         "fp32": None}[dtype]
        self.tokenizer = RnaTokenizer.from_pretrained(backbone, revision=revision)
        # The head reads last_hidden_state, so the pooler (not in the checkpoint) is not built.
        self.backbone = RiNALMoModel.from_pretrained(backbone, revision=revision, attn_implementation="sdpa",
                                                     add_pooling_layer=False).to(self.device).eval()
        self.head = load_head(weights, self.device)

    def _autocast(self):
        if self.dtype is None:
            return contextlib.nullcontext()
        return torch.autocast(device_type=self.device, dtype=self.dtype)

    def predict(self, sequences: list[str], batch_size: int = 16, per_sequence_head: bool = False) -> list[np.ndarray]:
        """APC-corrected dependency maps, one (L, L) float32 array per sequence.

        The backbone always runs on whole batches. The head runs once per batch when all sequences have the same
        length (as in the genome scan) and once per sequence otherwise, or always with `per_sequence_head=True`
        (as for the alignment scores of the paper); the two differ only at the level of bf16 rounding.
        """
        sequences = [s.upper().replace("T", "U") for s in sequences]
        low, high = TRAINED_LENGTHS
        if any(not low <= len(s) <= high for s in sequences):
            warnings.warn(f"the student was trained on {low}-{high} nt; other lengths are extrapolated", stacklevel=2)
        maps: list[np.ndarray] = []
        for start in range(0, len(sequences), batch_size):
            batch = sequences[start:start + batch_size]
            tokens = self.tokenizer(batch, return_tensors="pt", padding=True).to(self.device)
            with torch.inference_mode(), self._autocast():
                hidden = self.backbone(**tokens).last_hidden_state
                lengths = {len(s) for s in batch}
                if len(lengths) == 1 and not per_sequence_head:  # equal lengths (scan windows): one call
                    length = lengths.pop()
                    values = self.head(hidden[:, 1:length + 1].float()).float().cpu().numpy()
                    maps.extend(values.astype(np.float32))
                else:  # mixed lengths: drop padding per sequence
                    for i, sequence in enumerate(batch):
                        value = self.head(hidden[i:i + 1, 1:len(sequence) + 1].float())[0]
                        maps.append(value.float().cpu().numpy().astype(np.float32))
        return maps
