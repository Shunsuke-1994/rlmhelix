"""RLMHelix score.

Pair weights from the geometric mean of the dependency and the thermodynamic base-pair probability (BPP),
non-crossing Nussinov dynamic programming, removal of isolated pairs, and the sum of weights over stems of
at least four consecutive base pairs (minL4).

The arithmetic follows the implementation used for the paper (float32 weights and DP, per-stem float32 sums),
so that scores computed here match the published scan values.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - numba is a dependency; this keeps the module importable without it
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda function: function

MIN_LOOP_SIZE = 3
MIN_STEM = 4


def to_rna(sequence: str) -> str:
    return sequence.upper().replace("T", "U")


def base_pair_probabilities(sequence: str) -> np.ndarray:
    """Symmetric (L, L) float32 matrix of RNAfold base-pair probabilities (ViennaRNA, 37 °C, defaults)."""
    import RNA

    sequence = to_rna(sequence)
    compound = RNA.fold_compound(sequence)
    compound.pf()
    raw = np.asarray(compound.bpp())
    length = len(sequence)
    result = np.zeros((length, length), dtype=np.float32)
    upper = np.triu_indices(length, k=1)
    result[upper] = raw[1:, 1:][upper]
    result.T[upper] = result[upper]
    return result


def _encode(sequence: str) -> np.ndarray:
    """A=0, C=1, G=2, U/T=3, anything else -1."""
    table = {"A": 0, "C": 1, "G": 2, "U": 3, "T": 3}
    return np.asarray([table.get(base, -1) for base in sequence.upper()], dtype=np.int8)


@njit(cache=True)
def _nussinov(scores, seq_int, min_sep):
    """Maximum total weight of non-crossing AU/GC/GU pairs (ties keep the first option, as in the paper)."""
    length = scores.shape[0]
    dp = np.zeros((length, length), dtype=np.float32)
    choice = np.full((length, length), -2, dtype=np.int32)
    for span in range(1, length):
        for i in range(0, length - span):
            j = i + span
            best = dp[i + 1, j]
            best_choice = -1
            if dp[i, j - 1] > best:
                best = dp[i, j - 1]
                best_choice = -3
            if j - i >= min_sep:
                a = seq_int[i]
                b = seq_int[j]
                if ((a == 0 and b == 3) or (a == 3 and b == 0) or (a == 1 and b == 2) or (a == 2 and b == 1)
                        or (a == 2 and b == 3) or (a == 3 and b == 2)):
                    pair_total = scores[i, j]
                    if i + 1 <= j - 1:
                        pair_total += dp[i + 1, j - 1]
                    if pair_total > best:
                        best = pair_total
                        best_choice = j
            for k in range(i + 1, j):
                split_total = dp[i, k] + dp[k + 1, j]
                if split_total > best:
                    best = split_total
                    best_choice = k
            dp[i, j] = best
            choice[i, j] = best_choice
    return dp[0, length - 1], choice


def nussinov_pairs(weights: np.ndarray, sequence: str, min_loop_size: int = MIN_LOOP_SIZE) -> list[tuple[int, int]]:
    """Base pairs (0-based, i < j) of the non-crossing structure with the largest total weight."""
    scores = np.ascontiguousarray(np.asarray(weights, dtype=np.float32))
    length = scores.shape[0]
    if scores.shape != (length, length) or len(sequence) != length:
        raise ValueError("weights must be (L, L) for a sequence of length L")
    if length == 0:
        return []
    _, choice = _nussinov(scores, _encode(sequence), min_loop_size + 1)
    pairs: list[tuple[int, int]] = []
    stack = [(0, length - 1)]
    while stack:
        i, j = stack.pop()
        if i >= j:
            continue
        decision = int(choice[i, j])
        if decision == -1:
            stack.append((i + 1, j))
        elif decision == -3:
            stack.append((i, j - 1))
        elif decision == j:
            pairs.append((i, j))
            stack.append((i + 1, j - 1))
        elif decision >= i + 1:
            stack.append((decision + 1, j))
            stack.append((i, decision))
    return sorted(pairs)


def decode_pairs(weights: np.ndarray, sequence: str) -> set[tuple[int, int]]:
    """Nussinov decoding followed by removal of isolated pairs (pairs without a stacked neighbour)."""
    pairs = set(nussinov_pairs(weights, sequence))
    return {(i, j) for i, j in pairs if (i - 1, j + 1) in pairs or (i + 1, j - 1) in pairs}


def stems(pairs: set[tuple[int, int]]) -> list[list[tuple[int, int]]]:
    """Group pairs into stems of consecutive pairs (i, j), (i+1, j-1), ..., outermost pair first."""
    remaining = set(pairs)
    blocks: list[list[tuple[int, int]]] = []
    while remaining:
        i, j = min(remaining)
        while (i - 1, j + 1) in remaining:
            i, j = i - 1, j + 1
        block = []
        while (i, j) in remaining:
            block.append((i, j))
            remaining.remove((i, j))
            i, j = i + 1, j - 1
        blocks.append(block)
    return blocks


def dot_bracket(length: int, pairs) -> str:
    chars = ["."] * length
    for i, j in pairs:
        chars[i], chars[j] = "(", ")"
    return "".join(chars)


@dataclass
class HelixScore:
    """Result of the stem scoring of one weight matrix."""

    score: float                      # total weight of stems with >= min_stem pairs (the RLMHelix score)
    structure: str                    # dot-bracket of the scored stems
    pairs: list[tuple[int, int]] = field(repr=False)  # scored pairs (0-based)
    n_pairs: int = 0                  # decoded pairs after removing isolated pairs
    n_stems: int = 0                  # stems of any length
    n_stems_scored: int = 0           # stems with >= min_stem pairs
    longest_stem: int = 0


def helix_score(weights: np.ndarray, sequence: str, min_stem: int = MIN_STEM) -> HelixScore:
    """Decode `weights` and sum the weights of stems with at least `min_stem` consecutive pairs."""
    weights = np.asarray(weights, dtype=np.float32)
    decoded = decode_pairs(weights, sequence)
    blocks = stems(decoded)
    total, scored = 0.0, []
    for block in blocks:
        if len(block) >= min_stem:
            total += float(np.asarray([weights[i, j] for i, j in block], dtype=np.float32).sum())
            scored.extend(block)
    return HelixScore(score=total, structure=dot_bracket(len(sequence), scored), pairs=sorted(scored),
                      n_pairs=len(decoded), n_stems=len(blocks),
                      n_stems_scored=sum(len(b) >= min_stem for b in blocks),
                      longest_stem=max((len(b) for b in blocks), default=0))


def pair_weights(apc: np.ndarray, bpp: np.ndarray) -> np.ndarray:
    """W = sqrt(max(A, 0) * P), the weight used by the DP."""
    return np.sqrt(np.clip(np.asarray(apc, dtype=np.float32), 0, None) * np.asarray(bpp, dtype=np.float32)).astype(
        np.float32)


def rlmhelix_score(apc: np.ndarray, sequence: str, bpp: np.ndarray | None = None,
                   min_stem: int = MIN_STEM) -> HelixScore:
    """RLMHelix score of one sequence from its APC-corrected dependency map (student or teacher)."""
    sequence = to_rna(sequence)
    if bpp is None:
        bpp = base_pair_probabilities(sequence)
    return helix_score(pair_weights(apc, bpp), sequence, min_stem)
