"""MSA score (RLMHelix-MSA).

For each 200-column window (step 100) over the reference (RF) columns of an alignment, rows are restricted to
ungapped lengths of 30-300 nt and pruned to six sequences with the algorithm of RNAz's rnazSelectSeqs.pl, always
keeping the first row (the reference, e.g. human). The student's APC-corrected dependencies of each sequence are
mapped to alignment columns and averaged per column pair (pairs covered by fewer than half of the sequences are set
to 0), combined with the RNAalifold base-pair probabilities by the geometric mean, and scored with the same DP and
stem aggregation as single sequences. The MSA score is the maximum over windows.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .score import helix_score

GAPS = set("-.~_")
WINDOW, STEP = 200, 100
MAX_SEQUENCES = 6
LENGTHS = (30, 300)


@dataclass
class Alignment:
    names: list[str]
    rows: list[str]
    rf: str  # reference annotation; columns that are not gaps here are scored

    @property
    def rf_columns(self) -> list[int]:
        return [i for i, char in enumerate(self.rf) if char not in GAPS]


def read_stockholm(path: str | Path) -> Alignment:
    """Read one Stockholm alignment. Without a '#=GC RF' line, the residues of the first row define RF."""
    seqs: OrderedDict[str, str] = OrderedDict()
    rf = ""
    for raw in Path(path).open(encoding="latin-1"):
        line = raw.strip()
        if line == "//":
            break
        if line.startswith("#=GC RF"):
            rf += line.split(maxsplit=2)[2]
        elif line and not line.startswith("#"):
            name, seq = line.split()[:2]
            seqs[name] = seqs.get(name, "") + seq
    names, rows = list(seqs), list(seqs.values())
    if not rf:
        rf = "".join("." if char in GAPS else "x" for char in rows[0])
    return Alignment(names, rows, rf)


def _identity(rows: list[str]) -> np.ndarray:
    """Pairwise identity over columns where neither row has a gap."""
    n = len(rows)
    out = np.ones((n, n), dtype=np.float32)
    for i in range(n):
        for j in range(i + 1, n):
            compared = matches = 0
            for a, b in zip(rows[i], rows[j]):
                if a in "-.~" or b in "-.~":
                    continue
                compared += 1
                matches += a.upper() == b.upper()
            out[i, j] = out[j, i] = matches / compared if compared else 0.0
    return out


def prune(rows: list[str], n_target: int = MAX_SEQUENCES, max_id: float = 0.99, opt_id: float = 0.80) -> list[int]:
    """Indices kept by RNAz-style pruning; row 0 (the reference) is always kept."""
    identity = _identity(rows)
    kept = list(range(len(rows)))
    while len(kept) > 1:
        pairs = [(i, j) for i in kept for j in kept if i < j and identity[i, j] > max_id]
        if not pairs:
            break
        candidates = [x for x in pairs[0] if x != 0]
        means = {x: np.mean([identity[x, y] for y in kept if y != x]) for x in candidates}
        kept.remove(max(means, key=means.get))
    while len(kept) > n_target:
        costs = {i: sum((identity[i, j] - opt_id) ** 2 for j in kept if j != i) for i in kept if i != 0}
        if not costs:
            break
        kept.remove(max(costs, key=costs.get))
    return sorted(kept)


def window_starts(length: int, window: int = WINDOW, step: int = STEP) -> list[int]:
    """Window starts, with a final window anchored at the end."""
    if length <= window:
        return [0]
    starts = list(range(0, length - window + 1, step))
    if starts[-1] != length - window:
        starts.append(length - window)
    return starts


def _consensus(rows: list[str]) -> str:
    out = []
    for column in zip(*rows):
        bases = [c.upper().replace("T", "U") for c in column if c not in GAPS]
        counts = {base: bases.count(base) for base in "ACGU"}
        out.append(max(counts, key=counts.get))
    return "".join(out)


def _alifold(rows: list[str]) -> tuple[np.ndarray, float]:
    import RNA

    fold = RNA.fold_compound([r.upper().replace("T", "U").replace("~", "-").replace("_", "-") for r in rows])
    _, mfe = fold.mfe()
    fold.pf()
    packed = fold.bpp()
    length = len(rows[0])
    bpp = np.zeros((length, length), dtype=np.float32)
    for i in range(1, len(packed)):
        for j in range(i + 1, len(packed[i])):
            bpp[i - 1, j - 1] = bpp[j - 1, i - 1] = packed[i][j]
    return bpp, float(mfe)


def score_alignment(alignment: Alignment, student, batch_size: int = 16) -> dict[str, object]:
    """RLMHelix-MSA and its components for one alignment (maximum over W200/S100 windows of RF columns).

    Returns `score` (RLMHelix-MSA), `score_dependency_only`, `score_bpp_only`, `neg_alifold_mfe` and the
    structure of the best window (dot-bracket over that window's RF columns, with its start column).
    """
    rf_cols = alignment.rf_columns
    windows = []
    for start in window_starts(len(rf_cols)):
        columns = rf_cols[start:start + WINDOW]
        rows, sequences, positions = [], [], []
        for row in alignment.rows:
            window_row = "".join(row[c] for c in columns)
            sequence = "".join(x for x in window_row.upper().replace("T", "U") if x not in GAPS)
            if not LENGTHS[0] <= len(sequence) <= LENGTHS[1]:
                continue
            index, mapping = 0, []
            for char in window_row:
                mapping.append(index if char not in GAPS else -1)
                index += char not in GAPS
            rows.append(window_row), sequences.append(sequence), positions.append(mapping)
        if len(sequences) > MAX_SEQUENCES:
            keep = prune(rows)
            rows, sequences, positions = ([x[i] for i in keep] for x in (rows, sequences, positions))
        if len(sequences) < 2 or not columns:
            continue
        total = np.zeros((len(columns), len(columns)), dtype=np.float64)
        coverage = np.zeros_like(total)
        for dep, mapping in zip(student.predict(sequences, batch_size, per_sequence_head=True), positions):
            mapping = np.asarray(mapping)
            valid = np.flatnonzero(mapping >= 0)
            total[np.ix_(valid, valid)] += dep[np.ix_(mapping[valid], mapping[valid])]
            coverage[np.ix_(valid, valid)] += 1
        mean = np.divide(total, coverage, out=np.zeros_like(total), where=coverage > 0)
        mean[coverage < np.ceil(0.5 * len(sequences))] = 0
        mean = mean.astype(np.float32)
        consensus = _consensus(rows)
        bpp, mfe = _alifold(rows)
        combined = helix_score(np.sqrt(np.maximum(mean, 0) * bpp), consensus)
        windows.append({"start": start, "score": combined.score, "structure": combined.structure,
                        "score_dependency_only": helix_score(mean, consensus).score,
                        "score_bpp_only": helix_score(bpp, consensus).score,
                        "neg_alifold_mfe": -mfe, "n_sequences": len(sequences)})
    if not windows:
        return {"status": "insufficient", "score": float("nan")}
    best = max(windows, key=lambda w: w["score"])
    return {"status": "ok", "score": best["score"],
            "score_dependency_only": max(w["score_dependency_only"] for w in windows),
            "score_bpp_only": max(w["score_bpp_only"] for w in windows),
            "neg_alifold_mfe": max(w["neg_alifold_mfe"] for w in windows),
            "n_windows": len(windows), "best_window_start": best["start"],
            "best_window_structure": best["structure"], "n_sequences": best["n_sequences"]}
