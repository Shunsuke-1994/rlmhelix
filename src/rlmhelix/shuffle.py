"""Dinucleotide-preserving shuffles.

Uniform random sequences with the same dinucleotide counts (and the same first and last nucleotide) as the input,
by the Euler-path method of Altschul & Erickson (1985) as formulated by Kandel et al. (1996): a random last exit
edge is chosen for every nucleotide so that the last edges form a tree towards the final nucleotide, the other exit
edges are permuted, and the shuffled sequence is read off the resulting Euler path. The paper used uShuffle (k = 2),
which samples the same distribution.
"""
from __future__ import annotations

import random


def _reaches(last_edge: dict[str, str], end: str) -> bool:
    for start in last_edge:
        seen, node = set(), start
        while node != end:
            if node in seen:
                return False
            seen.add(node)
            node = last_edge[node]
    return True


def dinucleotide_shuffle(sequence: str, rng: random.Random | None = None) -> str:
    """One dinucleotide-preserving shuffle of `sequence`."""
    rng = rng or random.Random()
    if len(sequence) < 3:
        return sequence
    edges: dict[str, list[str]] = {}
    for a, b in zip(sequence, sequence[1:]):
        edges.setdefault(a, []).append(b)
    end = sequence[-1]
    while True:
        last_edge = {node: rng.choice(targets) for node, targets in edges.items() if node != end}
        if _reaches(last_edge, end):
            break
    order: dict[str, list[str]] = {}
    for node, targets in edges.items():
        rest = list(targets)
        if node != end:
            rest.remove(last_edge[node])
        rng.shuffle(rest)
        if node != end:
            rest.append(last_edge[node])
        order[node] = rest
    out, node = [sequence[0]], sequence[0]
    for _ in range(len(sequence) - 1):
        node = order[node].pop(0)
        out.append(node)
    return "".join(out)


def dinucleotide_shuffles(sequence: str, n: int, seed: int | None = None) -> list[str]:
    rng = random.Random(seed)
    return [dinucleotide_shuffle(sequence, rng) for _ in range(n)]
