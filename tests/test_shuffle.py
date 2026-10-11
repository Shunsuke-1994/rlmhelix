"""Dinucleotide-preserving shuffles."""
import collections
import itertools
import random

from rlmhelix.shuffle import dinucleotide_shuffle, dinucleotide_shuffles


def doublets(sequence):
    return collections.Counter(zip(sequence, sequence[1:]))


def test_counts_and_ends_are_kept():
    sequence = "GGCCGGUUAGCUCAGUUGGUUAGAGCGUGGUGCUAAUAACGCCAAGGUCGCGGGUUCGAUCCCCGUACGGGCCA"
    for shuffled in dinucleotide_shuffles(sequence, 100, seed=1):
        assert doublets(shuffled) == doublets(sequence)
        assert (shuffled[0], shuffled[-1]) == (sequence[0], sequence[-1])


def test_all_shuffles_are_reached_about_equally():
    sequence = "AACAGACA"  # 8! orderings to enumerate; 12 of them share the doublet counts
    space = {"".join(p) for p in itertools.permutations(sequence) if doublets("".join(p)) == doublets(sequence)}
    counts = collections.Counter(dinucleotide_shuffle(sequence, random.Random(i)) for i in range(24000))
    assert set(counts) == space
    assert max(counts.values()) / min(counts.values()) < 1.2
