"""RLMHelix: alignment-free discovery of conserved RNA structures from single sequences.

    from rlmhelix import Student, Calibrator, rlmhelix_score

    student = Student()                        # RiNALMo-giga + distilled head (GPU if available)
    apc = student.predict([sequence])[0]       # APC-corrected dependency map
    result = rlmhelix_score(apc, sequence)     # score, structure of the scored stems
    p = Calibrator.load().pvalues([result.score], [sequence])[0]
"""
from .calibrate import Calibrator, sequence_features
from .score import base_pair_probabilities, helix_score, rlmhelix_score
from .shuffle import dinucleotide_shuffles

__version__ = "0.1.0"


def __getattr__(name):  # torch-dependent classes are imported on first use
    if name == "Student":
        from .student import Student
        return Student
    if name in ("Alignment", "read_stockholm", "score_alignment"):
        from . import msa
        return getattr(msa, name)
    raise AttributeError(name)


__all__ = ["Student", "Calibrator", "rlmhelix_score", "helix_score", "base_pair_probabilities",
           "sequence_features", "dinucleotide_shuffles", "read_stockholm", "score_alignment", "Alignment"]
