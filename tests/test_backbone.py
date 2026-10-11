"""The backbone is loaded at the pinned Hugging Face revision of the paper, without the unused pooler."""
import re
import subprocess
import sys

import pytest

from rlmhelix.dependency import BACKBONE, BACKBONE_REVISION


def test_dependency_module_does_not_import_torch():
    code = "import sys, rlmhelix.dependency; assert 'torch' not in sys.modules"
    subprocess.run([sys.executable, "-c", code], check=True)


class _Recorder:
    """Stands in for a multimolecule class and records the from_pretrained calls."""

    calls: list[tuple[str, str, dict]] = []

    def __init__(self, name: str):
        self.name = name

    def from_pretrained(self, model_name: str, **kwargs):
        self.calls.append((self.name, model_name, kwargs))
        return self

    def to(self, device):
        return self

    def eval(self):
        return self

    def get_vocab(self) -> dict[str, int]:
        return {base: i for i, base in enumerate("ACGU", start=6)}


@pytest.fixture
def recorded(monkeypatch):
    pytest.importorskip("torch")
    multimolecule = pytest.importorskip("multimolecule")
    _Recorder.calls = []
    for name in ("RnaTokenizer", "RiNALMoModel", "RiNALMoForMaskedLM"):
        monkeypatch.setattr(multimolecule, name, _Recorder(name))
    return _Recorder.calls


def test_student_pins_revision_and_drops_pooler(recorded):
    from rlmhelix.student import Student

    assert re.fullmatch(r"[0-9a-f]{40}", BACKBONE_REVISION)
    Student(device="cpu")
    assert {name for name, _, _ in recorded} == {"RnaTokenizer", "RiNALMoModel"}
    for name, model_name, kwargs in recorded:
        assert (model_name, kwargs["revision"]) == (BACKBONE, BACKBONE_REVISION), name
    model_kwargs = next(kwargs for name, _, kwargs in recorded if name == "RiNALMoModel")
    assert model_kwargs["add_pooling_layer"] is False


def test_teacher_pins_revision(recorded):
    from rlmhelix.dependency import TeacherModel

    TeacherModel(device="cpu")
    assert {name for name, _, _ in recorded} == {"RnaTokenizer", "RiNALMoForMaskedLM"}
    assert all((model_name, kwargs["revision"]) == (BACKBONE, BACKBONE_REVISION)
               for _, model_name, kwargs in recorded)


@pytest.mark.model
def test_loaded_backbone_revision(student):
    assert student.backbone.config._commit_hash == BACKBONE_REVISION
    assert student.backbone.pooler is None
