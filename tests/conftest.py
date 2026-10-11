from pathlib import Path

import pandas as pd
import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def windows() -> pd.DataFrame:
    """Four genome-scan windows (chr22) with the score, stem count and P value of the paper's scan."""
    return pd.read_csv(DATA / "scan_windows.tsv", sep="\t")


@pytest.fixture(scope="session")
def student():
    torch = pytest.importorskip("torch")
    if not torch.cuda.is_available():
        pytest.skip("model tests compare bf16 GPU values with the paper; no GPU")
    from rlmhelix import Student

    return Student(device="cuda")
