# Environment

`environment.yml` describes the environment used for the paper (conda-forge and bioconda for ViennaRNA, pip for the
language-model stack and for pyarrow 26.0.0, which conda-forge does not have). It was rebuilt from scratch and tested on
2026-10-11 (Linux aarch64, NVIDIA GB10, CUDA 13): all tests passed, including the model tests (tests that need one of
the external programs below are skipped when it is not on `PATH`):

```bash
conda env create -f environments/environment.yml    # or mamba / micromamba
conda activate rlmhelix
pytest -m "not model"                                # tests without the language model
pytest -m model                                      # tests with RiNALMo-giga (downloaded on first use; GPU recommended)
```

`pyproject.toml` pins the language-model stack to the tested releases (`transformers>=5.5,<5.6`,
`multimolecule==0.0.9`): multimolecule 0.2.1 fails to import with transformers 5.19 and changes the RLMHelix-MSA score
with transformers 5.5.4, and multimolecule 0.0.9 fails with transformers 5.9 and later. PyTorch is left at
`torch>=2.11`; `pip install -e ".[test]"` in a clean Python 3.11 venv installed torch 2.14.1 (2026-10-11), which passed
all tests and gave the same scores as torch 2.11.0 for `examples/`. The backbone is loaded at the Hugging Face revision
used for the paper (`2a71f6f98fb41dd2e6542a5e131d3778111d1468`, `rlmhelix.dependency.BACKBONE_REVISION`).

openpyxl is needed only to read Supplementary Table S8 as .xlsx (`scripts/download_data.py --encode-table` and `rlmhelix catalogue build --eclip`). It is in
`environment.yml`; with pip, install the `data` extra (`pip install -e ".[test,data]"`).

The scoring package (`rlmhelix score`, `rlmhelix msa`) needs nothing else. The pipelines of the paper also call
external programs, which must be on `PATH` (or be given with the options of each command):

| Program | Version used | Needed by |
|---|---|---|
| tantan | 51 | `rlmhelix scan prepare` |
| Inverted Repeats Finder (IRF) | 3.09 | `rlmhelix scan prepare` (on aarch64 build with `make CFLAGS="-O2 -fsigned-char"`) |
| Riboseek | commit `eba085cb` (https://github.com/steineggerlab/riboseek) | `rlmhelix select search` |
| Infernal | 1.1.5 | `rlmhelix select seed`, `expand`, `families` |
| MAFFT | 7.526 | `rlmhelix catalogue covariation` |
| R-scape | 2.6.16 | `rlmhelix catalogue covariation` |

The comparisons of the paper also used RNAz 2.1, SISSIz 0.1, bedtools 2.31.1, ScanFold 2.0, RNAPhylo (CMfinder
0.4.1.18) and EvoFold (phy); see `reproduce/README.md`.
