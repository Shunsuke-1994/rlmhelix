# Tests

    pip install -e ".[test]"
    pytest                 # all tests; the model tests are skipped without a GPU
    pytest -m "not model"  # everything except the tests with the language model
    pytest -m model        # the tests with RiNALMo-giga (downloaded on first use)

The reference values come from the runs of the paper: four chr22 windows of the genome scan (`data/scan_windows.tsv`:
RLMHelix score, number of scored stems, calibrated P value), a student map of the first window with its score
(`data/window_apc.npz`), and the Rfam SECIS_1 alignment (RF00031) with its RLMHelix-MSA, dependency-only, BPP-only
and RNAalifold −MFE scores (`data/RF00031.sto`, `data/msa_expected.json`). The RF00031 scores come from the paper's
RLMHelix-MSA run with the settings of Fig. 3; RF00031 is not one of the 322 benchmark families of Fig. 3, so these
values are not in Fig. 3 or Supplementary Table S2.

The pipeline tests use small synthetic inputs and fake external programs; the few tests that run tantan, IRF or
R-scape are skipped when the program is not on `PATH`, and the tests that read or write parquet files need pyarrow.
