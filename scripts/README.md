# Scripts

- `download_data.py` downloads the public inputs of the pipelines (GENCODE 50, GRCh38, UCSC rmsk, the hg38 100-way
  tree, Rfam 15.1, the 55 Ensembl release-116 genomes and optional groups); see `data/README.md`. It needs pandas,
  and openpyxl to read Supplementary Table S8 as `.xlsx` (`pip install -e ".[test,data]"`).

The pipelines of the paper are sub-commands of the package's entry point `rlmhelix`:

| Pipeline | Command | Documentation |
|---|---|---|
| genome scan | `rlmhelix scan prepare \| score \| fdr` | `docs/pipeline_scan.md` |
| candidate selection | `rlmhelix select cores \| chunk \| search \| routes \| seed \| expand \| homologues \| families` | `docs/pipeline_select.md` |
| catalogue | `rlmhelix catalogue candidates \| icshape \| covariation \| build \| concordance` | `catalogue/README.md`, `reproduce/README.md` |
