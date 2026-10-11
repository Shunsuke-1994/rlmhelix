# Reproducing the paper

This file maps each figure and table of the paper to what regenerates it with this repository. The repository
reproduces the computational core of the paper: the RLMHelix and RLMHelix-MSA scores (`rlmhelix score`,
`rlmhelix msa`), the genome scan (`rlmhelix scan`), candidate selection and alignments (`rlmhelix select`) and the
RLM-Fam catalogue with its evidence tiers (`rlmhelix catalogue`). Several downstream analyses (benchmark baselines,
comparisons with earlier screens, eCLIP, the pull-down statistics) and the figure assembly are not part of the
repository; they are listed item by item below and collected under [Not included](#not-included). Each step's
arguments are listed by `rlmhelix <step> --help`.

## Hardware and run time

GPU steps were run on one NVIDIA GB10 with bf16 autocast and batches of 16 sequences (the defaults; other batch sizes
change scores by about 0.1%, see the top-level README). RNAfold and the dynamic programme run on CPU workers.

| Step | Measured time |
|---|---|
| Teacher dependency map (`rlmhelix.dependency`), 200 nt | 7.4 s per sequence |
| Student map (`rlmhelix.Student`), 200 nt | 0.037 s one at a time, 0.017 s in batches of 64 |
| DP and stem scoring (CPU), 200 nt | 0.033 s per sequence |
| `examples/` (4 RNAs, 100 shuffles each) | about 45 s on one GPU |
| Genome scan, 7,869,821 windows (`rlmhelix scan score`) | about 2 days on one GPU (an estimated 680 GPU-days with the teacher) |
| RLMHelix-MSA on the homologue alignments of 14,565 candidates and 19 shuffles each (Fig. 5c) | about 19 h on one GPU per alignment condition |

Run times of the homology searches and CM steps are in `docs/pipeline_select.md` (Run times in the paper), and those
of the scan steps in `docs/pipeline_scan.md`; those of R-scape and SISSIz were not recorded.

## External tools and data

The scoring package needs only the Python dependencies in `environments/environment.yml`. The pipelines also call
external programs, which must be on `PATH`:

| Tool | Version in the paper | Used for |
|---|---|---|
| ViennaRNA (RNAfold, RNAalifold) | 2.7 | all scores (Python bindings; installed with the package) |
| tantan, Inverted Repeats Finder (IRF) | 51, 3.09 | repeat filtering in `rlmhelix scan prepare` (on aarch64 build IRF with `make CFLAGS="-O2 -fsigned-char"`; see `docs/pipeline_scan.md`) |
| Riboseek | commit `eba085cb` | homology search in `rlmhelix select search` |
| Infernal | 1.1.5 | seed CMs, CM alignment and family scores in `rlmhelix select seed`, `select expand` and `select families` |
| MAFFT | 7.526 | structure-free alignments in `rlmhelix catalogue covariation` |
| R-scape | 2.6.16 | covariation in `rlmhelix catalogue covariation` |
| SISSIz | 0.1 | alignment shuffles (Fig. 3, Fig. 5c, Fig. 6f, Supplementary Fig. S3) |
| RNAz | 2.1 | RNAz baselines (Figs. 2, 3, S3) |
| bedtools | 2.31.1 | random placement (Fig. 4d, g, h; Supplementary Figs. S5, S7) |
| ScanFold 2.0 | published models, commit `5eec7c20` | ScanFold z baselines (Figs. 2, 4e, f) |
| RNAPhylo (CMfinder 0.4.1.18), EvoFold (EvoFoldV2 of the phy library) | commit `10ee240d` of phy | phylogenetic baselines (Fig. 3, Supplementary Table S2) |
| R2R | bundled with R-scape 2.6.16 | drawings of Fig. 5e–i |
| DIA-NN, proDA (R) | 2.0.2, 1.12.0 | pull-down quantification and binder calls (Fig. 6d–f, Supplementary Figs. S12, S13) |

Public inputs are fetched with `scripts/download_data.py` (see `data/README.md`; resumable, `--dry-run` lists the
URLs):

```bash
python scripts/download_data.py DATA --dry-run                # list the URLs without downloading
python scripts/download_data.py DATA                          # genome, rmsk, tree, rfam, ensembl (about 37 GB):
                                                              # GENCODE 50, GRCh38, UCSC rmsk, hg38.100way.nh,
                                                              # Rfam 15.1, the 55 Ensembl release-116 genomes
python scripts/download_data.py DATA --only icshape           # GSE145805 and GSE74353 tracks, GENCODE v26, Ensembl 74
python scripts/download_data.py DATA --only screens bprna polyasite   # earlier screens, bpRNA-new, PolyASite 2.0
python scripts/download_data.py DATA --only multiz taxonomy   # multiz100way (zarr, 39 GB) and NCBI taxonomy, only to
                                                              # rebuild the MSA benchmarks
python scripts/download_data.py DATA --only encode --encode-table Supplementary_Table_S8.xlsx \
    --encode-dataset "eCLIP peaks"                            # ENCODE files listed in Supplementary Table S8
python scripts/download_data.py DATA --list-manual            # EvoFam and the other inputs obtained by hand
```

Reading the `.xlsx` of Supplementary Table S8 (here and in `rlmhelix catalogue build --eclip`) needs openpyxl, which
is in the conda environment; with pip, install the `data` extra (`pip install -e ".[test,data]"`).

**Data to be deposited.** The RLM-Fam catalogue (candidate table, alignments, covariance models, provisional
families) and all 7,869,821 scan windows (coordinates, RLMHelix scores, calibrated P values, BH q values; about
280 MB) will be deposited on Zenodo (DOI to be added; public upon publication of the paper); `catalogue/README.md`
describes the files. The raw data and the protein identification and quantification of the pull-down mass
spectrometry will be deposited to ProteomeXchange via PRIDE and released upon publication (accession to be added).
Until then, the scan-dependent items below need the output of `rlmhelix scan`; once the deposit is public, they can
start from its window table instead of rerunning the two-day scan.

## Pipeline from genome to catalogue

```bash
# genome scan (Fig. 4; Supplementary Figs. S4, S8)
rlmhelix scan prepare ...    # 200-nt windows (step 100) of GENCODE v50 introns, UTRs and lncRNAs; RepeatMasker, tantan, IRF40 filters
rlmhelix scan score ...      # student + RNAfold scores and calibrated P values (one GPU, about 2 days)
rlmhelix scan fdr ...        # Benjamini-Hochberg q, Efron empirical null (mu0, sigma0, pi0), tail-area FDR and lfdr

# candidate selection (Fig. 4a)
rlmhelix select cores ...      # cores of the top 1% of windows (78,699 windows -> 71,196 representatives; GPU)
rlmhelix select chunk ...      # 10-kb chunks of the 55 genomes
rlmhelix select search ...     # Riboseek search in the 20-species panel (and later of the source strand and the 35 expansion genomes)
rlmhelix select routes ...     # four selection routes (15,913 cores)
rlmhelix select seed ...       # seed alignments and CMs (15,821 candidates)
rlmhelix select expand ...     # expansion to 55 species (final55 alignments and CMs) and the condition (iii) homologues
rlmhelix select families ...   # provisional families (14,037; RLMFam00001, ...)
# exact arguments of every step: docs/pipeline_scan.md and docs/pipeline_select.md

# catalogue (Fig. 5a, d, j; Supplementary Fig. S8; Supplementary Table S4); DATA as in download_data.py above
rlmhelix catalogue candidates select/expand/candidates.parquet --families select/families/families.tsv \
    -o candidates.tsv                                             # candidate table of the catalogue steps
I=DATA/icshape
rlmhelix catalogue icshape --candidates candidates.tsv --gtf $I/gencode26/gencode.v26.annotation.gtf.gz \
    --transcripts $I/gencode26/gencode.v26.transcripts.fa.gz \
    --track HEK293=$I/GSE145805/GSM4333256_293.out.txt.gz --track HeLa=$I/GSE145805/GSM4333258_HeLa.out.txt.gz \
    --track K562=$I/GSE145805/GSM4333260_K562.out.txt.gz --track HepG2=$I/GSE145805/GSM4333262_HepG2.out.txt.gz \
    --track H9=$I/GSE145805/GSM4333264_H9.out.txt.gz -o icshape.tsv
rlmhelix catalogue covariation --candidates candidates.tsv --homologues select/expand/homologues \
    --workdir rscape -o covariation_iii.tsv                       # condition (iii), the primary condition
rlmhelix catalogue covariation --candidates candidates.tsv --condition ii --alignments select/expand/candidates \
    --workdir rscape -o covariation_ii.tsv                        # same sequences realigned with MAFFT
rlmhelix catalogue covariation --candidates candidates.tsv --condition i --alignments select/expand/candidates \
    --workdir rscape -o covariation_i.tsv                         # the covariance-model alignments
rlmhelix catalogue build --candidates candidates.tsv --windows scan/eligible_windows.parquet --icshape icshape.tsv \
    --covariation covariation_iii.tsv covariation_ii.tsv covariation_i.tsv \
    [--eclip Supplementary_Table_S8.xlsx --eclip-dir DATA/encode] [--gtf DATA/gencode50/gencode.v50.annotation.gtf.gz] \
    [--rfam DATA/rfam/rfam15.1_human_GRCh38_hits.tsv] [--homologue-msa homologue_msa.tsv] \
    -o rlmfam_candidates.tsv --families rlmfam_families.tsv --tier-summary tier_summary.tsv \
    --evidence-summary evidence_summary.tsv [--lfdr-curve lfdr_curve.tsv]
rlmhelix catalogue concordance rlmfam_candidates.tsv -o concordance.tsv   # Supplementary Fig. S8c
```

`catalogue build` takes the lfdr from the window table of `rlmhelix scan fdr` (`scan/eligible_windows.parquet`), or
computes it with Lindsey's method under the paper's empirical null when the column is absent. `--lfdr-curve` writes the
curve of Supplementary Fig. S8a, b (columns `z`, `f`, `pi0_f0`, `lfdr`); the lfdr is then recomputed from the P values,
also when the table has an `lfdr` column. The window table must hold all eligible windows, because the lfdr, the BH q
values (when absent) and the top-1% cutoff of `n_top1pct_windows` are computed from it; for the paper's windows the
cutoff is 0.020876023786561414 (`--top1-p P` sets it, for example for a table that is not the full set). `catalogue
covariation` runs one candidate per process (`--jobs`), resumes from finished outputs, reports candidates without an
input file on standard error and stops when none has one.

## Main figures

**Fig. 1.** Panels a–c are the authors' concept diagram. Panels d and e (teacher versus student run time) can be
measured with the Python API: `rlmhelix.dependency.teacher_apc` (3L + 1 passes) against `rlmhelix.Student.predict`,
with `rlmhelix.score.rlmhelix_score` for the DP, on 20 chr22 scan windows per length (50, 100, 200, 300 nt). The
timing script is not included; the measured values are in the legend and in the table above.

**Fig. 2 (328 Rfam families, single sequences).** The human representative of each family is listed in
Supplementary Table S1 (sheet Rfam_328_families).
- RLMHelix (b–f): `rlmhelix score representatives.fa --shuffles 100 --seed 0 -o rfam328.tsv`. With one positive per
  family, the AUPRC of b is 1/(101 × `empirical_p`). For the pooled curves (c–e), generate the shuffles with
  `rlmhelix.dinucleotide_shuffles` and score them with `rlmhelix score`. The paper used uShuffle (k = 2), which samples
  the same distribution with different draws, so values agree within shuffle noise rather than digit for digit.
- Dependency-only and BPP-only ablations: `rlmhelix.helix_score` with the student map or
  `rlmhelix.base_pair_probabilities` as weights (Python API). RNAfold −MFE: ViennaRNA.
- ScanFold 2.0 z (mononucleotide and dinucleotide) and RNAz 2.1 z: external programs; the drivers are not included.
- Run time (estimate from the per-sequence times above): under an hour on one GPU for the 33,128 sequences.

**Fig. 3 (322 Rfam families, MSAs).** `rlmhelix msa family.sto ... -o msa.tsv` gives RLMHelix-MSA, dependency only, BPP
only and RNAalifold −MFE with the settings of the paper (human row first). `tests/` checks one Rfam alignment, RF00031
(SECIS_1), against its scores from the paper's RLMHelix-MSA run with the Fig. 3 settings; RF00031 is not one of the 322
benchmark families, so these values are not in Fig. 3 or Supplementary Table S2. Negatives are 100 SISSIz shuffles per
alignment (`SISSIz --simulate -d -n 100`). RNAz, RNAPhylo and EvoFold are external, and their drivers are not included;
their values are in Supplementary Table S2. The construction of the alignments (vertebrate rows of the Rfam 15.1 seed
alignments mapped to the assemblies of the UCSC hg38 100-way tree) is not included.

**Fig. 4 (genome scan and selection).**
- a: counts from the pipeline above (`rlmhelix scan`, `rlmhelix select`).
- b: GC content of the top windows, computable from `scan/eligible_windows.parquet` (column `GC`); the script is not
  included.
- c: `rlmhelix scan fdr` (empirical null and tail-area FDR by top fraction; Supplementary Fig. S4).
- d (Supplementary Fig. S5): overlap with Rfam 15.1 hits at each stage against random placement (bedtools); not
  included.
- e, f: RNAfold MFE and ScanFold 2.0 z on all windows with a block bootstrap; not included.
- g, h (Supplementary Fig. S7): overlap with RNAz 2005, EvoFold v1 and v2, CMfinder CRSs and PhastCons (lifted to
  hg38, bedtools shuffle); not included.

**Fig. 5 (structural support).**
- a: `rlmhelix catalogue icshape` (GSE145805, five cell lines; AUROC of reactivity for unpaired versus paired
  nucleotides with the circular-shift control). In the paper, 757 of the 71,196 top-1% cores from `rlmhelix select
  cores` were measured (`--candidates select/cores/representatives.parquet`), 474 of them RLM-Fam candidates
  (`--candidates candidates.tsv`). The Wilcoxon test, the RNAfold MFE and BPP-only structures and the grouping by
  earlier screens are not included.
- b: GSE74353 (HEK293T) mapped to Ensembl release 74 transcripts; not included.
- c: structure-free homologue alignments (written by `rlmhelix catalogue covariation` in its work directory) scored with
  `rlmhelix msa` with the human row removed, against 19 SISSIz shuffles each (about 19 h on one GPU). The driver
  (shuffling, human-row removal, Storey's estimate) and the 600 candidate–control pairs are not included.
- d: `rlmhelix catalogue covariation` (condition iii, the primary condition: MAFFT `--auto` alignment of the
  high-confidence homologues, predicted structure projected through the human row, R-scape `-s --nofigures`,
  E < 0.05). `--condition ii --alignments select/expand/candidates` realigns the sequences of the final55 alignments
  with MAFFT, and `--condition i --alignments select/expand/candidates` runs R-scape on the final55 alignments
  themselves (`<core_id>/final55.sto`; a folder of `<core_id>.sto` files, as planned for the Zenodo deposit, also
  works). The SISSIz background (dashed lines) is not included.
- e–i: R2R drawings; not included (coordinates in Supplementary Table S4).
- j: `rlmhelix catalogue build` (tiers A/B/C: 380, 4,631 and 10,810 candidates; tier and evidence summaries).

**Fig. 6 (protein binding).**
- a–c: eCLIP enrichment against matched within-gene controls, dose–response, candidate-specific partners and the
  genome-browser views; not included. `rlmhelix catalogue build --eclip Supplementary_Table_S8.xlsx --eclip-dir
  DATA/encode` reproduces only the per-candidate eCLIP overlap columns of Supplementary Table S4 (from the 250 rows of
  dataset `eCLIP peaks` of the Accessions sheet and their IDR peak files `<file_accession>.bed.gz`).
- d–f: bait scores with `rlmhelix score --shuffles 100` (single sequences) and `rlmhelix msa` against 100 SISSIz
  shuffles (alignments); for baits shorter than 30 nt, the teacher (`rlmhelix.dependency`). The bait sequences are in
  Supplementary Table S7; the mass-spectrometry data will be deposited to ProteomeXchange via PRIDE and released upon
  publication (accession to be added). Binder calls (DIA-NN, proDA), ROC, DeLong and correlation analyses are not
  included.

## Supplementary figures

| Figure | Reproduce with | Not included |
|---|---|---|
| S1 bpRNA base pairs | Teacher maps (`rlmhelix.dependency.teacher_apc`), student maps (`Student.predict`) and decoding (`rlmhelix.score.decode_pairs` with `pair_weights`) on the 5,313 bpRNA-new sequences of 30–300 nt | F1 evaluation and bootstrap |
| S2 RNA classes | As Fig. 2, grouped by Rfam family type | Class summaries; baselines as Fig. 2 |
| S3 EvoFam and RNAz 2005 MSAs | `rlmhelix msa` on the alignments, against 100 SISSIz shuffles (loci and regions in Supplementary Table S1) | Alignment construction from multiz100way; RNAz; evaluation |
| S4 Empirical null | `rlmhelix scan fdr` | Plotting |
| S5 Stage-wise Rfam enrichment | — | Whole analysis (as Fig. 4d) |
| S6 Enrichment by RNA class | — | Whole analysis (known structures of Supplementary Table S3, Rfam hits, same-gene background) |
| S7 Earlier screens | — | Whole analysis (as Fig. 4g, h) |
| S8 lfdr and covariation | a, b: lfdr curve (Lindsey's method under the empirical null) from `rlmhelix catalogue build --windows scan/eligible_windows.parquet --lfdr-curve lfdr_curve.tsv` (columns `z`, `f`, `pi0_f0`, `lfdr`; `scan fdr` writes the lfdr of each window but not the curve); c: `rlmhelix catalogue concordance` on the catalogue (rates by power class, Wilson intervals, Fisher's exact tests) | Plotting |
| S9 eCLIP design | — | Whole analysis (as Fig. 6a, b) |
| S10 Partner properties | — | Whole analysis (RBNS motifs, chance ratio) |
| S11 eCLIP and icSHAPE examples | Per-cell-line AUROC and shift control from `rlmhelix catalogue icshape` | GSE74353; genome-browser views |
| S12, S13 Pull-down | Bait scores as Fig. 6f | Binder calls, statistics and plots |

## Supplementary tables

| Table | Reproduce with | Not included |
|---|---|---|
| S1 Evaluation sets and denominators | Scan sheets from the scan windows (`rlmhelix scan`, or the Zenodo window table once public); the Rfam representatives, EvoFam loci and RNAz 2005 regions are listed in the table | Construction of the benchmark sets |
| S2 Phylogenetic comparison | RLMHelix-MSA, dependency only, BPP only and RNAalifold −MFE with `rlmhelix msa` | RNAz, RNAPhylo, EvoFold; evaluation and bootstraps |
| S3 Known functional structures | — (literature curation; elements and PMIDs listed in the table) | Recovery statistics and background intervals |
| S4 RLM-Fam catalogue | `rlmhelix catalogue build` writes the Catalogue and Families sheets and the tier and evidence summaries, from the outputs of `rlmhelix scan`, `rlmhelix select` (through `catalogue candidates`), `catalogue icshape` and `catalogue covariation`; optional columns need ENCODE eCLIP peaks (`--eclip`), the GENCODE v50 GTF (`--gtf`: gene annotation, pseudogene and splice-site columns, and the Families columns `genes` and `region_contexts`), Rfam 15.1 hits (`--rfam`) and the homologue-MSA table (`--homologue-msa`, SISSIz). To be deposited on Zenodo as `rlmfam_candidates.tsv` and `families.tsv` | Overlap with earlier screens (column `n_overlapping_earlier_screens` of the Families sheet); the homologue-MSA driver |
| S5 eCLIP per experiment | — | Whole analysis |
| S6 Implications | — (literature curation) | — |
| S7 Pull-down baits | Bait scores (P and z) with `rlmhelix score` and `rlmhelix msa` | Binder calls and labels |
| S8 Data accessions | `scripts/download_data.py` (the Genomes sheet lists the panel of `data/genome_panel_55.tsv`; GENCODE, Rfam, icSHAPE, PolyASite; ENCODE files with `--only encode --encode-table`) | — |

## Not included

The repository does not contain code for the following items. Their values are in the paper and its supplementary
tables as indicated.

- Student training and its synthetic training data, and fitting of the score calibrator (the trained head and the
  calibrator are included in `src/rlmhelix/models/`).
- Fig. 1d, e: the timing script (values in the legend of Fig. 1).
- Fig. 2, 3, Supplementary Figs. S2, S3: construction of the benchmark sets (Rfam 15.1 vertebrate alignments mapped to
  the hg38 100-way assemblies, EvoFam and RNAz 2005 alignments cut from multiz100way) and of their uShuffle and SISSIz
  negatives (sets listed in Supplementary Table S1).
- Fig. 2, 3, Supplementary Figs. S2, S3: drivers of the baselines (ScanFold 2.0 z, RNAz 2.1, RNAPhylo, EvoFold) and the
  evaluation code (AUPRC, pooled AP, paired bootstraps; values in Supplementary Table S2).
- Supplementary Fig. S1: base-pair F1 evaluation on bpRNA-new.
- Fig. 4b: GC content of the top windows.
- Fig. 4d and Supplementary Fig. S5: stage-wise Rfam enrichment.
- Fig. 4e, f: thermodynamic scanning of all windows (RNAfold MFE, ScanFold 2.0) and the block bootstrap.
- Fig. 4g, h, Fig. 5a (right), Supplementary Fig. S7 and the Families-sheet column `n_overlapping_earlier_screens` of
  Supplementary Table S4: comparison with the four earlier screens and PhastCons.
- Supplementary Fig. S6 and Supplementary Table S3: recovery of known functional structures and enrichment by RNA
  class (elements in Supplementary Table S3).
- Fig. 5a, b and Supplementary Fig. S11: icSHAPE comparisons with RNAfold MFE and BPP-only structures and the GSE74353
  analysis.
- Fig. 5c: shuffling and scoring driver for the homologue alignments, Storey's estimate and the 600 matched controls
  (per-candidate values in Supplementary Table S4).
- Fig. 5d: the SISSIz background of the covariation calls.
- Fig. 5e–i: R2R drawings.
- Fig. 6a–c, Supplementary Figs. S9–S11 and Supplementary Table S5: eCLIP enrichment, matched controls, dose–response,
  candidate-specific partners, partner properties and the genome-browser views (values in Supplementary Table S5;
  per-candidate binding in Supplementary Table S4).
- Fig. 6d–f, Supplementary Figs. S12, S13 and Supplementary Table S7: binder calls (proDA) and the pull-down statistics
  (values in Supplementary Table S7; the mass-spectrometry data will be deposited to ProteomeXchange via PRIDE and
  released upon publication, accession to be added).
- Supplementary Table S6: literature curation (no code).
- Assembly of the figures and of the supplementary workbooks.
