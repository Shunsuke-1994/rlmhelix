# Genome scan

This page describes how to rebuild the genome scan of the paper: 200-nt windows over the transcribed non-coding
regions of the human genome, repeat filters, RLMHelix scores and calibrated P values, the Benjamini–Hochberg (BH)
reference set and its loci, and the local false discovery rate (lfdr) from Efron's empirical null. The code is
`rlmhelix.scan`; each step is a sub-command of `rlmhelix scan`. The methods are described in the paper (Methods,
"Genome scan, calibration and false discovery"; Supplementary Methods, "Genome scan", "Score calibration and multiple
testing" and "Empirical null and local false discovery rates").

## Overview

| Step | Command | What it does | Hardware |
|---|---|---|---|
| 1 | `scan prepare` | windows, RepeatMasker, tantan and IRF filters, removal of intergenic windows | CPU, tantan, IRF |
| 2 | `scan score` | student dependency maps, RLMHelix score, calibrated P value | GPU (+ CPU workers) |
| 3 | `scan fdr` | genome-wide BH, loci, empirical null, tail-area Fdr, lfdr | CPU |

**Windows.** Window centres are the multiples of 100 (0-based genome coordinates) that lie in an intron or UTR of a
GENCODE mRNA (gene and transcript type `protein_coding`; introns are the gaps between consecutive exons, UTRs the exon
parts outside the span from the first to the last CDS base, so the stop codon falls in the 3′ UTR) or within a lncRNA
gene (gene type `lncRNA`, whole gene span); each window spans 100 nt on either side of its centre (200 nt, step
100 nt). Intervals are merged per strand and feature over all transcripts. A window is read 5′→3′ on the strand of
its annotation, so a position annotated on both strands gives two windows, and windows centred near an exon may
overlap coding sequence. Windows that run off the chromosome or contain a base other than A, C, G or T are dropped.
`features` lists the labels of a window (`intron`, `utr5`, `utr3`, `lncrna`). Windows lying entirely between genes
(any gene type) are labelled `intergenic`; they are not made unless `--keep-intergenic` is given and are always
excluded by `scan fdr`.

**Repeat filters**, in this order: more than 40 nt (20% of the window) overlapping RepeatMasker (UCSC `rmsk`, all
repeat classes, intervals merged); at least 20 nt masked by tantan (`tantan -f 3`, the union of the masked
intervals); and any inverted repeat reported by IRF with an alignment score of at least 40. IRF is run as
`irf <fasta> 2 7 7 80 10 40 10000 200 -h -ngs -t4 200 -t5 200 -t7 200` (match 2, mismatch 7, indel 7, match and
indel probabilities 80 and 10, minimum score 40, maximum length 10,000, maximum loop 200; a score of 40 corresponds to
a perfect stem of about 20 bp).

**Scores.** The student predicts the APC-corrected dependency map of each window (RiNALMo-giga backbone, bf16
autocast on the GPU, batches of 16 windows); the RLMHelix score combines it with the RNAfold base-pair probabilities
(`rlmhelix.score`), and the calibrated P value is the upper-tail probability of the conditional GEV model
(`rlmhelix.calibrate`).

**False discovery.** BH q values are computed across all eligible windows; windows with q ≤ 0.05 form the reference
set, and accepted windows that overlap on the same strand are merged into loci (windows that only touch are not).
For the empirical null, P values are turned into z = −Φ⁻¹(P) (P clipped to [10⁻³⁰⁰, 1 − 10⁻¹²]); a quadratic fitted
to the log counts of 60 bins between the 25th and 75th percentiles of z gives N(μ0, σ0²), and π0 is the observed
fraction in that range divided by the null mass in it. The tail-area Fdr of the top R windows is
π0 S0(t) N / R, with t the z of the R-th window. The lfdr is π0 f0(z) / f(z), where f is fitted by Lindsey's method:
counts in bins of 0.05 over the 10⁻⁶ to 1 − 10⁻⁶ quantile range of z, Poisson regression on a natural cubic spline
with knots at the two range ends and the 1, 5, 20, 50, 80, 95, 99, 99.5, 99.9 and 99.99% quantiles of z. The lfdr is
capped at 1 and made non-increasing in z to the right of μ0; windows with z < μ0 get lfdr 1 (the left tail of z is
also heavier than normal, but those windows have less structure than the null centre and are not discoveries).

## Software and data

| Input | Version | Source |
|---|---|---|
| Gene annotation | GENCODE release 50, `gencode.v50.annotation.gtf.gz` | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ |
| Genome | GENCODE GRCh38 primary assembly, `GRCh38.primary_assembly.genome.fa.gz` | same directory |
| RepeatMasker | UCSC hg38 `rmsk` table, `rmsk.txt.gz` (downloaded 2026-08-12 for the paper) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz |
| tantan | 51 (bioconda `tantan=51`) | https://gitlab.com/mcfrith/tantan |
| Inverted Repeats Finder (IRF) | 3.09 | https://github.com/Benson-Genomics-Lab/IRF |
| Student, calibration | as installed with this package | `rlmhelix.student`, `rlmhelix.calibrate` |

`python scripts/download_data.py data/raw --only genome rmsk` fetches the annotation, the genome and the RepeatMasker
table (to `data/raw/gencode50/` and `data/raw/ucsc_hg38/`). The genome can be given compressed (it is then read once,
in file order) or decompressed with a samtools index (`samtools faidx`, random access). tantan and IRF are taken
from `PATH` or given with `--tantan` and `--irf`. IRF is built from source with `cd src && make`; on ARM (aarch64)
use `make CFLAGS="-O2 -fsigned-char"`, without which the program is extremely slow (with the flag its output equals
that of the x86-64 release binary used for the paper).

## Commands

The commands below rebuild the scan of the paper. Paths are examples.

```bash
# 1. Windows and repeat filters (CPU; tantan and IRF)
gunzip -k data/raw/gencode50/GRCh38.primary_assembly.genome.fa.gz
samtools faidx data/raw/gencode50/GRCh38.primary_assembly.genome.fa      # optional, for random access
rlmhelix scan prepare \
    --gtf data/raw/gencode50/gencode.v50.annotation.gtf.gz \
    --genome data/raw/gencode50/GRCh38.primary_assembly.genome.fa \
    --rmsk data/raw/ucsc_hg38/rmsk.txt.gz \
    --irf /path/to/IRF/bin/irf.exe --threads 4 \
    -o scan/windows

# 2. Scores and calibrated P values (GPU; bf16, batches of 16 as in the paper)
rlmhelix scan score --windows scan/windows -o scan/scores --batch-size 16 --n-jobs 12 --resume

# 3. BH, loci, empirical null and lfdr (CPU; needs all chromosomes)
rlmhelix scan fdr --scores scan/scores -o scan
```

`--chrom chr21 chr22` restricts `prepare` and `score` to some chromosomes (useful for a test run); `fdr` must see all
24 chromosomes (chr1–chr22, chrX, chrY), because BH and the empirical null are genome-wide. `score --resume` keeps
finished chromosomes and finished shards of 10,000 windows (`--chunk`), so an interrupted run continues where it
stopped; `prepare --resume` skips chromosomes already written. `--skip-tantan` and `--skip-irf` drop a filter (not
the paper's setting). The same steps from Python:

```python
from rlmhelix import Student, scan

annotation = scan.read_annotation("gencode.v50.annotation.gtf.gz", ["chr22"])
repeats = scan.read_repeats("rmsk.txt.gz", ["chr22"])
(chrom, sequence), = scan.iter_chromosomes("GRCh38.primary_assembly.genome.fa", ["chr22"])
windows = scan.make_windows(chrom, sequence, annotation)
windows, counts = scan.filter_repeats(windows, repeats[chrom], irf="/path/to/irf.exe")
scored = scan.score_windows(windows, Student(), batch_size=16)
scored["q_value"] = scan.benjamini_hochberg(scored.p_value)     # genome-wide in the paper
loci = scan.merge_loci(scored[scored.q_value <= 0.05])
z = scan.z_scores(scored.p_value)
null = scan.empirical_null(z)
scored["lfdr"] = scan.local_fdr(z, null)
table = scan.tail_fdr(z, null, lfdr=scored.lfdr)
```

## Outputs

All coordinates are 0-based and half-open on the chromosome (GRCh38); `sequence` is DNA, 5′→3′ on the window's
strand.

**`scan/windows/<chrom>.parquet`** (step 1): the eligible windows of a chromosome, sorted by start and strand.

| Column | Content |
|---|---|
| `window_id` | `chrom:start-end:strand` |
| `chrom`, `start`, `end`, `center`, `strand` | position (`center` = start + 100) |
| `features` | labels joined by `;` (`intron`, `utr5`, `utr3`, `lncrna`; `intergenic` with `--keep-intergenic`) |
| `sequence` | 200 nt |
| `repeatmasker_bp`, `tantan_bp` | bases overlapping RepeatMasker; bases masked by tantan |

**`scan/windows/window_counts.tsv`**: per chromosome, `n_placed` (window positions), `n_non_acgt` (dropped for
other bases or the chromosome end), `n_windows` (windows tested), `n_repeatmasker`, `n_tantan`, `n_irf` (windows
removed by each filter, in this order) and `n_kept`.

**`scan/scores/<chrom>.parquet`** (step 2): the columns above and

| Column | Content |
|---|---|
| `score` | RLMHelix score (total weight of the stems of at least four pairs) |
| `n_pairs`, `n_stems`, `n_stems_scored`, `longest_stem` | decoded pairs after removing isolated pairs; stems; stems of ≥ 4 pairs; longest stem |
| `GC` | GC content |
| `p_value` | calibrated P value (0 above the upper end of the fitted GEV) |

**`scan/eligible_windows.parquet`** (step 3): every eligible window, in chromosome order (chr1–chr22, chrX, chrY)
and then by start and strand, with all columns of step 2 and

| Column | Content |
|---|---|
| `q_value` | BH q value across all eligible windows |
| `z` | −Φ⁻¹(`p_value`) |
| `lfdr` | local false discovery rate (1 for z < μ0) |
| `bh_accepted` | `q_value` ≤ 0.05 (`--q`) |

This is the input of `rlmhelix select cores` (candidate selection) and of the window statistics of `rlmhelix catalogue
build --windows`. In the paper's intermediate files `score`, `p_value` and `n_stems_scored` were called `score_min4`,
`p_min4_mixed` and `n_stems_ge4`; `rlmhelix select` accepts either name.

**Other tables of step 3** (tab-separated):

| File | Content |
|---|---|
| `bh_windows.tsv` | the accepted windows (columns of `eligible_windows.parquet`) |
| `bh_loci.tsv` | one row per locus: `locus` (`chrom:start-end:strand`), `chrom`, `start`, `end`, `strand`, `n_windows`, `peak_window` (highest score, then lowest P), `peak_score`, `min_p`, `min_q`, `features` (union), and `GC`, `n_stems_scored`, `longest_stem`, `sequence` of the peak window; sorted by `min_q`, then by decreasing `peak_score` |
| `empirical_null.tsv` | `mu0`, `sigma0`, `pi0`, `n` (windows), `central_low`, `central_high` (the z range of the fit) |
| `tail_fdr.tsv` | per set (top 2, 1, 0.5, 0.2, 0.1, 0.05 and 0.01% of windows, and the BH set): `n_windows`, `z_threshold`, `p_threshold` (theoretical-null P at the threshold), `fdr_theoretical_null` and `fdr_empirical_null` (tail-area Fdr under N(0, 1) and under the empirical null), `expected_null_windows`, `signal_windows` = R (1 − Fdr), `lfdr_at_boundary` and `mean_lfdr` (the mean lfdr of the set, another estimate of its Fdr) |
| `summary.tsv` | `n_windows`, `n_distinct_coordinates`, `n_p_le_0.05`, `n_p_le_0.01`, `n_bh_windows`, `n_bh_loci`, `n_lfdr_le_0.2`, `signal_windows_top_1pct` |

## Numbers of the paper

| Number in the paper | Where in the paper | Output of this pipeline |
|---|---|---|
| 7,869,821 windows (6,924,574 distinct coordinates) | Methods; Supplementary Methods ("Genome scan"); Results, Fig. 4a | `summary.tsv`: `n_windows`, `n_distinct_coordinates`; sum of `n_kept` in `window_counts.tsv` |
| 565 windows with BH q ≤ 0.05, merged into 481 loci | Results (genome scan); Methods | `summary.tsv`: `n_bh_windows`, `n_bh_loci`; `bh_windows.tsv`, `bh_loci.tsv` |
| μ0 = 0.673, σ0 = 0.548, π0 = 0.951 (fitted on z 0.22–1.01) | Methods; Supplementary Methods; Supplementary Fig. S4a | `empirical_null.tsv` |
| top 1% = 78,699 windows, about 30,000 of them signal; background fraction 0.613 directly and 0.612 from the lfdr | Results; Supplementary Methods; Fig. 4c | `tail_fdr.tsv`, row `top 1%`: `n_windows`, `signal_windows` (30,445), `fdr_empirical_null` (0.6131), `mean_lfdr` (0.6124) |
| tail Fdr of the top-fraction sets | Fig. 4c | `tail_fdr.tsv`: `fdr_empirical_null` (and `mean_lfdr`) |
| 4,379 RLM-Fam candidates with lfdr ≤ 0.2 | Results | the `lfdr` of `eligible_windows.parquet`, reduced to candidates by `rlmhelix catalogue build` (smallest lfdr of the windows overlapping the core on the same strand); at the window level 10,521 windows have lfdr ≤ 0.2 (`summary.tsv`: `n_lfdr_le_0.2`), from z ≈ 2.72 (P ≈ 0.0033) |

The top 1% of windows by calibrated P (78,699 windows, P ≤ 0.02088) is the starting set of the candidate selection
(`docs/pipeline_select.md`).

Window counts of step 1 for the whole genome (sums of `window_counts.tsv`): 22,773,867 window positions; 13,006
dropped for other bases; of the 22,760,861 tested, 13,448,163 removed by RepeatMasker, 1,440,638 by tantan and 2,239
by IRF, leaving 7,869,821.

## Run time and resources

| Step | Paper's run or measurement | Resources |
|---|---|---|
| `prepare` | 33 min for the whole genome (measured on the CPU of an NVIDIA GB10 system, one IRF process at a time; about 2 min of it reading the annotation and the RepeatMasker table) | 4.8 GB RAM; windows 0.7 GB |
| `score` | about 2 days for 7.87 million windows on one NVIDIA GB10 (bf16, batches of 16, 12 CPU workers for RNAfold and the DP). The paper's run scored 12,003,634 windows (the eligible windows plus intergenic and IRF-positive ones) in about 70 h of wall time, about 0.021 s per window; the Results quote the same estimate (about 2 days) from the runtime benchmark of Fig. 1. | about 5 GB of GPU memory in the runtime benchmark; scores about 0.9 GB |
| `fdr` | 40 s for the 7,869,821 windows | 2.7 GB RAM; `eligible_windows.parquet` 1.0 GB |

## Reproducibility

- With bf16 on a GPU and batches of 16 windows, scores and calibrated P values equal those of the paper. Other batch
  sizes change the bf16 rounding of the language model and shift scores by about 0.1%; which windows share a batch
  does not matter. `scan score` fills the last batch of each shard up to the batch size with copies of a window, so
  that every window is scored in a full batch whatever `--chunk` is. In the paper's run the last batch of each
  chromosome was shorter, so a few windows per chromosome (at most 15) can differ from the paper by about 0.1%.
- The paper's run also made the intergenic windows and scored all 12,003,634 windows that passed RepeatMasker and
  tantan; the 4,131,574 intergenic windows and the 2,239 IRF-positive transcribed windows were removed before BH.
  Scores do not depend on the other windows, so the eligible windows, their scores and everything computed from them
  are the same. `prepare --keep-intergenic` makes the larger set; `fdr` removes intergenic windows in any case.
- Checked against the paper's outputs: the 12,003,634 windows after RepeatMasker and tantan (all columns, all 24
  chromosomes) and their 3,584 IRF calls; the 7,869,821 windows written by `scan prepare`; student scores and stem
  counts (exact on 1,360 windows scored on a GPU); calibrated P values (to 10⁻¹³ relative); and the BH q values,
  loci, empirical null, tail Fdr and lfdr computed from the stored P values (BH and loci identical, null parameters
  to 10⁻¹⁴, lfdr to 10⁻⁹).
