# Candidate selection

This page describes how to go from the genome-scan windows to the RLM-Fam candidates: structural cores, homology
search and selection routes, covariance-model (CM) alignments of up to 55 species, and provisional families. The
code is `rlmhelix.select`; each step is a sub-command of `rlmhelix select`. The numbers in the right-hand column of
the tables are those of the paper.

## Overview

| Step | Command | What it does | Paper |
|---|---|---|---|
| 1 | `select cores` | top 1% of windows by calibrated P, loci, structural cores, representatives (GPU) | 78,699 windows (P ≤ 0.02088) → 71,459 loci → 71,196 representatives |
| 2 | `select chunk` | genomes cut into 10-kb chunks overlapping by 500 nt | 55 genomes |
| 3 | `select search` | Riboseek search of the representatives against the 20 seed genomes | 71,196 queries × 20 species |
| 4 | `select routes` | high-confidence hits, lineage counts, four selection routes | 15,913 selected |
| 5 | `select search` (×2) | source-strand search of the selected cores against the 20 seed and 35 expansion genomes | 71,196 queries (seed genomes; only the 15,913 selected are used) and 15,913 (expansion genomes) |
| 6 | `select seed` | hit loci, seed alignments and CMs, boundary search and seed CMs (Infernal) | 15,821 with a seed CM (92 without a high-confidence locus) |
| 7 | `select expand` | expansion loci, CM search, inclusion rule, final55 alignments and CMs (Infernal) | 15,821 candidates, median 34 species |
| 8 | `select homologues` | structure-free homologues per candidate, the input of covariation condition (iii) (also written by `select expand`) | 15,821 FASTA files, median 32 sequences |
| 9 | `select families` | HMM pairs, reciprocal CM scores, connected components (Infernal) | 14,037 provisional families |
| 10 | `catalogue candidates` | candidate table of the catalogue steps: `candidates.parquet` with the family of each candidate | 15,821 candidates in 14,037 families |

## Software

| Program | Version | Used in |
|---|---|---|
| Riboseek | commit `eba085cb` (https://github.com/steineggerlab/riboseek), CPU build | steps 3 and 5 |
| Infernal (`cmbuild`, `cmalign`, `cmcalibrate`, `cmsearch`) | 1.1.5 | steps 6, 7 and 9 |
| ViennaRNA (Python bindings) and the RLMHelix student | as in `environments/environment.yml` | step 1 |

Programs are taken from `PATH` (`shutil.which`); `--riboseek PATH` and `--infernal-bin DIR` point to other builds.
Each step looks up its programs before it reads any input and exits with status 1 and the name of a missing one:
`select search` needs Riboseek, `select seed` and `select expand` need `cmbuild`, `cmalign`, `cmcalibrate` and
`cmsearch`, and `select families` needs `cmsearch`.

## Data

- **Scan windows.** The eligible windows of the genome scan (7,869,821 windows of 200 nt in the paper; see the scan
  pipeline). Required columns: `window_id`, `chrom`, `start`, `end` (0-based, half-open), `strand`, `sequence` (in the
  scan's orientation), the RLMHelix score (`score_min4` or `score`) and the calibrated P value (`p_value` or
  `p_min4_mixed`; other names with `--p-column` and `--score-column`). `features`, `GC`, `n_stems_scored` (the paper's
  `n_stems_ge4`, carried under either name), `longest_stem`, `repeatmasker_bp` and `tantan_bp` are carried into the
  locus table when present. Ties in P are broken by the order of the table (chromosomes 1–22, X, Y in the paper).
- **Genomes.** Ensembl release 116 assemblies of the 55 species in `data/genome_panel_55.tsv` (unmasked `dna`; the
  primary assembly for human, mouse and zebrafish and the toplevel sequence otherwise). The panel gives the Ensembl
  slug, the lineage group (`primate`, `nonprimate_mammal`, `nonmammal_tetrapod`, `other_vertebrate`) and the search
  stage: 20 species for the seed search (9 non-human primates, 4 non-primate mammals, 3 non-mammalian tetrapods,
  4 other vertebrates) and 35 expansion genomes (34 vertebrates and the human genome, searched for paralogous loci).
  An identical copy is shipped with the package (`rlmhelix/data/genome_panel_55.tsv`) and is the default of
  `--panel` in `select search`, `routes`, `seed` and `expand`, so these work after a regular `pip install`;
  `--panel FILE` gives another panel. `scripts/download_data.py DEST --only ensembl` writes the genomes to
  `DEST/ensembl116/<slug>/<Slug>.<assembly>.dna.<kind>.fa.gz`. `--genomes DIR` finds `DIR/<slug>.fa.gz`,
  `DIR/<slug>/*.fa.gz` (this layout, with `--genomes DEST/ensembl116`) or Ensembl's `DIR/<Slug>.*.dna.*.fa.gz`.

## Commands

The commands below reproduce the selection of the paper from the scan output, with the genomes downloaded to
`data/raw/ensembl116/` by `scripts/download_data.py data/raw --only ensembl`. Other paths are examples.

```bash
# 1. Structural cores (GPU; RNAfold and the DP run on --n-jobs CPU workers)
rlmhelix select cores scan/eligible_windows.parquet -o select/cores --batch-size 16 --n-jobs 12

# 2. Genome chunks of the 55 genomes (the species is taken from each file name)
rlmhelix select chunk data/raw/ensembl116/*/*.fa.gz -o chunks

# 3. Search for the routes: the lexicographically normalized cores, as in the paper
rlmhelix select search select/cores/representatives_canonical.fa --chunks chunks --stage seed \
    -o select/search_routes --threads 20

# 4. Routes
rlmhelix select routes select/cores/representatives.parquet --hits select/search_routes \
    -o select/selected.parquet            # also writes select/selected.source.fa

# 5. Source-strand searches for the alignments
rlmhelix select search select/selected.source.fa --chunks chunks --stage seed -o select/search_seed
rlmhelix select search select/selected.source.fa --chunks chunks --stage expansion -o select/search_expansion

# 6. Seed alignments and CMs
rlmhelix select seed select/selected.parquet --hits select/search_seed --genomes data/raw/ensembl116 \
    -o select/seed --workers 20

# 7. Expansion to 55 species
rlmhelix select expand select/selected.parquet --seed select/seed --hits select/search_expansion \
    --genomes data/raw/ensembl116 -o select/expand --workers 20

# 8. Structure-free homologues for covariation condition (iii); only needed to rebuild them, because
#    select expand already writes the same files to select/expand/homologues
rlmhelix select homologues select/expand/candidates.parquet --seed select/seed --expand select/expand \
    -o select/homologues

# 9. Provisional families
rlmhelix select families select/expand/candidates.parquet --expand select/expand -o select/families --workers 15

# 10. Candidate table of the catalogue steps (catalogue/README.md): adds rlmfam_id and the family size
rlmhelix catalogue candidates select/expand/candidates.parquet --families select/families/families.tsv \
    -o candidates.tsv
```

The catalogue steps then read `candidates.tsv`, `select/expand/homologues` (`catalogue covariation --homologues`,
condition iii) and `select/expand/candidates` (`catalogue covariation --alignments`, conditions i and ii); see
`catalogue/README.md`.

`rlmhelix select search` skips species whose search finished (a `done` marker next to `hits.tsv`), so it can be
resumed. It records the SHA-256 checksum and path of the query FASTA in `<work>/queries.sha256` when it builds the
query database; a rerun into the same `-o` with another query FASTA (or into a work directory without the checksum)
is refused with status 1 unless `--force` is given, which removes the query database, the `done` markers and the
checksum and searches every species again. Use a separate `-o` for each of the three searches above. The Riboseek
commands per species are

```bash
riboseek createdb queries.fa query_db                        # once
riboseek createdb <species>.10kb_o500.fa.gz target_db
riboseek search query_db target_db aln tmp -a -s 7.5 -e 10 --max-seqs 100 --alt-ali 0 \
    --prefilter-mode 0 --gpu 0 --threads 20 --remove-tmp-files 1
riboseek convertalis query_db target_db aln hits.tsv \
    --format-output query,target,evalue,bits,qstart,qend,tstart,tend,alnlen,qcov,tcov,qaln,taln
```

and the Infernal commands per candidate are

```bash
cmbuild -F --hand initial.cm seed.sto                         # seed alignment, SS_cons = core structure, RF = all
cmalign --cpu 1 --cyk --notrunc --noprob -o refined.sto initial.cm sequences.fa
cmbuild -F --hand final.cm refined.sto
cmcalibrate -L 0.01 --cpu 1 calibrated.cm                     # copy of final.cm
cmsearch --max -g --cyk --acyk --notrunc --toponly --noali --cpu 1 --smxsize 2048 -T -1000 \
    --tblout cmsearch.tbl -o cmsearch.out calibrated.cm crops.fa
cmalign --cpu 1 --cyk --notrunc --noprob -o realigned_3l.sto calibrated.cm realign.fa
cmbuild -F --hand realigned_3l.cm realigned_3l.sto              # seed CM
# expansion: cmcalibrate and cmsearch as above with a copy of realigned_3l.cm on the expansion crops, then
cmalign --cpu 1 --cyk --notrunc --noprob --mapali realigned_3l.sto -o final55.sto realigned_3l.cm accepted.fa
cmbuild -F --hand final55.cm final55.sto
# families
cmsearch --hmmonly --toponly --noali --cpu 1 -E 10 --tblout hmm.tbl -o /dev/null seed20_3l.calibrated.cm human_cores.fa
cmsearch --max -g --cyk --acyk --notrunc --toponly --noali --cpu 1 --smxsize 2048 -T -1000 \
    --tblout rescore.tbl -o /dev/null seed20_3l.calibrated.cm targets.fa
```

**Failures and empty results.** `select seed` and `select expand` run one task per candidate. Failed tasks are
reported on standard error (`select seed: N of M tasks failed; first: <query_id>: <error>`; the errors of all tasks
are in the `status` and `error` columns of `seed_manifest.parquet` and `final55_manifest.parquet`), the run continues
with the successful candidates, and the command exits with status 1 when every task failed. `select families` stops
with an error when any of its Infernal searches fails, with the number of failures, the first error and the log
directory. Steps that would produce nothing exit with status 1 and a message: `select cores` (no windows, or no core),
`select routes` (no representative meets a route; empty outputs are written), `select seed` (no candidates, or no
high-confidence locus), `select expand` (no candidate with a seed CM) and `select homologues` and `select families`
(no candidates). `select expand` also finishes when a candidate has no expansion locus; its final55 alignment is then
the seed alignment. When no candidate pair is scored in `select families`, every candidate is its own family and
`edges.parquet` is empty.

## Steps, rules and outputs

### 1. Structural cores (`select cores`)

- **Top windows.** The ⌈0.01 n⌉ windows with the smallest calibrated P (stable sort).
  `top_windows.parquet`: the selected rows of the input, ordered by P.
- **Loci.** Windows on the same strand that overlap by at least 1 nt are merged (windows that only touch are not).
  The peak window has the smallest P, then the largest score. `loci.parquet`: `locus`, `chrom`, `start`, `end`,
  `strand`, `n_windows`, `peak_window`, `peak_start`, `peak_end`, `peak_score`, `peak_p` and the carried columns;
  ordered by peak P, peak score, chromosome and start.
- **Cores.** The peak window is scored again with the student (bf16, batches of 16 in locus order) and RNAfold;
  the weights √(d · p) are decoded by the Nussinov DP, isolated pairs and stems of fewer than four pairs are removed,
  and each outermost base-paired interval is a core, without flanks. The core score is the sum of the weights of its
  pairs, and the core with the largest score is the primary core of the locus. `cores.parquet`: `locus`,
  `peak_window`, `core_index`, `chrom`, `strand`, `features`, `window_score`, `core_start_in_window`,
  `core_end_in_window`, `core_length`, `core_genome_start`, `core_genome_end` (forward strand, 0-based half-open),
  `core_score`, `n_pairs`, `n_stems`, `longest_stem`, `touches_window_edge`, `core_sequence` and `core_dot_bracket`
  (scan orientation), `score_fraction` (share of the locus total) and `component_rank` (1 = primary).
- **Representatives.** Each primary core's (sequence, dot-bracket) and (reverse complement, reversed dot-bracket) are
  compared lexicographically and the smaller is its key; identical keys are collapsed, keeping the core with the
  smallest peak P, then the largest core score, then the peak window ID, and representatives are numbered in that
  order (`pvalue_top1_primary_0000001`, … ; `--id-prefix` changes the prefix). `representatives.parquet`: `query_id`,
  the core columns, `peak_p`, `GC`, `n_windows`, `sequence_u` and `core_dot_bracket` (source strand, RNA),
  `canonical_sequence` and `canonical_structure` (normalized orientation) and `dedup_key`; `membership.parquet`: every
  primary core with the `query_id` of its representative; `representatives_canonical.fa` and
  `representatives_source.fa` (DNA).

### 2–3. Chunks and search (`select chunk`, `select search`)

`select chunk` writes `<slug>.10kb_o500.fa.gz` with chunk IDs `<slug>|<contig>|<start>|<end>` (0-based, half-open;
chunks of only N are skipped) and a TSV of the chunks. `select search` writes `<work>/<species>/hits.tsv`
(columns `query_id`, `target_id`, `evalue`, `bits`, `qstart`, `qend`, `tstart`, `tend`, `alnlen`, `qcov`, `tcov`,
`qaln`, `taln`; Riboseek coordinates are 1-based inclusive within a chunk, and `qstart > qend` marks a hit on the
reverse strand of the target) and `search_status.tsv`.

### 4. Routes (`select routes`)

Hits with E ≤ 10⁻³ and query coverage ≥ 0.8 are high-confidence. For each representative, the seed species with a
high-confidence hit are counted per lineage group, and a representative is selected when it meets any route; it is
labelled by the first one met in this order:

| Route | Condition | Paper (of 15,913) | Paper (of the 15,821 candidates) |
|---|---|---|---|
| `vertebrate_wide` | ≥ 1 non-mammalian tetrapod and ≥ 1 other vertebrate | 2,110 | 2,110 |
| `mammal_conserved` | ≥ 3 non-primate mammals and ≥ 3 primates | 11,140 | 11,140 |
| `primate_focused` | ≥ 6 primates and peak P ≤ 0.001 | 1,531 | 1,531 |
| `high_score_rescue` | peak P ≤ 0.0005 | 1,132 | 1,040 |

Output (`-o`): the selected representatives with the four lineage counts, `route_<name>` flags, `selected` and
`route`, and `<output>.source.fa` (source strand, DNA).

### 5–6. Seed alignments (`select seed`)

Hits with E ≤ 0.1 of the source-strand search of the seed genomes are merged into loci per query, species, contig
and strand (hits from overlapping chunks overlap by at least 1 nt); the best hit (high-confidence first, then bits,
E value and coverage) represents the locus, and the locus is high-confidence when any of its hits is. Distinct loci in
one species (paralogues) are kept. For each selected core, its high-confidence loci (ordered by bits, E value and
coverage) are projected onto the human core through the pairwise alignments, one row per distinct target sequence,
and written with the core structure as SS_cons and every human column as RF (`seed.sto`). A CM is built, the
sequences realigned and the CM rebuilt (`final.cm`). Every high-confidence locus is cropped with a flank of one locus
length on each side, `[start − L, end + L)`, its boundaries are set by the best glocal CYK hit of the calibrated CM in
the crop, the distinct hit sequences (ordered by bits) are realigned with the human core, and the CM is rebuilt
(`realigned_3l.sto`, `realigned_3l.cm`: the seed alignment and CM). Cores without a high-confidence locus get no seed
CM (92 in the paper, all `high_score_rescue`).

Outputs: `seed_loci.parquet` (one row per locus: `query_id`, `species`, `contig`, `target_strand`, `locus_start`,
`locus_end`, `locus_id`, `n_chunk_hits`, `any_high_confidence` and the fields of the best hit),
`seed_locus_counts.parquet` (`n_species_e01`, `n_loci_e01`, `n_species_hc`, `n_loci_hc`), `seed_manifest.parquet`
(`status`, `seed_cm`, `n_hc_loci`, `n_hc_species`, `seed_q10_bits`, `threshold_bits`) and `candidates/<query_id>/`.

### 7. Expansion to 55 species (`select expand`)

Loci with E ≤ 0.1 in the 35 expansion genomes are cropped to three times the core length around the locus midpoint
and searched with the seed CM (calibrated with `-L 0.01`). The best hit per crop is accepted when its bit score is at
least the 10th percentile (linear interpolation) of the seed hits' scores minus 10 bits and it covers at least 80% of
the model. Accepted sequences are deduplicated (best bits first), those identical to a seed row are dropped, and the
rest are added with `cmalign --mapali`; the CM is rebuilt from `final55.sto`. The species and loci of a candidate are
those of its high-confidence seed loci plus its accepted expansion loci (the human genome counts as a species when it
has an accepted paralogous or self locus), and candidates with more than 100 high-confidence seed loci or final55
loci are flagged `high_copy`.

Outputs: `candidates.parquet` (one row per `query_id`; the selected columns; `threshold_bits`; seed
locus counts; `n_accepted_expansion_loci`, `n_accepted_expansion_species` and per lineage group;
`n_added_unique_sequences`; `n_final_alignment_sequences`; `n_final55_species`; `n_final55_loci`; `high_copy`;
`msa_path`; `cm_path`), `final55_manifest.parquet`, `candidates/<query_id>/` (`final55.sto`, `final55.cm`,
`expansion_hits.parquet`, `accepted_expansion_hits.parquet`, `crops.fa`) and `homologues/<query_id>.fa` (step 8).
`candidates.parquet` is the input of `select homologues`, `select families` and `rlmhelix catalogue candidates`, and
`candidates/` is the `--alignments` directory of `rlmhelix catalogue covariation` (conditions i and ii).

### 8. Structure-free homologues (`select homologues`)

For each candidate, one FASTA file `<core_id>.fa` (the core ID is the `query_id`) holds the sequences that covariation
condition (iii) aligns with MAFFT without structural information, built as for the paper:

1. the records of the seed `sequences.fa`: the human core (named `core_id`, always first) and one sequence per
   distinct target of the high-confidence seed loci, cut at the Riboseek pairwise alignment (`h00001_<species>`, …);
2. the expansion loci (35 genomes, including human paralogous and self loci) whose Riboseek hit has E ≤ 10⁻³ and
   query coverage ≥ 0.8, cut from the expansion crops at the merged locus interval (not at the CM hit), in the
   orientation of the hit, named `x<n>_<species>_<contig>_<locus_start>_<strand>` with `n` counting these loci in crop
   order (six digits) and characters other than letters, digits and `_.:-` (such as `+`) replaced by `_`.

Sequences are uppercase DNA without gaps, and exact duplicates (also of the human core) are removed, keeping the
first. The inclusion rule and CM scores are not used, so these homologues do not depend on the predicted structure.
Inputs: `candidates.parquet` and the output directories of `select seed` (`candidates/<id>/sequences.fa`) and
`select expand` (`candidates/<id>/crops.fa`, `expansion_hits.parquet`). Output: `<output>/<core_id>.fa`, read by
`rlmhelix catalogue covariation --homologues DIR` (in Python,
`catalogue.covariation_one(core_id, "iii", "<core_id>.fa", structure, workdir)`). In the paper the files have a median
of 32 sequences.

### 9. Provisional families (`select families`)

Each calibrated seed CM is searched in HMM-only mode against the human cores of all candidates (E ≤ 10); every pair
of distinct candidates with a hit in either direction is a candidate pair. Each candidate contributes up to four
held-out sequences from its added expansion sequences: one per species (the best by bits), chosen round-robin over
the lineage groups, each time the sequence most different in 4-mers (1 − Jaccard) from those already chosen. Each CM
scores the sequences of its own candidate and of its paired candidates (glocal CYK, all hits reported). With
b_{A→B} the median score of A's CM on B's held-out sequences and T_A A's inclusion threshold (10th percentile of the
seed hits minus 10 bits), r_{A→B} = (b_{A→B} − T_A)/T_A. Pairs scored in both directions whose mean of r_{A→B} and
r_{B→A} is at least −0.525 are joined, and families are the connected components (singletons included), numbered
`C01`, `C02`, … by decreasing size and then by their smallest `query_id`. The threshold was chosen as the most
stringent value at which no two candidates annotated with different Rfam families were joined while the pairs of the
14 Rfam families with several candidates were recovered best; it is neither an Rfam gathering threshold nor a P value.

Outputs: `families.tsv` (`query_id`, `family_id` (C01, …), `rlmfam_id` (RLMFam00001, …, the identifier used in the
paper), `family_size`, `family_rank`), `edges.parquet` (`unit_a`, `unit_b`, `b_ab`, `b_ba`, `r_ab`, `r_ba`,
`ratio_mean`; the pairs scored in both directions), `candidate_pairs.parquet`, `samples.parquet` and
`reciprocal_scores.parquet`. `families.tsv` holds only the family membership; `rlmhelix catalogue candidates` merges
it with `candidates.parquet` into the candidate table of the catalogue steps (step 10 of the commands).

## Numbers of the paper and where they come from

| Number | Value | Output of this package |
|---|---|---|
| Eligible windows | 7,869,821 | rows of `scan/eligible_windows.parquet` (`rlmhelix scan fdr`) |
| Top 1% windows; P cut-off | 78,699; 0.020876 | rows of `select/cores/top_windows.parquet`; `p_value` of its last row |
| Loci (primary cores) | 71,459 | rows of `loci.parquet`; rows of `membership.parquet` |
| Cores in the peak windows | 173,641 | rows of `cores.parquet` |
| Representatives (identical cores removed) | 71,196 (263 removed) | rows of `representatives.parquet` (and of `membership.parquet` minus these) |
| Selected by the routes | 15,913 (2,110 / 11,140 / 1,531 / 1,132) | rows of `select/selected.parquet`; its `route` column |
| High-confidence seed loci | 1,768,787 | rows of `select/seed/seed_loci.parquet` with `any_high_confidence` (sum of `n_loci_hc` in `seed_locus_counts.parquet`) |
| Without a high-confidence seed locus (no CM) | 92 | rows of `select/seed/seed_manifest.parquet` with `seed_cm` False |
| RLM-Fam candidates; median species | 15,821; 34 | rows of `select/expand/candidates.parquet`; median of `n_final55_species` |
| Routes among the candidates | 2,110 / 11,140 / 1,531 / 1,040 | `route` of `candidates.parquet` (`selection_route` of the catalogue table `rlmfam_candidates.tsv`) |
| Expansion loci searched; accepted | 5,829,188; 4,626,366 | rows of `select/expand/candidates/<query_id>/expansion_hits.parquet` (one per crop); sum of `n_accepted_expansion_loci` in `candidates.parquet` |
| High-copy candidates | 2,587 | `high_copy` of `candidates.parquet` |
| HMM candidate pairs; scored in both directions | 177,621; 177,600 | rows of `select/families/candidate_pairs.parquet`; rows of `edges.parquet` |
| Provisional families; largest | 14,037; 466 | distinct `family_id` of `families.tsv`; largest `family_size` (threshold −0.525) |

## Run times in the paper

Recorded wall-clock times (CPU searches with 20 threads):

| Step | Time |
|---|---|
| Structural cores | not recorded separately; the peak windows (71,459) are under 1% of the scan's 7.9 million windows, whose scoring took about two GPU-days |
| Search for the routes (71,196 queries × 20 genomes) | about 25 h for the species searched with all queries (part of the hits were reused from an earlier run) |
| Source-strand search of the 20 seed genomes (71,196 queries) | 29.1 h |
| Source-strand search of the 35 expansion genomes (15,913 queries) | 13.7 h |
| Seed alignments and CMs | 2.6 CPU-hours |
| Boundary search and seed CMs (3L crops) | 203 CPU-hours |
| Expansion CM search | 481 CPU-hours |
| final55 alignments and CMs | 2.2 CPU-hours |
| Families: HMM pairs; reciprocal CM scores | 1 min wall (15 workers); 2.2 h wall, 33 CPU-hours |

## Notes on reproducing the paper exactly

- The routes were computed from the search of the normalized cores (`representatives_canonical.fa`), and the
  alignments from separate searches of the source-strand cores. Riboseek searches both strands of the targets, but
  the reported hits of a query and its reverse complement are not identical, so the two searches are kept apart.
  For example, 2 of the 92 candidates without a high-confidence seed locus had high-confidence hits in the search for
  the routes (all 92 were selected by `high_score_rescue`).
- The paper searched all 71,196 source-strand representatives against the seed genomes; only the 15,913 selected
  were used afterwards, so searching the selected ones (step 5) is expected to give the same hits for them.
- Under bf16 the student's values depend slightly on how sequences are batched. The paper re-scored the peak windows
  in batches of 16 in locus order, but 39,999 of the 71,459 peak windows were taken from an earlier run with a
  different order, so a few cores near the stem-length or pairing thresholds can differ in a fresh run.
- Ties are broken as in the paper: by table order for P values, and by bits, E value, coverage and table order for
  hits, with species processed in alphabetical order of their Ensembl slugs.

## Differences from the manuscript text

The Methods describe the selection in brief. The code that produced the paper's results, which this package follows,
does the following where the text is shorter or worded differently:

- **Orientation of the searches.** The four routes were computed from a Riboseek search of the representatives in
  their lexicographically normalized orientation (the reverse complement of the source strand for 52,754 of 71,196).
  The seed alignments come from a second search of the 20 seed genomes with the source-strand cores, and the expansion
  from a source-strand search of the 35 other genomes.
- **Boundary refinement of the seed loci.** Seed loci are refined as well as expansion loci. Each high-confidence seed
  locus is re-cropped with one locus length of flank on either side (three times the locus span, not three times the
  model length), searched by glocal CYK with the calibrated CM of the first seed alignment, and the hit sequences are
  realigned to give the seed CM used for the expansion. The 10th percentile of the inclusion threshold is taken over
  these seed hits.
- **The 92 dropped candidates.** They have no high-confidence locus in the source-strand search of the seed genomes
  (2 of them had high-confidence hits in the search for the routes). A CM of the human core alone is built for them,
  but no seed CM, and they are not expanded; all 92 were selected by `high_score_rescue`.
- **Family scores.** Only candidate pairs proposed by an HMM-only search of each seed CM against all human cores
  (E ≤ 10; 177,621 pairs) are scored. "B's sequences" are up to four held-out expansion sequences of B (one per
  species, round-robin over lineage groups, chosen for 4-mer diversity), excluding the seed rows and the human core.
  The CMs are the calibrated seed CMs, not the final55 CMs. Pairs without held-out sequences in both directions
  (21) are left out.
- **Counts.** The species of a final55 alignment include the human genome when it has an accepted locus, which its own
  locus always provides. The high-copy flag counts merged loci (high-confidence seed loci, and seed plus accepted
  expansion loci), not raw hits.
- **Covariation condition (iii).** The "Riboseek high-confidence hits" are the seed sequences cut at the pairwise
  alignments plus the high-confidence expansion loci cut at the merged locus interval (step 8).
