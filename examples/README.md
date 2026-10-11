# Example

`known_rnas.fa` holds the human representative sequences of four Rfam 15.1 families (tRNA, the let-7 precursor, U2
snRNA and U3 snoRNA), built as for the Rfam benchmark of the paper; these four families are not among its 328
benchmark families. Score them and compare each with 100 of its own dinucleotide-preserving shuffles (about 45 s on one
GPU):

```bash
rlmhelix score examples/known_rnas.fa --shuffles 100 --seed 0 -o known_rnas_scores.tsv
```

The output has the columns `id`, `length`, `score`, `p_value`, `n_pairs`, `n_stems`, `n_stems_scored`, `sequence`
and `structure`, and with `--shuffles` also `shuffle_mean`, `shuffle_sd`, `z` and `empirical_p`. `n_pairs` counts the
decoded pairs after removal of isolated pairs (more than the scored pairs in `structure`), `n_stems` the stems of any
length and `n_stems_scored` the stems of at least four pairs, whose total weight is the score; the three columns mean
the same as in the output of `rlmhelix scan score`.

Expected output (`known_rnas_scores.tsv`; `shuffle_sd` and the sequence and structure columns are left out here):

| id | length | score | p_value | n_pairs | n_stems | n_stems_scored | shuffle_mean | z | empirical_p |
|---|---|---|---|---|---|---|---|---|---|
| RF00005_human (tRNA) | 74 | 23.2 | 0.029 | 23 | 6 | 4 | 7.3 | 3.5 | 0.020 |
| RF00027_human (let-7 precursor) | 88 | 32.7 | 0.0015 | 35 | 7 | 4 | 8.8 | 4.2 | 0.020 |
| RF00004_human (U2 snRNA) | 192 | 40.6 | 0.00082 | 64 | 15 | 9 | 15.0 | 6.2 | 0.0099 |
| RF00012_human (U3 snoRNA) | 217 | 75.7 | 0 | 76 | 14 | 11 | 15.9 | 11.3 | 0.0099 |

`p_value` is the calibrated P value of the genome scan, built for 200-nt windows. It is 0 when the score lies above
the upper end of the fitted null distribution (the GEV fitted to the shuffles has a bounded upper tail), as for the U3
snoRNA here and for the top windows of the scan. For short RNAs such as the tRNA and the let-7 precursor, `z` and
`empirical_p` are the appropriate comparison (with 100 shuffles the smallest possible `empirical_p` is 1/101). Scores
can differ in the last digits on other GPUs and with other batch sizes.
