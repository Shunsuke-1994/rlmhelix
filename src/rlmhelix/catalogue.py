"""RLM-Fam catalogue.

Statistical and structural evidence for each candidate, and the evidence tiers:

- **Statistics.** z = -Phi^-1(p) for the calibrated P values of all eligible scan windows. The local false discovery
  rate lfdr(z) = pi0 f0(z) / f(z) uses the empirical null N(mu0, sigma0^2) and pi0 of Efron's central matching and
  the mixture density f estimated by Lindsey's method (histogram of z with bins of 0.05, Poisson regression on a
  natural cubic spline with knots at the ends and at ten quantiles of z); it is capped at 1, made non-increasing
  to the right of mu0 and set to 1 to its left. Each candidate takes the smallest calibrated P value, lfdr and BH q
  among the windows that overlap its core on the same strand.
- **icSHAPE.** Cores lying in one same-strand GENCODE v26 exon whose sequence matches the transcript take the
  processed GSE145805 reactivities of the isoform with the most measured nucleotides in each cell line. Where at
  least 50% of the core and at least 10 paired and 10 unpaired nucleotides are measured, the AUROC of reactivity for
  unpaired versus paired nucleotides is compared with the AUROCs of the labels circularly shifted by every offset,
  P = (1 + #{shifts >= observed}) / (1 + #shifts). A candidate is supported when P <= 0.05 in at least half of its
  measured cell lines.
- **Covariation.** R-scape (`-s --nofigures`; E < 0.05 on the proposed pairs) under condition (iii): high-confidence
  homologues aligned with MAFFT (`--auto`) without structural information, with the predicted core structure
  projected through the human row. Conditions (i) (the covariance-model alignment) and (ii) (its sequences
  realigned with MAFFT) are descriptive. Absent covariation is never used against a candidate.
- **Annotation.** GENCODE v50 genes and transcripts overlapping each core (`gencode_annotation`) and the distance to
  the nearest same-strand splice site; descriptive only.
- **Tiers.** Statistical support: lfdr <= 0.2 or an overlapping window with BH q <= 0.05. Structural support:
  icSHAPE support or covariation under condition (iii). A, both; B, one of the two; C, neither.

`selection_candidates` turns the output of `rlmhelix select` into the candidate table, `build_catalogue` assembles one
row per candidate and `to_table` returns the columns of Supplementary Table S4. MAFFT and R-scape are external
programs; everything else uses NumPy, SciPy and pandas. Command line:

    rlmhelix catalogue candidates select/expand/candidates.parquet --families select/families/families.tsv \
        -o candidates.tsv
    rlmhelix catalogue icshape --candidates candidates.tsv --gtf gencode.v26.annotation.gtf.gz \
        --transcripts gencode.v26.transcripts.fa.gz --track HEK293=GSM4333256_293.out.txt.gz ... -o icshape.tsv
    rlmhelix catalogue covariation --candidates candidates.tsv --homologues select/expand/homologues \
        --workdir rscape -o covariation_iii.tsv
    rlmhelix catalogue covariation --candidates candidates.tsv --condition i --alignments select/expand/candidates \
        --workdir rscape -o covariation_i.tsv
    rlmhelix catalogue build --candidates candidates.tsv --windows scan/eligible_windows.parquet \
        --icshape icshape.tsv --covariation covariation_iii.tsv covariation_i.tsv --gtf gencode.v50.annotation.gtf.gz \
        -o rlmfam_candidates.tsv --families families.tsv
    rlmhelix catalogue concordance rlmfam_candidates.tsv -o concordance.tsv
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import json
import re
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy import stats

LFDR_SUPPORT = 0.2  # statistical support: lfdr <= 0.2 ...
BH_SUPPORT = 0.05  # ... or an overlapping window with BH q <= 0.05
TOP1_P = 0.020876023786561414  # calibrated P value at the top-1% boundary of the paper's scan windows
TOP1_FRACTION = 0.01
EMPIRICAL_NULL = (0.6729420119015076, 0.5481601926577551, 0.950765247091798)  # mu0, sigma0, pi0 of the paper's scan
LINDSEY_BIN = 0.05
LINDSEY_KNOTS = (0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.995, 0.999, 0.9999)  # quantiles of z (plus both ends)
ICSHAPE_COVERAGE = 0.5
ICSHAPE_MIN_CLASS = 10
ICSHAPE_P = 0.05
ICSHAPE_FRACTION = 0.5
ICSHAPE_TRACKS = {"HEK293": "GSM4333256_293.out.txt.gz", "HeLa": "GSM4333258_HeLa.out.txt.gz",
                  "K562": "GSM4333260_K562.out.txt.gz", "HepG2": "GSM4333262_HepG2.out.txt.gz",
                  "H9": "GSM4333264_H9.out.txt.gz"}  # GSE145805 in vivo tracks used in the paper
DSRBP = ("DGCR8", "DROSHA", "ILF3", "STAU2", "DHX30", "XPO5")
CONDITIONS = ("i", "ii", "iii")  # CM alignment; its sequences realigned with MAFFT; structure-free homologues
COV_STATUS = ("detected", "not evaluable", "not detected, low power (expected < 1 pair)",
              "not detected (expected >= 1 pair)")
TIER_ORDER = {"A": 0, "B": 1, "C": 2}
GAPS = "-."

# (column of Supplementary Table S4, column of the working table built by `build_catalogue`)
S4_COLUMNS = [
    ("rlmfam_id", "rlmfam_id"), ("core_id", "core_id"), ("family_size", "cluster_size"),
    ("evidence_tier", "tier"), ("tier_evidence", "tier_evidence"), ("chrom", "chrom"), ("strand", "strand"),
    ("core_start_0based", "core_start0"), ("core_end", "core_end0"), ("locus_hg38", "browser_locus"),
    ("core_length_nt", "core_length"), ("n_base_pairs", "n_pairs"), ("n_stems", "n_stems"),
    ("core_score", "core_score"), ("core_sequence", "human_core_rna_sequence"),
    ("core_structure", "human_core_predicted_ss"), ("scan_p", "scan_pvalue"), ("max_window_z", "max_window_z"),
    ("lfdr", "lfdr"), ("min_window_bh_q", "min_window_q_bh"), ("bh_window", "bh_window"),
    ("statistical_support", "stat_support"), ("n_top1pct_windows", "n_top1pct_windows"),
    ("best_window_id", "best_window_id"), ("selection_route", "lineage_route"),
    ("n_species_alignment", "msa_species"), ("n_sequences_alignment", "msa_sequences"),
    ("icshape_measured", "icshape_measured"), ("icshape_supported", "icshape_supported"),
    ("icshape_n_cell_lines", "icshape_n_contexts"), ("icshape_mean_auroc", "icshape_mean_auroc"),
    ("icshape_fraction_cell_lines_p_le_0.05", "icshape_frac_contexts_p05"), ("covariation_status", "cov_status"),
    ("covariation_iii_structure_free", "cov_detected"), ("covariation_iii_two_pairs", "cov_two_pairs"),
    ("covariation_iii_n_significant_pairs", "cov_n_significant_pairs"),
    ("covariation_iii_n_pairs_tested", "cov_n_proposed_pairs_tested"),
    ("covariation_iii_expected_pairs", "cov_expected_pairs_power"), ("covariation_iii_n_sequences", "cov_nseq"),
    ("covariation_iii_mean_identity_pct", "cov_avg_identity"),
    ("covariation_i_cm_alignment", "cov_detected_cm_aligned"),
    ("covariation_ii_same_sequences_mafft", "cov_detected_same_seqs_mafft"),
    ("structural_support", "structure_support"), ("structural_evidence_measured", "structure_evidence_measured"),
    ("homologue_msa_scored", "ms_scored"), ("homologue_rlmhelix_msa_no_human", "ms_wo_native"),
    ("homologue_rlmhelix_msa_no_human_z", "ms_wo_z"), ("homologue_rlmhelix_msa_no_human_p", "ms_wo_p"),
    ("homologue_above_all_19_shuffles", "ms_wo_pass"), ("homologue_rlmhelix_msa_with_human_z", "ms_wi_z"),
    ("homologue_above_all_19_with_human", "ms_wi_pass"), ("eclip_n_experiments", "eclip_n_experiments"),
    ("eclip_experiments", "eclip_experiments"), ("eclip_dsrbp", "eclip_dsrbp"),
    ("eclip_dsrbp_names", "eclip_dsrbp_names"), ("splice_site_distance_nt", "splice_site_distance_nt"),
    ("splice_site_class", "splice_site_class"), ("scan_region_types", "source_features"),
    ("annotation_status", "annotation_status"), ("same_strand_genes", "same_strand_gene_names"),
    ("same_strand_gene_types", "same_strand_gene_types"), ("same_strand_contexts", "same_strand_contexts"),
    ("representative_gene", "representative_gene_name"), ("representative_gene_type", "representative_gene_type"),
    ("representative_transcript", "representative_transcript_name"),
    ("representative_context", "representative_context"), ("high_copy", "high_copy"),
    ("pseudogene_context", "pseudogene_context"), ("pseudogene_genes", "pseudogene_gene_names"),
    ("rfam_accessions", "rfam_accessions"), ("rfam_names", "rfam_names"),
]
FAMILY_COLUMNS = ["rlmfam_id", "n_candidates", "candidate_loci", "representative_locus_hg38", "representative_strand",
                  "genes", "region_contexts", "best_tier", "min_lfdr", "n_statistical_support",
                  "n_structural_support", "n_overlapping_earlier_screens", "rfam"]
# working name <- names of the same column in other tables (`select expand`, `select families`), in order of preference
_CANDIDATE_ALIASES = {
    "core_id": ("query_id",), "core_start0": ("core_genome_start",), "core_end0": ("core_genome_end",),
    "human_core_rna_sequence": ("sequence_u", "sequence"), "human_core_predicted_ss": ("core_dot_bracket", "structure"),
    "lineage_route": ("route", "primary_tier", "selection_route"), "msa_species": ("n_final55_species",),
    "msa_sequences": ("n_final_alignment_sequences",), "source_features": ("features",), "scan_pvalue": ("peak_p",),
    "provisional_cluster_id": ("family_id",), "cluster_size": ("family_size",)}
# candidate table written by `catalogue candidates` (working names)
CANDIDATE_COLUMNS = ["core_id", "rlmfam_id", "provisional_cluster_id", "cluster_size", "chrom", "strand", "core_start0",
                     "core_end0", "browser_locus", "core_length", "n_pairs", "n_stems", "core_score",
                     "human_core_rna_sequence", "human_core_predicted_ss", "scan_pvalue", "source_features",
                     "lineage_route", "msa_species", "msa_sequences", "high_copy", "msa_path", "cm_path"]
# GENCODE gene annotation of the cores (`gencode_annotation`)
GENCODE_COLUMNS = ["annotation_status", "same_strand_gene_names", "same_strand_gene_types", "same_strand_contexts",
                   "representative_gene_name", "representative_gene_type", "representative_transcript_name",
                   "representative_context", "pseudogene_gene_names"]


# ----------------------------------------------------------------------------------------------------------------
# identifiers, files and intervals
# ----------------------------------------------------------------------------------------------------------------
def family_id(cluster_id: str | int) -> str:
    """RLM-Fam family ID of a provisional family, as used in the text: 'C01' -> 'RLMFam00001'."""
    text = str(cluster_id)
    if text.startswith("RLMFam"):
        return text
    return f"RLMFam{int(text.lstrip('C')):05d}"


def browser_locus(chrom: str, start0: int, end: int) -> str:
    """Genome-browser locus (1-based inclusive) of a BED interval."""
    return f"{chrom}:{int(start0) + 1}-{int(end)}"


def read_table(path: str | Path) -> pd.DataFrame:
    """Parquet, CSV or TSV (also gzipped) by file name."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    sep = "," if ".csv" in path.suffixes else "\t"
    return pd.read_csv(path, sep=sep, low_memory=False)


def _open(path: str | Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """(identifier, sequence) pairs; the identifier is the first word of the header."""
    records, name, chunks = [], None, []
    with _open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                words = line[1:].split()
                name, chunks = (words[0] if words else ""), []
            elif line.strip():
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _interval_index(frame: pd.DataFrame, start: str = "start", end: str = "end") -> dict:
    """(chrom, strand) -> (starts, ends, running maximum of ends, row labels), sorted by start."""
    index = {}
    for key, group in frame.groupby(["chrom", "strand"], sort=False):
        group = group.sort_values(start, kind="mergesort")
        ends = group[end].to_numpy()
        index[key] = (group[start].to_numpy(), ends, np.maximum.accumulate(ends), group.index.to_numpy())
    return index


def _overlapping(index: dict, chrom: str, strand: str, start: int, end: int) -> np.ndarray:
    """Row labels of the indexed intervals overlapping [start, end) on the same strand, in order of start."""
    if (chrom, strand) not in index:
        return np.array([], dtype=int)
    starts, ends, running, rows = index[(chrom, strand)]
    hi = np.searchsorted(starts, end, side="left")
    lo = np.searchsorted(running[:hi], start, side="right")
    selected = np.arange(lo, hi)
    return rows[selected[ends[selected] > start]]


def normalize_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Candidate table of the selection step with the column names used here.

    Required: core_id (or query_id), chrom, strand, core_start0, core_end0 (BED convention), and for the structural
    evidence human_core_rna_sequence and human_core_predicted_ss. provisional_cluster_id (or rlmfam_id) gives the
    family. Annotation columns (source_features, annotation_status, same_strand_*, representative_*, high_copy,
    pseudogene_gene_names, msa_species, msa_sequences, core_score, n_stems, lineage_route) are carried over. The
    column names of `rlmhelix select expand` (core_genome_start, core_genome_end, sequence_u, core_dot_bracket, route,
    n_final55_species, n_final_alignment_sequences, features, peak_p) and `select families` (family_id, family_size)
    are accepted.
    """
    renames = {}
    for target, aliases in _CANDIDATE_ALIASES.items():
        present = [a for a in aliases if a in candidates.columns]
        if target not in candidates.columns and present:
            renames[present[0]] = target
    frame = candidates.rename(columns=renames).copy()
    missing = {"core_id", "chrom", "strand", "core_start0", "core_end0"} - set(frame.columns)
    if missing:
        raise ValueError(f"candidate table lacks {sorted(missing)}")
    if not frame.core_id.is_unique:
        raise ValueError("core_id must be unique")
    if "rlmfam_id" not in frame.columns and "provisional_cluster_id" in frame.columns:
        frame["rlmfam_id"] = frame.provisional_cluster_id.map(family_id)
    if "rlmfam_id" in frame.columns and "cluster_size" not in frame.columns:
        frame["cluster_size"] = frame.groupby("rlmfam_id").core_id.transform("size")
    if "browser_locus" not in frame.columns:
        frame["browser_locus"] = [browser_locus(c, s, e) for c, s, e in
                                  zip(frame.chrom, frame.core_start0, frame.core_end0)]
    if "core_length" not in frame.columns:
        frame["core_length"] = frame.core_end0 - frame.core_start0
    if "n_pairs" not in frame.columns and "human_core_predicted_ss" in frame.columns:
        frame["n_pairs"] = frame.human_core_predicted_ss.map(lambda s: len(structure_pairs(s)))
    return frame.reset_index(drop=True)


def selection_candidates(candidates: pd.DataFrame, families: pd.DataFrame) -> pd.DataFrame:
    """The candidate table of the catalogue steps from the output of the candidate selection.

    candidates: candidates.parquet of `rlmhelix select expand` (one row per candidate with a final55 alignment);
    families: families.tsv of `rlmhelix select families` (query_id, family_id, rlmfam_id, family_size). Returns the
    columns of CANDIDATE_COLUMNS that are available, with the working names, sorted by family and core_id. Every
    candidate must have a family.
    """
    if candidates.empty:
        raise ValueError("the candidate table is empty")
    left = candidates.rename(columns={"query_id": "core_id"}) if "core_id" not in candidates else candidates
    right = families.rename(columns={"query_id": "core_id"}) if "core_id" not in families else families
    keep = [c for c in ("core_id", "rlmfam_id", "family_id", "family_size") if c in right.columns]
    if "core_id" not in keep or "rlmfam_id" not in keep:
        raise ValueError("the family table needs query_id (or core_id) and rlmfam_id")
    table = left.drop(columns=[c for c in keep[1:] + ["provisional_cluster_id", "cluster_size"] if c in left]).merge(
        right[keep], on="core_id", how="left", validate="one_to_one")
    lacking = table.core_id[table.rlmfam_id.isna()]
    if len(lacking):
        raise ValueError(f"{len(lacking):,} candidates are not in the family table, e.g. {lacking.iloc[0]}")
    frame = normalize_candidates(table)
    frame = frame[[c for c in CANDIDATE_COLUMNS if c in frame.columns]]
    return frame.sort_values(["rlmfam_id", "core_id"], kind="mergesort").reset_index(drop=True)


# ----------------------------------------------------------------------------------------------------------------
# statistics: BH, lfdr (Lindsey's method under Efron's empirical null), assignment to candidates
# ----------------------------------------------------------------------------------------------------------------
def bh_qvalues(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg q values."""
    p = np.asarray(pvalues, dtype=np.float64)
    order = np.argsort(p, kind="mergesort")
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1)
    return out


def z_scores(pvalues: Iterable[float]) -> np.ndarray:
    """z = -Phi^-1(p) (P clipped to [1e-300, 1 - 1e-12])."""
    return -stats.norm.ppf(np.clip(np.asarray(pvalues, dtype=np.float64), 1e-300, 1 - 1e-12))


def _ns_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Natural cubic spline basis (Hastie et al., ESL eq. 5.4-5.5): 1, x, d_k - d_{K-1}."""
    knots = np.asarray(knots, float)

    def d(j):
        return (np.clip(x - knots[j], 0, None) ** 3 - np.clip(x - knots[-1], 0, None) ** 3) / (knots[-1] - knots[j])

    return np.column_stack([np.ones_like(x), x] + [d(j) - d(len(knots) - 2) for j in range(len(knots) - 2)])


def _poisson_irls(X: np.ndarray, y: np.ndarray, n_iter: int = 100) -> np.ndarray:
    beta = np.linalg.lstsq(X, np.log(y + 0.5), rcond=None)[0]
    deviance_old = np.inf
    for _ in range(n_iter):
        eta = np.clip(X @ beta, -50, 50)
        mu = np.exp(eta)
        working = eta + (y - mu) / mu
        weighted = X.T * mu
        beta = np.linalg.solve(weighted @ X + 1e-10 * np.eye(X.shape[1]), weighted @ working)
        mu = np.exp(np.clip(X @ beta, -50, 50))
        deviance = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0) - (y - mu))
        if abs(deviance_old - deviance) < 1e-8 * max(1.0, deviance):
            break
        deviance_old = deviance
    return beta


def window_lfdr(pvalues: Iterable[float], null: tuple[float, float, float] = EMPIRICAL_NULL,
                return_curve: bool = False):
    """Local false discovery rate of each window from the calibrated P values of all eligible windows.

    null: (mu0, sigma0, pi0) of the empirical null of z (Efron's central matching; the default is the paper's scan).
    The lfdr is a right-tail quantity: windows with z < mu0 get 1 (the paper reported the top 2% of windows only).
    With return_curve, also returns a table of z, the mixture density f, pi0 f0 and the lfdr on a 0.01 grid.
    """
    mu0, sigma0, pi0 = null
    z = z_scores(pvalues)
    lo, hi = np.quantile(z, [1e-6, 1 - 1e-6])
    lo, hi = np.floor(lo * 20) / 20, np.ceil(hi * 20) / 20
    knots = (np.r_[lo, np.quantile(z, LINDSEY_KNOTS), hi] - lo) / (hi - lo)
    edges = np.arange(lo, hi + LINDSEY_BIN, LINDSEY_BIN)
    counts, _ = np.histogram(z, edges)
    beta = _poisson_irls(_ns_basis(((edges[:-1] + edges[1:]) / 2 - lo) / (hi - lo), knots), counts.astype(float))

    def density(values: np.ndarray) -> np.ndarray:
        x = _ns_basis((np.asarray(values, float) - lo) / (hi - lo), knots)
        return np.exp(np.clip(x @ beta, -700, 700)) / (len(z) * LINDSEY_BIN)

    grid = np.arange(z.min() - 0.001, z.max() + 0.002, 0.001)
    curve = np.minimum(1.0, pi0 * stats.norm.pdf(grid, mu0, sigma0) / density(grid))
    right = grid >= mu0
    curve[right] = np.minimum.accumulate(curve[right])
    curve[~right] = 1.0  # the heavy left tail (scores below expectation) carries no evidence for structure
    lfdr = np.interp(z, grid, curve)
    if not return_curve:
        return lfdr
    zz = np.arange(np.floor(lo), np.ceil(z.max()) + 0.01, 0.01)
    table = pd.DataFrame({"z": zz, "f": density(zz), "pi0_f0": pi0 * stats.norm.pdf(zz, mu0, sigma0),
                          "lfdr": np.interp(zz, grid, curve)})
    return lfdr, table


def prepare_windows(windows: pd.DataFrame, null: tuple[float, float, float] = EMPIRICAL_NULL) -> pd.DataFrame:
    """Scan windows with p, z, q_bh, lfdr and window_id.

    The P value column may be named p, p_value or pvalue, and the BH q value q_bh, q_value or bh_q (the output of
    `rlmhelix scan fdr` works as it is). When q_bh or lfdr is absent it is computed here, which requires the table to
    hold all eligible windows of the scan (7,869,821 in the paper). Windows with z < mu0 get lfdr 1 (see
    `window_lfdr`), so that a given lfdr column is used for the right tail only, as in the paper.
    """
    frame = windows.copy()
    for target, aliases in (("p", ("p_value", "pvalue")), ("q_bh", ("q_value", "bh_q"))):
        for alias in aliases:
            if target not in frame.columns and alias in frame.columns:
                frame = frame.rename(columns={alias: target})
    missing = {"chrom", "start", "end", "strand", "p"} - set(frame.columns)
    if missing:
        raise ValueError(f"window table lacks {sorted(missing)}")
    if frame.empty:
        raise ValueError("the window table is empty")
    if "window_id" not in frame.columns:
        frame["window_id"] = (frame.chrom + ":" + frame.start.astype(str) + "-" + frame.end.astype(str) + ":"
                              + frame.strand)
    if "z" not in frame.columns:
        frame["z"] = z_scores(frame.p)
    if "q_bh" not in frame.columns:
        frame["q_bh"] = bh_qvalues(frame.p)
    if "lfdr" not in frame.columns:
        frame["lfdr"] = window_lfdr(frame.p, null)
    frame.loc[frame.z < null[0], "lfdr"] = 1.0
    return frame


def top1_cutoff(pvalues: Iterable[float], fraction: float = TOP1_FRACTION) -> float:
    """Calibrated P value at the top-1% boundary: the P value of the ceil(fraction * n)-th window in a stable sort,
    as in `rlmhelix select cores`. Needs all eligible windows of the scan; for the paper's 7,869,821 windows it is
    TOP1_P."""
    p = np.sort(np.asarray(pvalues, dtype=np.float64), kind="mergesort")
    if len(p) == 0:
        raise ValueError("no windows")
    return float(p[int(np.ceil(len(p) * fraction)) - 1])


def window_statistics(candidates: pd.DataFrame, windows: pd.DataFrame, top1_p: float = TOP1_P) -> pd.DataFrame:
    """Scan statistics of each candidate from the windows overlapping its core on the same strand.

    windows: output of `prepare_windows`. Returns core_id, min_window_p (smallest calibrated P), max_window_z, lfdr
    (smallest), min_window_q_bh, bh_window (q <= 0.05), n_top1pct_windows (P <= top1_p; the default is the paper's
    cutoff, `top1_cutoff` gives it for another scan) and best_window_id (the window with the smallest P).
    """
    candidates = normalize_candidates(candidates)
    windows = windows.reset_index(drop=True)
    index = _interval_index(windows)
    p, z, q, lfdr, wid = (windows[c].to_numpy() for c in ("p", "z", "q_bh", "lfdr", "window_id"))
    rows = []
    for core, chrom, strand, start, end in zip(candidates.core_id, candidates.chrom, candidates.strand,
                                               candidates.core_start0, candidates.core_end0):
        hit = _overlapping(index, chrom, strand, start, end)
        if len(hit) == 0:
            rows.append({"core_id": core, "n_top1pct_windows": 0})
            continue
        best = hit[np.argmin(p[hit])]
        rows.append({"core_id": core, "n_top1pct_windows": int((p[hit] <= top1_p).sum()),
                     "best_window_id": wid[best], "min_window_p": float(p[hit].min()),
                     "max_window_z": float(z[hit].max()), "min_window_q_bh": float(q[hit].min()),
                     "lfdr": float(lfdr[hit].min())})
    out = pd.DataFrame(rows, columns=["core_id", "n_top1pct_windows", "best_window_id", "min_window_p",
                                      "max_window_z", "min_window_q_bh", "lfdr"])
    out["bh_window"] = pd.to_numeric(out.min_window_q_bh) <= BH_SUPPORT
    return out


# ----------------------------------------------------------------------------------------------------------------
# icSHAPE
# ----------------------------------------------------------------------------------------------------------------
_ATTRIBUTE = re.compile(r'(\S+) "([^"]+)"')
_BIN = 10_000


def read_gencode_transcripts(gtf: str | Path, transcripts_fasta: str | Path) -> tuple[dict, dict]:
    """Exons (in transcript order) and sequences (DNA alphabet) of GENCODE transcripts (release 26 for icSHAPE)."""
    records: dict[str, dict] = {}
    with _open(gtf) as handle:
        for line in handle:
            if line.startswith("#") or "\texon\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            attributes = dict(_ATTRIBUTE.findall(fields[8]))
            transcript = attributes.get("transcript_id")
            if not transcript or fields[6] not in "+-":
                continue
            record = records.setdefault(transcript, {"chrom": fields[0], "strand": fields[6], "exons": [],
                                                     "gene_name": attributes.get("gene_name", "")})
            record["exons"].append((int(fields[3]) - 1, int(fields[4])))
    for record in records.values():
        record["exons"] = sorted(record["exons"], reverse=record["strand"] == "-")
        record["length"] = sum(end - start for start, end in record["exons"])
    sequences = {name.split("|")[0]: sequence.upper().replace("U", "T")
                 for name, sequence in read_fasta(transcripts_fasta)}
    return records, sequences


def map_cores_to_transcripts(candidates: pd.DataFrame, records: dict, sequences: dict) -> pd.DataFrame:
    """Transcripts with one same-strand exon containing the whole core and an exactly matching sequence."""
    candidates = normalize_candidates(candidates)
    bins: dict[tuple, list] = defaultdict(list)
    for transcript, record in records.items():
        offset = 0
        for start, end in record["exons"]:
            for b in range(start // _BIN, (end - 1) // _BIN + 1):
                bins[(record["chrom"], record["strand"], b)].append((start, end, transcript, offset))
            offset += end - start
    rows = []
    for core, chrom, strand, start, end, sequence in zip(
            candidates.core_id, candidates.chrom, candidates.strand, candidates.core_start0, candidates.core_end0,
            candidates.human_core_rna_sequence):
        query = str(sequence).upper().replace("U", "T")
        if len(query) != end - start:
            raise ValueError(f"{core}: sequence length differs from the core interval")
        seen: set[str] = set()
        for exon_start, exon_end, transcript, exon_offset in bins.get((chrom, strand, start // _BIN), []):
            if transcript in seen or exon_start > start or exon_end < end:
                continue
            seen.add(transcript)
            offset = exon_offset + (start - exon_start if strand == "+" else exon_end - end)
            if sequences.get(transcript, "")[offset:offset + len(query)] == query:
                rows.append({"core_id": core, "transcript_id": transcript, "transcript_start": offset,
                             "transcript_end": offset + len(query), "transcript_length": records[transcript]["length"],
                             "gene_name": records[transcript]["gene_name"]})
    return pd.DataFrame(rows, columns=["core_id", "transcript_id", "transcript_start", "transcript_end",
                                       "transcript_length", "gene_name"])


def read_icshape_track(path: str | Path, transcripts: set[str]) -> dict[str, np.ndarray]:
    """Processed icSHAPE-pipe reactivities (GSE145805 '.out.txt': transcript, length, '*', one value per nucleotide;
    NULL = not measured) of the given transcripts."""
    out = {}
    with _open(path) as handle:
        for line in handle:
            name = line.split("\t", 1)[0]
            if name not in transcripts:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != int(fields[1]) + 3 or fields[2] != "*":
                raise ValueError(f"{path}: unexpected row for {name}")
            out[name] = np.array([np.nan if v == "NULL" else float(v) for v in fields[3:]])
    return out


def extract_reactivities(mapping: pd.DataFrame, tracks: Mapping[str, str | Path]) -> pd.DataFrame:
    """Reactivities of each core in each cell line, from the mapped isoform with the most measured nucleotides
    (ties: smallest transcript ID; reactivity values are not consulted).

    tracks: cell line -> GSE145805 '.out.txt.gz' file. Returns core_id, cell_line, transcript_id, n_measured,
    coverage and reactivity (array over the core, 5'->3').
    """
    wanted = set(mapping.transcript_id)
    rows = []
    for cell, path in tracks.items():
        values = read_icshape_track(path, wanted)
        options = []
        for m in mapping[mapping.transcript_id.isin(values)].itertuples(index=False):
            track = values[m.transcript_id]
            if len(track) != m.transcript_length:
                raise ValueError(f"{path}: length of {m.transcript_id} differs from the annotation")
            segment = track[m.transcript_start:m.transcript_end]
            options.append({"core_id": m.core_id, "cell_line": cell, "transcript_id": m.transcript_id,
                            "n_measured": int(np.isfinite(segment).sum()),
                            "coverage": float(np.isfinite(segment).mean()), "reactivity": segment})
        if options:
            chosen = pd.DataFrame(options).sort_values(["core_id", "n_measured", "transcript_id"],
                                                       ascending=[True, False, True]).drop_duplicates("core_id")
            rows.extend(chosen.to_dict("records"))
    return pd.DataFrame(rows, columns=["core_id", "cell_line", "transcript_id", "n_measured", "coverage", "reactivity"])


def auroc(labels: np.ndarray, values: np.ndarray) -> float:
    """Area under the ROC curve of `values` for boolean `labels` (ties count one half).

    Computed along the same floating-point path as scikit-learn's roc_auc_score (trapezoids over the ROC points),
    so that ties between observed and shifted AUROCs in `icshape_metrics` resolve as in the paper.
    """
    labels, values = np.asarray(labels, bool), np.asarray(values, float)
    order = np.argsort(values, kind="mergesort")[::-1]
    y, s = labels[order].astype(float), values[order]
    ends = np.r_[np.where(np.diff(s))[0], y.size - 1]
    tps = np.cumsum(y, dtype=np.float64)[ends]
    fps = 1 + ends - tps
    if len(fps) > 2:
        keep = np.where(np.r_[True, np.logical_or(np.diff(fps, 2), np.diff(tps, 2)), True])[0]
        fps, tps = fps[keep], tps[keep]
    fpr, tpr = np.r_[0, fps] / fps[-1], np.r_[0, tps] / tps[-1]
    return float(np.sum(np.diff(fpr) * (tpr[1:] + tpr[:-1]) / 2.0))


def icshape_metrics(structure: str, reactivity: np.ndarray, min_class: int = ICSHAPE_MIN_CLASS) -> dict:
    """AUROC of reactivity for unpaired versus paired nucleotides and its circular-shift control.

    NaN reactivities are not measured. The AUROC is computed when at least `min_class` paired and unpaired
    nucleotides are measured; each nonzero circular shift of the paired/unpaired labels that keeps both classes at
    that size gives a control AUROC, and P = (1 + #{control >= observed}) / (1 + #controls).
    """
    unpaired = np.fromiter((c == "." for c in structure), dtype=bool)
    reactivity = np.asarray(reactivity, float)
    if len(unpaired) != len(reactivity):
        raise ValueError("structure and reactivity differ in length")
    finite = np.isfinite(reactivity)
    labels, observed = unpaired[finite], reactivity[finite]
    out = {"n_measured": int(finite.sum()), "n_measured_unpaired": int(labels.sum()),
           "n_measured_paired": int((~labels).sum()), "reactivity_auroc": np.nan, "circular_shift_mean_auroc": np.nan,
           "circular_shift_empirical_p": np.nan, "n_valid_circular_shifts": 0}
    if out["n_measured_unpaired"] < min_class or out["n_measured_paired"] < min_class:
        return out
    auc = auroc(labels, observed)
    null = []
    for shift in range(1, len(structure)):
        shifted = np.roll(unpaired, shift)[finite]
        if shifted.sum() >= min_class and (~shifted).sum() >= min_class:
            null.append(auroc(shifted, observed))
    null = np.asarray(null)
    out.update(reactivity_auroc=auc, circular_shift_mean_auroc=float(null.mean()) if len(null) else np.nan,
               circular_shift_empirical_p=float((1 + (null >= auc).sum()) / (1 + len(null))),
               n_valid_circular_shifts=len(null))
    return out


def icshape_support(reactivities: pd.DataFrame, structures: Mapping[str, str],
                    coverage: float = ICSHAPE_COVERAGE) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per cell line metrics and per candidate summary.

    reactivities: output of `extract_reactivities`; structures: core_id -> dot-bracket. Cell lines with at least
    `coverage` of the core measured and a defined AUROC count as measured. The summary has icshape_n_contexts (cell
    lines), icshape_mean_auroc, icshape_frac_contexts_p05 (fraction with P <= 0.05), icshape_mean_shift_auroc and
    icshape_supported (that fraction >= 0.5), for the measured candidates.
    """
    rows = [{"core_id": r.core_id, "cell_line": r.cell_line, "transcript_id": r.transcript_id,
             "coverage": r.coverage, **icshape_metrics(structures[r.core_id], r.reactivity)}
            for r in reactivities.itertuples(index=False)]
    detail = pd.DataFrame(rows, columns=["core_id", "cell_line", "transcript_id", "coverage", "n_measured",
                                         "n_measured_unpaired", "n_measured_paired", "reactivity_auroc",
                                         "circular_shift_mean_auroc", "circular_shift_empirical_p",
                                         "n_valid_circular_shifts"])
    valid = detail[(detail.coverage >= coverage) & detail.reactivity_auroc.notna()]
    summary = valid.groupby("core_id", as_index=False).agg(
        icshape_n_contexts=("cell_line", "nunique"), icshape_mean_auroc=("reactivity_auroc", "mean"),
        icshape_frac_contexts_p05=("circular_shift_empirical_p", lambda x: x.le(ICSHAPE_P).mean()),
        icshape_mean_shift_auroc=("circular_shift_mean_auroc", "mean"))
    summary["icshape_supported"] = summary.icshape_frac_contexts_p05 >= ICSHAPE_FRACTION
    return detail, summary


def icshape_evidence(candidates: pd.DataFrame, gtf: str | Path, transcripts_fasta: str | Path,
                     tracks: Mapping[str, str | Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """icSHAPE detail and summary for the candidates (GENCODE v26 annotation and GSE145805 tracks)."""
    candidates = normalize_candidates(candidates)
    mapping = map_cores_to_transcripts(candidates, *read_gencode_transcripts(gtf, transcripts_fasta))
    reactivities = extract_reactivities(mapping, tracks)
    return icshape_support(reactivities, dict(zip(candidates.core_id, candidates.human_core_predicted_ss)))


# ----------------------------------------------------------------------------------------------------------------
# covariation (R-scape)
# ----------------------------------------------------------------------------------------------------------------
_OPEN, _CLOSE = "(<[{" + "ABCDEFGHIJKLMNOPQRSTUVWXYZ", ")>]}" + "abcdefghijklmnopqrstuvwxyz"
_ID = re.compile(r"[^A-Za-z0-9_.:-]+")
_MSA = re.compile(r"^# MSA \S+ nseq (\d+) \((\d+)\) alen (\d+) \((\d+)\) avgid ([0-9.]+) \(([0-9.]+)\) "
                  r"nbpairs (\d+) \((\d+)\)")
_GTP = re.compile(r"^# GTp\s+\S+\s+\[[^]]+\]\s+\[(\d+) \| (\d+) (\d+) (\d+) \|")
_EXPECTED = re.compile(r"^# BPAIRS expected to covary ([0-9.]+) \+/- ([0-9.]+)")
_HELIX = re.compile(r"^# RM_HELIX .* nbp = (\d+) nbp_cov = (\d+)")
_UNINFORMATIVE = "covariation scores are almost constant"


def structure_pairs(structure: str) -> list[tuple[int, int]]:
    """Base pairs (i < j, 0-based) of a dot-bracket structure; (), <>, [], {} and Aa... brackets."""
    partner = dict(zip(_CLOSE, _OPEN))
    stacks: dict[str, list[int]] = {c: [] for c in _OPEN}
    pairs = []
    for k, c in enumerate(structure):
        if c in stacks:
            stacks[c].append(k)
        elif c in partner:
            if not stacks[partner[c]]:
                raise ValueError(f"unbalanced structure at position {k}")
            pairs.append((stacks[partner[c]].pop(), k))
    if any(stacks.values()):
        raise ValueError("unbalanced structure")
    return sorted(pairs)


def _ungapped(sequence: str) -> str:
    return sequence.upper().replace("U", "T").replace("-", "").replace(".", "")


def unique_sequences(records: list[tuple[str, str]], human_id: str) -> list[tuple[str, str]]:
    """The human row first, then the other sequences without exact duplicates (of each other or of the human)."""
    human = [s for n, s in records if n == human_id]
    if len(human) != 1:
        raise ValueError(f"expected one row named {human_id}")
    seen = {_ungapped(human[0])}
    out = [(human_id, _ungapped(human[0]))]
    for name, sequence in records:
        sequence = _ungapped(sequence)
        if name != human_id and sequence not in seen:
            out.append((_ID.sub("_", name), sequence))
            seen.add(sequence)
    return out


def projected_stockholm(aligned: list[tuple[str, str]], human_id: str, pairs: list[tuple[int, int]]) -> str:
    """Stockholm text of an alignment with the human structure (pairs over ungapped human positions) as SS_cons and
    the human residues as RF."""
    human = [s for n, s in aligned if n == human_id]
    if len(human) != 1:
        raise ValueError(f"expected one row named {human_id}")
    columns = [k for k, c in enumerate(human[0]) if c not in GAPS]
    if not pairs or max(max(p) for p in pairs) >= len(columns):
        raise ValueError(f"{human_id}: structure does not fit the human row")
    ss = ["."] * len(human[0])
    for i, j in pairs:
        ss[columns[i]], ss[columns[j]] = "<", ">"
    rf = "".join("x" if c not in GAPS else "." for c in human[0])
    width = max(max(len(n) for n, _ in aligned), len("#=GC SS_cons")) + 2
    lines = ["# STOCKHOLM 1.0"] + [f"{n:<{width}}{s.upper().replace('T', 'U')}" for n, s in aligned]
    lines += [f"{'#=GC SS_cons':<{width}}{''.join(ss)}", f"{'#=GC RF':<{width}}{rf}", "//", ""]
    return "\n".join(lines)


def read_stockholm_rows(path: str | Path) -> tuple[list[tuple[str, str]], str]:
    """Rows and SS_cons of a single-block or interleaved Stockholm alignment."""
    rows: dict[str, str] = {}
    ss = ""
    for line in open(path, encoding="latin-1"):
        line = line.rstrip("\n")
        if line.startswith("//"):
            break
        if line.startswith("#=GC SS_cons"):
            ss += line.split(maxsplit=2)[2]
        elif line.strip() and not line.startswith("#"):
            name, sequence = line.split()[:2]
            rows[name] = rows.get(name, "") + sequence
    return list(rows.items()), ss


def parse_rscape(log: str, status: str = "success") -> dict:
    """Counts from R-scape's standard output (`-s`): sequences and mean identity used, proposed pairs tested and
    significant (E < 0.05), covarying pairs expected from the power analysis, and helices with >= 2 covarying pairs."""
    out: dict = {"status": status}
    helices = []
    for line in log.splitlines():
        if m := _MSA.match(line):
            out.update(nseq=int(m.group(1)), avg_identity=float(m.group(5)), n_pairs_input=int(m.group(8)))
        if m := _GTP.match(line):
            out.update(n_significant_other=int(m.group(1)), n_significant_pairs=int(m.group(2)),
                       n_pairs_tested=int(m.group(3)))
        if m := _EXPECTED.match(line):
            out["expected_pairs"] = float(m.group(1))
        if m := _HELIX.match(line):
            helices.append(int(m.group(2)))
    if "n_significant_pairs" not in out:
        if status != "uninformative":
            raise ValueError("no GTp line in the R-scape output")
        out.update(n_significant_other=0, n_significant_pairs=0, n_pairs_tested=out.get("n_pairs_input", 0))
    out["n_helices_with_2_cov_pairs"] = sum(n >= 2 for n in helices)
    return out


def run_rscape(stockholm: str | Path, outdir: str | Path, rscape: str = "R-scape") -> str:
    """Run R-scape 2.6.16 (`-s --nofigures`); keeps its standard output as rscape.log. Returns the status.

    R-scape runs in `outdir` and is given it as an absolute path (it resolves a relative --outdir from its working
    directory)."""
    outdir = Path(outdir).resolve()
    result = subprocess.run([rscape, "-s", "--nofigures", "--outdir", str(outdir), "--outname", "rscape",
                             str(Path(stockholm).resolve())], text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, cwd=outdir)
    (outdir / "rscape.log").write_text(result.stdout)
    if result.returncode == 0 and (outdir / "rscape.cov").exists():
        return "success"
    return "uninformative" if _UNINFORMATIVE in result.stdout else "failed"


def covariation_one(core_id: str, condition: str, source: str | Path, structure: str | None, workdir: str | Path,
                    mafft: str = "mafft", rscape: str = "R-scape", resume: bool = True) -> dict:
    """Covariation of one candidate under condition i, ii or iii.

    source: for (iii), a FASTA of the human core (named core_id) and its high-confidence homologues; for (i) and
    (ii), the covariance-model alignment (Stockholm, human row named core_id). (i) runs R-scape on the alignment as
    it is; (ii) and (iii) align the unique sequences with MAFFT (`--auto`) and project `structure` (the predicted
    core structure) through the human row.
    """
    outdir = Path(workdir) / condition / core_id
    outdir.mkdir(parents=True, exist_ok=True)
    summary_path = outdir / "summary.json"
    if resume and summary_path.exists():
        previous = json.loads(summary_path.read_text())
        if previous.get("status") in ("success", "uninformative"):
            return previous
    if condition == "i":
        stockholm = Path(source)
    else:
        records = read_fasta(source) if condition == "iii" else read_stockholm_rows(source)[0]
        records = unique_sequences(records, core_id)
        fasta = outdir / "input.fa"
        fasta.write_text("".join(f">{n}\n{s}\n" for n, s in records))
        with open(outdir / "mafft.fa", "w") as stdout, open(outdir / "mafft.log", "w") as stderr:
            aligned = subprocess.run([mafft, "--auto", "--thread", "1", "--quiet", str(fasta.resolve())],
                                     stdout=stdout, stderr=stderr, cwd=outdir)
        if aligned.returncode != 0:
            summary = {"core_id": core_id, "condition": condition, "status": "failed"}
            summary_path.write_text(json.dumps(summary) + "\n")
            return summary
        stockholm = outdir / "mafft.sto"
        stockholm.write_text(projected_stockholm(read_fasta(outdir / "mafft.fa"), core_id, structure_pairs(structure)))
    status = run_rscape(stockholm, outdir, rscape)
    summary = {"core_id": core_id, "condition": condition, "status": status}
    if status != "failed":
        summary.update(parse_rscape((outdir / "rscape.log").read_text(), status))
    summary_path.write_text(json.dumps(summary) + "\n")
    return summary


def read_covariation(workdir: str | Path, condition: str, core_ids: Iterable[str]) -> pd.DataFrame:
    """Parse existing R-scape outputs (<workdir>/<condition>/<core_id>/rscape.log) without running anything."""
    rows = []
    for core in core_ids:
        outdir = Path(workdir) / condition / core
        log = outdir / "rscape.log"
        if not log.exists():
            rows.append({"core_id": core, "condition": condition, "status": "missing"})
            continue
        text = log.read_text()
        status = ("success" if (outdir / "rscape.cov").exists() else
                  "uninformative" if _UNINFORMATIVE in text else "failed")
        row = {"core_id": core, "condition": condition, "status": status}
        if status != "failed":
            row.update(parse_rscape(text, status))
        rows.append(row)
    return covariation_table(rows)


def covariation_table(rows: list[dict]) -> pd.DataFrame:
    """Covariation results with the calls: detected (>= 1 significant proposed pair), two_pairs (>= 2) and
    stringent (>= 3 with 2 in one helix); calls are missing where R-scape did not run or failed."""
    frame = pd.DataFrame(rows)
    for column in ("core_id", "condition", "status", "nseq", "avg_identity", "n_pairs_tested", "n_significant_pairs",
                   "expected_pairs", "n_helices_with_2_cov_pairs"):
        if column not in frame.columns:
            frame[column] = np.nan
    ok = frame.status.isin(["success", "uninformative"])
    n = frame.n_significant_pairs
    frame["detected"] = (n >= 1).where(ok)
    frame["two_pairs"] = (n >= 2).where(ok)
    frame["stringent"] = ((n >= 3) & (frame.n_helices_with_2_cov_pairs >= 1)).where(ok)
    return frame[["core_id", "condition", "status", "nseq", "avg_identity", "n_pairs_tested", "n_significant_pairs",
                  "expected_pairs", "n_helices_with_2_cov_pairs", "detected", "two_pairs", "stringent"]]


def _covariation_task(task: tuple) -> dict:
    return covariation_one(*task)


def covariation_inputs(directory: str | Path, core_ids: Iterable[str], condition: str) -> dict[str, Path]:
    """Input file of each candidate in `directory`, for the candidates that have one.

    Condition (iii): <core_id>.fa (the `homologues/` directory of `rlmhelix select expand`). Conditions (i) and (ii):
    <core_id>.sto, or <core_id>/final55.sto as written by `select expand` (its `candidates/` directory).
    """
    directory = Path(directory)
    names = ("{}.fa",) if condition == "iii" else ("{}.sto", "{}/final55.sto")
    found = {}
    for core in core_ids:
        for name in names:
            path = directory / name.format(core)
            if path.is_file():
                found[core] = path
                break
    return found


# ----------------------------------------------------------------------------------------------------------------
# eCLIP, splice sites, Rfam and homologue scores (descriptive columns)
# ----------------------------------------------------------------------------------------------------------------
def eclip_overlaps(candidates: pd.DataFrame, experiments: pd.DataFrame) -> pd.DataFrame:
    """Same-strand overlap (>= 1 nt) of each core with the IDR peaks of each eCLIP experiment.

    experiments: one row per experiment with rbp, cell and path (BED with chrom, start, end, name, score, strand).
    """
    candidates = normalize_candidates(candidates)
    start, end = candidates.core_start0.to_numpy(), candidates.core_end0.to_numpy()
    bound = np.zeros((len(candidates), len(experiments)), bool)
    for k, experiment in enumerate(experiments.itertuples(index=False)):
        peaks = pd.read_csv(experiment.path, sep="\t", header=None, usecols=range(6),
                            names=["chrom", "start", "end", "name", "score", "strand"])
        index = _interval_index(peaks)
        for row, (chrom, strand) in enumerate(zip(candidates.chrom, candidates.strand)):
            bound[row, k] = len(_overlapping(index, chrom, strand, start[row], end[row])) > 0
    rbp, cell = experiments.rbp.to_numpy(), experiments.cell.to_numpy()
    rbps = [sorted({rbp[k] for k in np.nonzero(row)[0]}) for row in bound]
    return pd.DataFrame({
        "core_id": candidates.core_id, "eclip_n_experiments": bound.sum(1),
        "eclip_n_rbps": [len(r) for r in rbps], "eclip_rbps": [";".join(r) for r in rbps],
        "eclip_experiments": [";".join(f"{rbp[k]} ({cell[k]})" for k in np.nonzero(row)[0]) for row in bound],
        "eclip_dsrbp": [any(x in DSRBP for x in r) for r in rbps],
        "eclip_dsrbp_names": [";".join(x for x in r if x in DSRBP) for r in rbps]})


def eclip_experiments(table: pd.DataFrame) -> pd.DataFrame:
    """eCLIP experiments for `eclip_overlaps`, in table order (the order of eclip_experiments).

    table: a manifest with rbp, cell and path or file_accession, or the Accessions sheet of Supplementary Table S8
    (target and cell_line; only the rows of dataset 'eCLIP peaks' are kept).
    """
    frame = table[table.dataset == "eCLIP peaks"] if "dataset" in table.columns else table
    frame = frame.rename(columns={alias: target for target, alias in (("rbp", "target"), ("cell", "cell_line"))
                                  if target not in frame.columns and alias in frame.columns})
    if {"rbp", "cell"} - set(frame.columns) or not {"path", "file_accession"} & set(frame.columns):
        raise ValueError("eCLIP table needs rbp (or target), cell (or cell_line) and path or file_accession")
    if frame.empty:
        raise ValueError("no eCLIP experiments in the table")
    return frame.reset_index(drop=True)


def splice_sites(gtf: str | Path) -> dict[tuple[str, str], np.ndarray]:
    """Splice sites (exon ends and next exon starts, 0-based) of all annotated transcripts, by chrom and strand."""
    exons = defaultdict(list)
    with _open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 9)
            if fields[2] == "exon":
                transcript = re.search('transcript_id "([^"]+)"', fields[8]).group(1)
                exons[(fields[0], fields[6], transcript)].append((int(fields[3]) - 1, int(fields[4])))
    sites = defaultdict(set)
    for (chrom, strand, _), intervals in exons.items():
        intervals.sort()
        for (_, left_end), (right_start, _) in zip(intervals, intervals[1:]):
            sites[(chrom, strand)].update((left_end, right_start))
    return {key: np.array(sorted(values)) for key, values in sites.items()}


def splice_site_distance(candidates: pd.DataFrame, sites: dict[tuple[str, str], np.ndarray]) -> pd.DataFrame:
    """Distance from each core to the nearest same-strand splice site (0 if the core spans one) and its class."""
    candidates = normalize_candidates(candidates)
    distance = np.full(len(candidates), np.inf)
    for k, (chrom, strand, start, end) in enumerate(zip(candidates.chrom, candidates.strand, candidates.core_start0,
                                                        candidates.core_end0)):
        a = sites.get((chrom, strand))
        if a is None:
            continue
        j = np.searchsorted(a, start)
        if j < len(a) and a[j] <= end:
            distance[k] = 0
            continue
        options = ([a[j] - end] if j < len(a) else []) + ([start - a[j - 1]] if j > 0 else [])
        distance[k] = min(options) if options else np.inf
    return pd.DataFrame({"core_id": candidates.core_id, "splice_site_distance_nt": distance,
                         "splice_site_class": np.select([distance == 0, distance <= 100, np.isfinite(distance)],
                                                        ["spans splice site", "within 100 nt", "farther"],
                                                        "no same-strand splice site")})


_GTF_ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')
_GTF_TAG = re.compile(r'tag "([^"]*)"')
_GENCODE_FEATURES = ("gene", "transcript", "exon", "CDS", "UTR")


def _union_length(spans: Iterable[tuple[int, int]]) -> int:
    """Number of positions covered by half-open intervals."""
    total, right = 0, -1
    for start, end in sorted(spans):
        total += max(0, end - max(start, right))
        right = max(right, end)
    return total


def gencode_annotation(candidates: pd.DataFrame, gtf: str | Path) -> pd.DataFrame:
    """GENCODE genes and transcripts overlapping each core (release 50 in the paper; Supplementary Table S4).

    Gene, transcript, exon, CDS and UTR records overlapping the core by >= 1 nt are collected on both strands; the
    overlaps of each feature type of a transcript are united. The context of a transcript lists CDS, 5UTR and 3UTR
    (UTR segments before and after the CDS of that transcript, by strand; UTR_unresolved otherwise),
    exon_non_CDS_UTR (exonic overlap without CDS or UTR) and intron (transcript overlap beyond its exons).

    Returns core_id and GENCODE_COLUMNS: annotation_status (same_strand_exonic, else same_strand_intronic,
    same_strand_gene_only, opposite_strand_only or no_gene_overlap; combined over transcripts),
    same_strand_gene_names and same_strand_gene_types (genes on the strand of the core), same_strand_contexts (union
    of the contexts of the same-strand transcripts), pseudogene_gene_names (same-strand genes whose type contains
    'pseudogene'), and the representative same-strand transcript (largest exonic overlap, then MANE_Select, then
    largest transcript overlap, then smallest transcript ID; a display choice, not the expressed isoform) with its
    name, gene name, gene type and context. Lists are unique, sorted and ';'-separated; empty (None) where nothing
    overlaps.
    """
    candidates = normalize_candidates(candidates)
    cores: dict[str, tuple[list[int], list[tuple]]] = {}
    for chrom, part in candidates.groupby("chrom", sort=False):
        part = part.sort_values("core_start0", kind="mergesort")
        cores[chrom] = (part.core_start0.tolist(), list(zip(part.index, part.core_start0, part.core_end0, part.strand)))
    longest = int((candidates.core_end0 - candidates.core_start0).max()) if len(candidates) else 0

    def overlapping(chrom: str, start: int, end: int) -> list[tuple]:
        starts, rows = cores.get(chrom, ([], []))
        return [r for r in rows[bisect.bisect_right(starts, start - longest):bisect.bisect_left(starts, end)]
                if r[2] > start]

    genes: dict[int, list[tuple[str, str, bool]]] = defaultdict(list)  # core -> (gene name, gene type, same strand)
    transcripts: dict[str, dict] = {}
    cds: dict[str, tuple[int, int]] = {}  # CDS bounds of the transcripts overlapping a core
    hits: dict[tuple[int, str], dict] = {}
    with _open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in _GENCODE_FEATURES:
                continue
            feature, strand, start, end = fields[2], fields[6], int(fields[3]) - 1, int(fields[4])
            attributes = None
            if feature == "CDS":
                attributes = dict(_GTF_ATTRIBUTE.findall(fields[8]))
                transcript = attributes.get("transcript_id", "")
                if transcript in transcripts:
                    lo, hi = cds.get(transcript, (start, end))
                    cds[transcript] = (min(lo, start), max(hi, end))
            rows = overlapping(fields[0], start, end)
            if not rows:
                continue
            attributes = attributes or dict(_GTF_ATTRIBUTE.findall(fields[8]))
            if feature == "gene":
                for row in rows:
                    genes[row[0]].append((attributes.get("gene_name", ""), attributes.get("gene_type", ""),
                                          row[3] == strand))
                continue
            transcript = attributes.get("transcript_id", "")
            if feature == "transcript":
                transcripts[transcript] = {
                    "transcript_name": attributes.get("transcript_name", ""),
                    "gene_name": attributes.get("gene_name", ""), "gene_type": attributes.get("gene_type", ""),
                    "strand": strand, "mane": "MANE_Select" in ";".join(_GTF_TAG.findall(fields[8]))}
            for row, core_start, core_end, core_strand in rows:
                hit = hits.setdefault((row, transcript), {"same": core_strand == strand, "spans": defaultdict(list)})
                hit["spans"][feature].append((max(start, core_start), min(end, core_end)))

    same_strand: dict[int, list[dict]] = defaultdict(list)
    for (row, transcript), hit in hits.items():
        if transcript not in transcripts:
            raise ValueError(f"{gtf}: no transcript record for {transcript}")
        if not hit["same"]:
            continue
        meta, spans = transcripts[transcript], hit["spans"]
        bp = {k: _union_length(spans[k]) for k in ("transcript", "exon", "CDS", "UTR")}
        labels, five, three = (["CDS"] if bp["CDS"] else []), [], []
        if bp["UTR"]:
            if transcript in cds:
                lo, hi = cds[transcript]
                for a, b in spans["UTR"]:
                    if b <= lo:
                        (five if meta["strand"] == "+" else three).append((a, b))
                    elif a >= hi:
                        (three if meta["strand"] == "+" else five).append((a, b))
            labels += (["5UTR"] if five else []) + (["3UTR"] if three else []) + ([] if five or three
                                                                                 else ["UTR_unresolved"])
        if bp["exon"] and not bp["CDS"] and not bp["UTR"]:
            labels.append("exon_non_CDS_UTR")
        if bp["transcript"] > bp["exon"]:
            labels.append("intron")
        same_strand[row].append({"transcript_id": transcript, **meta, "context": ";".join(labels),
                                 "exon_bp": bp["exon"], "transcript_bp": bp["transcript"]})

    def join(values: Iterable[str]) -> str | None:
        return ";".join(sorted({v for v in values if v})) or None

    rows = []
    for row, core in enumerate(candidates.core_id):
        same, near = same_strand.get(row, []), genes.get(row, [])
        same_genes = [g for g in near if g[2]]
        best = min(same, key=lambda t: (-t["exon_bp"], not t["mane"], -t["transcript_bp"], t["transcript_id"]),
                   default={})
        status = ("same_strand_exonic" if any(t["exon_bp"] > 0 for t in same) else "same_strand_intronic" if same
                  else "same_strand_gene_only" if same_genes else "opposite_strand_only" if near
                  else "no_gene_overlap")
        rows.append({"core_id": core, "annotation_status": status,
                     "same_strand_gene_names": join(g[0] for g in same_genes),
                     "same_strand_gene_types": join(g[1] for g in same_genes),
                     "same_strand_contexts": join(c for t in same for c in t["context"].split(";")),
                     "representative_gene_name": best.get("gene_name") or None,
                     "representative_gene_type": best.get("gene_type") or None,
                     "representative_transcript_name": best.get("transcript_name") or None,
                     "representative_context": best.get("context") or None,
                     "pseudogene_gene_names": join(g[0] for g in same_genes if "pseudogene" in g[1])})
    return pd.DataFrame(rows, columns=["core_id", *GENCODE_COLUMNS])


def rfam_overlaps(candidates: pd.DataFrame, hits: pd.DataFrame) -> pd.DataFrame:
    """Rfam families with a significant hit overlapping each core on the same strand by >= 1 nt.

    hits: chrom, start (0-based), end, strand, rfam_acc, rfam_id (e.g. the significant 'full' hits of Rfam 15.1 on
    the GRCh38 chromosomes). Descriptive only; Rfam overlap is not used for the tiers.
    """
    candidates = normalize_candidates(candidates)
    hits = hits.reset_index(drop=True)
    index = _interval_index(hits)
    accessions, names = [], []
    for chrom, strand, start, end in zip(candidates.chrom, candidates.strand, candidates.core_start0,
                                         candidates.core_end0):
        found = hits.loc[_overlapping(index, chrom, strand, start, end)]
        found = found.drop_duplicates("rfam_acc").sort_values("rfam_acc")
        accessions.append(";".join(found.rfam_acc) if len(found) else None)
        names.append(";".join(found.rfam_id.astype(str)) if len(found) else None)
    return pd.DataFrame({"core_id": candidates.core_id, "rfam_accessions": accessions, "rfam_names": names})


def homologue_columns(scores: pd.DataFrame) -> pd.DataFrame:
    """Structural signal in homologues: RLMHelix-MSA of the condition (iii) alignment against 19 SISSIz shuffles.

    scores: core_id (or query_id), human_row ('without' or 'with'), native, z, empirical_p, passes, and optionally
    score (rows other than 'RLMHelix-MSA' are ignored) and n_valid_null (19 for a complete set of shuffles).
    """
    scores = scores.rename(columns={"query_id": "core_id"})
    if "score" in scores.columns:
        scores = scores[scores.score == "RLMHelix-MSA"]
    without = scores[scores.human_row == "without"].set_index("core_id")
    with_human = scores[scores.human_row == "with"].set_index("core_id")
    if "n_valid_null" in without.columns:
        without = without[without.n_valid_null == 19]
    out = pd.DataFrame({"ms_wo_native": without.native, "ms_wo_z": without.z, "ms_wo_p": without.empirical_p,
                        "ms_wo_pass": without.passes.astype(bool)})
    out = out.join(pd.DataFrame({"ms_wi_z": with_human.z, "ms_wi_pass": with_human.passes.astype(bool)}),
                   how="outer")
    out["ms_scored"] = out.index.isin(without.index)
    return out.rename_axis("core_id").reset_index()


# ----------------------------------------------------------------------------------------------------------------
# tiers and the catalogue table
# ----------------------------------------------------------------------------------------------------------------
def evidence_tiers(catalogue: pd.DataFrame) -> pd.DataFrame:
    """Add the covariation status, statistical and structural support, the tier and the evidence behind it.

    Needs lfdr, bh_window, cov_detected (missing where not evaluated), cov_nseq, cov_expected_pairs_power and the
    icSHAPE summary columns (missing for unmeasured candidates).
    """
    cat = catalogue.copy()
    detected = cat.cov_detected.astype("boolean")
    power = pd.to_numeric(cat.cov_expected_pairs_power)
    cat["cov_status"] = np.select(
        [detected.fillna(False).to_numpy(bool),
         (detected.isna() | cat.cov_nseq.isna() | power.isna()).to_numpy(bool),
         (power < 1).to_numpy(bool)], COV_STATUS[:3], COV_STATUS[3])
    cat["cov_detected"] = detected.fillna(False).astype(bool)
    cat["icshape_measured"] = cat.icshape_n_contexts.notna()
    cat["icshape_supported"] = (pd.to_numeric(cat.icshape_frac_contexts_p05) >= ICSHAPE_FRACTION).fillna(False)
    cat["bh_window"] = cat.bh_window.astype("boolean").fillna(False).astype(bool)
    cat["stat_support"] = cat.bh_window | (cat.lfdr <= LFDR_SUPPORT).fillna(False)
    cat["structure_support"] = cat.icshape_supported | cat.cov_detected
    k = cat.stat_support.astype(int) + cat.structure_support.astype(int)
    cat["tier"] = np.select([k == 2, k == 1], ["A", "B"], "C")
    cat["tier_evidence"] = [
        "+".join(name for name, flag in (("statistics", s), ("icSHAPE", i), ("covariation", c)) if flag) or "none"
        for s, i, c in zip(cat.stat_support, cat.icshape_supported, cat.cov_detected)]
    cat["structure_evidence_measured"] = cat.icshape_measured | cat.cov_status.isin([COV_STATUS[0], COV_STATUS[3]])
    return cat


def _covariation_columns(covariation: pd.DataFrame) -> pd.DataFrame:
    """Working columns of the catalogue from a covariation table with conditions i, ii and/or iii."""
    covariation = covariation.rename(columns={"query_id": "core_id"})
    out = pd.DataFrame({"core_id": pd.unique(covariation.core_id)}).set_index("core_id")
    main = covariation[covariation.condition.astype(str) == "iii"].set_index("core_id")
    for target, source in (("cov_detected", "detected"), ("cov_two_pairs", "two_pairs"),
                           ("cov_n_significant_pairs", "n_significant_pairs"),
                           ("cov_n_proposed_pairs_tested", "n_pairs_tested"),
                           ("cov_expected_pairs_power", "expected_pairs"), ("cov_nseq", "nseq"),
                           ("cov_avg_identity", "avg_identity")):
        out[target] = main[source]
    for condition, target in (("i", "cov_detected_cm_aligned"), ("ii", "cov_detected_same_seqs_mafft")):
        other = covariation[covariation.condition.astype(str) == condition].set_index("core_id")
        out[target] = other.detected
    return out.reset_index()


def build_catalogue(candidates: pd.DataFrame, windows: pd.DataFrame, icshape: pd.DataFrame | None = None,
                    covariation: pd.DataFrame | None = None, eclip: pd.DataFrame | None = None,
                    splice: pd.DataFrame | None = None, rfam: pd.DataFrame | None = None,
                    homologues: pd.DataFrame | None = None, annotation: pd.DataFrame | None = None,
                    top1_p: float = TOP1_P) -> pd.DataFrame:
    """One row per candidate with its scan statistics, evidence and tier (working column names; see `to_table`).

    candidates: `normalize_candidates`; windows: `prepare_windows`; icshape: summary of `icshape_support`;
    covariation: `covariation_table` (conditions i, ii, iii); eclip: `eclip_overlaps`; splice:
    `splice_site_distance`; rfam: `rfam_overlaps` (replaces the candidate table's Rfam columns); homologues:
    `homologue_columns`; annotation: `gencode_annotation` (replaces the candidate table's GENCODE columns); top1_p:
    cutoff of n_top1pct_windows (default: the paper's, TOP1_P; see `top1_cutoff`). Missing evidence is treated as
    not measured, never as evidence against a candidate. Sorted by tier, lfdr and the smallest window P value.
    """
    cat = normalize_candidates(candidates).merge(window_statistics(candidates, windows, top1_p), on="core_id",
                                                 how="left")
    if "scan_pvalue" not in cat.columns:
        cat["scan_pvalue"] = cat.min_window_p
    cov = _covariation_columns(covariation) if covariation is not None else pd.DataFrame(
        columns=["core_id", "cov_detected", "cov_nseq", "cov_expected_pairs_power"])
    cat = cat.merge(cov, on="core_id", how="left")
    ic = icshape if icshape is not None else pd.DataFrame(
        columns=["core_id", "icshape_n_contexts", "icshape_frac_contexts_p05"])
    ic = ic.rename(columns={"query_id": "core_id"}).drop(columns=["icshape_supported"], errors="ignore")
    cat = cat.merge(ic, on="core_id", how="left")
    if rfam is not None:
        cat = cat.drop(columns=["rfam_accessions", "rfam_names"], errors="ignore")
    if annotation is not None:
        cat = cat.drop(columns=GENCODE_COLUMNS, errors="ignore")
    for extra in (eclip, splice, rfam, homologues, annotation):
        if extra is not None:
            cat = cat.merge(extra, on="core_id", how="left")
    if "ms_scored" in cat.columns:
        cat["ms_scored"] = cat.ms_scored.astype("boolean").fillna(False).astype(bool)
    if "pseudogene_gene_names" in cat.columns:
        names = cat.pseudogene_gene_names
        cat["pseudogene_context"] = names.notna() & (names.astype(str).str.len() > 0)
    cat = evidence_tiers(cat)
    return cat.sort_values(["tier", "lfdr", "min_window_p"], kind="mergesort").reset_index(drop=True)


def to_table(catalogue: pd.DataFrame) -> pd.DataFrame:
    """The catalogue with the columns of Supplementary Table S4 (Catalogue sheet), in that order. Columns whose
    evidence was not supplied are left empty."""
    return pd.DataFrame({public: catalogue[working] if working in catalogue.columns else pd.NA
                         for public, working in S4_COLUMNS})


def families_table(catalogue: pd.DataFrame) -> pd.DataFrame:
    """One row per RLM-Fam family (Families sheet of Supplementary Table S4), sorted by best tier and lowest lfdr.

    n_overlapping_earlier_screens counts candidates whose `overlaps_earlier_screens` column is true (left empty
    when the column is absent).
    """
    f = catalogue.assign(tier_rank=catalogue.tier.map(TIER_ORDER),
                         member=catalogue.browser_locus + "(" + catalogue.strand + ")")
    f = f.sort_values(["rlmfam_id", "tier_rank", "lfdr", "scan_pvalue"], kind="mergesort")

    def unique_join(values: pd.Series) -> str:
        out = []
        for value in values.dropna().astype(str):
            for item in value.split(";"):
                if item.strip() and item.strip() not in out:
                    out.append(item.strip())
        return ";".join(out)

    def column(name: str) -> pd.Series:
        return f[name] if name in f.columns else pd.Series(pd.NA, index=f.index, dtype="object")

    g = f.assign(_genes=column("representative_gene_name"), _contexts=column("representative_context"),
                 _rfam=column("rfam_names")).groupby("rlmfam_id", sort=False)
    out = pd.DataFrame({
        "rlmfam_id": list(g.groups), "n_candidates": g.size().to_numpy(),
        "candidate_loci": g.member.agg(";".join).to_numpy(),
        "representative_locus_hg38": g.browser_locus.first().to_numpy(),
        "representative_strand": g.strand.first().to_numpy(), "genes": g._genes.agg(unique_join).to_numpy(),
        "region_contexts": g._contexts.agg(unique_join).to_numpy(), "best_tier": g.tier.first().to_numpy(),
        "min_lfdr": g.lfdr.min().to_numpy(), "n_statistical_support": g.stat_support.sum().to_numpy(),
        "n_structural_support": g.structure_support.sum().to_numpy(),
        "n_overlapping_earlier_screens": (g.overlaps_earlier_screens.sum().to_numpy()
                                          if "overlaps_earlier_screens" in f.columns else pd.NA),
        "rfam": g._rfam.agg(unique_join).to_numpy()})
    out = out.assign(_rank=out.best_tier.map(TIER_ORDER)).sort_values(["_rank", "min_lfdr", "rlmfam_id"])
    return out.drop(columns="_rank").reset_index(drop=True)[FAMILY_COLUMNS]


def tier_summary(catalogue: pd.DataFrame) -> pd.DataFrame:
    """Counts of evidence per tier (Tier_summary sheet of Supplementary Table S4)."""
    rows = []
    for tier, g in [("all", catalogue)] + list(catalogue.groupby("tier")):
        def count(name: str):
            return int(g[name].fillna(False).astype(bool).sum()) if name in g.columns else pd.NA
        rows.append({"tier": tier, "n_candidates": len(g), "fraction_of_catalogue": len(g) / len(catalogue),
                     "statistical_support": int(g.stat_support.sum()), "bh_window": int(g.bh_window.sum()),
                     "lfdr_le_0.2": int((g.lfdr <= LFDR_SUPPORT).sum()),
                     "icshape_measured": int(g.icshape_measured.sum()),
                     "icshape_supported": int(g.icshape_supported.sum()), "covariation_iii": int(g.cov_detected.sum()),
                     "eclip_dsrbp": count("eclip_dsrbp"),
                     "eclip_any": int((g.eclip_n_experiments > 0).sum()) if "eclip_n_experiments" in g else pd.NA,
                     "high_copy": count("high_copy"), "median_lfdr": float(g.lfdr.median())})
    return pd.DataFrame(rows)


def evidence_summary(catalogue: pd.DataFrame) -> pd.DataFrame:
    """Numbers of candidates with each kind of evidence (Evidence_summary sheet of Supplementary Table S4)."""
    cat, n = catalogue, len(catalogue)
    evaluable = cat.cov_status != COV_STATUS[1]
    rows = [("Statistical support (BH window or lfdr <= 0.2)", n, int(cat.stat_support.sum())),
            ("  BH q <= 0.05 window", n, int(cat.bh_window.sum())),
            ("  lfdr <= 0.2", n, int((cat.lfdr <= LFDR_SUPPORT).sum())),
            ("icSHAPE supported (of measured)", int(cat.icshape_measured.sum()), int(cat.icshape_supported.sum())),
            ("Covariation detected, structure-free alignment (of evaluable)", int(evaluable.sum()),
             int(cat.cov_detected.sum())),
            ("  evaluable with expected >= 1 covarying pair", int(evaluable.sum()),
             int(cat.cov_status.isin([COV_STATUS[0], COV_STATUS[3]]).sum()))]
    if "eclip_dsrbp" in cat.columns:
        rows += [("dsRBP eCLIP peak on the core (6 dsRBPs)", n, int(cat.eclip_dsrbp.fillna(False).sum())),
                 ("Any eCLIP peak on the core", n, int((cat.eclip_n_experiments > 0).sum()))]
    if "splice_site_class" in cat.columns:
        rows.append(("Core spans a splice site", n, int((cat.splice_site_class == "spans splice site").sum())))
    out = pd.DataFrame(rows, columns=["evidence", "denominator", "n_with_evidence"])
    return out.assign(fraction=out.n_with_evidence / out.denominator)


def working_columns(table: pd.DataFrame) -> pd.DataFrame:
    """A table with the columns of Supplementary Table S4 (e.g. the candidate table of the Zenodo deposit, once public) renamed to the
    working names of `build_catalogue`."""
    return table.rename(columns={public: working for public, working in S4_COLUMNS
                                 if public in table.columns and working not in table.columns})


def _wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return np.nan, np.nan
    p = k / n
    centre = (p + z * z / (2 * n)) / (1 + z * z / n)
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return centre - half, centre + half


def evidence_concordance(catalogue: pd.DataFrame) -> pd.DataFrame:
    """Independent evidence at candidates with and without statistical support (Supplementary Fig. S8c).

    Covariation under condition (iii) among evaluable candidates, overall and by R-scape power class (expected
    covarying pairs in (-inf, 1], (1, 3], (3, 10] and above 10, labelled '<1', '1-3', '3-10', '>=10' as in the
    paper's tables), icSHAPE support among measured candidates and eCLIP binding (when present): numbers, rates with
    Wilson 95% intervals, odds ratios and two-sided Fisher's exact P values. Accepts the working table or the
    columns of Supplementary Table S4.
    """
    cat = working_columns(catalogue)
    rows = []

    def compare(evidence: str, stratum: str, sub: pd.DataFrame, flags: pd.Series) -> None:
        support = sub.stat_support.astype(bool).to_numpy()
        flags = flags.fillna(False).astype(bool).to_numpy()
        ka, na = int(flags[support].sum()), int(support.sum())
        kb, nb = int(flags[~support].sum()), int((~support).sum())
        odds, p = stats.fisher_exact([[ka, na - ka], [kb, nb - kb]])
        (la, ha), (lb, hb) = _wilson(ka, na), _wilson(kb, nb)
        rows.append({"evidence": evidence, "stratum": stratum, "n_stat": na, "k_stat": ka,
                     "rate_stat": ka / na if na else np.nan, "ci_low_stat": la, "ci_high_stat": ha, "n_other": nb,
                     "k_other": kb, "rate_other": kb / nb if nb else np.nan, "ci_low_other": lb, "ci_high_other": hb,
                     "odds_ratio": odds, "p_fisher": p})

    evaluable = cat[cat.cov_status != COV_STATUS[1]]
    power = pd.cut(pd.to_numeric(evaluable.cov_expected_pairs_power), [-np.inf, 1, 3, 10, np.inf],
                   labels=["<1", "1-3", "3-10", ">=10"])
    compare("covariation, structure-free arm", "all evaluable", evaluable, evaluable.cov_detected)
    for label, sub in evaluable.groupby(power, observed=True):
        compare("covariation, structure-free arm", f"expected pairs {label}", sub, sub.cov_detected)
    measured = cat[cat.icshape_measured.astype(bool)]
    compare("icSHAPE supported", "measured", measured, measured.icshape_supported)
    if "eclip_dsrbp" in cat.columns and cat.eclip_dsrbp.notna().any():
        compare("dsRBP eCLIP peak", "all", cat, cat.eclip_dsrbp)
        compare("any eCLIP peak", "all", cat, pd.to_numeric(cat.eclip_n_experiments) > 0)
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------------------
# command line
# ----------------------------------------------------------------------------------------------------------------
def _write(frame: pd.DataFrame, path: str | None) -> None:
    """Parquet by file name, else TSV; floats are written in full (lfdr is not rounded)."""
    if path and path.endswith(".parquet"):
        frame.to_parquet(path, index=False)
    else:
        frame.to_csv(path or sys.stdout, sep="\t", index=False)


def _read_candidates(path: str) -> pd.DataFrame:
    candidates = read_table(path)
    if candidates.empty:
        raise SystemExit(f"{path}: the candidate table is empty")
    return candidates


def _candidates_command(args: argparse.Namespace) -> None:
    try:
        table = selection_candidates(_read_candidates(args.candidates), read_table(args.families))
    except ValueError as error:
        raise SystemExit(str(error)) from None
    _write(table, args.output)
    print(f"candidates: {len(table):,} in {table.rlmfam_id.nunique():,} families", file=sys.stderr)


def _icshape_command(args: argparse.Namespace) -> None:
    tracks = dict(item.split("=", 1) for item in args.track)
    detail, summary = icshape_evidence(_read_candidates(args.candidates), args.gtf, args.transcripts, tracks)
    if args.detail:
        _write(detail, args.detail)
    _write(summary, args.output)
    print(f"icSHAPE: {len(summary)} candidates measured, {int(summary.icshape_supported.sum())} supported",
          file=sys.stderr)


def _covariation_command(args: argparse.Namespace) -> None:
    candidates = normalize_candidates(_read_candidates(args.candidates))
    if args.parse_only:
        table = read_covariation(args.workdir, args.condition, candidates.core_id)
    else:
        directory = args.homologues if args.condition == "iii" else args.alignments
        if not directory or not Path(directory).is_dir():
            raise SystemExit("give --homologues (condition iii) or --alignments (conditions i and ii) as a directory")
        if args.condition != "i" and "human_core_predicted_ss" not in candidates.columns:
            raise SystemExit("the candidate table has no core structure (human_core_predicted_ss or core_dot_bracket)")
        inputs = covariation_inputs(directory, candidates.core_id, args.condition)
        expected = "<core_id>.fa" if args.condition == "iii" else "<core_id>.sto or <core_id>/final55.sto"
        if not inputs:
            raise SystemExit(f"no input in {directory} for any of the {len(candidates):,} candidates "
                             f"(expected {expected})")
        if len(inputs) < len(candidates):
            print(f"covariation ({args.condition}): no {expected} in {directory} for "
                  f"{len(candidates) - len(inputs):,} of {len(candidates):,} candidates (left out)", file=sys.stderr)
        structures = candidates.set_index("core_id").get("human_core_predicted_ss", pd.Series(dtype=object))
        tasks = [(core, args.condition, path, structures.get(core), args.workdir, args.mafft, args.rscape,
                  not args.no_resume) for core, path in inputs.items()]
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            rows = list(pool.map(_covariation_task, tasks, chunksize=8))
        table = covariation_table(rows)
    _write(table, args.output)
    print(f"covariation ({args.condition}): {table.status.value_counts().to_dict()}, "
          f"{int(table.detected.fillna(False).astype(bool).sum())} detected", file=sys.stderr)


def _eclip_command_input(path: str, folder: str | None) -> pd.DataFrame:
    table = pd.read_excel(path, sheet_name="Accessions") if path.endswith(".xlsx") else read_table(path)
    try:
        experiments = eclip_experiments(table)
    except ValueError as error:
        raise SystemExit(f"{path}: {error}") from None
    if "path" not in experiments.columns:
        folder = Path(folder or Path(path).parent)
        experiments["path"] = [str(folder / f"{a}.bed.gz") for a in experiments.file_accession]
    missing = [p for p in experiments.path if not Path(p).is_file()]
    if missing:
        raise SystemExit(f"{len(missing)} of {len(experiments)} eCLIP peak files are missing, e.g. {missing[0]}")
    return experiments


def _build_command(args: argparse.Namespace) -> None:
    candidates = _read_candidates(args.candidates)
    windows, null = read_table(args.windows), tuple(args.null)
    if args.lfdr_curve:
        if "lfdr" in windows.columns:
            print("--lfdr-curve: the lfdr of the window table is recomputed from its P values", file=sys.stderr)
            windows = windows.drop(columns="lfdr")
        windows = prepare_windows(windows.assign(lfdr=np.nan), null)  # p, z, q_bh and window_id
        windows["lfdr"], curve = window_lfdr(windows.p, null, return_curve=True)
        _write(curve, args.lfdr_curve)
    windows = prepare_windows(windows, null)
    top1_p = args.top1_p if args.top1_p is not None else top1_cutoff(windows.p)
    print(f"n_top1pct_windows: P <= {top1_p!r} ("
          + ("--top1-p" if args.top1_p is not None else f"top 1% of the {len(windows):,} windows") + ")",
          file=sys.stderr)
    covariation = (pd.concat([read_table(p) for p in args.covariation], ignore_index=True)
                   if args.covariation else None)
    eclip = eclip_overlaps(candidates, _eclip_command_input(args.eclip, args.eclip_dir)) if args.eclip else None
    splice = splice_site_distance(candidates, splice_sites(args.gtf)) if args.gtf else None
    annotation = gencode_annotation(candidates, args.gtf) if args.gtf else None
    cat = build_catalogue(
        candidates, windows, icshape=read_table(args.icshape) if args.icshape else None, covariation=covariation,
        eclip=eclip, splice=splice, rfam=rfam_overlaps(candidates, read_table(args.rfam)) if args.rfam else None,
        homologues=homologue_columns(read_table(args.homologue_msa)) if args.homologue_msa else None,
        annotation=annotation, top1_p=top1_p)
    _write(to_table(cat), args.output)
    if args.families:
        _write(families_table(cat), args.families)
    if args.tier_summary:
        _write(tier_summary(cat), args.tier_summary)
    if args.evidence_summary:
        _write(evidence_summary(cat), args.evidence_summary)
    print("tiers: " + ", ".join(f"{t} {n:,}" for t, n in cat.tier.value_counts().sort_index().items())
          + f"; statistical support {int(cat.stat_support.sum()):,}, structural support "
          f"{int(cat.structure_support.sum()):,}", file=sys.stderr)


def register(subparsers) -> None:
    """Add the `catalogue` command (subcommands candidates, icshape, covariation, build, concordance) to argparse
    subparsers. Each subcommand sets `func`, which takes the parsed arguments and writes the outputs."""
    parser = subparsers.add_parser("catalogue", help="RLM-Fam catalogue: icSHAPE, covariation, lfdr and tiers")
    sub = parser.add_subparsers(dest="catalogue_command", required=True)

    cand = sub.add_parser("candidates", help="candidate table of the catalogue steps from the selection output")
    cand.add_argument("candidates", help="candidates.parquet of 'select expand'")
    cand.add_argument("--families", required=True, help="families.tsv of 'select families'")
    cand.add_argument("-o", "--output", help="candidate table (TSV, or parquet by name; default: standard output)")
    cand.set_defaults(func=_candidates_command)

    ic = sub.add_parser("icshape", help="icSHAPE agreement of the candidate cores (GSE145805)")
    ic.add_argument("--candidates", required=True, help="candidate table (core_id, coordinates, sequence, structure)")
    ic.add_argument("--gtf", required=True, help="GENCODE v26 annotation (gencode.v26.annotation.gtf.gz)")
    ic.add_argument("--transcripts", required=True, help="GENCODE v26 transcript FASTA (gencode.v26.transcripts.fa.gz)")
    ic.add_argument("--track", action="append", required=True, metavar="CELL=FILE",
                    help="GSE145805 reactivity file of one cell line (repeat; paper: "
                         + ", ".join(f"{c}={f}" for c, f in ICSHAPE_TRACKS.items()) + ")")
    ic.add_argument("--detail", help="also write the per cell line metrics")
    ic.add_argument("-o", "--output", help="per candidate summary (default: standard output)")
    ic.set_defaults(func=_icshape_command)

    cov = sub.add_parser("covariation", help="R-scape covariation of the predicted structures (needs MAFFT, R-scape)")
    cov.add_argument("--candidates", required=True)
    cov.add_argument("--condition", choices=CONDITIONS, default="iii",
                     help="iii: structure-free homologue alignment (primary); i: CM alignment; ii: its sequences "
                          "realigned with MAFFT")
    cov.add_argument("--homologues", help="directory of <core_id>.fa: human core and high-confidence homologues (iii; "
                                          "the homologues/ directory of 'select expand')")
    cov.add_argument("--alignments", help="directory of <core_id>.sto or <core_id>/final55.sto covariance-model "
                                          "alignments (i, ii; the candidates/ directory of 'select expand')")
    cov.add_argument("--workdir", required=True, help="R-scape and MAFFT outputs, <workdir>/<condition>/<core_id>/")
    cov.add_argument("--mafft", default="mafft")
    cov.add_argument("--rscape", default="R-scape")
    cov.add_argument("--jobs", type=int, default=8)
    cov.add_argument("--no-resume", action="store_true", help="rerun candidates with finished outputs")
    cov.add_argument("--parse-only", action="store_true", help="only read existing R-scape outputs in --workdir")
    cov.add_argument("-o", "--output")
    cov.set_defaults(func=_covariation_command)

    build = sub.add_parser("build", help="assemble the catalogue: lfdr, evidence, tiers (Supplementary Table S4)")
    build.add_argument("--candidates", required=True, help="candidate table ('catalogue candidates')")
    build.add_argument("--windows", required=True,
                       help="all eligible scan windows (chrom, start, end, strand, p; q_bh and lfdr are computed if "
                            "absent)")
    build.add_argument("--null", nargs=3, type=float, default=list(EMPIRICAL_NULL), metavar=("MU0", "SIGMA0", "PI0"),
                       help="empirical null of z for the lfdr (default: the paper's scan)")
    build.add_argument("--icshape", help="output of 'catalogue icshape'")
    build.add_argument("--covariation", nargs="+", help="outputs of 'catalogue covariation' (conditions i, ii, iii)")
    build.add_argument("--eclip", help="eCLIP experiments (rbp, cell, and path or file_accession), or the Accessions "
                                       "sheet of Supplementary Table S8 (.xlsx, or the sheet as TSV)")
    build.add_argument("--eclip-dir", help="directory of <file_accession>.bed.gz (default: next to --eclip)")
    build.add_argument("--gtf", help="GENCODE v50 annotation: gene annotation columns and distance to splice sites")
    build.add_argument("--rfam", help="Rfam hits (chrom, start, end, strand, rfam_acc, rfam_id)")
    build.add_argument("--homologue-msa", help="RLMHelix-MSA of homologue alignments against 19 SISSIz shuffles")
    build.add_argument("-o", "--output", help="catalogue (default: standard output)")
    build.add_argument("--families", help="also write one row per family")
    build.add_argument("--tier-summary", help="also write the counts of evidence per tier")
    build.add_argument("--evidence-summary", help="also write the numbers of candidates with each kind of evidence")
    build.add_argument("--lfdr-curve", help="also write the lfdr curve (z, f, pi0_f0, lfdr); the lfdr is then "
                                           "computed here from the P values, also when the table has an lfdr column")
    build.add_argument("--top1-p", type=float, metavar="P",
                       help="P cutoff of n_top1pct_windows (default: the P value at the top-1%% boundary of the window "
                            f"table, {TOP1_P} for the paper's 7,869,821 windows)")
    build.set_defaults(func=_build_command)

    conc = sub.add_parser("concordance", help="evidence at candidates with and without statistical support (Fig. S8c)")
    conc.add_argument("catalogue", help="catalogue written by 'catalogue build' (or the candidate table of the Zenodo deposit, once public)")
    conc.add_argument("-o", "--output")
    conc.set_defaults(func=lambda args: _write(evidence_concordance(read_table(args.catalogue)), args.output))


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="rlmhelix")
    register(parser.add_subparsers(dest="command", required=True))
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
