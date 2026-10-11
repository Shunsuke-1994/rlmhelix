"""Catalogue step on small synthetic tables: identifiers, window statistics and lfdr, icSHAPE, R-scape parsing,
evidence tiers and the Supplementary Table S4 layout."""
import argparse

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from rlmhelix import catalogue as rc

STRUCTURE = "((((((((((" + "." * 20 + "))))))))))"  # 20 paired, 20 unpaired


@pytest.fixture
def candidates() -> pd.DataFrame:
    return pd.DataFrame({
        "query_id": ["c1", "c2", "c3", "c4"], "provisional_cluster_id": ["C01", "C01", "C02", "C14037"],
        "chrom": ["chr1", "chr1", "chr2", "chr3"], "strand": ["+", "-", "+", "+"],
        "core_start0": [1000, 1000, 5000, 100], "core_end0": [1040, 1040, 5040, 140],
        "human_core_rna_sequence": ["GGGGGGGGGG" + "A" * 20 + "CCCCCCCCCC"] * 4,
        "human_core_predicted_ss": [STRUCTURE] * 4, "representative_gene_name": ["G1", "G1", "G2", "G3"]})


@pytest.fixture
def windows() -> pd.DataFrame:
    return pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1", "chr1", "chr2"], "start": [900, 1000, 1040, 950, 4900],
        "end": [1100, 1200, 1240, 1150, 5100], "strand": ["+", "+", "+", "-", "+"],
        "p": [1e-4, 0.01, 1e-6, 0.03, 0.001], "q_bh": [0.01, 0.3, 1e-4, 0.5, 0.2],
        "lfdr": [0.05, 0.8, 0.001, 0.95, 0.15]})


@pytest.fixture
def evidence() -> tuple[pd.DataFrame, pd.DataFrame]:
    icshape = pd.DataFrame({"core_id": ["c2", "c4"], "icshape_n_contexts": [2, 3], "icshape_mean_auroc": [0.8, 0.55],
                            "icshape_frac_contexts_p05": [1.0, 0.0]})
    covariation = rc.covariation_table([
        {"core_id": "c1", "condition": "iii", "status": "success", "nseq": 30, "avg_identity": 80.0,
         "n_pairs_tested": 10, "n_significant_pairs": 0, "expected_pairs": 5.0, "n_helices_with_2_cov_pairs": 0},
        {"core_id": "c3", "condition": "iii", "status": "success", "nseq": 40, "avg_identity": 70.0,
         "n_pairs_tested": 10, "n_significant_pairs": 3, "expected_pairs": 4.0, "n_helices_with_2_cov_pairs": 1},
        {"core_id": "c4", "condition": "iii", "status": "success", "nseq": 4, "avg_identity": 95.0,
         "n_pairs_tested": 10, "n_significant_pairs": 0, "expected_pairs": 0.5, "n_helices_with_2_cov_pairs": 0},
        {"core_id": "c1", "condition": "i", "status": "success", "nseq": 30, "avg_identity": 80.0,
         "n_pairs_tested": 10, "n_significant_pairs": 2, "expected_pairs": 5.0, "n_helices_with_2_cov_pairs": 1},
        {"core_id": "c2", "condition": "iii", "status": "failed"}])
    return icshape, covariation


def test_family_id():
    assert rc.family_id("C01") == "RLMFam00001"
    assert rc.family_id("C14037") == "RLMFam14037"
    assert rc.family_id("RLMFam00092") == "RLMFam00092"
    assert rc.browser_locus("chr7", 149071457, 149071600) == "chr7:149071458-149071600"


def test_window_statistics(candidates, windows):
    out = rc.window_statistics(candidates, rc.prepare_windows(windows)).set_index("core_id")
    c1 = out.loc["c1"]  # w1 and w2 overlap; w3 starts at the core end and does not
    assert c1.min_window_p == 1e-4 and c1.lfdr == 0.05 and c1.min_window_q_bh == 0.01 and bool(c1.bh_window)
    assert c1.best_window_id == "chr1:900-1100:+" and c1.n_top1pct_windows == 2
    assert c1.max_window_z == pytest.approx(stats.norm.isf(1e-4))
    assert out.loc["c2", "min_window_p"] == 0.03 and out.loc["c2", "n_top1pct_windows"] == 0  # minus strand only
    assert not out.loc["c3", "bh_window"] and out.loc["c3", "lfdr"] == 0.15
    assert np.isnan(out.loc["c4", "min_window_p"]) and not out.loc["c4", "bh_window"]


def test_bh_qvalues():
    np.testing.assert_allclose(rc.bh_qvalues([0.01, 0.04, 0.03, 0.5]), [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 0.5])


def test_window_lfdr_shape():
    rng = np.random.default_rng(1)
    null = (0.67, 0.55, 0.95)
    z = np.r_[rng.normal(0.67, 0.55, 190_000), rng.normal(3.0, 1.0, 10_000)]
    lfdr, curve = rc.window_lfdr(stats.norm.sf(z), null, return_curve=True)
    assert ((lfdr >= 0) & (lfdr <= 1)).all() and (lfdr[z < null[0] - 0.01] == 1).all()
    right = np.argsort(z)[np.sort(z) >= null[0]]
    assert (np.diff(lfdr[right]) <= 1e-12).all()  # non-increasing to the right of mu0
    assert lfdr[z > 5].max() < 0.05 and np.median(lfdr[np.abs(z - 0.67) < 0.2]) > 0.95
    assert list(curve.columns) == ["z", "f", "pi0_f0", "lfdr"]


def test_auroc_and_icshape_metrics():
    rng = np.random.default_rng(2)
    labels, values = rng.random(80) < 0.4, np.round(rng.random(80), 1)  # many ties
    u = stats.mannwhitneyu(values[labels], values[~labels]).statistic
    assert rc.auroc(labels, values) == pytest.approx(u / (labels.sum() * (~labels).sum()), abs=1e-12)
    reactivity = np.where(np.array(list(STRUCTURE)) == ".", 0.9, 0.1)
    reactivity[[0, 15]] = np.nan
    out = rc.icshape_metrics(STRUCTURE, reactivity)
    assert out["reactivity_auroc"] == 1.0 and out["n_measured_paired"] == 19 and out["n_measured_unpaired"] == 19
    assert out["circular_shift_empirical_p"] == 1 / (1 + out["n_valid_circular_shifts"])
    assert np.isnan(rc.icshape_metrics(STRUCTURE, np.where(reactivity > 0.5, np.nan, 0.1))["reactivity_auroc"])


def test_icshape_support():
    good = np.where(np.array(list(STRUCTURE)) == ".", 0.9, 0.1)
    sparse = good.copy()
    sparse[:25] = np.nan  # 37.5% measured: below the 50% coverage
    reactivities = pd.DataFrame({"core_id": ["a", "a", "b"], "cell_line": ["HeLa", "K562", "HeLa"],
                                 "transcript_id": ["t1", "t1", "t2"], "coverage": [1.0, 0.375, 1.0],
                                 "reactivity": [good, sparse, good[::-1].copy()]})
    detail, summary = rc.icshape_support(reactivities, {"a": STRUCTURE, "b": STRUCTURE})
    assert len(detail) == 3
    summary = summary.set_index("core_id")
    assert summary.loc["a", "icshape_n_contexts"] == 1 and summary.loc["a", "icshape_supported"]


def test_map_cores_and_extract(tmp_path, candidates):
    records = {"t1": {"chrom": "chr1", "strand": "+", "exons": [(900, 1100)], "length": 200, "gene_name": "G1"},
               "t2": {"chrom": "chr1", "strand": "-", "exons": [(1200, 1300), (990, 1050)], "length": 160,
                      "gene_name": "G1"}}
    core = candidates.human_core_rna_sequence[0].replace("U", "T")
    minus = "N" * 100 + "N" * 10 + core + "N" * 10  # second exon of t2 on the minus strand: 1050 -> 990
    sequences = {"t1": "N" * 100 + core + "N" * 60, "t2": minus}
    mapping = rc.map_cores_to_transcripts(candidates, records, sequences)
    assert set(zip(mapping.core_id, mapping.transcript_id)) == {("c1", "t1"), ("c2", "t2")}
    assert mapping.set_index("core_id").loc["c2", "transcript_start"] == 110
    track = tmp_path / "cell.out.txt"
    values = ["NULL"] * 100 + ["0.5"] * 40 + ["NULL"] * 60
    track.write_text("t1\t200\t*\t" + "\t".join(values) + "\n")
    reactivities = rc.extract_reactivities(mapping[mapping.core_id == "c1"], {"HeLa": track})
    assert reactivities.coverage.iloc[0] == 1.0 and reactivities.reactivity.iloc[0].shape == (40,)


def test_structure_projection_and_rscape_parsing():
    assert rc.structure_pairs("((..[[..))..]]") == [(0, 9), (1, 8), (4, 13), (5, 12)]
    records = rc.unique_sequences([("x", "ACGU"), ("h", "AC-GUACGU"), ("y", "ACGT"), ("z", "acgu-acgg")], "h")
    assert [n for n, _ in records] == ["h", "x", "z"]
    text = rc.projected_stockholm([("h", "AC-GU-ACGU"), ("x", "ACGGUAAC-U")], "h", rc.structure_pairs("((....))"))
    lines = {line.split()[0] + (line.split()[1] if line.startswith("#=GC") else ""): line.split()[-1]
             for line in text.splitlines() if line and not line.startswith(("# ", "//"))}
    assert lines["#=GCSS_cons"] == "<<......>>" and lines["#=GCRF"] == "xx.xx.xxxx" and lines["x"] == "ACGGUAAC-U"
    log = "\n".join([
        "# MSA rscape nseq 120 (125) alen 175 (1393) avgid 80.57 (79.36) nbpairs 47 (49)",
        "# GTp    0.05         [-9.95,926.98]     [1 | 3 47 4 | 0.06 0.75 0.12] ",
        "# BPAIRS expected to covary 6.6 +/- 0.6", "# RM_HELIX NESTED 9-35 174-308, nbp = 9 nbp_cov = 2"])
    out = rc.parse_rscape(log)
    assert (out["nseq"], out["avg_identity"], out["n_pairs_tested"], out["n_significant_pairs"],
            out["expected_pairs"], out["n_helices_with_2_cov_pairs"]) == (120, 80.57, 47, 3, 6.6, 1)
    table = rc.covariation_table([{"core_id": "a", "condition": "iii", **out}])
    assert table.detected.iloc[0] and table.two_pairs.iloc[0] and table.stringent.iloc[0]
    assert rc.parse_rscape("covariation scores are almost constant", "uninformative")["n_significant_pairs"] == 0


def test_tiers_and_table(candidates, windows, evidence):
    icshape, covariation = evidence
    cat = rc.build_catalogue(candidates, rc.prepare_windows(windows), icshape=icshape, covariation=covariation)
    c = cat.set_index("core_id")
    assert c.tier.to_dict() == {"c1": "B", "c2": "B", "c3": "A", "c4": "C"}
    assert c.tier_evidence.to_dict() == {"c1": "statistics", "c2": "icSHAPE", "c3": "statistics+covariation",
                                         "c4": "none"}
    assert c.cov_status.to_dict() == {"c1": "not detected (expected >= 1 pair)", "c2": "not evaluable",
                                      "c3": "detected", "c4": "not detected, low power (expected < 1 pair)"}
    assert c.loc["c1", "cov_detected_cm_aligned"] and c.loc["c4", "icshape_measured"]
    assert not c.loc["c4", "icshape_supported"] and c.structure_evidence_measured.to_dict() == {
        "c1": True, "c2": True, "c3": True, "c4": True}
    assert cat.core_id.tolist() == ["c3", "c1", "c2", "c4"]  # tier, then lfdr
    table = rc.to_table(cat)
    assert list(table.columns) == [public for public, _ in rc.S4_COLUMNS]
    assert table.rlmfam_id.tolist() == ["RLMFam00002", "RLMFam00001", "RLMFam00001", "RLMFam14037"]
    assert table.family_size.tolist() == [1, 2, 2, 1] and table.eclip_n_experiments.isna().all()
    families = rc.families_table(cat)
    assert list(families.columns) == rc.FAMILY_COLUMNS
    assert families.rlmfam_id.tolist() == ["RLMFam00002", "RLMFam00001", "RLMFam14037"]
    first = families.set_index("rlmfam_id").loc["RLMFam00001"]
    assert first.candidate_loci == "chr1:1001-1040(+);chr1:1001-1040(-)" and first.best_tier == "B"
    assert (first.n_statistical_support, first.n_structural_support, first.min_lfdr) == (1, 1, 0.05)
    summary = rc.tier_summary(cat).set_index("tier")
    assert summary.loc["all", "n_candidates"] == 4 and summary.loc["B", "statistical_support"] == 1
    assert rc.evidence_summary(cat).n_with_evidence.tolist()[:4] == [2, 1, 2, 1]
    concordance = rc.evidence_concordance(table)  # also from the published column names
    assert concordance.stratum.iloc[0] == "all evaluable" and concordance.n_stat.iloc[0] == 2


def test_cli_build(tmp_path, candidates, windows, evidence):
    icshape, covariation = evidence
    paths = {}
    for name, frame in (("candidates", candidates), ("windows", windows), ("icshape", icshape),
                        ("covariation", covariation)):
        paths[name] = tmp_path / f"{name}.tsv"
        frame.to_csv(paths[name], sep="\t", index=False)
    out, families = tmp_path / "catalogue.tsv", tmp_path / "families.tsv"
    rc.main(["catalogue", "build", "--candidates", str(paths["candidates"]), "--windows", str(paths["windows"]),
             "--icshape", str(paths["icshape"]), "--covariation", str(paths["covariation"]), "-o", str(out),
             "--families", str(families), "--tier-summary", str(tmp_path / "tiers.tsv")])
    table = pd.read_csv(out, sep="\t")
    assert list(table.columns) == [public for public, _ in rc.S4_COLUMNS]
    assert table.evidence_tier.tolist() == ["A", "B", "B", "C"]
    assert pd.read_csv(families, sep="\t").shape == (3, len(rc.FAMILY_COLUMNS))
    rc.main(["catalogue", "concordance", str(out), "-o", str(tmp_path / "concordance.tsv")])
    assert (tmp_path / "concordance.tsv").exists()
    parser = argparse.ArgumentParser()
    rc.register(parser.add_subparsers(dest="command"))
    assert callable(parser.parse_args(["catalogue", "concordance", str(out)]).func)


# ----------------------------------------------------------------------------------------------------------------
# from the selection output to the catalogue
# ----------------------------------------------------------------------------------------------------------------
def _hairpin_weights(length: int, pairs: list[tuple[int, int]]) -> np.ndarray:
    weights = np.zeros((length, length), dtype=np.float32)
    for i, j in pairs:
        weights[i, j] = weights[j, i] = 0.5
    return weights


@pytest.fixture
def selection(tmp_path) -> tuple:
    """candidates.parquet of `select expand` and families.tsv of `select families`, made with rlmhelix.select."""
    from rlmhelix import select as sel

    pairs = [(5, 18), (6, 17), (7, 16), (8, 15)]
    records = [
        {"sequence": "AAAAAGGGGAAAAAACCCCAA", "strand": "+", "chrom": "chr1", "peak_start": 1000, "peak_end": 1021,
         "locus": "L1", "peak_window": "chr1:1000-1200:+", "peak_score": 9.0, "features": "intron"},
        {"sequence": "AAAAACCCCAAAAAAGGGGAA", "strand": "-", "chrom": "chr2", "peak_start": 5000, "peak_end": 5021,
         "locus": "L2", "peak_window": "chr2:5000-5200:-", "peak_score": 8.0, "features": "lncrna"}]
    cores = sel.rank_cores(pd.DataFrame([row for r in records for row in sel.window_cores(r, _hairpin_weights(21, pairs))]))
    loci = pd.DataFrame({"peak_window": [r["peak_window"] for r in records], "peak_p": [1e-5, 1e-3], "GC": [0.4, 0.4],
                         "n_windows": [2, 1]})
    reps, _ = sel.collapse_cores(cores, loci)
    counts = pd.DataFrame({"primate": [6, 6], "nonprimate_mammal": [3, 0], "nonmammal_tetrapod": [1, 0],
                           "other_vertebrate": [1, 0]}, index=reps.query_id)
    selected = sel.assign_routes(reps, counts)
    selected = selected.loc[selected.selected].assign(threshold_bits=[20.0, 25.0])
    loci_counts = pd.DataFrame({"n_species_e01": [10, 8], "n_loci_e01": [12, 9], "n_species_hc": [9, 7],
                                "n_loci_hc": [10, 8]}, index=reps.query_id)
    final = pd.DataFrame({"query_id": reps.query_id, "n_accepted_expansion_loci": [30, 4],
                          "n_accepted_expansion_species": [25, 3], "n_final_alignment_sequences": [40, 11]})
    table = sel.candidate_table(selected, loci_counts, final)
    table["msa_path"] = [str(tmp_path / "expand/candidates" / q / "final55.sto") for q in table.query_id]
    families = sel.provisional_families(table.query_id, pd.DataFrame(columns=["unit_a", "unit_b", "ratio_mean"]))
    table.to_parquet(tmp_path / "candidates.parquet", index=False)
    families.to_csv(tmp_path / "families.tsv", sep="\t", index=False)
    return tmp_path / "candidates.parquet", tmp_path / "families.tsv", table


def test_candidates_from_selection(tmp_path, selection):
    candidates_path, families_path, table = selection
    out = tmp_path / "candidates.tsv"
    rc.main(["catalogue", "candidates", str(candidates_path), "--families", str(families_path), "-o", str(out)])
    cand = rc.read_table(out)
    assert list(cand.columns) == [c for c in rc.CANDIDATE_COLUMNS if c != "cm_path"]
    first = cand.set_index("core_id").loc[table.query_id.iloc[0]]
    assert (first.core_start0, first.core_end0, first.browser_locus) == (1005, 1019, "chr1:1006-1019")
    assert first.human_core_rna_sequence == "GGGGAAAAAACCCC" and first.human_core_predicted_ss == "((((......))))"
    assert (first.lineage_route, first.msa_species, first.msa_sequences, first.scan_pvalue) == (
        "vertebrate_wide", 34, 40, 1e-5)
    assert first.source_features == "intron" and first.rlmfam_id.startswith("RLMFam") and first.cluster_size == 1
    # the select table itself is also accepted by the catalogue steps
    assert rc.normalize_candidates(table).human_core_predicted_ss.tolist() == ["((((......))))"] * 2
    windows = pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [1000, 5000], "end": [1200, 5200],
                            "strand": ["+", "-"], "p": [1e-5, 1e-3], "q_bh": [0.01, 0.5], "lfdr": [0.1, 0.9]})
    windows.to_csv(tmp_path / "windows.tsv", sep="\t", index=False)
    rc.main(["catalogue", "build", "--candidates", str(out), "--windows", str(tmp_path / "windows.tsv"),
             "-o", str(tmp_path / "catalogue.tsv"), "--families", str(tmp_path / "fam.tsv")])
    s4 = pd.read_csv(tmp_path / "catalogue.tsv", sep="\t")
    for column in ("rlmfam_id", "core_sequence", "core_structure", "selection_route", "n_species_alignment",
                   "n_sequences_alignment", "scan_region_types", "scan_p"):
        assert s4[column].notna().all(), column
    assert pd.read_csv(tmp_path / "fam.tsv", sep="\t").shape[0] == 2
    pd.read_csv(families_path, sep="\t").head(1).to_csv(tmp_path / "one.tsv", sep="\t", index=False)
    with pytest.raises(SystemExit, match="not in the family table"):
        rc.main(["catalogue", "candidates", str(candidates_path), "--families", str(tmp_path / "one.tsv")])


FAKE_RSCAPE = """#!{python}
import os, sys
args = sys.argv[1:]
outdir = args[args.index("--outdir") + 1]
if not os.path.isdir(outdir):  # as R-scape 2.6.16, relative to its working directory
    print(f"Outdir {{outdir}} does not seem to exist")
    sys.exit(1)
open(os.path.join(outdir, "rscape.cov"), "w").close()
print("# MSA rscape nseq 4 (4) alen 14 (14) avgid 80.00 (80.00) nbpairs 4 (4)")
print("# GTp    0.05         [-9.95,926.98]     [0 | 1 4 4 | 0.06 0.75 0.12] ")
print("# BPAIRS expected to covary 1.5 +/- 0.6")
"""


def test_covariation_inputs_relative_workdir_and_no_inputs(tmp_path, monkeypatch, selection):
    import sys

    candidates_path, _, table = selection
    first, second = table.query_id
    (tmp_path / "expand/candidates" / first).mkdir(parents=True)
    alignment = "# STOCKHOLM 1.0\n{} GGGGAAAAAACCCC\nx GGGAAAAAAACCCC\n#=GC SS_cons <<<<......>>>>\n//\n"
    (tmp_path / "expand/candidates" / first / "final55.sto").write_text(alignment.format(first))
    (tmp_path / "expand/candidates" / f"{second}.sto").write_text(alignment.format(second))
    assert set(rc.covariation_inputs(tmp_path / "expand/candidates", table.query_id, "i")) == {first, second}
    assert rc.covariation_inputs(tmp_path / "expand/candidates", table.query_id, "iii") == {}
    rscape = tmp_path / "fake_rscape"
    rscape.write_text(FAKE_RSCAPE.format(python=sys.executable))
    rscape.chmod(0o755)
    monkeypatch.chdir(tmp_path)  # relative --workdir, as in the documented commands
    rc.main(["catalogue", "covariation", "--candidates", str(candidates_path), "--condition", "i", "--alignments",
             "expand/candidates", "--workdir", "rscape", "--rscape", str(rscape), "--jobs", "1", "-o", "cov_i.tsv"])
    cov = pd.read_csv("cov_i.tsv", sep="\t")
    assert cov.status.tolist() == ["success", "success"] and cov.detected.all() and cov.n_pairs_tested.eq(4).all()
    (tmp_path / "empty").mkdir()
    with pytest.raises(SystemExit, match="no input"):
        rc.main(["catalogue", "covariation", "--candidates", str(candidates_path), "--homologues", "empty",
                 "--workdir", "rscape"])
    empty = rc.covariation_table([])
    assert empty.empty and "detected" in empty.columns


def test_rscape_relative_workdir(tmp_path, monkeypatch):
    """R-scape itself with a relative work directory (skipped without R-scape on PATH)."""
    import shutil
    from pathlib import Path

    if shutil.which("R-scape") is None:
        pytest.skip("R-scape not on PATH")
    alignment = Path(__file__).parent / "data" / "RF00031.sto"
    monkeypatch.chdir(tmp_path)
    out = rc.covariation_one("hg38/689-759", "i", alignment, None, "rel", resume=False)
    assert out["status"] == "success" and out["nseq"] == 9


def test_lfdr_curve_top1_cutoff_and_precision(tmp_path):
    rng = np.random.default_rng(3)
    z = np.r_[rng.normal(0.67, 0.55, 19_000), rng.normal(3.0, 1.0, 1_000)]
    p = stats.norm.sf(z)
    n = len(p)
    windows = pd.DataFrame({"chrom": "chr1", "start": np.arange(n) * 100, "end": np.arange(n) * 100 + 200,
                            "strand": "+", "p_value": p, "lfdr": 0.5})  # an lfdr column, as from `scan fdr`
    windows.to_parquet(tmp_path / "windows.parquet", index=False)
    cutoff = rc.top1_cutoff(p)
    assert cutoff == np.sort(p)[int(np.ceil(0.01 * n)) - 1] and (p <= cutoff).sum() == 200
    candidates = pd.DataFrame({"core_id": ["a", "b"], "chrom": "chr1", "strand": "+", "core_start0": [1000, 5050],
                               "core_end0": [1150, 5060], "rlmfam_id": ["RLMFam00001", "RLMFam00002"]})
    candidates.to_csv(tmp_path / "cand.tsv", sep="\t", index=False)
    args = ["catalogue", "build", "--candidates", str(tmp_path / "cand.tsv"), "--windows",
            str(tmp_path / "windows.parquet"), "-o", str(tmp_path / "cat.tsv")]
    rc.main(args + ["--lfdr-curve", str(tmp_path / "curve.tsv")])
    curve = pd.read_csv(tmp_path / "curve.tsv", sep="\t")
    assert list(curve.columns) == ["z", "f", "pi0_f0", "lfdr"]
    table = pd.read_csv(tmp_path / "cat.tsv", sep="\t").set_index("core_id")
    hits = [np.arange(9, 12), np.arange(49, 51)]  # windows overlapping a and b
    assert table.n_top1pct_windows.tolist() == [int((p[h] <= cutoff).sum()) for h in hits]
    lfdr = rc.window_lfdr(p)
    lfdr[z < rc.EMPIRICAL_NULL[0]] = 1.0
    assert table.loc["a", "lfdr"] == pytest.approx(lfdr[hits[0]].min(), rel=1e-15)  # recomputed, not rounded
    rc.main(args + ["--top1-p", "1"])
    table = pd.read_csv(tmp_path / "cat.tsv", sep="\t").set_index("core_id")
    assert table.n_top1pct_windows.tolist() == [3, 2] and (table.lfdr == 0.5).all()


def test_lfdr_written_in_full(tmp_path, candidates, windows):
    """lfdr is not rounded in the written tables, so that re-thresholding it agrees with statistical_support."""
    windows.assign(lfdr=[0.05, 0.8, 0.001, 0.95, 0.2000437]).to_csv(tmp_path / "w.tsv", sep="\t", index=False)
    candidates.iloc[[2]].to_csv(tmp_path / "c.tsv", sep="\t", index=False)  # c3: one window, BH q 0.2
    rc.main(["catalogue", "build", "--candidates", str(tmp_path / "c.tsv"), "--windows", str(tmp_path / "w.tsv"),
             "-o", str(tmp_path / "cat.tsv"), "--families", str(tmp_path / "fam.tsv")])
    table = pd.read_csv(tmp_path / "cat.tsv", sep="\t")
    assert table.lfdr.iloc[0] == 0.2000437 and not table.statistical_support.iloc[0]
    assert pd.read_csv(tmp_path / "fam.tsv", sep="\t").min_lfdr.iloc[0] == 0.2000437


def test_eclip_from_s8_accessions(tmp_path, candidates, windows):
    import gzip

    s8 = pd.DataFrame({"dataset": ["eCLIP peaks", "eCLIP signal tracks", "eCLIP peaks"],
                       "target": ["DGCR8", "DGCR8", "HNRNPC"], "cell_line": ["K562", "K562", "HepG2"],
                       "file_accession": ["ENCFF000AAA", "ENCFF000BBB", "ENCFF000CCC"]})
    experiments = rc.eclip_experiments(s8)
    assert experiments[["rbp", "cell", "file_accession"]].values.tolist() == [["DGCR8", "K562", "ENCFF000AAA"],
                                                                            ["HNRNPC", "HepG2", "ENCFF000CCC"]]
    folder = tmp_path / "encode"
    folder.mkdir()
    for accession, line in (("ENCFF000AAA", "chr1\t1010\t1020\tp\t0\t+\n"), ("ENCFF000CCC", "chr2\t1\t5\tp\t0\t+\n")):
        with gzip.open(folder / f"{accession}.bed.gz", "wt") as handle:
            handle.write(line)
    s8.to_csv(tmp_path / "S8_accessions.tsv", sep="\t", index=False)
    candidates.to_csv(tmp_path / "c.tsv", sep="\t", index=False)
    windows.to_csv(tmp_path / "w.tsv", sep="\t", index=False)
    rc.main(["catalogue", "build", "--candidates", str(tmp_path / "c.tsv"), "--windows", str(tmp_path / "w.tsv"),
             "--eclip", str(tmp_path / "S8_accessions.tsv"), "--eclip-dir", str(folder), "-o", str(tmp_path / "x.tsv")])
    table = pd.read_csv(tmp_path / "x.tsv", sep="\t").set_index("core_id")
    assert table.loc["c1", "eclip_experiments"] == "DGCR8 (K562)" and table.loc["c1", "eclip_dsrbp"]
    assert table.loc["c3", "eclip_n_experiments"] == 0


GTF = [  # chrom, feature, start (1-based), end, strand, attributes
    ("chr1", "gene", 1001, 3000, "+", 'gene_id "G1"; gene_type "protein_coding"; gene_name "GENE1";'),
    ("chr1", "transcript", 1001, 3000, "+", 'gene_id "G1"; transcript_id "T1"; gene_type "protein_coding"; '
                                             'gene_name "GENE1"; transcript_name "GENE1-201"; tag "basic";'),
    ("chr1", "exon", 1001, 1100, "+", 'gene_id "G1"; transcript_id "T1";'),
    ("chr1", "exon", 2001, 3000, "+", 'gene_id "G1"; transcript_id "T1";'),
    ("chr1", "CDS", 2051, 2900, "+", 'gene_id "G1"; transcript_id "T1";'),
    ("chr1", "UTR", 1001, 1100, "+", 'gene_id "G1"; transcript_id "T1";'),
    ("chr1", "UTR", 2001, 2050, "+", 'gene_id "G1"; transcript_id "T1";'),
    ("chr1", "UTR", 2901, 3000, "+", 'gene_id "G1"; transcript_id "T1";'),
    ("chr1", "transcript", 2001, 3000, "+", 'gene_id "G1"; transcript_id "T2"; gene_type "protein_coding"; '
                                             'gene_name "GENE1"; transcript_name "GENE1-202"; tag "MANE_Select";'),
    ("chr1", "exon", 2001, 3000, "+", 'gene_id "G1"; transcript_id "T2";'),
    ("chr1", "gene", 5001, 6000, "-", 'gene_id "G2"; gene_type "processed_pseudogene"; gene_name "PSG2";'),
    ("chr1", "gene", 7001, 8000, "+", 'gene_id "G3"; gene_type "unprocessed_pseudogene"; gene_name "PSG3";'),
    ("chr1", "gene", 10001, 11000, "-", 'gene_id "G4"; gene_type "protein_coding"; gene_name "GENE4";'),
    ("chr1", "transcript", 10001, 11000, "-", 'gene_id "G4"; transcript_id "T4"; gene_type "protein_coding"; '
                                               'gene_name "GENE4"; transcript_name "GENE4-201";'),
    ("chr1", "exon", 10001, 11000, "-", 'gene_id "G4"; transcript_id "T4";'),
    ("chr1", "CDS", 10201, 10800, "-", 'gene_id "G4"; transcript_id "T4";'),
    ("chr1", "UTR", 10001, 10200, "-", 'gene_id "G4"; transcript_id "T4";'),
    ("chr1", "UTR", 10801, 11000, "-", 'gene_id "G4"; transcript_id "T4";'),
]


def test_gencode_annotation(tmp_path):
    import gzip

    with gzip.open(tmp_path / "gencode.gtf.gz", "wt") as handle:
        handle.write("##description: test\n")
        for chrom, feature, start, end, strand, attributes in GTF:
            handle.write(f"{chrom}\tHAVANA\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
    cores = {"utr5": ("+", 1010), "intron": ("+", 1500), "cds": ("+", 2100), "utr3": ("+", 2950),
             "opposite": ("+", 5100), "gene_only": ("+", 7500), "minus_utr3": ("-", 10050), "none": ("+", 20000)}
    candidates = pd.DataFrame({"core_id": list(cores), "chrom": "chr1", "strand": [s for s, _ in cores.values()],
                               "core_start0": [b for _, b in cores.values()],
                               "core_end0": [b + 40 for _, b in cores.values()]})
    out = rc.gencode_annotation(candidates, tmp_path / "gencode.gtf.gz").set_index("core_id")
    assert list(out.columns) == rc.GENCODE_COLUMNS
    assert out.annotation_status.to_dict() == {
        "utr5": "same_strand_exonic", "intron": "same_strand_intronic", "cds": "same_strand_exonic",
        "utr3": "same_strand_exonic", "opposite": "opposite_strand_only", "gene_only": "same_strand_gene_only",
        "minus_utr3": "same_strand_exonic", "none": "no_gene_overlap"}
    assert out.representative_context.fillna("").to_dict() == {
        "utr5": "5UTR", "intron": "intron", "cds": "exon_non_CDS_UTR", "utr3": "exon_non_CDS_UTR",
        "opposite": "", "gene_only": "", "minus_utr3": "3UTR", "none": ""}
    # same exonic overlap: the MANE_Select transcript T2 is the representative
    assert out.loc["cds", "representative_transcript_name"] == "GENE1-202"
    assert out.loc["cds", "same_strand_contexts"] == "CDS;exon_non_CDS_UTR"
    assert out.loc["utr3", "same_strand_contexts"] == "3UTR;exon_non_CDS_UTR"
    assert out.loc["utr5", "same_strand_gene_names"] == "GENE1" and pd.isna(out.loc["opposite", "same_strand_gene_names"])
    assert out.loc["gene_only", "pseudogene_gene_names"] == "PSG3" and pd.isna(out.loc["opposite", "pseudogene_gene_names"])
    # `catalogue build --gtf` fills the annotation and splice-site columns and the Families genes/region_contexts
    candidates.assign(rlmfam_id=["RLMFam00001"] * 4 + ["RLMFam00002"] * 4).to_csv(tmp_path / "c.tsv", sep="\t",
                                                                                index=False)
    pd.DataFrame({"chrom": "chr1", "start": [1000, 2000], "end": [1200, 2200], "strand": "+", "p": [1e-4, 0.01],
                  "q_bh": [0.01, 0.5], "lfdr": [0.1, 0.8]}).to_csv(tmp_path / "w.tsv", sep="\t", index=False)
    rc.main(["catalogue", "build", "--candidates", str(tmp_path / "c.tsv"), "--windows", str(tmp_path / "w.tsv"),
             "--gtf", str(tmp_path / "gencode.gtf.gz"), "-o", str(tmp_path / "cat.tsv"), "--families",
             str(tmp_path / "fam.tsv")])
    table = pd.read_csv(tmp_path / "cat.tsv", sep="\t").set_index("core_id")
    assert table.loc["cds", "representative_gene"] == "GENE1" and table.loc["gene_only", "pseudogene_context"]
    assert table.loc["utr5", "splice_site_class"] == "within 100 nt" and table.loc["intron", "splice_site_distance_nt"] == 400
    families = pd.read_csv(tmp_path / "fam.tsv", sep="\t").set_index("rlmfam_id")
    assert families.loc["RLMFam00001", "genes"] == "GENE1"
    assert set(families.loc["RLMFam00001", "region_contexts"].split(";")) == {"5UTR", "intron", "exon_non_CDS_UTR"}


def test_empty_inputs(tmp_path, candidates, windows):
    with pytest.raises(ValueError, match="empty"):
        rc.prepare_windows(windows.head(0))
    far = candidates.assign(chrom="chrY")
    out = rc.window_statistics(far, rc.prepare_windows(windows))
    assert out.min_window_p.isna().all() and not out.bh_window.any()
    candidates.head(0).to_csv(tmp_path / "c.tsv", sep="\t", index=False)
    windows.to_csv(tmp_path / "w.tsv", sep="\t", index=False)
    with pytest.raises(SystemExit, match="empty"):
        rc.main(["catalogue", "build", "--candidates", str(tmp_path / "c.tsv"), "--windows", str(tmp_path / "w.tsv")])
