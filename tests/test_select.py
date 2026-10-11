"""Candidate selection on small synthetic inputs (no GPU, network or external tools)."""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rlmhelix import select as sel


def windows_table(rows):
    return pd.DataFrame(rows, columns=["window_id", "chrom", "start", "end", "strand", "sequence", "score", "p_value"])


def test_top_windows_keeps_table_order_for_ties():
    frame = windows_table([(f"w{k}", "chr1", 100 * k, 100 * k + 200, "+", "A", 1.0, p)
                           for k, p in enumerate([0.5, 0.01, 0.01, 0.2, 0.01])])
    top = sel.top_windows(frame, fraction=0.5)                # ceil(2.5) = 3 windows
    assert list(top.window_id) == ["w1", "w2", "w4"]


def test_merge_loci_overlap_strand_and_peak():
    frame = windows_table([
        ("a", "chr1", 0, 200, "+", "A", 10.0, 0.010),
        ("b", "chr1", 100, 300, "+", "A", 12.0, 0.010),      # overlaps a; same P, larger score -> peak
        ("c", "chr1", 300, 500, "+", "A", 50.0, 0.001),      # starts at the end of b: a new locus
        ("d", "chr1", 100, 300, "-", "A", 5.0, 0.020),       # other strand
    ])
    loci = sel.merge_loci(frame)
    assert list(loci.locus) == ["chr1:300-500:+", "chr1:0-300:+", "chr1:100-300:-"]
    merged = loci.set_index("locus").loc["chr1:0-300:+"]
    assert merged.n_windows == 2 and merged.peak_window == "b" and merged.peak_score == 12.0


def test_structure_cores_drop_short_stems_and_split_at_outermost_pairs():
    structure = "(((((((...)))((((...)))))))).(((...)))..((((....))))"
    cores = sel.structure_cores(sel.parse_dot_bracket(structure), len(structure))
    # The three-pair stems at 4-12 and 29-37 are removed; a multi-branch core and a hairpin remain.
    assert [(c.start, c.end) for c in cores] == [(0, 28), (40, 52)]
    assert len(cores[0].pairs) == 8 and (4, 12) not in cores[0].pairs and (0, 27) in cores[0].pairs
    assert sel.structure_cores(sel.parse_dot_bracket("((...))"), 7) == []


def test_structure_cores_reject_crossing_pairs():
    with pytest.raises(ValueError):
        sel.structure_cores([(0, 10), (5, 15)], 20)
    with pytest.raises(ValueError):
        sel.structure_cores([(0, 10), (0, 12)], 20)


def hairpin_weights(length, pairs, value=1.0):
    weights = np.zeros((length, length), dtype=np.float32)
    for i, j in pairs:
        weights[i, j] = weights[j, i] = value
    return weights


@pytest.mark.parametrize("strand,expected", [("+", (1005, 1019)), ("-", (1002, 1016))])
def test_window_cores_coordinates_and_score(strand, expected):
    sequence = "AAAAAGGGGAAAAAACCCCAA"                                # 21 nt; hairpin at 5-18
    pairs = [(5, 18), (6, 17), (7, 16), (8, 15)]
    record = {"sequence": sequence, "strand": strand, "peak_start": 1000, "peak_end": 1021, "locus": "L",
              "peak_window": "W", "chrom": "chr1", "peak_score": 3.0}
    rows = sel.window_cores(record, hairpin_weights(21, pairs, 0.5))
    assert len(rows) == 1
    core = rows[0]
    assert (core["core_genome_start"], core["core_genome_end"]) == expected
    assert core["core_sequence"] == "GGGGAAAAAACCCC" and core["core_dot_bracket"] == "((((......))))"
    assert core["core_score"] == pytest.approx(2.0) and core["n_stems"] == 1 and not core["touches_window_edge"]


def test_canonical_key_is_strand_independent():
    sequence, structure = "GGGACUUCGGUCCCAA", "((((....))))...."
    reverse = sequence.translate(str.maketrans("ACGU", "UGCA"))[::-1]
    flipped = structure[::-1].translate(str.maketrans("()", ")("))
    assert sel.canonical_key(sequence, structure) == sel.canonical_key(reverse, flipped)
    assert sel.canonical_key(sequence.replace("U", "T"), structure) == sel.canonical_key(sequence, structure)


def test_collapse_cores_keeps_best_and_source_strand():
    sequence, structure = "GGGACTTCGGTCCCAA", "((((....))))...."
    reverse = "TTGGGACCGAAGTCCC"
    flipped = "....((((....))))"
    cores = pd.DataFrame({
        "locus": ["L1", "L1", "L2", "L3"], "peak_window": ["W1", "W1", "W2", "W3"],
        "core_sequence": [sequence, "GGGGAAACCCC", reverse, "GCGCAAAAGCGC"],
        "core_dot_bracket": [structure, "((((...))))", flipped, "((((....))))"],
        "core_score": [5.0, 1.0, 9.0, 2.0], "component_rank": [1, 2, 1, 1]})
    loci = pd.DataFrame({"peak_window": ["W1", "W2", "W3"], "peak_p": [0.01, 0.01, 0.001]})
    reps, members = sel.collapse_cores(cores, loci, id_prefix="q")
    assert list(reps.query_id) == ["q0000001", "q0000002"]
    assert list(reps.peak_window) == ["W3", "W2"]               # smallest P, then the larger core score wins
    assert reps.loc[1, "sequence_u"] == reverse.replace("T", "U") and reps.loc[1, "core_dot_bracket"] == flipped
    assert members.set_index("peak_window").query_id.to_dict() == {"W3": "q0000001", "W2": "q0000002",
                                                                   "W1": "q0000002"}


def test_genome_chunks_overlap_and_all_n():
    chunks = list(sel.genome_chunks("ACGT" * 7 + "N" * 7, chunk_size=10, overlap=2))
    assert [(s, e) for s, e, _ in chunks] == [(0, 10), (8, 18), (16, 26), (24, 34)]   # (32, 35) is all N
    assert all(len(c) == e - s for s, e, c in chunks)


def test_chunk_genome_writes_ids(tmp_path):
    genome = tmp_path / "g.fa"
    genome.write_text(">chr1 description\n" + "ACGT" * 5000 + "\n>chrN\n" + "N" * 50 + "\n")
    out = sel.chunk_genome(genome, "toy_species", tmp_path / "chunks")
    with gzip.open(out, "rt") as handle:
        names = [line[1:].strip() for line in handle if line.startswith(">")]
    assert out.name == "toy_species.10kb_o500.fa.gz"
    assert names == ["toy_species|chr1|0|10000", "toy_species|chr1|9500|19500", "toy_species|chr1|19000|20000"]


def test_riboseek_commands_use_paper_settings():
    commands = sel.riboseek_commands("riboseek", "q_db", "t.fa.gz", "work", threads=8)
    search = commands[1]
    assert search[:2] == ["riboseek", "search"]
    for flag, value in [("-s", "7.5"), ("-e", "10"), ("--max-seqs", "100"), ("--alt-ali", "0"),
                        ("--prefilter-mode", "0"), ("--threads", "8")]:
        assert search[search.index(flag) + 1] == value
    assert "-a" in search and commands[2][-1] == sel.HIT_FORMAT


def hits_table(rows):
    frame = pd.DataFrame(rows, columns=["query_id", "target_id", "evalue", "bits", "qstart", "qend", "tstart", "tend",
                                        "qcov", "qaln", "taln"])
    frame["alnlen"], frame["tcov"] = frame.qaln.str.len(), 0.0
    return frame[sel.HIT_COLUMNS]


def test_high_confidence_boundaries():
    hits = pd.DataFrame({"evalue": [1e-3, 1.1e-3, 1e-5, 1e-5], "qcov": [0.8, 0.9, 0.79, 1.0]})
    assert list(sel.high_confidence(hits)) == [True, False, False, True]


def test_locate_and_merge_hits_from_overlapping_chunks():
    hits = hits_table([
        # the same locus seen in two overlapping chunks (genome 9,600-9,620)
        ("q1", "sp|chr1|0|10000", 1e-5, 50.0, 1, 20, 9601, 9620, 1.0, "A" * 20, "A" * 20),
        ("q1", "sp|chr1|9500|19500", 1e-4, 40.0, 1, 20, 101, 120, 1.0, "A" * 20, "A" * 20),
        # a paralogue far away, minus strand (qstart > qend)
        ("q1", "sp|chr1|0|10000", 0.05, 20.0, 20, 1, 3000, 3019, 0.5, "ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"),
    ])
    located = sel.locate_hits(hits, "sp")
    assert list(located.genome_start) == [9600, 9600, 2999] and list(located.target_strand) == ["+", "+", "-"]
    assert located.oriented_qaln[2] == sel.reverse_complement("ACGTACGTACGTACGTACGT")
    loci = sel.merge_hit_loci(located)
    assert len(loci) == 2
    plus = loci.loc[loci.target_strand.eq("+")].iloc[0]
    assert (plus.locus_start, plus.locus_end, plus.n_chunk_hits, plus.bits) == (9600, 9620, 2, 50.0)
    assert bool(plus.any_high_confidence) and plus.locus_id == "sp|chr1:9600-9620:+"
    counts = sel.locus_counts(loci)
    assert counts.loc["q1"].tolist() == [1, 2, 1, 1]


def test_project_hit_places_target_on_query_columns():
    # query positions 3-8 of a 10-nt query; one deletion in the target and one insertion (dropped)
    assert sel.project_hit("GGA-UUC", "GG-AUUC", 3, 8, 10) == "--GG-TTC--"
    with pytest.raises(ValueError):
        sel.project_hit("GGAUUC", "GGAUUC", 3, 9, 10)


def test_routes_counts_and_priority():
    lineage = {"p1": "primate", "p2": "primate", "p3": "primate", "p4": "primate", "p5": "primate", "p6": "primate",
               "m1": "nonprimate_mammal", "m2": "nonprimate_mammal", "m3": "nonprimate_mammal",
               "t1": "nonmammal_tetrapod", "o1": "other_vertebrate"}
    hc = pd.DataFrame([("vw", s) for s in ("p1", "t1", "o1", "m1", "m2", "m3", "p2", "p3")]
                      + [("mc", s) for s in ("p1", "p2", "p3", "m1", "m2", "m3")]
                      + [("pf", s) for s in ("p1", "p2", "p3", "p4", "p5", "p6", "p6")]
                      + [("pf_weak", s) for s in ("p1", "p2", "p3", "p4", "p5", "p6")]
                      + [("none", "t1")], columns=["query_id", "species"])
    counts = sel.lineage_counts(hc, lineage)
    assert counts.loc["pf"].tolist() == [6, 0, 0, 0]
    candidates = pd.DataFrame({"query_id": ["vw", "mc", "pf", "pf_weak", "none", "rescue"],
                               "peak_p": [0.01, 0.01, 0.001, 0.0011, 0.01, 0.0005]})
    routed = sel.assign_routes(candidates, counts).set_index("query_id")
    assert routed.route.to_dict() == {"vw": "vertebrate_wide", "mc": "mammal_conserved", "pf": "primate_focused",
                                      "pf_weak": "not_selected", "none": "not_selected", "rescue": "high_score_rescue"}
    assert bool(routed.loc["vw", "route_mammal_conserved"])                # labelled by the first route only
    assert routed.loc["rescue", sel.LINEAGES[0]] == 0


def test_crops_seed_and_expansion(tmp_path):
    assert sel.seed_crop(100, 150) == (50, 200) and sel.seed_crop(20, 50) == (0, 80)
    assert sel.expansion_crop(100, 150, 40) == (65, 185)
    genome = tmp_path / "g.fa.gz"
    contig = "AAAAACCCCCGGGGGTTTTT"
    with gzip.open(genome, "wt") as handle:
        handle.write(f">c1\n{contig}\n")
    loci = pd.DataFrame({"query_id": ["q", "q"], "species": ["sp", "sp"], "contig": ["c1", "c1"],
                         "target_strand": ["+", "-"], "locus_id": ["a", "b"], "evalue": [1e-5, 1e-5],
                         "bits": [30.0, 30.0], "qcov": [1.0, 1.0], "locus_start": [5, 5], "locus_end": [10, 10]})
    crops = sel.extract_crops(genome, loci, "seed").set_index("strand")
    assert crops.loc["+", "crop_sequence"] == contig[0:15]
    assert crops.loc["-", "crop_sequence"] == sel.reverse_complement(contig[0:15])
    assert (crops.loc["-", "locus_offset_start"], crops.loc["-", "locus_offset_end"]) == (5, 10)
    assert crops.loc["+", "locus_offset_start"] == 5


TBLOUT = """#target name accession query name accession mdl mdl from mdl to seq from seq to strand trunc pass gc bias score E-value inc description
c000001 - m - cm 1 20 31 50 + no 1 0.50 0.0 25.0 1e-5 ! -
c000001 - m - cm 1 10 5 14 + no 1 0.50 0.0 3.0 1.0 ? -
c000002 - m - cm 3 18 20 5 - no 1 0.50 0.0 12.0 1e-2 ? -
"""


def test_parse_tblout_and_best_crop_hits():
    hits = sel.parse_tblout(TBLOUT)
    assert len(hits) == 3 and hits.cm_bits.tolist() == [25.0, 3.0, 12.0] and hits.included.tolist()[0]
    crops = pd.DataFrame({"crop_id": ["c000001", "c000002"], "crop_sequence": ["A" * 30 + "C" * 20, "G" * 4 + "ACGT" * 4 + "G" * 6],
                          "model_length": [20, 20], "locus_offset_start": [30, 4], "locus_offset_end": [50, 20]})
    best = sel.best_crop_hits(hits, crops).set_index("crop_id")
    assert best.loc["c000001", "cm_hit_sequence"] == "C" * 20 and best.loc["c000001", "cm_model_coverage"] == 1.0
    assert best.loc["c000002", "cm_hit_sequence"] == sel.reverse_complement(("G" * 4 + "ACGT" * 4)[4:20])
    assert best.loc["c000002", "cm_model_coverage"] == pytest.approx(0.8)
    assert best.loc["c000002", "left_boundary_shift"] == 0


def test_inclusion_rule_and_final_additions():
    seed_bits = [30.0, 40.0, 50.0, 60.0, 70.0]                       # 10th percentile 34 -> threshold 24
    assert sel.inclusion_threshold(seed_bits) == pytest.approx(24.0)
    hits = pd.DataFrame({"species": ["a", "b", "c", "d", "e"], "cm_bits": [24.0, 23.9, 30.0, 31.0, 50.0],
                         "cm_model_coverage": [0.8, 1.0, 0.79, 1.0, 1.0], "riboseek_evalue": [1e-3] * 5,
                         "cm_hit_sequence": ["ACGU", "ACGA", "ACGC", "AAAA", "ACGU"]})
    hits["cm_hit_sequence"] = hits.cm_hit_sequence.str.replace("U", "T")
    judged = sel.apply_inclusion(hits, seed_bits)
    assert judged.accepted.tolist() == [True, False, False, True, True]
    added = sel.final_additions(judged, seed_sequences=["AA-A.A"])      # AAAA is already a seed row
    assert added.cm_hit_sequence.tolist() == ["ACGT"] and added.final_sequence_id.tolist() == ["e00001_e"]
    assert added.represented_loci.tolist() == [2] and added.represented_species.tolist() == [2]


def test_high_copy_flag():
    flags = sel.high_copy(pd.Series([100, 101, 20]), pd.Series([100, 101, 150]))
    assert flags.tolist() == [False, True, True]


def test_sample_heldout_round_robin_one_per_species():
    hits = pd.DataFrame({
        "species": ["p1", "p1", "p2", "m1", "t1", "o1", "o2"],
        "lineage": ["primate", "primate", "primate", "nonprimate_mammal", "nonmammal_tetrapod",
                    "other_vertebrate", "other_vertebrate"],
        "cm_bits": [50.0, 60.0, 55.0, 40.0, 30.0, 28.0, 25.0],
        "cm_hit_sequence": ["ACGTACGT", "ACGTACGA", "ACGTACGT", "TTTTACGT", "GGGGCCCC", "ACGTACGG", "CCCCGGGG"]})
    chosen = sel.sample_heldout(hits, k=4)
    assert chosen.species.tolist() == ["p1", "m1", "t1", "o2"]         # o2: lower bits, but more 4-mer-diverse
    assert chosen.cm_bits.tolist()[0] == 60.0                        # best hit of the species
    assert len(sel.sample_heldout(hits, k=10)) == 6


def test_relative_scores_and_families():
    scores = pd.DataFrame({
        "model_query_id": ["A", "A", "B", "B", "B", "C", "A"],
        "target_query_id": ["B", "B", "A", "A", "C", "B", "A"],
        "sample_kind": ["heldout", "heldout", "heldout", "heldout", "heldout", "heldout", "human"],
        "cm_bits": [10.0, 14.0, 20.0, 30.0, 1.0, 2.0, 99.0]})
    medians = sel.direction_medians(scores)
    assert medians.loc[("A", "B"), "b"] == 12.0 and medians.loc[("B", "A"), "n_heldout"] == 2
    pairs = pd.DataFrame({"unit_a": ["A", "B", "A"], "unit_b": ["B", "C", "D"]})
    edges = sel.relative_scores(pairs, medians, {"A": 20.0, "B": 20.0, "C": 10.0, "D": 10.0})
    assert len(edges) == 2                                           # A-D has no scores
    ab = edges.set_index(["unit_a", "unit_b"]).loc[("A", "B")]
    assert ab.r_ab == pytest.approx(-0.4) and ab.r_ba == pytest.approx(0.25) and ab.ratio_mean == pytest.approx(-0.075)
    families = sel.provisional_families(["D", "C", "B", "A", "E"], edges).set_index("query_id")
    assert families.family_id.to_dict() == {"A": "C01", "B": "C01", "C": "C02", "D": "C03", "E": "C04"}
    assert families.loc["A", "family_size"] == 2
    joined = sel.provisional_families(["A", "B", "C"], edges, threshold=-1.0)
    assert joined.family_id.nunique() == 1


def test_candidate_pairs_are_undirected():
    directed = pd.DataFrame({"source_query_id": ["A", "B", "A", "C"], "target_query_id": ["B", "A", "A", "A"]})
    assert sel.candidate_pairs(directed).values.tolist() == [["A", "B"], ["A", "C"]]


def test_genome_panel_groups():
    panel = sel.load_panel()
    seed, expansion = sel.seed_species(panel), sel.expansion_species(panel)
    assert len(seed) == 20 and len(expansion) == 35 and not set(seed) & set(expansion)
    assert "homo_sapiens" in expansion
    groups = pd.Series(sel.lineage_map(panel)).loc[seed].value_counts().to_dict()
    assert groups == {"primate": 9, "nonprimate_mammal": 4, "nonmammal_tetrapod": 3, "other_vertebrate": 4}


def test_infernal_commands_use_paper_flags():
    search = sel.cmsearch_command("cmsearch", "m.cm", "x.fa", "x.tbl")
    assert search[1:8] == ["--max", "-g", "--cyk", "--acyk", "--notrunc", "--toponly", "--noali"]
    assert search[search.index("-T") + 1] == "-1000"
    assert sel.cmalign_command("cmalign", "m.cm", "x.fa", "o.sto", mapali="s.sto")[1:9] == [
        "--cpu", "1", "--cyk", "--notrunc", "--noprob", "--mapali", "s.sto", "-o"]
    assert sel.cmcalibrate_command("cmcalibrate", "m.cm")[1:3] == ["-L", "0.01"]
    assert sel.cmbuild_command("cmbuild", "m.cm", "s.sto")[1:3] == ["-F", "--hand"]


def expansion_crops():
    return pd.DataFrame({
        "crop_sequence": ["aaGGGGCCCCtt", "AAGGGGCCCCTT", "TTACGUACGTAA", "CCCCCCCCCCCC", "GGAAAAGG"],
        "locus_offset_start": [2, 2, 2, 0, 2], "locus_offset_end": [10, 10, 10, 6, 6],
        "riboseek_evalue": [1e-3, 1e-5, 1e-4, 1e-2, 1e-6], "riboseek_qcov": [0.8, 1.0, 0.95, 1.0, 0.79],
        "species": ["sp_a", "sp_b", "homo_sapiens", "sp_c", "sp_d"], "contig": ["c1", "KQ1.1", "1", "c9", "c2"],
        "locus_start": [100, 200, 300, 400, 500], "strand": ["+", "-", "+", "+", "-"]})


def test_homologue_records_condition_iii():
    seed = [("h00001_sp_x", "GGGG-CCCC"), ("core1", "ggggaaaacccc"), ("h00002_sp_y", "GGGGAAAACCCC")]
    records = sel.homologue_records("core1", seed, expansion_crops())
    names = [name for name, _ in records]
    # human first; the seed copy of the human core and the duplicate expansion locus are removed;
    # E > 1e-3 or coverage < 0.8 loci are left out; '+' and other symbols in names become '_'
    assert names == ["core1", "h00001_sp_x", "x000003_homo_sapiens_1_300__"]
    assert dict(records) == {"core1": "GGGGAAAACCCC", "h00001_sp_x": "GGGGCCCC", "x000003_homo_sapiens_1_300__": "ACGTACGT"}
    with pytest.raises(ValueError):
        sel.homologue_records("missing", seed)


def test_write_homologues_from_pipeline_files(tmp_path):
    pytest.importorskip("pyarrow")
    seed_dir, expand_dir = tmp_path / "seed", tmp_path / "expand"
    sel.write_fasta([("core1", "GGGGAAAACCCC"), ("h00001_sp_x", "GGGAACCC")], seed_dir / "sequences.fa")
    crops = expansion_crops()
    crops.insert(0, "crop_id", [f"c{i:06d}" for i in range(1, len(crops) + 1)])
    sel.write_fasta(zip(crops.crop_id, crops.crop_sequence), expand_dir / "crops.fa", dna=False)
    expand_dir.mkdir(exist_ok=True)
    crops.drop(columns="crop_sequence").to_parquet(expand_dir / "expansion_hits.parquet", index=False)
    assert sel.write_homologues("core1", seed_dir, expand_dir, tmp_path / "out" / "core1.fa") == 4
    written = list(sel.iter_fasta(tmp_path / "out" / "core1.fa"))
    assert [n for n, _ in written] == ["core1", "h00001_sp_x", "x000001_sp_a_c1_100__", "x000003_homo_sapiens_1_300__"]
    assert written[2][1] == "GGGGCCCC" and sel.write_homologues("core1", seed_dir, tmp_path / "none", tmp_path / "s.fa") == 2


def test_cli_registers_steps():
    parser = argparse.ArgumentParser()
    sel.register(parser.add_subparsers(dest="command"))
    for argv in (["select", "cores", "w.parquet", "-o", "out"], ["select", "chunk", "g.fa.gz", "-o", "c"],
                 ["select", "search", "q.fa", "--chunks", "c", "-o", "w"],
                 ["select", "routes", "r.parquet", "--hits", "w", "-o", "s.parquet"],
                 ["select", "seed", "s.parquet", "--hits", "w", "--genomes", "g", "-o", "seed"],
                 ["select", "expand", "s.parquet", "--seed", "seed", "--hits", "w", "--genomes", "g", "-o", "e"],
                 ["select", "homologues", "c.parquet", "--seed", "seed", "--expand", "e", "-o", "h"],
                 ["select", "families", "c.parquet", "--expand", "e", "-o", "f"]):
        assert callable(parser.parse_args(argv).func)


# ---------------------------------------------------------------------------------------------------------------
# Package data, column names, empty inputs and external programs (fake Riboseek/Infernal scripts; Linux /bin/sh)


def test_genome_panel_is_package_data():
    package = Path(sel.__file__).resolve().parent
    assert sel.DEFAULT_PANEL.resolve() == package / "data" / "genome_panel_55.tsv" and sel.DEFAULT_PANEL.is_file()
    checkout = package.parents[1] / "data" / "genome_panel_55.tsv"
    if checkout.is_file():                                    # the copy of the repository must not drift
        assert checkout.read_bytes() == sel.DEFAULT_PANEL.read_bytes()


@pytest.mark.parametrize("stem_column", ["n_stems_scored", "n_stems_ge4"])
def test_merge_loci_carries_scan_and_paper_stem_counts(stem_column):
    frame = windows_table([("a", "chr1", 0, 200, "+", "A", 10.0, 0.01), ("b", "chr1", 100, 300, "+", "A", 12.0, 0.01)])
    frame[stem_column] = [2, 3]
    loci = sel.merge_loci(frame)
    assert loci[stem_column].tolist() == [3]                  # the peak window's count


def write_hits(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    hits_table(rows).to_csv(path, sep="\t", header=False, index=False)


def test_hit_loci_of_species_without_hits_keep_dtypes(tmp_path):
    hc = (1e-5, 50.0, 1, 20, 11, 30, 1.0, "A" * 20, "A" * 20)
    write_hits(tmp_path / "sp1" / "hits.tsv", [("q1", "sp1|c1|0|1000", *hc), ("q2", "sp1|c1|0|1000", *hc),
                                               ("q2", "sp1|c2|0|1000", *hc)])
    weak = (5.0, 10.0, 1, 20, 11, 30, 1.0, "A" * 20, "A" * 20)                 # E = 5: no locus
    write_hits(tmp_path / "sp2" / "hits.tsv", [("q1", "sp2|c1|0|1000", *weak)])
    (tmp_path / "sp3").mkdir()
    (tmp_path / "sp3" / "hits.tsv").write_text("")         # Riboseek writes an empty table without hits
    empty = sel.read_hits(tmp_path / "sp3" / "hits.tsv")
    assert empty.empty and empty.evalue.dtype == np.float64 and list(empty) == sel.HIT_COLUMNS
    loci = sel._seed_loci(tmp_path, ["sp1", "sp2", "sp3"], {"q1", "q2"})
    assert loci.any_high_confidence.dtype == bool and loci.evalue.dtype == np.float64
    counts = sel.locus_counts(loci)
    assert counts.n_loci_hc.tolist() == [1, 2] and counts.n_loci_hc.dtype == np.int64
    none = sel._seed_loci(tmp_path, ["sp2", "sp3"], {"q1", "q2"})
    assert none.empty and {"locus_id", "any_high_confidence", "species", "target_sequence"} <= set(none)
    assert sel.locus_counts(none).empty


class ZeroStudent:
    """Stand-in for the student: maps without any dependency, so that no stem is decoded."""

    def predict(self, sequences, batch_size=16):
        return [np.zeros((len(s), len(s)), dtype=np.float32) for s in sequences]


def test_extract_cores_without_cores_gives_typed_empty_table():
    pytest.importorskip("RNA")
    loci = pd.DataFrame({"locus": ["L"], "peak_window": ["W"], "chrom": ["chr1"], "strand": ["+"], "peak_start": [0],
                         "peak_end": [24], "peak_score": [1.0], "sequence": ["GGGGAAAACCCCAAAAGGGGAAAA"]})
    cores = sel.extract_cores(loci, ZeroStudent())
    assert cores.empty and list(cores) == [*sel.CORE_DTYPES, "score_fraction", "component_rank"]
    assert cores.core_score.dtype == np.float64


def fake_program(directory, name, script):
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / name
    path.write_text("#!/bin/sh\n" + script)
    path.chmod(0o755)
    return str(path)


FAKE_RIBOSEEK = """echo "$1" >> "$(dirname "$0")/calls.log"
case "$1" in
  createdb) cp "$2" "$3" ;;
  search) : > "$4" ;;
  convertalis) : > "$5" ;;
  rmdb) rm -f "$2" ;;
esac
"""


def test_search_refuses_a_work_directory_of_other_queries(tmp_path):
    riboseek = fake_program(tmp_path / "bin", "riboseek", FAKE_RIBOSEEK)
    queries, other = tmp_path / "q.fa", tmp_path / "other.fa"
    queries.write_text(">q1\nGGGGAAAACCCC\n")
    other.write_text(">q2\nGCGCAAAAGCGC\n")
    (tmp_path / "chunks").mkdir()
    (tmp_path / "chunks" / "sp.10kb_o500.fa.gz").write_bytes(b"")
    chunks, work = {"sp": tmp_path / "chunks" / "sp.10kb_o500.fa.gz"}, tmp_path / "work"
    assert sel.run_riboseek(queries, chunks, work, riboseek).status.tolist() == ["success"]
    assert sel.run_riboseek(queries, chunks, work, riboseek).status.tolist() == ["reused"]
    assert (work / "queries.sha256").read_text().split()[0] == sel.file_sha256(queries)
    command = ["search", str(other), "--chunks", str(tmp_path / "chunks"), "--species", "sp", "--riboseek", riboseek,
               "-o", str(work)]
    with pytest.raises(SystemExit, match="other queries.*--force"):
        sel.main(command)
    sel.main(command + ["--force"])                          # rebuilds query_db and searches the species again
    assert (work / "query_db").read_text() == other.read_text()
    assert (work / "queries.sha256").read_text().split()[0] == sel.file_sha256(other)
    assert pd.read_csv(work / "search_status.tsv", sep="\t").status.tolist() == ["success"]
    assert (tmp_path / "bin" / "calls.log").read_text().split().count("createdb") == 4   # query and target, twice
    (work / "queries.sha256").unlink()                       # unrecorded queries are refused as well
    with pytest.raises(FileExistsError):
        sel.run_riboseek(other, chunks, work, riboseek)


def test_steps_name_a_missing_program_before_any_task(tmp_path):
    for argv, program in ((["seed", "s.tsv", "--hits", "h", "--genomes", "g"], "cmbuild"),
                          (["expand", "s.tsv", "--seed", "x", "--hits", "h", "--genomes", "g"], "cmbuild"),
                          (["families", "c.tsv", "--expand", "e"], "cmsearch")):
        with pytest.raises(SystemExit, match=f"{program} not found in"):
            sel.main([*argv, "--infernal-bin", str(tmp_path / "none"), "-o", str(tmp_path / "out")])
    with pytest.raises(SystemExit, match="not found or not executable"):
        sel.main(["search", "q.fa", "--chunks", "c", "--riboseek", str(tmp_path / "none"), "-o", str(tmp_path / "w")])


def test_check_manifest_counts_failed_tasks(tmp_path, capsys):
    manifest = pd.DataFrame({"query_id": ["a", "b"], "status": ["success", "failed"],
                             "error": ["", "RuntimeError('x')"]})
    sel._check_manifest(manifest, "select seed", tmp_path / "m.parquet")
    assert "select seed: 1 of 2 tasks failed; first: b: RuntimeError('x')" in capsys.readouterr().err
    with pytest.raises(SystemExit, match="2 of 2 tasks failed; first: a"):
        sel._check_manifest(manifest.assign(status="failed"), "select seed", tmp_path / "m.parquet")
    sel._check_manifest(manifest.assign(status="success"), "select seed", tmp_path / "m.parquet")
    assert capsys.readouterr().err == ""
    assert list(sel._run_parallel([], 2)) == ["query_id", "status", "error"]


def small_panel(path):
    pd.DataFrame({"ensembl_slug": ["sp1", "sp2", "sp3"], "lineage_group": ["primate", "primate", "other_vertebrate"],
                  "search_stage": ["seed search (20 species)"] * 2 + ["expansion (35 genomes)"]}).to_csv(
        path, sep="\t", index=False)
    return str(path)


def test_seed_exits_non_zero_when_every_task_fails(tmp_path):
    pytest.importorskip("pyarrow")
    infernal = tmp_path / "infernal"
    for name in sel.INFERNAL_PROGRAMS:
        fake_program(infernal, name, "exit 1\n" if name == "cmbuild" else "exit 0\n")
    hc = (1e-5, 50.0, 1, 20, 11, 30, 1.0, "A" * 20, "A" * 20)
    write_hits(tmp_path / "hits" / "sp1" / "hits.tsv", [("q1", "sp1|c1|0|1000", *hc), ("q2", "sp1|c1|0|1000", *hc),
                                                        ("q2", "sp1|c2|0|1000", *hc)])
    write_hits(tmp_path / "hits" / "sp2" / "hits.tsv", [("q1", "sp2|c1|0|1000", 5.0, 10.0, 1, 20, 11, 30, 1.0,
                                                         "A" * 20, "A" * 20)])
    (tmp_path / "genomes").mkdir()
    (tmp_path / "genomes" / "sp1.fa").write_text(">c1\n" + "A" * 100 + "\n>c2\n" + "A" * 100 + "\n")
    pd.DataFrame({"query_id": ["q1", "q2"], "sequence_u": ["A" * 20] * 2, "core_dot_bracket": ["." * 20] * 2}).to_csv(
        tmp_path / "selected.tsv", sep="\t", index=False)
    out = tmp_path / "seed"
    with pytest.raises(SystemExit, match="select seed: 2 of 2 tasks failed; first: q1: RuntimeError"):
        sel.main(["seed", str(tmp_path / "selected.tsv"), "--hits", str(tmp_path / "hits"), "--genomes",
                  str(tmp_path / "genomes"), "--panel", small_panel(tmp_path / "panel.tsv"), "--infernal-bin",
                  str(infernal), "--workers", "1", "-o", str(out)])
    counts = pd.read_parquet(out / "seed_locus_counts.parquet").set_index("query_id")
    assert counts.n_loci_hc.to_dict() == {"q1": 1, "q2": 2}
    assert pd.read_parquet(out / "seed_manifest.parquet").status.eq("failed").all()
    # 'select expand' on that output stops with a message instead of a missing-column error
    with pytest.raises(SystemExit, match="no candidate with a seed CM"):
        sel.main(["expand", str(tmp_path / "selected.tsv"), "--seed", str(out), "--hits", str(tmp_path / "hits"),
                  "--genomes", str(tmp_path / "genomes"), "--panel", str(tmp_path / "panel.tsv"),
                  "--infernal-bin", str(infernal), "-o", str(tmp_path / "expand")])


FAKE_CMSEARCH = """while [ $# -gt 0 ]; do
  if [ "$1" = "--tblout" ]; then printf '# no hits\\n' > "$2"; fi
  shift
done
"""


def test_families_without_pairs_are_singletons(tmp_path):
    pytest.importorskip("pyarrow")
    fake_program(tmp_path / "infernal", "cmsearch", FAKE_CMSEARCH)
    candidates = pd.DataFrame({"query_id": ["q2", "q1"], "sequence_u": ["GGGGAAAACCCC", "GCGCAAAAGCGC"],
                               "threshold_bits": [20.0, 15.0]})
    families, edges = sel.score_families(candidates, tmp_path / "expand", tmp_path / "out", tmp_path / "infernal", 2)
    assert families.family_id.tolist() == ["C01", "C02"] and families.query_id.tolist() == ["q1", "q2"]
    assert edges.empty and list(edges) == list(sel.EDGE_DTYPES)
    fake_program(tmp_path / "failing", "cmsearch", "exit 1\n")
    with pytest.raises(RuntimeError, match="HMM searches: 2 of 2 tasks failed; first: q1: RuntimeError"):
        sel.score_families(candidates, tmp_path / "expand", tmp_path / "out2", tmp_path / "failing", 2)


def test_expand_without_expansion_loci_keeps_the_seed_alignment(tmp_path):
    pytest.importorskip("pyarrow")
    infernal = tmp_path / "infernal"
    for name in sel.INFERNAL_PROGRAMS:
        fake_program(infernal, name, 'touch "$3"\n' if name == "cmbuild" else "exit 0\n")   # cmbuild -F --hand <cm>
    seed, query = tmp_path / "seed", "GGGGAAAACCCC"
    rows = [("q1", query), ("h00001_sp1", "GGGGAAAACCCA")]
    (seed / "candidates" / "q1").mkdir(parents=True)
    sel.write_stockholm(seed / "candidates" / "q1" / "realigned_3l.sto", "q1", rows, "((((....))))")
    sel.write_fasta(rows, seed / "candidates" / "q1" / "sequences.fa")
    pd.DataFrame({"query_id": ["q1"], "status": ["success"], "error": [""], "seed_cm": [True],
                  "threshold_bits": [20.0]}).to_parquet(seed / "seed_manifest.parquet", index=False)
    pd.DataFrame({"query_id": ["q1"], "n_species_e01": [1], "n_loci_e01": [1], "n_species_hc": [1],
                  "n_loci_hc": [1]}).to_parquet(seed / "seed_locus_counts.parquet", index=False)
    pd.DataFrame({"query_id": ["q1"], "sequence_u": [query], "core_length": [12]}).to_csv(
        tmp_path / "selected.tsv", sep="\t", index=False)
    write_hits(tmp_path / "hits" / "sp3" / "hits.tsv", [("q1", "sp3|c1|0|1000", 5.0, 10.0, 1, 12, 11, 22, 1.0,
                                                         "A" * 12, "A" * 12)])     # E = 5: no expansion locus
    (tmp_path / "genomes").mkdir()
    sel.main(["expand", str(tmp_path / "selected.tsv"), "--seed", str(seed), "--hits", str(tmp_path / "hits"),
              "--genomes", str(tmp_path / "genomes"), "--panel", small_panel(tmp_path / "panel.tsv"),
              "--infernal-bin", str(infernal), "--workers", "1", "-o", str(tmp_path / "expand")])
    table = pd.read_parquet(tmp_path / "expand" / "candidates.parquet")
    assert table.n_accepted_expansion_loci.tolist() == [0] and table.n_final55_species.tolist() == [1]
    assert table.n_final_alignment_sequences.tolist() == [2]
    assert [n for n, _ in sel.iter_fasta(tmp_path / "expand" / "homologues" / "q1.fa")] == ["q1", "h00001_sp1"]
