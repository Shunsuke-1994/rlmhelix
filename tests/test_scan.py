"""Genome-scan steps on small synthetic inputs (no GPU, network or external tools; tool tests skip without them)."""
import argparse
import random
import shutil

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from rlmhelix import scan

CHROM_LENGTH = 2000


def gtf_line(start, end, strand, feature, attributes):
    """GTF row for a 0-based half-open interval."""
    text = " ".join(f'{key} "{value}";' for key, value in attributes.items())
    return "\t".join(["chr1", "TEST", feature, str(start + 1), str(end), ".", strand, ".", text])


@pytest.fixture
def toy(tmp_path):
    """chr1 of 2,000 nt: one protein-coding gene (+), lncRNAs on both strands, an N at 650, soft-masked bases."""
    pc = {"gene_id": "G1", "gene_type": "protein_coding"}
    t1 = {**pc, "transcript_id": "T1", "transcript_type": "protein_coding"}
    t2 = {**pc, "transcript_id": "T2", "transcript_type": "retained_intron"}  # not a protein-coding transcript
    lines = [gtf_line(300, 1300, "+", "gene", pc)]
    lines += [gtf_line(s, e, "+", "exon", t1) for s, e in [(300, 500), (800, 1000), (1200, 1300)]]
    lines += [gtf_line(s, e, "+", "CDS", t1) for s, e in [(400, 500), (800, 1000), (1200, 1250)]]
    lines += [gtf_line(450, 900, "+", "exon", t2)]
    for gene_id, (s, e, strand) in {"L1": (0, 150, "+"), "L2": (950, 1150, "+"), "L3": (1500, 1750, "-")}.items():
        lines.append(gtf_line(s, e, strand, "gene", {"gene_id": gene_id, "gene_type": "lncRNA"}))
    gtf = tmp_path / "toy.gtf"
    gtf.write_text("##description: toy\n" + "\n".join(lines) + "\n")
    rng = random.Random(0)
    sequence = [rng.choice("ACGT") for _ in range(CHROM_LENGTH)]
    sequence[650] = "N"
    sequence = "".join(sequence)
    sequence = sequence[:100] + sequence[100:180].lower() + sequence[180:]
    fasta = tmp_path / "toy.fa"
    body = "\n".join(sequence[k:k + 60] for k in range(0, len(sequence), 60))
    fasta.write_text(">chrM\nACGT\n>chr1 toy\n" + body + "\n")
    return gtf, fasta, sequence.upper()


def test_mrna_regions_name_utrs_by_strand():
    exons, cds = [(300, 400), (100, 200)], [(150, 200), (300, 350)]
    assert sorted(scan.mrna_regions(exons, cds, "+")) == [(100, 150, "utr5"), (200, 300, "intron"), (350, 400, "utr3")]
    assert sorted(scan.mrna_regions(exons, cds, "-")) == [(100, 150, "utr3"), (200, 300, "intron"), (350, 400, "utr5")]
    assert scan.mrna_regions(exons, [], "+") == []


def test_annotation_regions(toy):
    gtf, _, _ = toy
    annotation = scan.read_annotation(gtf)
    regions = [tuple(r) for r in annotation.regions[["start", "end", "strand", "feature"]].itertuples(index=False)]
    assert sorted(regions) == [(300, 400, "+", "utr5"), (500, 800, "+", "intron"), (1000, 1200, "+", "intron"),
                               (1250, 1300, "+", "utr3")]
    assert annotation.genes["chr1"] == [(0, 150), (300, 1300), (1500, 1750)]
    assert annotation.lncrna[("chr1", "-")] == [(1500, 1750)]


def test_window_placement(toy):
    gtf, fasta, sequence = toy
    annotation = scan.read_annotation(gtf)
    (name, chromosome), = scan.iter_chromosomes(fasta, ["chr1"])
    assert name == "chr1" and chromosome == sequence
    windows = scan.make_windows("chr1", chromosome, annotation)
    # centres at multiples of 100 inside the regions; 500 and 600 contain the N; the lncRNA centre 0 starts at -100
    expected = [(0, "+", "lncrna"), (200, "+", "utr5"), (400, "+", "intron"), (900, "+", "intron;lncrna"),
                (1000, "+", "intron;lncrna"), (1400, "-", "lncrna"), (1500, "-", "lncrna"), (1600, "-", "lncrna")]
    assert list(zip(windows.start, windows.strand, windows.features)) == expected
    assert windows.attrs["n_placed"] == len(expected) + 2
    assert (windows.end - windows.start == 200).all() and (windows.center == windows.start + 100).all()
    first, last = windows.iloc[0], windows.iloc[-1]
    assert first.window_id == "chr1:0-200:+" and first.sequence == sequence[0:200]
    assert last.sequence == scan.reverse_complement(sequence[1600:1800])
    with_intergenic = scan.make_windows("chr1", chromosome, annotation, intergenic=True)
    intergenic = with_intergenic[scan.is_intergenic(with_intergenic.features)]
    assert list(zip(intergenic.start, intergenic.strand)) == [(1300, "+"), (1300, "-"), (1800, "+"), (1800, "-")]
    assert with_intergenic[~scan.is_intergenic(with_intergenic.features)].window_id.tolist() == windows.window_id.tolist()


def test_fasta_index_and_stream_agree(toy):
    _, fasta, sequence = toy
    text = fasta.read_text()
    offset = text.index("\n", text.index(">chr1")) + 1
    fasta.with_name("toy.fa.fai").write_text(f"chrM\t4\t6\t4\t5\nchr1\t{len(sequence)}\t{offset}\t60\t61\n")
    indexed = dict(scan.iter_chromosomes(fasta, ["chr1", "chrM"]))
    assert indexed == {"chr1": sequence, "chrM": "ACGT"}


def test_repeat_overlap_matches_brute_force():
    rng = np.random.default_rng(1)
    intervals = scan.merge_intervals((int(s), int(s + rng.integers(1, 80))) for s in rng.integers(0, 5000, 120))
    starts, ends = np.array(intervals).T
    window_starts = rng.integers(0, 5000, 500)
    covered = np.zeros(6000, dtype=bool)
    for s, e in intervals:
        covered[s:e] = True
    expected = [covered[s:s + 200].sum() for s in window_starts]
    assert scan.repeat_overlap(window_starts, window_starts + 200, starts, ends).tolist() == expected
    assert scan.merge_intervals([(5, 10), (10, 12), (0, 3)]) == [(0, 3), (5, 12)]


def test_filter_repeats_without_tools(toy):
    gtf, fasta, sequence = toy
    windows = scan.make_windows("chr1", sequence, scan.read_annotation(gtf))
    repeats = (np.array([159, 1050]), np.array([200, 1090]))  # 41 nt in window 0-200; 40 nt in 900-1100, 1000-1200
    kept, counts = scan.filter_repeats(windows, repeats, use_tantan=False, use_irf=False)
    assert counts == {"n_windows": 8, "n_repeatmasker": 1, "n_kept": 7}
    assert kept.start.tolist() == [200, 400, 900, 1000, 1400, 1500, 1600]
    assert kept.repeatmasker_bp.tolist() == [0, 0, 40, 40, 0, 0, 0]


def test_prepare_step_without_tools(toy, tmp_path):
    pytest.importorskip("pyarrow")
    gtf, fasta, _ = toy
    bed = tmp_path / "repeats.bed"
    bed.write_text("chr1\t159\t200\nchr2\t0\t100\n")
    counts = scan.run_prepare(gtf, fasta, bed, tmp_path / "windows", ["chr1"], use_tantan=False, use_irf=False)
    row = counts.iloc[0]
    assert (row.n_placed, row.n_non_acgt, row.n_repeatmasker, row.n_kept) == (10, 2, 1, 7)
    assert len(pd.read_parquet(tmp_path / "windows" / "chr1.parquet")) == 7


class RecordingStudent:
    """Stands in for `Student`: zero maps, and records the number of sequences per call."""

    def __init__(self):
        self.calls = []

    def predict(self, sequences, batch_size=16):
        self.calls.append(len(sequences))
        return [np.zeros((len(s), len(s)), dtype=np.float32) for s in sequences]


def test_score_windows_fills_the_last_batch():
    pytest.importorskip("RNA")
    rng = random.Random(6)
    sequences = ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(21)]
    windows = pd.DataFrame({"window_id": [f"w{k}" for k in range(21)], "sequence": sequences})
    student = RecordingStudent()
    scored = scan.score_windows(windows, student, batch_size=4, n_jobs=1, chunk=10)
    assert student.calls == [12, 12, 4]  # chunks of 10, 10 and 1 windows, each filled up to a multiple of 4
    assert scored.window_id.tolist() == windows.window_id.tolist() and (scored.score == 0).all()
    assert scored.p_value.between(0, 1).all() and scored.GC.between(0, 1).all()


def test_benjamini_hochberg():
    p = np.array([0.01, 0.04, 0.03, 0.005])
    # ranks 2, 4, 3, 1: p*n/rank = 0.02, 0.04, 0.04, 0.02, then the running minimum from the largest rank
    np.testing.assert_allclose(scan.benjamini_hochberg(p), [0.02, 0.04, 0.04, 0.02])
    rng = np.random.default_rng(2)
    p = np.concatenate([rng.uniform(size=300), rng.uniform(0, 1e-3, size=30), [0.5, 0.5]])
    n = len(p)
    ranked = np.sort(p)
    reference = [min(1.0, min(ranked[k] * n / (k + 1) for k in range(np.searchsorted(ranked, x), n))) for x in p]
    np.testing.assert_allclose(scan.benjamini_hochberg(p), reference)
    assert len(scan.benjamini_hochberg([])) == 0


def test_merge_loci():
    windows = pd.DataFrame({
        "window_id": ["a", "b", "c", "d", "e", "f"],
        "chrom": ["chr1"] * 5 + ["chr2"],
        "start": [100, 200, 400, 300, 200, 0],
        "end": [300, 400, 600, 500, 400, 200],
        "strand": ["+", "+", "+", "+", "-", "+"],
        "features": ["intron", "intron", "utr3", "intron;lncrna", "lncrna", "utr5"],
        "score": [10.0, 30.0, 5.0, 30.0, 8.0, 7.0],
        "p_value": [1e-3, 1e-6, 1e-2, 1e-5, 1e-4, 1e-3],
        "q_value": [0.04, 0.01, 0.05, 0.02, 0.03, 0.04]})
    loci = merge_loci_by_id(scan.merge_loci(windows))
    assert set(loci) == {"chr1:100-600:+", "chr1:200-400:-", "chr2:0-200:+"}
    plus = loci["chr1:100-600:+"]
    assert (plus.n_windows, plus.peak_window, plus.peak_score, plus.min_p, plus.min_q) == (4, "b", 30.0, 1e-6, 0.01)
    assert plus.features == "intron;lncrna;utr3"
    # windows that only touch are not merged
    touching = windows.iloc[[0]].assign(window_id="g", start=300, end=500)
    assert len(scan.merge_loci(pd.concat([windows.iloc[[0]], touching]))) == 2
    assert scan.merge_loci(windows).locus.iloc[0] == "chr1:100-600:+"  # smallest q first


def merge_loci_by_id(loci):
    return {row.locus: row for row in loci.itertuples(index=False)}


def mixture(seed, n=2_000_000):
    rng = np.random.default_rng(seed)
    null = rng.normal(0.7, 0.55, int(0.95 * n))
    return np.concatenate([null, rng.normal(3.0, 1.0, n - len(null))])


def test_empirical_null_recovers_the_centre():
    z = mixture(3)
    null = scan.empirical_null(z)
    assert abs(null.mu0 - 0.7) < 0.02 and abs(null.sigma0 - 0.55) < 0.02 and abs(null.pi0 - 0.95) < 0.02
    np.testing.assert_allclose(scan.z_scores(stats.norm.sf(z[:1000])), z[:1000], atol=1e-6)
    p_empirical = scan.empirical_pvalues([null.mu0], null)
    assert p_empirical[0] == pytest.approx(0.5)


def test_local_fdr_and_tail_fdr():
    z = mixture(4)
    null = scan.empirical_null(z)
    lfdr = scan.local_fdr(z, null)
    np.testing.assert_array_equal(lfdr, scan.local_fdr(z, (null.mu0, null.sigma0, null.pi0)))
    f0, f1 = 0.95 * stats.norm.pdf(z, 0.7, 0.55), 0.05 * stats.norm.pdf(z, 3.0, 1.0)
    middle = (z > 1.5) & (z < 4)
    assert np.abs(lfdr - f0 / (f0 + f1))[middle].max() < 0.15
    assert lfdr.min() >= 0 and lfdr.max() <= 1
    right = np.argsort(z)[z[np.argsort(z)] >= null.mu0]
    assert (np.diff(lfdr[right]) <= 1e-12).all()  # non-increasing to the right of mu0
    table = scan.tail_fdr(z, null, (0.05, 0.01), {"top 1000": 1000}, lfdr)
    assert table.n_windows.tolist() == [100_000, 20_000, 1000]
    top5 = table.iloc[0]
    assert abs(top5.fdr_empirical_null - top5.mean_lfdr) < 0.03  # the tail Fdr is the mean lfdr above the threshold
    assert top5.signal_windows == pytest.approx(top5.n_windows * (1 - top5.fdr_empirical_null))


def test_fdr_step(tmp_path):
    pytest.importorskip("pyarrow")
    rng = np.random.default_rng(5)
    for k, chrom in enumerate(["chr1", "chr2"]):
        n = 300_000
        p = np.concatenate([rng.uniform(size=n - 200), rng.uniform(0, 1e-7, 200)])
        start = np.arange(n) * 100
        pd.DataFrame({"window_id": [f"{chrom}:{s}-{s + 200}:+" for s in start], "chrom": chrom, "start": start,
                      "end": start + 200, "strand": "+", "features": "intron", "score": -np.log(p),
                      "p_value": p}).to_parquet(tmp_path / f"{chrom}.parquet")
    result = scan.run_fdr([tmp_path / "chr1.parquet", tmp_path / "chr2.parquet"], tmp_path / "fdr")
    windows = result["windows"]
    assert len(windows) == 600_000
    assert len(result["bh_windows"]) == int((scan.benjamini_hochberg(windows.p_value) <= 0.05).sum()) >= 400
    assert result["summary"].n_bh_loci.iloc[0] == len(result["bh_loci"])
    written = pd.read_parquet(tmp_path / "fdr" / "eligible_windows.parquet")
    assert written.window_id.tolist() == windows.window_id.tolist()
    np.testing.assert_array_equal(written.lfdr, windows.lfdr)
    assert written.bh_accepted.sum() == len(result["bh_windows"]) and written.bh_accepted.dtype == bool
    for name in ["eligible_windows.parquet", "bh_windows.tsv", "bh_loci.tsv", "empirical_null.tsv", "tail_fdr.tsv"]:
        assert (tmp_path / "fdr" / name).exists()


def test_parse_irf_output():
    text = ("@r000000000\n117 153 37 164 200 37 10 97.3 0.0 65 82.4 17.6 83.3 16.7 0.0 317 317 TCAT ATGA\n"
            "@r000000001\n@r000000002\n1 36 36 79 114 36 42 91.7 0.0 45 73.6 26.4 75.8 24.2 0.0 115 115 TTCC GGAA\n")
    assert scan.parse_irf_ngs(text) == {"r000000000", "r000000002"}


@pytest.mark.skipif(shutil.which("tantan") is None, reason="tantan not on PATH")
def test_tantan_masks_tandem_repeats():
    rng = random.Random(3)
    plain = "".join(rng.choice("ACGT") for _ in range(200))
    repeat = plain[:80] + "CA" * 20 + plain[120:]
    masked = scan.tantan_masked_bp([plain, repeat])
    assert masked[1] >= 20 and masked[0] < masked[1]


@pytest.mark.skipif(not any(shutil.which(n) for n in scan.IRF_NAMES), reason="IRF not on PATH")
def test_irf_finds_a_hairpin():
    rng = random.Random(4)
    arm = "".join(rng.choice("ACGT") for _ in range(30))
    loop = "".join(rng.choice("ACGT") for _ in range(20))
    hairpin = "A" * 60 + arm + loop + scan.reverse_complement(arm) + "T" * 60
    plain = "".join(rng.choice("ACGT") for _ in range(200))
    assert scan.irf_positive([plain, hairpin, plain], shard=1, threads=2).tolist() == [False, True, False]


def test_cli_registration():
    parser = argparse.ArgumentParser()
    scan.register(parser.add_subparsers(dest="command"))
    args = parser.parse_args(["scan", "prepare", "--gtf", "a.gtf", "--genome", "g.fa", "--rmsk", "r.txt", "-o", "w",
                              "--chrom", "chr21", "chr22"])
    assert args.func is scan.run and args.scan_step == "prepare" and args.chrom == ["chr21", "chr22"]
    args = parser.parse_args(["scan", "score", "--windows", "w", "-o", "s"])
    assert args.batch_size == 16 and args.chrom == list(scan.CHROMOSOMES)
    assert parser.parse_args(["scan", "fdr", "--scores", "s", "-o", "f"]).q == 0.05
