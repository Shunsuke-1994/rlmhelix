"""Genome scanning.

The scan of the paper in four steps (`rlmhelix scan prepare | score | fdr`, or the functions below):

1. Windows (`read_annotation`, `make_windows`). 200-nt windows with a 100-nt step, whose centres (multiples of
   100 in 0-based genome coordinates) lie in the introns or UTRs of GENCODE protein-coding mRNAs (gene and transcript
   type protein_coding; UTRs are the exon parts outside the CDS span) or within lncRNA genes. Each window is read
   5'->3' on the strand of its annotation, so a position covered on both strands gives two windows. Windows with a
   base other than A, C, G or T are dropped. Windows lying entirely between genes (intergenic) were scored for the
   paper but excluded before multiple testing; they are made only on request.
2. Repeat filters (`filter_repeats`). Windows are removed when more than 40 nt overlap RepeatMasker (UCSC rmsk, all
   classes), when tantan (`tantan -f 3`) masks at least 20 nt, or when IRF 3.09 finds an inverted repeat with an
   alignment score of at least 40 (match 2, mismatch 7, indel 7, PM 80, PI 10, maximum length 10,000, maximum loop 200).
3. Scoring (`score_windows`). Student APC maps in batches of 16 (bf16 on a GPU), the RLMHelix score and the
   calibrated P value of every window.
4. False discovery (`benjamini_hochberg`, `merge_loci`, `empirical_null`, `tail_fdr`, `local_fdr`). Benjamini-Hochberg
   across all eligible windows; windows with q <= 0.05 that overlap on the same strand are merged into loci. Efron's
   empirical null N(mu0, sigma0^2) is fitted to the centre of z = -Phi^-1(P) by central matching, and the local false
   discovery rate lfdr(z) = pi0 f0(z) / f(z) uses the mixture density f from Lindsey's method.

Steps 1-2 and 4 need no GPU; tantan and IRF are external programs (found on PATH or given as paths).
"""
from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import numpy as np
import pandas as pd

WINDOW = 200
STEP = 100
CHROMOSOMES = tuple([f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"])
RMSK_MAX_BP = 40                 # windows with more RepeatMasker overlap are removed (20% of 200 nt)
TANTAN_MIN_BP = 20               # windows with at least this many tantan-masked bases are removed
IRF_PARAMETERS = ("2", "7", "7", "80", "10", "40", "10000", "200",   # match mismatch indel PM PI minscore maxlen maxloop
                  "-h", "-ngs", "-t4", "200", "-t5", "200", "-t7", "200")
IRF_NAMES = ("irf", "irf.exe", "irf309.linux.exe", "irf309")
IRF_SHARD = 50_000
BH_LEVEL = 0.05
TOP_FRACTIONS = (0.02, 0.01, 0.005, 0.002, 0.001, 0.0005, 0.0001)
LINDSEY_KNOT_QUANTILES = (0.01, 0.05, 0.2, 0.5, 0.8, 0.95, 0.99, 0.995, 0.999, 0.9999)
WINDOW_COLUMNS = ["window_id", "chrom", "start", "end", "center", "strand", "features", "sequence"]

_ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')
_GENE_TYPE = re.compile(r'gene_type "([^"]+)"')
_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_ACGT = frozenset("ACGT")


# ----------------------------------------------------------------------------------------------- input files


def _open_text(path: str | Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_table(path: str | Path, columns: list[str] | None = None) -> pd.DataFrame:
    """Read a table written by this module (Parquet, or tab-separated text for any other suffix)."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path, columns=columns)
    return pd.read_csv(path, sep="\t", usecols=columns)


def write_table(frame: pd.DataFrame, path: str | Path) -> None:
    """Write Parquet (needs pyarrow) or tab-separated text, atomically."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    if path.suffix == ".parquet":
        frame.to_parquet(temporary, index=False)
    else:
        frame.to_csv(temporary, sep="\t", index=False)
    temporary.replace(path)


def iter_chromosomes(fasta: str | Path, names: Iterable[str]) -> Iterator[tuple[str, str]]:
    """Yield (name, upper-case sequence) for the requested chromosomes of a genome FASTA.

    An uncompressed FASTA with a samtools `.fai` index is read by random access in the order of `names`; otherwise
    the file (plain or gzip) is streamed once and chromosomes are yielded in file order.
    """
    wanted = list(names)
    fai = Path(str(fasta) + ".fai")
    if fai.exists() and Path(fasta).suffix != ".gz":
        rows = (line.split("\t") for line in fai.read_text().splitlines())
        index = {f[0]: tuple(map(int, f[1:5])) for f in rows if len(f) >= 5}
        with open(fasta, "rb") as handle:
            for name in wanted:
                if name not in index:
                    raise KeyError(f"{name} not in {fai}")
                length, offset, line_bases, line_width = index[name]
                handle.seek(offset)
                raw = handle.read(length + (length // line_bases + 1) * (line_width - line_bases))
                yield name, raw.replace(b"\n", b"").replace(b"\r", b"")[:length].decode().upper()
        return
    remaining = set(wanted)
    name, chunks = None, []
    with _open_text(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                if name in remaining:
                    remaining.discard(name)
                    yield name, "".join(chunks).upper()
                    if not remaining:
                        return
                fields = line[1:].split()
                name, chunks = (fields[0] if fields else ""), []
            elif name in remaining:
                chunks.append(line.strip())
    if name in remaining:
        remaining.discard(name)
        yield name, "".join(chunks).upper()
    if remaining:
        raise KeyError(f"not in {fasta}: {sorted(remaining)}")


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge half-open intervals that overlap or touch."""
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def read_repeats(path: str | Path, chromosomes: Iterable[str] = CHROMOSOMES) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Merged repeat intervals per chromosome, from UCSC `rmsk.txt(.gz)` (columns genoName, genoStart, genoEnd; all
    repeat classes) or, for a `.bed`/`.bed.gz` file, from the first three BED columns."""
    columns = [0, 1, 2] if ".bed" in Path(path).suffixes else [5, 6, 7]
    table = pd.read_csv(path, sep="\t", header=None, usecols=columns, dtype=str)
    table.columns = ["chrom", "start", "end"]
    table = table[table.chrom.isin(set(chromosomes))]
    table = table.assign(start=pd.to_numeric(table.start, errors="coerce"),
                         end=pd.to_numeric(table.end, errors="coerce")).dropna()
    out = {}
    for chrom, group in table.groupby("chrom", sort=False):
        merged = np.asarray(merge_intervals(zip(group.start.astype(np.int64), group.end.astype(np.int64))),
                            dtype=np.int64)
        out[chrom] = (merged[:, 0], merged[:, 1])
    return out


@dataclass
class Annotation:
    """Scan targets from a GENCODE GTF.

    regions: intron, utr5 and utr3 intervals of protein-coding mRNAs, merged per chromosome, strand and feature
        (0-based, half-open; columns chrom, start, end, strand, feature).
    genes: merged spans of all genes (any type and strand) per chromosome, which delimit intergenic windows.
    lncrna: merged lncRNA gene spans per (chromosome, strand).
    """

    regions: pd.DataFrame
    genes: dict[str, list[tuple[int, int]]]
    lncrna: dict[tuple[str, str], list[tuple[int, int]]]


def _subtract(intervals: list[tuple[int, int]], cut_start: int, cut_end: int) -> list[tuple[int, int]]:
    out = []
    for start, end in intervals:
        if start < cut_start:
            out.append((start, min(end, cut_start)))
        if end > cut_end:
            out.append((max(start, cut_end), end))
    return [(start, end) for start, end in out if end > start]


def mrna_regions(exons: list[tuple[int, int]], cds: list[tuple[int, int]], strand: str) -> list[tuple[int, int, str]]:
    """Introns and UTRs of one transcript: gaps between consecutive exons, and exon parts outside the CDS span
    (the stop codon, which GENCODE annotates outside CDS, falls in the 3' UTR)."""
    exons, cds = sorted(set(exons)), sorted(set(cds))
    if not exons or not cds:
        return []
    cds_start, cds_end = min(s for s, _ in cds), max(e for _, e in cds)
    rows = [(left[1], right[0], "intron") for left, right in zip(exons, exons[1:]) if left[1] < right[0]]
    for start, end in _subtract(exons, cds_start, cds_end):
        upstream = end <= cds_start
        rows.append((start, end, ("utr5" if upstream else "utr3") if strand == "+" else ("utr3" if upstream else "utr5")))
    return rows


def read_annotation(gtf: str | Path, chromosomes: Iterable[str] = CHROMOSOMES) -> Annotation:
    """Parse a GENCODE GTF (plain or gzip) once into the scan targets (see `Annotation`)."""
    chromosomes = set(chromosomes)
    transcripts: dict[str, dict] = defaultdict(lambda: {"exons": [], "cds": [], "chrom": "", "strand": ""})
    genes: dict[str, list[tuple[int, int]]] = defaultdict(list)
    lncrna: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with _open_text(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] not in chromosomes:
                continue
            chrom, feature, strand = fields[0], fields[2], fields[6]
            interval = (int(fields[3]) - 1, int(fields[4]))  # GTF is 1-based inclusive
            if feature == "gene":
                genes[chrom].append(interval)
                match = _GENE_TYPE.search(fields[8])
                if match and match.group(1) == "lncRNA" and strand in ("+", "-"):
                    lncrna[(chrom, strand)].append(interval)
            elif feature in ("exon", "CDS"):
                attributes = dict(_ATTRIBUTE.findall(fields[8]))
                if attributes.get("gene_type") != "protein_coding" or attributes.get("transcript_type") != "protein_coding":
                    continue
                record = transcripts[attributes["transcript_id"]]
                record["chrom"], record["strand"] = chrom, strand
                record["exons" if feature == "exon" else "cds"].append(interval)
    pieces: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    for record in transcripts.values():
        for start, end, feature in mrna_regions(record["exons"], record["cds"], record["strand"]):
            pieces[(record["chrom"], record["strand"], feature)].append((start, end))
    rows = [(chrom, start, end, strand, feature) for (chrom, strand, feature), values in sorted(pieces.items())
            for start, end in merge_intervals(values)]
    regions = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "feature"])
    return Annotation(regions, {c: merge_intervals(v) for c, v in genes.items()},
                      {k: merge_intervals(v) for k, v in lncrna.items()})


# ----------------------------------------------------------------------------------------------- windows


def _centres(start: int, end: int, step: int) -> range:
    """Multiples of `step` in [start, end)."""
    return range(((start + step - 1) // step) * step, end, step)


def window_labels(chrom: str, chrom_length: int, annotation: Annotation, window: int = WINDOW, step: int = STEP,
                  intergenic: bool = False) -> dict[tuple[int, str], set[str]]:
    """Labels (intron, utr5, utr3, lncrna, intergenic) of every window start and strand on one chromosome."""
    half = window // 2
    labels: dict[tuple[int, str], set[str]] = defaultdict(set)

    def add(span_start: int, span_end: int, strand: str, label: str) -> None:
        for centre in _centres(span_start, span_end, step):
            start = centre - half
            if start >= 0 and start + window <= chrom_length:
                labels[(start, strand)].add(label)

    regions = annotation.regions[annotation.regions.chrom == chrom]
    for start, end, strand, feature in zip(regions.start, regions.end, regions.strand, regions.feature):
        add(int(start), int(end), str(strand), str(feature))
    for strand in "+-":
        for start, end in annotation.lncrna.get((chrom, strand), []):
            add(start, end, strand, "lncrna")
    if intergenic:  # windows lying entirely between merged gene spans, both strands
        previous_end = 0
        for gene_start, gene_end in [*annotation.genes.get(chrom, []), (chrom_length, chrom_length)]:
            for centre in _centres(previous_end + half, gene_start - half + 1, step):
                for strand in "+-":
                    labels[(centre - half, strand)].add("intergenic")
            previous_end = max(previous_end, gene_end)
    return labels


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def make_windows(chrom: str, sequence: str, annotation: Annotation, window: int = WINDOW, step: int = STEP,
                 intergenic: bool = False) -> pd.DataFrame:
    """Scan windows of one chromosome (DNA, 5'->3' on the window's strand), sorted by start and strand.

    Columns: window_id (`chrom:start-end:strand`), chrom, start, end, center, strand, features (labels joined by
    ';'), sequence. Windows with a base other than A/C/G/T are dropped; `frame.attrs["n_placed"]` counts them in.
    """
    sequence = sequence.upper()
    labels = window_labels(chrom, len(sequence), annotation, window, step, intergenic)
    rows = []
    for (start, strand), features in sorted(labels.items()):
        piece = sequence[start:start + window]
        if len(piece) != window or not set(piece) <= _ACGT:
            continue
        rows.append((f"{chrom}:{start}-{start + window}:{strand}", chrom, start, start + window, start + window // 2,
                     strand, ";".join(sorted(features)), reverse_complement(piece) if strand == "-" else piece))
    frame = pd.DataFrame(rows, columns=WINDOW_COLUMNS)
    frame.attrs["n_placed"] = len(labels)
    return frame


def is_intergenic(features: pd.Series) -> np.ndarray:
    """True for windows carrying the `intergenic` label."""
    return features.fillna("").str.split(";").apply(lambda labels: "intergenic" in labels).to_numpy(bool)


# ----------------------------------------------------------------------------------------------- repeat filters


def repeat_overlap(starts: np.ndarray, ends: np.ndarray, repeat_starts: np.ndarray,
                   repeat_ends: np.ndarray) -> np.ndarray:
    """Bases of each [start, end) covered by merged, sorted, non-overlapping repeat intervals."""
    starts, ends = np.asarray(starts, np.int64), np.asarray(ends, np.int64)
    if len(repeat_starts) == 0:
        return np.zeros(len(starts), dtype=np.int64)
    lengths = repeat_ends - repeat_starts
    before = np.concatenate([[0], np.cumsum(lengths)])  # bases covered by intervals 0..k-1

    def covered(x: np.ndarray) -> np.ndarray:  # covered bases in [0, x)
        k = np.searchsorted(repeat_starts, x, side="right") - 1
        inside = np.where(k >= 0, np.clip(x - repeat_starts[np.maximum(k, 0)], 0, lengths[np.maximum(k, 0)]), 0)
        return np.where(k >= 0, before[np.maximum(k, 0)], 0) + inside

    return covered(ends) - covered(starts)


def find_tool(given: str | None, names: Sequence[str]) -> str:
    """Path of an external program: `given` if set, otherwise the first of `names` on PATH."""
    for candidate in ([given] if given else list(names)):
        path = shutil.which(candidate)
        if path:
            return path
    raise FileNotFoundError(f"{given or names[0]} not found; install it or pass its path")


def _write_fasta(path: Path, sequences: Sequence[str], offset: int = 0) -> None:
    with path.open("w") as handle:
        for index, sequence in enumerate(sequences, start=offset):
            handle.write(f">r{index:09d}\n{sequence}\n")


def tantan_masked_bp(sequences: Sequence[str], tantan: str | None = None) -> np.ndarray:
    """Bases of each sequence masked as tandem repeats by `tantan -f 3` (union of the reported intervals)."""
    executable = find_tool(tantan, ["tantan"])
    sequences = list(sequences)
    masked: dict[int, list[tuple[int, int]]] = defaultdict(list)
    with tempfile.TemporaryDirectory(prefix="rlmhelix_tantan_") as tmp:
        fasta = Path(tmp) / "windows.fa"
        _write_fasta(fasta, sequences)
        result = subprocess.run([executable, "-f", "3", str(fasta)], check=True, capture_output=True, text=True)
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 3:
            masked[int(fields[0][1:])].append((int(fields[1]), int(fields[2])))
    counts = np.zeros(len(sequences), dtype=np.int64)
    for index, intervals in masked.items():
        counts[index] = sum(end - start for start, end in merge_intervals(intervals))
    return counts


def parse_irf_ngs(text: str) -> set[str]:
    """Names of sequences with at least one repeat in IRF `-ngs` output (`@name` lines followed by repeat lines)."""
    positive: set[str] = set()
    current = None
    for line in text.splitlines():
        if line.startswith("@"):
            current = line[1:].strip()
        elif current is not None and line[:1].isdigit():
            positive.add(current)
    return positive


def _irf_shard(executable: str, sequences: Sequence[str], offset: int, parameters: Sequence[str]) -> set[str]:
    with tempfile.TemporaryDirectory(prefix="rlmhelix_irf_") as tmp:
        _write_fasta(Path(tmp) / "windows.fa", sequences, offset=offset)
        result = subprocess.run([executable, "windows.fa", *parameters], cwd=tmp, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"IRF failed ({result.returncode}): {result.stderr[-2000:]}")
    return parse_irf_ngs(result.stdout)


def irf_positive(sequences: Sequence[str], irf: str | None = None, shard: int = IRF_SHARD,
                 parameters: Sequence[str] = IRF_PARAMETERS, threads: int = 1) -> np.ndarray:
    """True for sequences in which IRF reports an inverted repeat (score >= 40 with the default parameters).

    IRF runs on shards of `shard` sequences, `threads` at a time, each in a temporary directory (about 12 s per shard
    of 50,000 windows). On ARM (aarch64), build IRF with `make CFLAGS="-O2 -fsigned-char"`; without it the program
    is very slow.
    """
    from concurrent.futures import ThreadPoolExecutor

    executable = find_tool(irf, IRF_NAMES)
    sequences = list(sequences)
    hits = np.zeros(len(sequences), dtype=bool)
    begins = range(0, len(sequences), shard)
    with ThreadPoolExecutor(max(1, threads)) as pool:
        for names in pool.map(lambda b: _irf_shard(executable, sequences[b:b + shard], b, parameters), begins):
            hits[[int(name[1:]) for name in names]] = True
    return hits


def filter_repeats(windows: pd.DataFrame, repeats: tuple[np.ndarray, np.ndarray] | None = None,
                   tantan: str | None = None, irf: str | None = None, use_tantan: bool = True,
                   use_irf: bool = True, max_repeatmasker_bp: int = RMSK_MAX_BP,
                   min_tantan_bp: int = TANTAN_MIN_BP, threads: int = 1) -> tuple[pd.DataFrame, dict[str, int]]:
    """Remove repeat windows: RepeatMasker overlap > 40 nt, then tantan >= 20 nt, then an IRF inverted repeat.

    `repeats` is the (starts, ends) pair of merged repeats of the windows' chromosome (`read_repeats`). The first two
    tests add the columns repeatmasker_bp and tantan_bp. A test is skipped when its input is missing (`repeats=None`,
    `use_tantan=False`, `use_irf=False`), which does not reproduce the paper. Returns the kept windows and the counts
    (n_windows, then the windows removed by each test: n_repeatmasker, n_tantan, n_irf, and n_kept).
    """
    counts = {"n_windows": len(windows)}
    frame = windows.reset_index(drop=True)
    if repeats is not None:
        frame["repeatmasker_bp"] = repeat_overlap(frame.start, frame.end, *repeats)
        keep = frame.repeatmasker_bp.to_numpy() <= max_repeatmasker_bp
        counts["n_repeatmasker"] = int((~keep).sum())
        frame = frame[keep].reset_index(drop=True)
    if use_tantan and len(frame):
        frame["tantan_bp"] = tantan_masked_bp(list(frame.sequence), tantan)
        keep = frame.tantan_bp.to_numpy() < min_tantan_bp
        counts["n_tantan"] = int((~keep).sum())
        frame = frame[keep].reset_index(drop=True)
    if use_irf and len(frame):
        positive = irf_positive(list(frame.sequence), irf, threads=threads)
        counts["n_irf"] = int(positive.sum())
        frame = frame[~positive].reset_index(drop=True)
    counts["n_kept"] = len(frame)
    return frame, counts


# ----------------------------------------------------------------------------------------------- scoring


def _score_window(sequence: str, apc: np.ndarray) -> dict[str, float]:
    from .score import rlmhelix_score

    result = rlmhelix_score(apc, sequence)
    return {"score": result.score, "n_pairs": result.n_pairs, "n_stems": result.n_stems,
            "n_stems_scored": result.n_stems_scored, "longest_stem": result.longest_stem}


def iter_scores(windows: pd.DataFrame, student, calibrator=None, batch_size: int = 16, n_jobs: int = 4,
                chunk: int = 10_000) -> Iterator[pd.DataFrame]:
    """Score windows `chunk` at a time; yields the windows of each chunk with score columns and the P value.

    The last batch of a chunk is filled up to `batch_size` with copies of its last window (whose extra maps are
    discarded), so that every window is scored in a batch of the same shape: under bf16 the result of a window
    depends on the batch size but not on the other windows in the batch.
    """
    from joblib import Parallel, delayed

    from .calibrate import Calibrator

    calibrator = calibrator or Calibrator.load()
    for begin in range(0, len(windows), chunk):
        part = windows.iloc[begin:begin + chunk].reset_index(drop=True)
        sequences = list(part.sequence)
        padding = [sequences[-1]] * (-len(sequences) % batch_size)
        maps = student.predict(sequences + padding, batch_size)[:len(sequences)]
        rows = Parallel(n_jobs=n_jobs, batch_size=4)(delayed(_score_window)(s, m) for s, m in zip(sequences, maps))
        scored = pd.concat([part, pd.DataFrame(rows)], axis=1)
        scored["GC"] = [(s.count("G") + s.count("C")) / len(s) for s in sequences]
        scored["p_value"] = calibrator.pvalues(scored.score, sequences)
        yield scored


def score_windows(windows: pd.DataFrame, student, calibrator=None, batch_size: int = 16, n_jobs: int = 4,
                  chunk: int = 10_000) -> pd.DataFrame:
    """RLMHelix score, stem counts, GC and calibrated P value of every window (`student` is a `Student`).

    The paper scored 200-nt windows in batches of 16 on a GPU (bf16); other batch sizes change the bf16 rounding
    (about 0.1% of the score).
    """
    parts = list(iter_scores(windows, student, calibrator, batch_size, n_jobs, chunk))
    return pd.concat(parts, ignore_index=True) if parts else windows.iloc[:0].copy()


# ----------------------------------------------------------------------------------------------- multiple testing


def benjamini_hochberg(pvalues) -> np.ndarray:
    """Benjamini-Hochberg adjusted P values (q values)."""
    p = np.asarray(pvalues, dtype=np.float64)
    n = len(p)
    if n == 0:
        return p.copy()
    order = np.argsort(p, kind="mergesort")
    adjusted = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(adjusted, 0, 1)
    return q


def merge_loci(windows: pd.DataFrame, score: str = "score", p: str = "p_value", q: str = "q_value") -> pd.DataFrame:
    """Merge windows that overlap on the same chromosome and strand into loci (touching windows stay apart).

    One row per locus: locus id, span, number of windows, the peak window (highest score, then lowest P) with its
    score, the smallest P and q, the union of the window labels, and the GC, stem counts and sequence of the peak
    window where present. Sorted by q, then by decreasing peak score.
    """
    extra = [c for c in ("GC", "n_stems_scored", "longest_stem", "sequence") if c in windows.columns]
    rows = []
    for (chrom, strand), group in windows.groupby(["chrom", "strand"], sort=False):
        group = group.sort_values(["start", "end"], kind="mergesort")
        starts, ends = group.start.to_numpy(), group.end.to_numpy()
        locus_index = np.zeros(len(group), dtype=int)
        current_end = -1
        for k in range(len(group)):
            new = k > 0 and starts[k] >= current_end
            locus_index[k] = locus_index[k - 1] + new if k else 0
            current_end = ends[k] if new else max(current_end, ends[k])
        for _, frame in group.groupby(locus_index, sort=True):
            peak = frame.sort_values([score, p], ascending=[False, True], kind="mergesort").iloc[0]
            start, end = int(frame.start.min()), int(frame.end.max())
            rows.append({"locus": f"{chrom}:{start}-{end}:{strand}", "chrom": chrom, "start": start, "end": end,
                         "strand": strand, "n_windows": len(frame), "peak_window": peak.window_id,
                         "peak_score": peak[score], "min_p": frame[p].min(), "min_q": frame[q].min(),
                         "features": ";".join(sorted(set(";".join(frame.features).split(";")))),
                         **{c: peak[c] for c in extra}})
    if not rows:
        return pd.DataFrame(columns=["locus", "chrom", "start", "end", "strand", "n_windows", "peak_window",
                                     "peak_score", "min_p", "min_q", "features", *extra])
    return pd.DataFrame(rows).sort_values(["min_q", "peak_score"], ascending=[True, False]).reset_index(drop=True)


def z_scores(pvalues) -> np.ndarray:
    """z = -Phi^-1(P), with P clipped to [1e-300, 1 - 1e-12]; large z means more structure than expected."""
    from scipy import stats

    return -stats.norm.ppf(np.clip(np.asarray(pvalues, dtype=np.float64), 1e-300, 1 - 1e-12))


@dataclass
class EmpiricalNull:
    """Efron's empirical null N(mu0, sigma0^2) of z, with null proportion pi0 (central matching)."""

    mu0: float
    sigma0: float
    pi0: float
    n: int
    central_low: float
    central_high: float

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame([asdict(self)])


def empirical_null(z, central: tuple[float, float] = (0.25, 0.75), bins: int = 60) -> EmpiricalNull:
    """Central matching: a quadratic fitted to the log counts of `bins` histogram bins between the given quantiles
    of z gives the normal null; pi0 is the observed fraction in that range over the null mass in it (capped at 1)."""
    from scipy import stats

    z = np.asarray(z, dtype=np.float64)
    low, high = np.quantile(z, central)
    edges = np.linspace(low, high, bins + 1)
    counts, _ = np.histogram(z, edges)
    middle = (edges[:-1] + edges[1:]) / 2
    used = counts > 0
    b2, b1, _ = np.polyfit(middle[used], np.log(counts[used]), 2)
    if b2 >= 0:
        raise ValueError("central matching failed: the log counts are not concave")
    sigma0 = float(np.sqrt(-1 / (2 * b2)))
    mu0 = float(b1 * sigma0 ** 2)
    null_mass = stats.norm.cdf(high, mu0, sigma0) - stats.norm.cdf(low, mu0, sigma0)
    pi0 = float(min(1.0, (counts.sum() / len(z)) / null_mass))
    return EmpiricalNull(mu0, sigma0, pi0, len(z), float(low), float(high))


def _null_parameters(null) -> tuple[float, float, float]:
    """(mu0, sigma0, pi0) of an `EmpiricalNull` or of such a tuple."""
    return (null.mu0, null.sigma0, null.pi0) if isinstance(null, EmpiricalNull) else tuple(map(float, null))


def empirical_pvalues(z, null: EmpiricalNull | tuple[float, float, float]) -> np.ndarray:
    """Upper-tail P values of z under the empirical null."""
    from scipy import stats

    mu0, sigma0, _ = _null_parameters(null)
    return stats.norm.sf((np.asarray(z, dtype=np.float64) - mu0) / sigma0)


def tail_fdr(z, null: EmpiricalNull, top_fractions: Sequence[float] = TOP_FRACTIONS,
             sets: dict[str, int] | None = None, lfdr=None) -> pd.DataFrame:
    """Efron's tail-area Fdr of the top windows: Fdr(t) = pi0 S0(t) / S(t), t the z of the last window in the set.

    Rows for each top fraction (ceil(fraction * N) windows) and for named sets of the top n windows (`sets`, e.g.
    {"BH q<=0.05": 565}). Columns: z and P thresholds, Fdr under the theoretical N(0, 1) and the empirical null,
    expected null windows, estimated signal windows R (1 - Fdr), and with `lfdr` the lfdr at the boundary and the
    mean lfdr of the set (which estimates the same Fdr).
    """
    from scipy import stats

    z = np.asarray(z, dtype=np.float64)
    n = len(z)
    order = np.argsort(-z, kind="mergesort")
    descending = z[order]
    named = [(f"top {100 * f:g}%", f, int(np.ceil(f * n))) for f in top_fractions]
    named += [(label, np.nan, int(r)) for label, r in (sets or {}).items() if r > 0]
    rows = []
    for label, fraction, r in named:
        t = descending[r - 1]
        null_tail = null.pi0 * stats.norm.sf(t, null.mu0, null.sigma0)
        fdr = min(1.0, null_tail * n / r)
        row = {"set": label, "top_fraction": fraction, "n_windows": r, "z_threshold": t,
               "p_threshold": stats.norm.sf(t), "fdr_theoretical_null": min(1.0, stats.norm.sf(t) * n / r),
               "fdr_empirical_null": fdr, "expected_null_windows": null_tail * n, "signal_windows": r * (1 - fdr)}
        if lfdr is not None:
            lfdr = np.asarray(lfdr, dtype=np.float64)
            row["lfdr_at_boundary"] = lfdr[order[r - 1]]
            row["mean_lfdr"] = float(lfdr[z >= t].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def _natural_spline_basis(x: np.ndarray, knots: np.ndarray) -> np.ndarray:
    """Natural cubic spline basis (Hastie et al., ESL eq. 5.4-5.5): 1, x and K - 2 truncated-power terms."""
    k = np.asarray(knots, dtype=float)

    def d(j: int) -> np.ndarray:
        return (np.clip(x - k[j], 0, None) ** 3 - np.clip(x - k[-1], 0, None) ** 3) / (k[-1] - k[j])

    return np.column_stack([np.ones_like(x), x] + [d(j) - d(len(k) - 2) for j in range(len(k) - 2)])


def _poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    return float(2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0) - (y - mu)))


def _poisson_irls(x: np.ndarray, y: np.ndarray, n_iter: int = 100) -> np.ndarray:
    """Poisson regression (log link) by iteratively reweighted least squares."""
    beta = np.linalg.lstsq(x, np.log(y + 0.5), rcond=None)[0]
    previous = np.inf
    for _ in range(n_iter):
        eta = np.clip(x @ beta, -50, 50)
        mu = np.exp(eta)
        weighted = x.T * mu
        beta = np.linalg.solve(weighted @ x + 1e-10 * np.eye(x.shape[1]), weighted @ (eta + (y - mu) / mu))
        deviance = _poisson_deviance(y, np.exp(np.clip(x @ beta, -50, 50)))
        if abs(previous - deviance) < 1e-8 * max(1.0, deviance):
            break
        previous = deviance
    return beta


def lindsey_density(z, bin_width: float = 0.05, knot_quantiles: Sequence[float] = LINDSEY_KNOT_QUANTILES):
    """Mixture density f of z by Lindsey's method: histogram counts (bins of `bin_width` over the 1e-6 to 1 - 1e-6
    quantile range, rounded outwards to 0.05) fitted by Poisson regression on a natural cubic spline with knots at
    the range ends and the given quantiles of z. Returns a function of z."""
    z = np.asarray(z, dtype=np.float64)
    low, high = np.quantile(z, [1e-6, 1 - 1e-6])
    low, high = np.floor(low * 20) / 20, np.ceil(high * 20) / 20
    knots = (np.r_[low, np.quantile(z, list(knot_quantiles)), high] - low) / (high - low)
    edges = np.arange(low, high + bin_width, bin_width)
    counts, _ = np.histogram(z, edges)
    middle = (edges[:-1] + edges[1:]) / 2
    beta = _poisson_irls(_natural_spline_basis((middle - low) / (high - low), knots), counts.astype(float))
    n = len(z)

    def density(values) -> np.ndarray:
        x = _natural_spline_basis((np.asarray(values, dtype=float) - low) / (high - low), knots)
        return np.exp(np.clip(x @ beta, -700, 700)) / (n * bin_width)

    return density


def local_fdr(z, null: EmpiricalNull | tuple[float, float, float], bin_width: float = 0.05,
              knot_quantiles: Sequence[float] = LINDSEY_KNOT_QUANTILES, one_sided: bool = True) -> np.ndarray:
    """lfdr(z) = pi0 f0(z) / f(z) per window, f from `lindsey_density`; capped at 1 and made non-increasing in z to
    the right of mu0 (cumulative minimum on a 0.001 grid, interpolated to each z).

    `z` holds all eligible windows (f is fitted to them); `null` is an `EmpiricalNull` or (mu0, sigma0, pi0), e.g.
    the paper's (0.6729420119015076, 0.5481601926577551, 0.950765247091798).

    With `one_sided` (default), windows with z < mu0 (less structure than the null centre) get 1: the left tail of z
    is also heavier than normal, so pi0 f0 / f falls there too, but those windows are not discoveries. The paper used
    the lfdr of the top 2% of windows only, which this does not change.
    """
    from scipy import stats

    mu0, sigma0, pi0 = _null_parameters(null)
    z = np.asarray(z, dtype=np.float64)
    density = lindsey_density(z, bin_width, knot_quantiles)
    grid = np.arange(z.min() - 0.001, z.max() + 0.002, 0.001)
    values = np.minimum(1.0, pi0 * stats.norm.pdf(grid, mu0, sigma0) / density(grid))
    right = grid >= mu0
    values[right] = np.minimum.accumulate(values[right])
    lfdr = np.interp(z, grid, values)
    if one_sided:
        lfdr[z < mu0] = 1.0
    return lfdr


# ----------------------------------------------------------------------------------------------- pipeline steps


def _log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def run_prepare(gtf: str | Path, genome: str | Path, rmsk: str | Path, outdir: str | Path,
                chromosomes: Sequence[str] = CHROMOSOMES, tantan: str | None = None, irf: str | None = None,
                use_tantan: bool = True, use_irf: bool = True, intergenic: bool = False,
                resume: bool = False, threads: int = 1) -> pd.DataFrame:
    """Windows and repeat filters for each chromosome: writes `<outdir>/<chrom>.parquet` (kept windows) and
    `<outdir>/window_counts.tsv`; returns the counts."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pending = [c for c in chromosomes if not (resume and (outdir / f"{c}.parquet").exists())]
    if use_tantan:
        find_tool(tantan, ["tantan"])
    if use_irf:
        find_tool(irf, IRF_NAMES)
    counts_path = outdir / "window_counts.tsv"
    previous = read_table(counts_path) if resume and counts_path.exists() else pd.DataFrame()
    rows = [r for r in previous.to_dict("records") if r["chrom"] not in pending] if len(previous) else []
    if pending:
        _log("reading annotation and repeats")
        annotation = read_annotation(gtf, pending)
        repeats = read_repeats(rmsk, pending)
        for chrom, sequence in iter_chromosomes(genome, pending):
            windows = make_windows(chrom, sequence, annotation, intergenic=intergenic)
            empty = (np.zeros(0, np.int64), np.zeros(0, np.int64))
            kept, counts = filter_repeats(windows, repeats.get(chrom, empty), tantan, irf, use_tantan, use_irf,
                                          threads=threads)
            rows.append({"chrom": chrom, "n_placed": windows.attrs["n_placed"],
                         "n_non_acgt": windows.attrs["n_placed"] - len(windows), **counts})
            write_table(kept, outdir / f"{chrom}.parquet")
            _log(" ".join(f"{k}={v}" for k, v in rows[-1].items()))
    order = {c: k for k, c in enumerate(chromosomes)}
    table = pd.DataFrame(rows).sort_values("chrom", key=lambda s: s.map(order))
    write_table(table, counts_path)
    return table


def run_score(windows_dir: str | Path, outdir: str | Path, chromosomes: Sequence[str] = CHROMOSOMES,
              device: str | None = None, dtype: str = "bf16", batch_size: int = 16, n_jobs: int = 4,
              chunk: int = 10_000, resume: bool = False) -> None:
    """Score the windows of each chromosome: `<outdir>/<chrom>.parquet`, written in shards of `chunk` windows so
    that an interrupted run continues with `resume=True`."""
    from .calibrate import Calibrator
    from .student import Student

    outdir = Path(outdir)
    calibrator = Calibrator.load()
    student = None
    for chrom in chromosomes:
        source, destination = Path(windows_dir) / f"{chrom}.parquet", outdir / f"{chrom}.parquet"
        if resume and destination.exists():
            continue
        windows = read_table(source)
        student = student or Student(device=device, dtype=dtype)
        shard_dir = outdir / "shards" / chrom
        shard_dir.mkdir(parents=True, exist_ok=True)
        for index, begin in enumerate(range(0, len(windows), chunk)):
            shard = shard_dir / f"shard_{index:05d}.parquet"
            if resume and shard.exists():
                continue
            part = windows.iloc[begin:begin + chunk]
            write_table(score_windows(part, student, calibrator, batch_size, n_jobs, chunk), shard)
            _log(f"{chrom} {min(begin + chunk, len(windows)):,}/{len(windows):,}")
        shards = [shard_dir / f"shard_{k:05d}.parquet" for k in range(-(-len(windows) // chunk))]
        scored = pd.concat([read_table(s) for s in shards], ignore_index=True) if shards else windows
        if list(scored.window_id) != list(windows.window_id):
            raise ValueError(f"{chrom}: shards do not match the windows; rerun without --resume")
        write_table(scored, destination)
        shutil.rmtree(shard_dir)
    if (outdir / "shards").is_dir() and not any((outdir / "shards").iterdir()):
        (outdir / "shards").rmdir()


def run_fdr(score_files: Sequence[str | Path], outdir: str | Path, level: float = BH_LEVEL,
            top_fractions: Sequence[float] = TOP_FRACTIONS) -> dict[str, pd.DataFrame]:
    """Benjamini-Hochberg, loci, empirical null and lfdr over all scored windows (intergenic windows excluded).

    Writes to `outdir`: eligible_windows.parquet (every eligible window with all score columns, q_value, z, lfdr and
    bh_accepted = q_value <= level, in the order of `score_files`), bh_windows.tsv and bh_loci.tsv (accepted windows
    and their loci), empirical_null.tsv, tail_fdr.tsv and summary.tsv. Returns the tables, with the eligible windows
    without sequences.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    columns = ["window_id", "chrom", "start", "end", "strand", "features", "score", "p_value"]
    parts = [read_table(path, columns) for path in score_files]
    masks = [~is_intergenic(part.features) for part in parts]
    windows = pd.concat([part[mask] for part, mask in zip(parts, masks)], ignore_index=True).drop(columns="features")
    windows["q_value"] = benjamini_hochberg(windows.p_value)
    windows["z"] = z_scores(windows.p_value)
    null = empirical_null(windows.z)
    windows["lfdr"] = local_fdr(windows.z, null)
    windows["bh_accepted"] = windows.q_value.to_numpy() <= level

    destination = outdir / "eligible_windows.parquet"  # written one chromosome at a time
    temporary = destination.with_name(destination.name + ".tmp")
    writer, offset, bh_parts = None, 0, []
    for path, mask in zip(score_files, masks):
        part = read_table(path)[mask].reset_index(drop=True)
        added = windows.iloc[offset:offset + len(part)]
        offset += len(part)
        if list(added.window_id) != list(part.window_id):
            raise ValueError(f"{path} changed while running")
        part = part.assign(**{c: added[c].to_numpy() for c in ("q_value", "z", "lfdr", "bh_accepted")})
        bh_parts.append(part[part.bh_accepted])
        table = pa.Table.from_pandas(part, preserve_index=False)
        writer = writer or pq.ParquetWriter(temporary, table.schema)
        writer.write_table(table.cast(writer.schema))
    if writer is not None:
        writer.close()
        temporary.replace(destination)
    bh = pd.concat(bh_parts, ignore_index=True)
    loci = merge_loci(bh)
    tails = tail_fdr(windows.z, null, top_fractions, {f"BH q<={level:g}": len(bh)}, windows.lfdr)
    summary = pd.DataFrame([{
        "n_windows": len(windows),
        "n_distinct_coordinates": len(windows[["chrom", "start", "end"]].drop_duplicates()),
        "n_p_le_0.05": int((windows.p_value <= 0.05).sum()), "n_p_le_0.01": int((windows.p_value <= 0.01).sum()),
        "n_bh_windows": len(bh), "n_bh_loci": len(loci), "n_lfdr_le_0.2": int((windows.lfdr <= 0.2).sum()),
        "signal_windows_top_1pct": float(tails.loc[tails.top_fraction == 0.01, "signal_windows"].sum())}])
    for name, frame in [("bh_windows", bh), ("bh_loci", loci), ("empirical_null", null.to_frame()),
                        ("tail_fdr", tails), ("summary", summary)]:
        write_table(frame, outdir / f"{name}.tsv")
    return {"windows": windows, "bh_windows": bh, "bh_loci": loci, "empirical_null": null.to_frame(),
            "tail_fdr": tails, "summary": summary}


# ----------------------------------------------------------------------------------------------- command line


def _chromosome_files(directory: str | Path, chromosomes: Sequence[str]) -> list[Path]:
    paths = [Path(directory) / f"{c}.parquet" for c in chromosomes]
    missing = [p for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} chromosome files missing, first {missing[0]}")
    return paths


def run(args) -> None:
    """Run the `rlmhelix scan` step selected on the command line."""
    if args.scan_step == "prepare":
        run_prepare(args.gtf, args.genome, args.rmsk, args.outdir, args.chrom, args.tantan, args.irf,
                    not args.skip_tantan, not args.skip_irf, args.keep_intergenic, args.resume, args.threads)
    elif args.scan_step == "score":
        run_score(args.windows, args.outdir, args.chrom, args.device, args.dtype, args.batch_size, args.n_jobs,
              args.chunk, args.resume)
    else:
        result = run_fdr(_chromosome_files(args.scores, args.chrom), args.outdir, args.q)
        for name, value in result["summary"].iloc[0].items():
            print(f"{name}\t{value:,.0f}")


def register(subparsers) -> None:
    """Add `rlmhelix scan prepare | score | fdr` to the subparsers of the main command."""
    scan = subparsers.add_parser("scan", help="genome scan: windows, repeat filters, scores, FDR",
                                 description="Genome scan of the paper (see docs/pipeline_scan.md).")
    steps = scan.add_subparsers(dest="scan_step", required=True)
    chrom = {"nargs": "+", "default": list(CHROMOSOMES), "metavar": "CHROM",
             "help": "chromosomes (default: chr1-chr22, chrX, chrY)"}

    step = steps.add_parser("prepare", help="place windows and remove repeats (CPU, tantan and IRF)")
    step.add_argument("--gtf", required=True, help="GENCODE annotation GTF (plain or .gz)")
    step.add_argument("--genome", required=True, help="genome FASTA (plain with .fai for random access, or .gz)")
    step.add_argument("--rmsk", required=True, help="UCSC rmsk.txt(.gz), or a BED of repeats")
    step.add_argument("--tantan", help="tantan executable (default: tantan on PATH)")
    step.add_argument("--irf", help="IRF executable (default: irf, irf.exe or irf309.linux.exe on PATH)")
    step.add_argument("--skip-tantan", action="store_true", help="no tandem-repeat filter (not the paper's setting)")
    step.add_argument("--skip-irf", action="store_true", help="no inverted-repeat filter (not the paper's setting)")
    step.add_argument("--keep-intergenic", action="store_true",
                      help="also make intergenic windows (scored in the paper's run, excluded before testing)")
    step.add_argument("--chrom", **chrom)
    step.add_argument("--threads", type=int, default=1, help="IRF processes run in parallel")
    step.add_argument("--resume", action="store_true", help="skip chromosomes already written")
    step.add_argument("-o", "--outdir", required=True, help="directory for <chrom>.parquet and window_counts.tsv")

    step = steps.add_parser("score", help="score windows with the student and calibrate (GPU recommended)")
    step.add_argument("--windows", required=True, help="output directory of `scan prepare`")
    step.add_argument("--chrom", **chrom)
    step.add_argument("--device", help="torch device (default: cuda if available, else cpu)")
    step.add_argument("--dtype", choices=("bf16", "fp16", "fp32"), default="bf16", help="autocast type on GPU")
    step.add_argument("--batch-size", type=int, default=16, help="windows per language-model batch (paper: 16)")
    step.add_argument("--n-jobs", type=int, default=4, help="CPU workers for RNAfold and the DP")
    step.add_argument("--chunk", type=int, default=10_000, help="windows per shard")
    step.add_argument("--resume", action="store_true", help="keep finished chromosomes and shards")
    step.add_argument("-o", "--outdir", required=True, help="directory for scored <chrom>.parquet")

    step = steps.add_parser("fdr", help="Benjamini-Hochberg, loci, empirical null and local FDR (CPU)")
    step.add_argument("--scores", required=True, help="output directory of `scan score`")
    step.add_argument("--chrom", **chrom)
    step.add_argument("--q", type=float, default=BH_LEVEL, help="BH level for the window set and loci")
    step.add_argument("-o", "--outdir", required=True,
                      help="directory for eligible_windows.parquet and the BH, loci and FDR tables")
    scan.set_defaults(func=run)
