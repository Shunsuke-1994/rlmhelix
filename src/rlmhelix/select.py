"""Candidate selection.

Stage 1, structural cores. The top 1% of the eligible scan windows by calibrated P value (a stable sort, so ties keep
the order of the scan table) are merged into loci: windows on the same strand that overlap by at least 1 nt. The peak
window of each locus (smallest P value, then largest score) is scored again with the student. In its DP decoding,
isolated pairs and stems of fewer than four pairs are removed, and every outermost base-paired interval, over which
the bracket depth rises from and returns to 0, becomes a core (as in EvoFold), scored by the total weight
sqrt(d * p) of its pairs. The core with the largest score is the primary core of the locus. Primary cores whose
(sequence, dot-bracket) pairs are identical after lexicographic normalization with the reverse complement and the
reversed dot-bracket are collapsed to one representative, the one with the smallest peak P value and then the largest
core score. Representatives keep both orientations: the normalized one, which was searched for the selection
routes in the paper, and the source strand of the scan, which is used for the alignments.

Stage 2, homology search and routes. Genomes are cut into 10-kb chunks overlapping by 500 nt and searched with
Riboseek (sensitivity 7.5, E <= 10, at most 100 target sequences). Hits with E <= 1e-3 and query coverage >= 0.8
are high-confidence. For each core, the species of the 20-species seed panel with a high-confidence hit are counted
per lineage group, and cores meeting any of four routes are selected; each is labelled by the first route it meets:

    vertebrate_wide     >= 1 non-mammalian tetrapod and >= 1 other vertebrate
    mammal_conserved    >= 3 non-primate mammals and >= 3 primates
    primate_focused     >= 6 primates and peak P <= 0.001
    high_score_rescue   peak P <= 0.0005

Alignments. The selected cores are searched again in their source orientation. Hits with E <= 0.1 from overlapping
chunks are merged into loci per query, species, contig and strand; a locus is high-confidence when one of its hits
is. The high-confidence loci of the seed panel, projected onto the human core through the pairwise alignments, form
the seed alignment with the core structure as SS_cons and the human columns as RF; a CM is built (cmbuild --hand),
the sequences are realigned (cmalign --cyk --notrunc --noprob) and the CM is rebuilt. Each high-confidence locus is
then cropped again with one locus length of flank on either side, its boundaries are set by a glocal CYK search with
the calibrated CM, and the sequences are realigned to give the seed CM. Cores without a high-confidence locus get no
seed CM and are dropped. In the 35 expansion genomes (34 vertebrates and the human genome, for paralogous loci),
regions of three times the core length around loci with E <= 0.1 are searched with the seed CM, and hits are
accepted with a bit score of at least the 10th percentile of the seed hits minus 10 bits and a model coverage of at
least 0.8. Accepted sequences not identical to a seed row are added with cmalign --mapali (final55), and the CM is
rebuilt.

Structure-free homologues. The human core, the distinct seed sequences (cut at the Riboseek alignments) and the
expansion loci whose Riboseek hit has E <= 1e-3 and query coverage >= 0.8 (cut at the merged locus) form the input of
covariation condition (iii), which is aligned without structural information.

Provisional families. Each seed CM is searched in HMM-only mode against the human cores of all candidates
(E <= 10) to propose pairs. For each pair, each CM scores up to four held-out expansion sequences of the other
candidate (one per species, chosen round-robin over lineage groups for 4-mer diversity). With b_{A->B} the median
score of A's CM on B's held-out sequences and T_A the inclusion threshold of A,
r_{A->B} = (b_{A->B} - T_A) / T_A, and pairs whose mean of r_{A->B} and r_{B->A} is at least -0.525 are joined into
connected components. Candidates with more than 100 high-confidence or final55 loci are flagged as high-copy.

The pure-Python steps (windows, loci, cores, normalization, hit parsing and merging, routes, crops, inclusion rule,
families) need no external tools; Riboseek and Infernal are located with `shutil.which` or given as options.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np
import pandas as pd

from .score import decode_pairs, dot_bracket, stems

# Stage 1
TOP_FRACTION = 0.01
MIN_STEM = 4
ID_PREFIX = "pvalue_top1_primary_"   # representative IDs of the paper (query_id in the catalogue)
P_COLUMNS = ("p_value", "p_min4_mixed")
SCORE_COLUMNS = ("score_min4", "score")
# carried from the peak window; n_stems_scored is the scan's name of the paper's n_stems_ge4
LOCUS_CARRY = ("features", "sequence", "GC", "n_stems_scored", "n_stems_ge4", "longest_stem", "repeatmasker_bp",
               "tantan_bp")

CORE_DTYPES = {"locus": "str", "peak_window": "str", "core_index": "int64", "chrom": "str", "strand": "str",
               "features": "object", "window_score": "float64", "core_start_in_window": "int64",
               "core_end_in_window": "int64", "core_length": "int64", "core_genome_start": "int64",
               "core_genome_end": "int64", "core_score": "float64", "n_pairs": "int64", "n_stems": "int64",
               "longest_stem": "int64", "touches_window_edge": "bool", "core_sequence": "str",
               "core_dot_bracket": "str"}

# Stage 2
CHUNK_SIZE, CHUNK_OVERLAP = 10_000, 500
HIT_FORMAT = "query,target,evalue,bits,qstart,qend,tstart,tend,alnlen,qcov,tcov,qaln,taln"
HIT_COLUMNS = ["query_id", "target_id", "evalue", "bits", "qstart", "qend", "tstart", "tend", "alnlen", "qcov",
               "tcov", "qaln", "taln"]
HIT_DTYPES = {"query_id": "string", "target_id": "string", "evalue": "float64", "bits": "float64", "qstart": "int64",
              "qend": "int64", "tstart": "int64", "tend": "int64", "alnlen": "int64", "qcov": "float64",
              "tcov": "float64", "qaln": "string", "taln": "string"}
RIBOSEEK_SEARCH = ("-a", "-s", "7.5", "-e", "10", "--max-seqs", "100", "--alt-ali", "0", "--prefilter-mode", "0",
                   "--gpu", "0")
HC_EVALUE, HC_QCOV = 1e-3, 0.8
LOCUS_EVALUE = 0.1
LINEAGES = ("primate", "nonprimate_mammal", "nonmammal_tetrapod", "other_vertebrate")
ROUTES = ("vertebrate_wide", "mammal_conserved", "primate_focused", "high_score_rescue")
PRIMATE_FOCUSED_P, RESCUE_P = 1e-3, 5e-4

# Alignments
INCLUSION_QUANTILE, INCLUSION_MARGIN, INCLUSION_COVERAGE = 0.10, -10.0, 0.8
HIGH_COPY_LOCI = 100

# Families
HMM_EVALUE = 10.0
HELDOUT_PER_CANDIDATE = 4
FAMILY_THRESHOLD = -0.525
EDGE_DTYPES = {"unit_a": "str", "unit_b": "str", "b_ab": "float64", "b_ba": "float64", "r_ab": "float64",
               "r_ba": "float64", "ratio_mean": "float64"}

# the 55-species panel, shipped with the package (a copy of data/genome_panel_55.tsv of the repository)
DEFAULT_PANEL = Path(str(resources.files("rlmhelix").joinpath("data/genome_panel_55.tsv")))
INFERNAL_PROGRAMS = ("cmbuild", "cmalign", "cmcalibrate", "cmsearch")
_COMPLEMENT = str.maketrans("ACGTUNacgtun-.", "TGCAANtgcaan-.")
_RNA_COMPLEMENT = str.maketrans("ACGTU", "UGCAA")


# ---------------------------------------------------------------------------------------------------------------
# Small helpers


def read_table(path: str | Path) -> pd.DataFrame:
    """Read a parquet, TSV or CSV table (by suffix)."""
    path = Path(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="," if path.suffix == ".csv" else "\t")


def write_table(frame: pd.DataFrame, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".parquet":
        frame.to_parquet(path, index=False)
    else:
        frame.to_csv(path, sep="," if path.suffix == ".csv" else "\t", index=False)


def write_fasta(records: Iterable[tuple[str, str]], path: str | Path, dna: bool = True) -> None:
    """Write (name, sequence) records, as DNA (U -> T) by default."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence.replace('U', 'T') if dna else sequence}\n")


def iter_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    """(first word of the header, sequence) for each record of a plain or gzipped FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                words = line[1:].split()
                name, chunks = (words[0] if words else ""), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def reverse_complement(sequence: str) -> str:
    """DNA reverse complement (keeps N and gap symbols)."""
    return sequence.translate(_COMPLEMENT)[::-1]


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", value)


def _column(frame: pd.DataFrame, explicit: str | None, choices: Sequence[str], what: str) -> str:
    if explicit:
        return explicit
    for choice in choices:
        if choice in frame:
            return choice
    raise KeyError(f"no {what} column (tried {', '.join(choices)}); pass it explicitly")


def _empty_frame(dtypes: Mapping[str, str]) -> pd.DataFrame:
    """An empty table with typed columns (so that concatenation and aggregation keep the dtypes)."""
    return pd.DataFrame({column: pd.Series(dtype=dtype) for column, dtype in dtypes.items()})


def file_sha256(path: str | Path) -> str:
    """SHA-256 of a file's content (hex)."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tool(name: str, directory: str | Path | None = None) -> str:
    """Path of an external program: `directory/name` when a directory is given, otherwise from PATH."""
    if directory is not None:
        path = Path(directory) / name
        if not path.is_file():
            raise FileNotFoundError(f"{name} not found in {directory}")
        return str(path)
    found = shutil.which(name)
    if found is None:
        raise FileNotFoundError(f"{name} not found on PATH; install it or pass its location")
    return found


def require_tools(names: Iterable[str], directory: str | Path | None = None) -> dict[str, str]:
    """Look up external programs before any task starts; exit with the name of the first one that is missing."""
    try:
        return {name: tool(name, directory) for name in names}
    except FileNotFoundError as exc:
        raise SystemExit(str(exc)) from None


def run_logged(command: Sequence[str], log: str | Path | None = None) -> None:
    """Run a command, writing its output to `log`; raise on failure."""
    if log is None:
        subprocess.run(list(command), check=True, stdout=subprocess.DEVNULL)
        return
    with Path(log).open("w") as handle:
        completed = subprocess.run(list(command), stdout=handle, stderr=subprocess.STDOUT, text=True)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): {' '.join(command)}; see {log}")


# ---------------------------------------------------------------------------------------------------------------
# Genome panel


def load_panel(path: str | Path | None = None) -> pd.DataFrame:
    """The 55-species genome panel: Ensembl slug, lineage group and search stage.

    By default the copy shipped with the package (`rlmhelix/data/genome_panel_55.tsv`, identical to
    `data/genome_panel_55.tsv` of the repository); `path` overrides it.
    """
    path = Path(path) if path is not None else DEFAULT_PANEL
    if not path.exists():
        raise FileNotFoundError(f"{path}; pass the genome panel TSV with --panel")
    return pd.read_csv(path, sep="\t")


def seed_species(panel: pd.DataFrame) -> list[str]:
    """The 20 species of the seed search, whose high-confidence hits define the routes and seed alignments."""
    return sorted(panel.loc[panel.search_stage.str.startswith("seed search"), "ensembl_slug"])


def expansion_species(panel: pd.DataFrame) -> list[str]:
    """The 35 expansion genomes (34 vertebrates and the human genome)."""
    return sorted(panel.loc[panel.search_stage.str.contains("expansion"), "ensembl_slug"])


def lineage_map(panel: pd.DataFrame) -> dict[str, str]:
    """Ensembl slug -> lineage group."""
    return dict(zip(panel.ensembl_slug, panel.lineage_group))


def genome_path(directory: str | Path, species: str) -> Path:
    """Genome FASTA of a species: `<dir>/<slug>.fa(.gz)`, `<dir>/<slug>/*.fa.gz` or Ensembl's `<Slug>.*.dna.*.fa.gz`."""
    directory = Path(directory)
    candidates = [directory / f"{species}.fa.gz", directory / f"{species}.fa"]
    candidates += sorted((directory / species).glob("*.fa.gz")) if (directory / species).is_dir() else []
    candidates += sorted(directory.glob(f"{species.capitalize()}.*.dna.*.fa.gz"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"no genome FASTA for {species} under {directory}")


# ---------------------------------------------------------------------------------------------------------------
# Stage 1: top windows, loci, cores, representatives


def top_windows(windows: pd.DataFrame, fraction: float = TOP_FRACTION, p_column: str | None = None) -> pd.DataFrame:
    """The ceil(fraction * n) windows with the smallest calibrated P values (stable: ties keep the table order)."""
    p_column = _column(windows, p_column, P_COLUMNS, "P value")
    n_select = int(np.ceil(len(windows) * fraction))
    order = np.argsort(windows[p_column].to_numpy(dtype=np.float64), kind="mergesort")
    return windows.iloc[order[:n_select]].reset_index(drop=True)


def _locus_record(records: list[dict], p_column: str, score_column: str) -> dict[str, object]:
    peak = min(records, key=lambda row: (row[p_column], -row[score_column]))
    start = min(int(row["start"]) for row in records)
    end = max(int(row["end"]) for row in records)
    out = {"locus": f"{peak['chrom']}:{start}-{end}:{peak['strand']}", "chrom": peak["chrom"], "start": start,
           "end": end, "strand": peak["strand"], "n_windows": len(records), "peak_window": peak["window_id"],
           "peak_start": int(peak["start"]), "peak_end": int(peak["end"]), "peak_score": float(peak[score_column]),
           "peak_p": float(peak[p_column])}
    out.update({column: peak[column] for column in LOCUS_CARRY if column in peak})
    return out


def merge_loci(windows: pd.DataFrame, p_column: str | None = None, score_column: str | None = None) -> pd.DataFrame:
    """Merge windows that overlap on the same strand into loci; the peak window has the smallest P, then largest score.

    Loci are ordered by ascending peak P, descending peak score, chromosome and start.
    """
    p_column = _column(windows, p_column, P_COLUMNS, "P value")
    score_column = _column(windows, score_column, SCORE_COLUMNS, "score")
    windows = windows.copy()
    if "window_id" not in windows:
        windows["window_id"] = (windows.chrom + ":" + windows.start.astype(str) + "-" + windows.end.astype(str)
                                + ":" + windows.strand)
    rows: list[dict[str, object]] = []
    for _, group in windows.groupby(["chrom", "strand"], sort=False):
        current: list[dict] = []
        current_end = -1
        for record in group.sort_values(["start", "end"]).to_dict("records"):
            if current and int(record["start"]) >= current_end:
                rows.append(_locus_record(current, p_column, score_column))
                current, current_end = [], -1
            current.append(record)
            current_end = max(current_end, int(record["end"]))
        if current:
            rows.append(_locus_record(current, p_column, score_column))
    return pd.DataFrame(rows).sort_values(["peak_p", "peak_score", "chrom", "start"],
                                          ascending=[True, False, True, True]).reset_index(drop=True)


@dataclass(frozen=True)
class Core:
    """An outermost base-paired interval [start, end) of a structure and its retained pairs (0-based, i < j)."""

    start: int
    end: int
    pairs: tuple[tuple[int, int], ...]


def parse_dot_bracket(structure: str) -> list[tuple[int, int]]:
    """Pairs (i < j) of a dot-bracket string with '(' and ')'."""
    stack, pairs = [], []
    for index, symbol in enumerate(structure):
        if symbol == "(":
            stack.append(index)
        elif symbol == ")":
            if not stack:
                raise ValueError("unbalanced dot-bracket")
            pairs.append((stack.pop(), index))
    if stack:
        raise ValueError("unbalanced dot-bracket")
    return sorted(pairs)


def structure_cores(pairs: Iterable[tuple[int, int]], length: int, min_stem: int = MIN_STEM) -> list[Core]:
    """EvoFold-style cores: drop stems shorter than `min_stem`, then split at the outermost pairs (no flanks)."""
    canonical = {(min(int(i), int(j)), max(int(i), int(j))) for i, j in pairs}
    used = [x for pair in canonical for x in pair]
    if any(i == j for i, j in canonical) or len(used) != len(set(used)):
        raise ValueError("each nucleotide may occur in at most one pair")
    if any(j >= length or i < 0 for i, j in canonical):
        raise ValueError("pair outside the sequence")
    if set(parse_dot_bracket(dot_bracket(length, canonical))) != canonical:
        raise ValueError("crossing pairs cannot be split into cores")
    retained = {pair for block in stems(canonical) if len(block) >= min_stem for pair in block}
    cores, depth, begin = [], 0, 0
    for index, symbol in enumerate(dot_bracket(length, retained)):
        if symbol == "(":
            begin = index if depth == 0 else begin
            depth += 1
        elif symbol == ")":
            depth -= 1
            if depth == 0:
                inside = tuple(sorted(p for p in retained if begin <= p[0] < p[1] < index + 1))
                cores.append(Core(begin, index + 1, inside))
    return cores


def window_cores(record: Mapping[str, object], weights: np.ndarray, min_stem: int = MIN_STEM) -> list[dict]:
    """Cores of one peak window from its pair weights (student APC and BPP, sqrt(d * p)).

    `record` needs `sequence` (scan orientation), `strand`, `peak_start`, `peak_end`, `locus`, `peak_window`,
    `chrom` and `peak_score`. Genomic coordinates are 0-based half-open on the forward strand.
    """
    sequence = str(record["sequence"])
    weights = np.asarray(weights, dtype=np.float32)
    rows = []
    for index, core in enumerate(structure_cores(decode_pairs(weights, sequence), len(sequence), min_stem)):
        blocks = stems(set(core.pairs))
        if record["strand"] == "+":
            genome = (int(record["peak_start"]) + core.start, int(record["peak_start"]) + core.end)
        else:
            genome = (int(record["peak_end"]) - core.end, int(record["peak_end"]) - core.start)
        symbols = ["."] * (core.end - core.start)
        for left, right in core.pairs:
            symbols[left - core.start], symbols[right - core.start] = "(", ")"
        rows.append({
            "locus": record["locus"], "peak_window": record["peak_window"], "core_index": index,
            "chrom": record["chrom"], "strand": record["strand"], "features": record.get("features"),
            "window_score": float(record["peak_score"]), "core_start_in_window": core.start,
            "core_end_in_window": core.end, "core_length": core.end - core.start,
            "core_genome_start": genome[0], "core_genome_end": genome[1],
            "core_score": float(sum(weights[i, j] for i, j in core.pairs)),   # float32 accumulation, as in the paper
            "n_pairs": len(core.pairs), "n_stems": len(blocks), "longest_stem": max(map(len, blocks), default=0),
            "touches_window_edge": core.start == 0 or core.end == len(sequence),
            "core_sequence": sequence[core.start:core.end], "core_dot_bracket": "".join(symbols),
        })
    return rows


def rank_cores(cores: pd.DataFrame) -> pd.DataFrame:
    """Add each core's share of its locus score and its rank in the locus (1 = primary core)."""
    cores = cores.copy()
    cores["score_fraction"] = cores.core_score / cores.groupby("locus").core_score.transform("sum")
    cores["component_rank"] = cores.groupby("locus").core_score.rank(method="first", ascending=False).astype(int)
    return cores.sort_values(["window_score", "locus", "component_rank"],
                             ascending=[False, True, True]).reset_index(drop=True)


def _peak_window_cores(record: dict, apc: np.ndarray) -> list[dict]:
    from .score import base_pair_probabilities, pair_weights, to_rna

    return window_cores(record, pair_weights(apc, base_pair_probabilities(to_rna(str(record["sequence"])))))


def extract_cores(loci: pd.DataFrame, student, batch_size: int = 16, n_jobs: int = 1, chunk: int = 1024) -> pd.DataFrame:
    """Score the peak window of every locus with the student and extract its cores (GPU recommended).

    Loci are processed in table order in batches of `batch_size` (16 in the paper); under bf16 the batch composition
    shifts weights at the level of rounding. Without any core the result is an empty table with the core columns.
    """
    from joblib import Parallel, delayed

    records = loci.to_dict("records")
    rows: list[dict] = []
    for start in range(0, len(records), chunk):
        part = records[start:start + chunk]
        maps = student.predict([str(r["sequence"]) for r in part], batch_size)
        nested = Parallel(n_jobs=n_jobs)(delayed(_peak_window_cores)(r, m) for r, m in zip(part, maps))
        rows.extend(row for found in nested for row in found)
    return rank_cores(pd.DataFrame(rows) if rows else _empty_frame(CORE_DTYPES))


def canonical_key(sequence: str, structure: str) -> tuple[str, str]:
    """Lexicographically smaller of (sequence, structure) and (reverse complement, reversed structure), in RNA."""
    sequence = sequence.upper().replace("T", "U")
    reverse = (sequence.translate(_RNA_COMPLEMENT)[::-1], structure[::-1].translate(str.maketrans("()", ")(")))
    return min((sequence, structure), reverse)


def collapse_cores(cores: pd.DataFrame, loci: pd.DataFrame, id_prefix: str = ID_PREFIX
                   ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collapse identical primary cores (either strand); returns (representatives, membership of every primary core).

    Representatives are chosen by ascending peak P, descending core score and peak window, and numbered in that
    order. `sequence_u` and `core_dot_bracket` are the source strand; `canonical_sequence` and `canonical_structure`
    the normalized orientation.
    """
    columns = [c for c in ("peak_window", "peak_p", "GC", "n_windows") if c in loci]
    primary = cores.loc[cores.component_rank.eq(1)].merge(loci[columns], on="peak_window", how="left",
                                                          validate="one_to_one")
    if primary.peak_p.isna().any():
        raise ValueError("a primary core lacks its calibrated P value")
    primary["sequence_u"] = primary.core_sequence.str.upper().str.replace("T", "U")
    if (~primary.core_dot_bracket.str.fullmatch(r"[().]+", na=False)).any():
        raise ValueError("non-dot-bracket symbols found")
    if not (primary.sequence_u.str.len() == primary.core_dot_bracket.str.len()).all():
        raise ValueError("sequence/structure length mismatch")
    keys = [canonical_key(s, d) for s, d in zip(primary.sequence_u, primary.core_dot_bracket)]
    primary["canonical_sequence"] = [key[0] for key in keys]
    primary["canonical_structure"] = [key[1] for key in keys]
    primary["dedup_key"] = primary.canonical_sequence + "|" + primary.canonical_structure
    primary = primary.sort_values(["peak_p", "core_score", "peak_window"],
                                  ascending=[True, False, True]).reset_index(drop=True)
    representatives = primary.drop_duplicates("dedup_key").copy()
    representatives.insert(0, "query_id", [f"{id_prefix}{i:07d}" for i in range(1, len(representatives) + 1)])
    membership = primary.merge(representatives[["dedup_key", "query_id"]], on="dedup_key", how="left",
                               validate="many_to_one")
    return representatives.reset_index(drop=True), membership


# ---------------------------------------------------------------------------------------------------------------
# Stage 2: genome chunks, Riboseek, hits, routes


def genome_chunks(sequence: str, chunk_size: int = CHUNK_SIZE, overlap: int = CHUNK_OVERLAP
                  ) -> Iterator[tuple[int, int, str]]:
    """(start, end, chunk) of overlapping chunks; chunks made only of N are skipped; the last chunk ends the contig."""
    if not 0 <= overlap < chunk_size:
        raise ValueError("overlap must satisfy 0 <= overlap < chunk_size")
    sequence = sequence.upper().replace("U", "T")
    for start in range(0, len(sequence), chunk_size - overlap):
        end = min(start + chunk_size, len(sequence))
        chunk = sequence[start:end]
        if chunk and not set(chunk) <= {"N"}:
            yield start, end, chunk
        if end == len(sequence):
            break


def chunk_genome(genome: str | Path, species: str, out_dir: str | Path, chunk_size: int = CHUNK_SIZE,
                 overlap: int = CHUNK_OVERLAP) -> Path:
    """Write `<species>.10kb_o500.fa.gz` (chunk IDs `species|contig|start|end`) and a TSV of the chunks."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    label = f"{chunk_size // 1000}kb" if chunk_size % 1000 == 0 else f"{chunk_size}nt"
    output = out_dir / f"{species}.{label}_o{overlap}.fa.gz"
    rows = []
    with gzip.open(output.with_suffix(".gz.part"), "wt", compresslevel=1) as handle:
        for contig, sequence in iter_fasta(genome):
            for start, end, chunk in genome_chunks(sequence, chunk_size, overlap):
                chunk_id = f"{species}|{contig}|{start}|{end}"
                handle.write(f">{chunk_id}\n{chunk}\n")
                rows.append((chunk_id, species, contig, start, end, end - start))
    output.with_suffix(".gz.part").replace(output)
    pd.DataFrame(rows, columns=["chunk_id", "species", "contig", "start", "end", "length"]).to_csv(
        out_dir / f"{species}.{label}_o{overlap}.tsv", sep="\t", index=False)
    return output


def riboseek_commands(riboseek: str, query_db: str | Path, target_fasta: str | Path, work_dir: str | Path,
                      threads: int = 20) -> list[list[str]]:
    """createdb, search and convertalis commands of one species, with the settings of the paper."""
    work_dir = Path(work_dir)
    target_db, aln, tmp = work_dir / "target_db", work_dir / "aln", work_dir / "tmp"
    return [
        [riboseek, "createdb", str(target_fasta), str(target_db)],
        [riboseek, "search", str(query_db), str(target_db), str(aln), str(tmp), *RIBOSEEK_SEARCH,
         "--threads", str(threads), "--remove-tmp-files", "1"],
        [riboseek, "convertalis", str(query_db), str(target_db), str(aln), str(work_dir / "hits.tsv"),
         "--format-output", HIT_FORMAT],
    ]


def _clear_search(work_root: Path) -> None:
    """Remove the query database and the `done` markers of a work directory (hits are then searched again)."""
    for path in [*work_root.glob("query_db*"), *work_root.glob("*/done"), work_root / "queries.sha256"]:
        if path.is_dir() and not path.is_symlink():
            shutil.rmtree(path)
        elif path.exists() or path.is_symlink():
            path.unlink()


def run_riboseek(queries: str | Path, chunk_fastas: Mapping[str, str | Path], work_root: str | Path,
                 riboseek: str | None = None, threads: int = 20, cleanup: bool = True,
                 force: bool = False) -> pd.DataFrame:
    """Search the queries against each species' chunk database; writes `<work_root>/<species>/hits.tsv`.

    Completed species (with `hits.tsv` and a `done` marker) are skipped, so a run can be resumed. The SHA-256 of the
    query FASTA is kept in `<work_root>/queries.sha256` beside `query_db`; a work directory whose query database or
    finished species come from other (or unrecorded) queries is refused with FileExistsError unless `force`, which
    removes the query database and the `done` markers so that every species is searched again.
    """
    binary = riboseek or tool("riboseek")
    work_root = Path(work_root)
    work_root.mkdir(parents=True, exist_ok=True)
    query_db, stamp = work_root / "query_db", work_root / "queries.sha256"
    digest = file_sha256(queries)
    if force:
        _clear_search(work_root)
    elif query_db.exists() or any(work_root.glob("*/done")):
        recorded = stamp.read_text().split()[:1] if stamp.exists() else []
        if recorded != [digest]:
            raise FileExistsError(
                f"{work_root} holds a query database or hits of other queries than {queries} "
                f"({'checksum differs from' if recorded else 'no checksum in'} {stamp}); use another output "
                "directory, or --force to rebuild the query database and search every species again")
    if not query_db.exists():
        run_logged([binary, "createdb", str(Path(queries).resolve()), str(query_db)])
        stamp.write_text(f"{digest}  {Path(queries).resolve()}\n")
    rows = []
    for species, fasta in sorted(chunk_fastas.items()):
        species_dir = work_root / species
        species_dir.mkdir(parents=True, exist_ok=True)
        if (species_dir / "done").exists() and (species_dir / "hits.tsv").exists():
            rows.append({"species": species, "status": "reused"})
            continue
        started = time.monotonic()
        for command in riboseek_commands(binary, query_db, Path(fasta).resolve(), species_dir, threads):
            run_logged(command)
        if cleanup:
            for database in ("aln", "target_db"):
                run_logged([binary, "rmdb", str(species_dir / database)])
            for path in species_dir.glob("target_db_h*"):
                path.unlink()
            shutil.rmtree(species_dir / "tmp", ignore_errors=True)
        (species_dir / "done").write_text("riboseek search complete\n")
        rows.append({"species": species, "status": "success", "wall_seconds": time.monotonic() - started})
    return pd.DataFrame(rows)


def read_hits(path: str | Path, query_ids: set[str] | None = None, max_evalue: float | None = None,
              columns: Sequence[str] | None = None, chunksize: int = 500_000) -> pd.DataFrame:
    """Read a Riboseek `convertalis` table (HIT_FORMAT), optionally keeping some queries and hits with E <= max_evalue.

    Without any kept hit (or for an empty file) the result is an empty table with typed columns.
    """
    names = HIT_COLUMNS if columns is None else [c for c in HIT_COLUMNS if c in columns]
    dtype = {c: HIT_DTYPES[c] for c in ("query_id", "target_id", "qaln", "taln") if c in names}
    frames = []
    for chunk in pd.read_csv(path, sep="\t", header=None, names=HIT_COLUMNS, usecols=names, dtype=dtype,
                             chunksize=chunksize):
        keep = np.ones(len(chunk), dtype=bool)
        if query_ids is not None:
            keep &= chunk.query_id.isin(query_ids).to_numpy()
        if max_evalue is not None:
            keep &= (chunk.evalue <= max_evalue).to_numpy()
        if keep.any():
            frames.append(chunk.loc[keep])
    if not frames:
        return _empty_frame({c: HIT_DTYPES[c] for c in names})
    return pd.concat(frames, ignore_index=True)


def high_confidence(hits: pd.DataFrame) -> pd.Series:
    """E <= 1e-3 and query coverage >= 0.8."""
    return (hits.evalue <= HC_EVALUE) & (hits.qcov >= HC_QCOV)


def high_confidence_queries(path: str | Path) -> set[str]:
    """Queries with at least one high-confidence hit in one species' hit table."""
    hits = read_hits(path, columns=["query_id", "evalue", "qcov"], max_evalue=HC_EVALUE)
    return set(hits.loc[high_confidence(hits), "query_id"])


def lineage_counts(hc_species: pd.DataFrame, lineage: Mapping[str, str]) -> pd.DataFrame:
    """Number of species with a high-confidence hit per query and lineage group.

    `hc_species` has one row per (query_id, species) with a high-confidence hit (duplicates are ignored).
    """
    pairs = hc_species[["query_id", "species"]].drop_duplicates().copy()
    pairs["lineage"] = pairs.species.map(lineage)
    if pairs.lineage.isna().any():
        raise ValueError(f"species without lineage group: {sorted(set(pairs.species[pairs.lineage.isna()]))}")
    counts = pairs.groupby(["query_id", "lineage"]).size().unstack(fill_value=0)
    for group in LINEAGES:
        if group not in counts:
            counts[group] = 0
    return counts[list(LINEAGES)].astype(int)


def assign_routes(candidates: pd.DataFrame, counts: pd.DataFrame | None = None) -> pd.DataFrame:
    """Add `route_<name>` flags, `selected` and `route` (first route met, or 'not_selected').

    `candidates` needs `query_id` and `peak_p`, and either the four lineage columns or `counts` (from
    `lineage_counts`, indexed by query_id); queries without high-confidence hits count 0.
    """
    out = candidates.copy()
    if counts is not None:
        out = out.drop(columns=[c for c in LINEAGES if c in out]).merge(
            counts, left_on="query_id", right_index=True, how="left")
    out[list(LINEAGES)] = out[list(LINEAGES)].fillna(0).astype(int)
    flags = {
        "vertebrate_wide": (out.nonmammal_tetrapod >= 1) & (out.other_vertebrate >= 1),
        "mammal_conserved": (out.nonprimate_mammal >= 3) & (out.primate >= 3),
        "primate_focused": (out.primate >= 6) & (out.peak_p <= PRIMATE_FOCUSED_P),
        "high_score_rescue": out.peak_p <= RESCUE_P,
    }
    for name in ROUTES:
        out[f"route_{name}"] = flags[name]
    out["selected"] = out[[f"route_{name}" for name in ROUTES]].any(axis=1)
    out["route"] = "not_selected"
    for name in reversed(ROUTES):
        out.loc[out[f"route_{name}"], "route"] = name
    return out


# ---------------------------------------------------------------------------------------------------------------
# Hit loci and seed alignments


def locate_hits(hits: pd.DataFrame, species: str | None = None) -> pd.DataFrame:
    """Genomic coordinates, target strand and query-oriented alignment strings of Riboseek hits on chunks.

    Riboseek coordinates are 1-based inclusive within a chunk; loci are 0-based half-open. A hit is on the minus
    strand when qstart > qend; its alignment strings are reverse-complemented to the query orientation.
    """
    if hits.empty:   # an empty table with every output column
        target = pd.DataFrame(columns=range(4), dtype="str")
    else:
        target = hits.target_id.astype(str).str.rsplit("|", n=3, expand=True)
    if target.shape[1] != 4:
        raise ValueError("malformed target identifiers (expected species|contig|start|end)")
    target.columns = ["target_species", "contig", "chunk_start", "chunk_end"]
    if species is not None and len(target) and set(target.target_species.dropna()) != {species}:
        raise ValueError(f"{species}: target species mismatch")
    target["chunk_start"] = pd.to_numeric(target.chunk_start, errors="raise")
    target["chunk_end"] = pd.to_numeric(target.chunk_end, errors="raise")
    hits = pd.concat([hits.reset_index(drop=True), target.reset_index(drop=True)], axis=1)
    hits["genome_start"] = hits.chunk_start + hits[["tstart", "tend"]].min(axis=1) - 1
    hits["genome_end"] = hits.chunk_start + hits[["tstart", "tend"]].max(axis=1)
    if ((hits.genome_start < hits.chunk_start) | (hits.genome_end > hits.chunk_end)).any():
        raise ValueError("hits leave their source chunk")
    reverse = hits.qstart > hits.qend
    hits["target_strand"] = reverse.map({False: "+", True: "-"})
    hits["oriented_qaln"] = hits.qaln.where(~reverse, hits.qaln.map(reverse_complement))
    hits["oriented_taln"] = hits.taln.where(~reverse, hits.taln.map(reverse_complement))
    hits["target_sequence"] = (hits.oriented_taln.str.replace("-", "", regex=False).str.replace(".", "", regex=False)
                               .str.upper().str.replace("U", "T", regex=False))
    hits["high_confidence"] = high_confidence(hits)
    return hits


def merge_hit_loci(hits: pd.DataFrame) -> pd.DataFrame:
    """Merge located hits that overlap per query, contig and strand into loci (one representative row each).

    The representative is the best hit (high-confidence first, then bits, E value and coverage); `locus_start`,
    `locus_end`, `n_chunk_hits` and `any_high_confidence` describe the merged locus. Distinct loci (paralogues) are kept.
    """
    keys = ["query_id", "contig", "target_strand"]
    ordered = hits.sort_values(keys + ["genome_start", "genome_end"]).copy()
    previous_end = ordered.groupby(keys, sort=False).genome_end.cummax().groupby(
        [ordered[c] for c in keys], sort=False).shift()
    new_locus = previous_end.isna() | (ordered.genome_start >= previous_end)
    ordered["locus_component"] = new_locus.groupby([ordered[c] for c in keys], sort=False).cumsum().astype("int32")
    locus_keys = keys + ["locus_component"]
    bounds = ordered.groupby(locus_keys, sort=False).agg(
        locus_start=("genome_start", "min"), locus_end=("genome_end", "max"), n_chunk_hits=("target_id", "size"),
        any_high_confidence=("high_confidence", "any"))
    loci = (ordered.sort_values(locus_keys + ["high_confidence", "bits", "evalue", "qcov"],
                                ascending=[True, True, True, True, False, False, True, False])
            .drop_duplicates(locus_keys).set_index(locus_keys).join(bounds).reset_index())
    loci["species"] = loci.target_species.astype(str)
    loci["locus_id"] = (loci.species + "|" + loci.contig.astype(str) + ":" + loci.locus_start.astype(str) + "-"
                        + loci.locus_end.astype(str) + ":" + loci.target_strand.astype(str))
    return loci


def hit_loci(path: str | Path, species: str, query_ids: set[str] | None = None,
             max_evalue: float = LOCUS_EVALUE) -> pd.DataFrame:
    """Loci of one species' hit table (hits with E <= 0.1)."""
    return merge_hit_loci(locate_hits(read_hits(path, query_ids, max_evalue), species))


def locus_counts(loci: pd.DataFrame) -> pd.DataFrame:
    """Per query: species and loci with E <= 0.1 and with a high-confidence hit."""
    per_species = loci.assign(_hc=loci.any_high_confidence.astype(bool).astype("int64")).groupby(
        ["query_id", "species"]).agg(n_loci=("locus_id", "nunique"), n_hc=("_hc", "sum")).reset_index()
    return per_species.groupby("query_id").agg(
        n_species_e01=("species", "nunique"), n_loci_e01=("n_loci", "sum"),
        n_species_hc=("n_hc", lambda values: int((values > 0).sum())), n_loci_hc=("n_hc", "sum"))


def project_hit(oriented_qaln: str, oriented_taln: str, qstart: int, qend: int, query_length: int) -> str:
    """Target residues placed on the query columns of a pairwise alignment (gaps elsewhere, insertions dropped)."""
    query_alignment = oriented_qaln.upper().replace("U", "T")
    target_alignment = oriented_taln.upper().replace("U", "T")
    if len(query_alignment) != len(target_alignment):
        raise ValueError("pairwise alignment width mismatch")
    position = min(int(qstart), int(qend)) - 1
    projected = ["-"] * query_length
    for query_base, target_base in zip(query_alignment, target_alignment):
        if query_base in "-.":
            continue
        if position >= query_length:
            raise ValueError("pairwise alignment exceeds query length")
        projected[position] = "-" if target_base in "-." else target_base
        position += 1
    if position != max(int(qstart), int(qend)):
        raise ValueError("pairwise alignment disagrees with the reported query coordinates")
    return "".join(projected)


def seed_rows(query_id: str, query_sequence: str, hc_loci: pd.DataFrame) -> tuple[list[tuple[str, str]], pd.DataFrame, int]:
    """Seed alignment rows: the human core, then one row per distinct target sequence of its high-confidence loci.

    Loci are ordered by bits, E value and coverage; returns (rows, the kept loci with their sequence IDs, number of
    pairwise alignments that could not be projected).
    """
    sequence = query_sequence.upper().replace("U", "T")
    loci = hc_loci.sort_values(["bits", "evalue", "qcov"], ascending=[False, True, False]).copy()
    projected, invalid = [], 0
    for row in loci.itertuples(index=False):
        try:
            projected.append(project_hit(row.oriented_qaln, row.oriented_taln, row.qstart, row.qend, len(sequence)))
        except ValueError:
            projected.append("")
            invalid += 1
    loci["projected"] = projected
    unique = loci.loc[loci.projected.ne("")].drop_duplicates("target_sequence").copy()
    unique["sequence_id"] = [safe_name(f"h{i:05d}_{s}") for i, s in enumerate(unique.species, start=1)]
    rows = [(safe_name(query_id), sequence)] + list(zip(unique.sequence_id, unique.projected))
    return rows, unique, invalid


def write_stockholm(path: str | Path, name: str, rows: Sequence[tuple[str, str]], structure: str) -> None:
    """Stockholm alignment with the structure as SS_cons and every column as RF ('x')."""
    width = max(max(len(n) for n, _ in rows), len("#=GC SS_cons")) + 2
    with Path(path).open("w") as handle:
        handle.write(f"# STOCKHOLM 1.0\n#=GF ID {safe_name(name)}\n")
        for row_name, sequence in rows:
            handle.write(f"{row_name:<{width}}{sequence}\n")
        handle.write(f"{'#=GC SS_cons':<{width}}{structure}\n{'#=GC RF':<{width}}{'x' * len(structure)}\n//\n")


# ---------------------------------------------------------------------------------------------------------------
# Crops, Infernal and the inclusion rule


def seed_crop(locus_start: int, locus_end: int) -> tuple[int, int]:
    """Crop of a seed locus: one locus length of flank on either side, [start - L, end + L) (end clipped later)."""
    span = locus_end - locus_start
    if span <= 0:
        raise ValueError("non-positive locus span")
    return max(0, locus_start - span), locus_end + span


def expansion_crop(locus_start: int, locus_end: int, model_length: int) -> tuple[int, int]:
    """Crop of an expansion locus: 3 x model length centred on the locus midpoint (end clipped later)."""
    length = 3 * int(model_length)
    start = max(0, (locus_start + locus_end) // 2 - length // 2)
    return start, start + length


def extract_crops(genome: str | Path, loci: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Cut crops around loci from one species' genome, in the orientation of the locus strand.

    `mode` is 'seed' (flank of one locus length) or 'expansion' (3 x `model_length`, a column of `loci`).
    """
    requested: dict[str, list[dict]] = defaultdict(list)
    for row in loci.itertuples(index=False):
        if mode == "seed":
            start, end = seed_crop(int(row.locus_start), int(row.locus_end))
        elif mode == "expansion":
            start, end = expansion_crop(int(row.locus_start), int(row.locus_end), int(row.model_length))
        else:
            raise ValueError("mode must be 'seed' or 'expansion'")
        record = {"query_id": str(row.query_id), "species": str(row.species), "contig": str(row.contig),
                  "strand": str(row.target_strand), "locus_id": str(row.locus_id), "riboseek_evalue": float(row.evalue),
                  "riboseek_bits": float(row.bits), "riboseek_qcov": float(row.qcov),
                  "locus_start": int(row.locus_start), "locus_end": int(row.locus_end),
                  "crop_start": start, "crop_end_requested": end}
        if mode == "expansion":
            record["model_length"] = int(row.model_length)
        requested[str(row.contig)].append(record)
    records, seen = [], set()
    for contig, sequence in iter_fasta(genome):
        if contig not in requested:
            continue
        seen.add(contig)
        sequence = sequence.upper().replace("U", "T")
        for record in requested[contig]:
            crop_end = min(len(sequence), record["crop_end_requested"])
            crop = sequence[record["crop_start"]:crop_end]
            if record["strand"] == "-":
                crop = reverse_complement(crop)
                offsets = (crop_end - record["locus_end"], crop_end - record["locus_start"])
            else:
                offsets = (record["locus_start"] - record["crop_start"], record["locus_end"] - record["crop_start"])
            records.append({**record, "crop_end": crop_end, "crop_length": len(crop), "locus_offset_start": offsets[0],
                            "locus_offset_end": offsets[1], "crop_sequence": crop})
    missing = set(requested) - seen
    if missing:
        raise RuntimeError(f"{len(missing)} requested contigs absent from {genome}")
    if not records:
        return pd.DataFrame()
    return pd.DataFrame(records).sort_values(["query_id", "contig", "crop_start", "strand"]).reset_index(drop=True)


def cmbuild_command(cmbuild: str, cm: str | Path, alignment: str | Path) -> list[str]:
    return [cmbuild, "-F", "--hand", str(cm), str(alignment)]


def cmalign_command(cmalign: str, cm: str | Path, fasta: str | Path, output: str | Path,
                    mapali: str | Path | None = None) -> list[str]:
    extra = ["--mapali", str(mapali)] if mapali is not None else []
    return [cmalign, "--cpu", "1", "--cyk", "--notrunc", "--noprob", *extra, "-o", str(output), str(cm), str(fasta)]


def cmcalibrate_command(cmcalibrate: str, cm: str | Path) -> list[str]:
    return [cmcalibrate, "-L", "0.01", "--cpu", "1", str(cm)]


def cmsearch_command(cmsearch: str, cm: str | Path, fasta: str | Path, tblout: str | Path,
                     output: str | Path = "/dev/null") -> list[str]:
    """Glocal CYK search reporting every hit (boundaries of crops and reciprocal scores)."""
    return [cmsearch, "--max", "-g", "--cyk", "--acyk", "--notrunc", "--toponly", "--noali", "--cpu", "1",
            "--smxsize", "2048", "-T", "-1000", "--tblout", str(tblout), "-o", str(output), str(cm), str(fasta)]


def hmmsearch_command(cmsearch: str, cm: str | Path, fasta: str | Path, tblout: str | Path) -> list[str]:
    """HMM-only filter search (E <= 10) proposing candidate pairs for the families."""
    return [cmsearch, "--hmmonly", "--toponly", "--noali", "--cpu", "1", "-E", str(HMM_EVALUE), "--tblout",
            str(tblout), "-o", "/dev/null", str(cm), str(fasta)]


def parse_tblout(source: str | Path) -> pd.DataFrame:
    """Hits of an Infernal `--tblout` file (a Path or file name) or of its text."""
    is_file = isinstance(source, Path) or ("\n" not in source and Path(source).is_file())
    text = Path(source).read_text() if is_file else str(source)
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=17)
        rows.append({"target": fields[0], "model_from": int(fields[5]), "model_to": int(fields[6]),
                     "sequence_from": int(fields[7]), "sequence_to": int(fields[8]), "cm_strand": fields[9],
                     "cm_trunc": fields[10], "cm_gc": float(fields[12]), "cm_bias": float(fields[13]),
                     "cm_bits": float(fields[14]), "cm_evalue": float(fields[15]), "included": fields[16] == "!"})
    return pd.DataFrame(rows, columns=["target", "model_from", "model_to", "sequence_from", "sequence_to",
                                       "cm_strand", "cm_trunc", "cm_gc", "cm_bias", "cm_bits", "cm_evalue",
                                       "included"])


def best_crop_hits(hits: pd.DataFrame, crops: pd.DataFrame) -> pd.DataFrame:
    """Best CM hit per crop (by bits), its sequence in the CM orientation and, with `model_length`, model coverage."""
    best = hits.rename(columns={"target": "crop_id"}).sort_values(["crop_id", "cm_bits"], ascending=[True, False])
    best = crops.merge(best.drop_duplicates("crop_id"), on="crop_id", how="left", validate="one_to_one")
    if best.cm_bits.isna().any():
        raise RuntimeError(f"cmsearch missed {int(best.cm_bits.isna().sum())} crops")
    left = best[["sequence_from", "sequence_to"]].min(axis=1).astype(int) - 1
    right = best[["sequence_from", "sequence_to"]].max(axis=1).astype(int)
    sequences = [crop[a:b] for crop, a, b in zip(best.crop_sequence, left, right)]
    best["cm_hit_sequence"] = [reverse_complement(s) if strand == "-" else s
                               for s, strand in zip(sequences, best.cm_strand)]
    best["cm_hit_length"] = best.cm_hit_sequence.str.len()
    if "model_length" in best:
        best["cm_model_coverage"] = (best.model_to - best.model_from + 1) / best.model_length
    best["cm_hit_offset_start"], best["cm_hit_offset_end"] = left, right
    best["left_boundary_shift"] = best.cm_hit_offset_start - best.locus_offset_start
    best["right_boundary_shift"] = best.cm_hit_offset_end - best.locus_offset_end
    return best


def inclusion_threshold(seed_bits: Iterable[float], quantile: float = INCLUSION_QUANTILE,
                        margin: float = INCLUSION_MARGIN) -> float:
    """10th percentile (linear interpolation) of the seed hits' bit scores minus 10 bits."""
    values = pd.Series(list(seed_bits), dtype=float)
    if values.empty or not np.isfinite(values).all():
        raise ValueError("seed bit scores must be finite and non-empty")
    return float(values.quantile(quantile)) + margin


def apply_inclusion(hits: pd.DataFrame, seed_bits: Iterable[float],
                    min_coverage: float = INCLUSION_COVERAGE) -> pd.DataFrame:
    """Flag expansion hits with bits >= the inclusion threshold and model coverage >= 0.8 as `accepted`."""
    out = hits.copy()
    out["inclusion_threshold_bits"] = inclusion_threshold(seed_bits)
    out["accepted"] = out.cm_bits.ge(out.inclusion_threshold_bits) & out.cm_model_coverage.ge(min_coverage)
    return out


def final_additions(hits: pd.DataFrame, seed_sequences: Iterable[str]) -> pd.DataFrame:
    """Accepted hits added to the final alignment: one per distinct sequence (best bits, then E value), minus seed rows."""
    accepted = hits.loc[hits.accepted].sort_values(["cm_bits", "riboseek_evalue"], ascending=[False, True]).copy()
    accepted["sequence_group"] = pd.factorize(accepted.cm_hit_sequence, sort=False)[0] + 1
    added = accepted.drop_duplicates("cm_hit_sequence").copy()
    added["represented_loci"] = added.sequence_group.map(accepted.sequence_group.value_counts())
    added["represented_species"] = added.sequence_group.map(accepted.groupby("sequence_group").species.nunique())
    seed = {s.replace("-", "").replace(".", "").upper().replace("U", "T") for s in seed_sequences}
    added = added.loc[~added.cm_hit_sequence.isin(seed)].copy()
    added["final_sequence_id"] = [f"e{i:05d}_{s}" for i, s in enumerate(added.species, start=1)]
    return added


def high_copy(n_loci_hc: pd.Series | int, n_final55_loci: pd.Series | int, limit: int = HIGH_COPY_LOCI):
    """More than 100 high-confidence seed loci or final55 loci."""
    return (n_loci_hc > limit) | (n_final55_loci > limit)


# ---------------------------------------------------------------------------------------------------------------
# Per-candidate Infernal pipelines


def _stockholm_sequences(path: str | Path) -> list[str]:
    from .msa import read_stockholm

    return [row.replace("-", "").replace(".", "").upper().replace("U", "T") for row in read_stockholm(path).rows]


def build_seed(query_id: str, query_sequence: str, structure: str, hc_loci: pd.DataFrame, crops: pd.DataFrame,
               out_dir: str | Path, infernal_bin: str | Path | None = None) -> dict[str, object]:
    """Seed alignment and CM of one candidate, then the boundary search in seed crops and the seed CM.

    Writes seed.sto, sequences.fa, final.cm (alignment of the high-confidence loci) and, when there are crops,
    calibrated.cm, cmsearch_best_hits.parquet, realigned_3l.sto and realigned_3l.cm (the seed CM).
    """
    cmbuild, cmalign = tool("cmbuild", infernal_bin), tool("cmalign", infernal_bin)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    sequence = query_sequence.upper().replace("U", "T")
    if len(sequence) != len(structure):
        raise ValueError(f"{query_id}: query/structure length mismatch")
    rows, unique, invalid = seed_rows(query_id, sequence, hc_loci)
    write_stockholm(out / "seed.sto", query_id, rows, structure)
    write_fasta([rows[0]] + list(zip(unique.sequence_id, unique.target_sequence)), out / "sequences.fa")
    unique.to_parquet(out / "sequence_mapping.parquet", index=False)
    run_logged(cmbuild_command(cmbuild, out / "initial.cm", out / "seed.sto"), out / "cmbuild_initial.log")
    run_logged(cmalign_command(cmalign, out / "initial.cm", out / "sequences.fa", out / "refined.sto"),
               out / "cmalign_initial.log")
    run_logged(cmbuild_command(cmbuild, out / "final.cm", out / "refined.sto"), out / "cmbuild_final.log")
    summary = {"query_id": query_id, "n_hc_loci": int(hc_loci.locus_id.nunique()) if len(hc_loci) else 0,
               "n_hc_species": int(hc_loci.species.nunique()) if len(hc_loci) else 0,
               "n_invalid_pairwise_alignments": invalid, "n_unique_target_sequences": len(unique),
               "seed_cm": False}
    if crops.empty:
        return summary
    crops = crops.reset_index(drop=True).copy()
    crops["crop_id"] = [f"c{i:06d}" for i in range(1, len(crops) + 1)]
    write_fasta(zip(crops.crop_id, crops.crop_sequence), out / "crops.fa")
    shutil.copy2(out / "final.cm", out / "calibrated.cm")
    run_logged(cmcalibrate_command(tool("cmcalibrate", infernal_bin), out / "calibrated.cm"), out / "cmcalibrate.log")
    run_logged(cmsearch_command(tool("cmsearch", infernal_bin), out / "calibrated.cm", out / "crops.fa",
                                out / "cmsearch.tbl", out / "cmsearch.out"), out / "cmsearch.log")
    best = best_crop_hits(parse_tblout(out / "cmsearch.tbl"), crops)
    best.drop(columns="crop_sequence").to_parquet(out / "cmsearch_best_hits.parquet", index=False)
    realign = best.sort_values(["cm_bits", "riboseek_evalue"], ascending=[False, True]).drop_duplicates(
        "cm_hit_sequence")
    names = [f"h{i:05d}_{s}" for i, s in enumerate(realign.species, start=1)]
    write_fasta([(query_id, sequence)] + list(zip(names, realign.cm_hit_sequence)), out / "realign.fa")
    run_logged(cmalign_command(cmalign, out / "calibrated.cm", out / "realign.fa", out / "realigned_3l.sto"),
               out / "cmalign.log")
    run_logged(cmbuild_command(cmbuild, out / "realigned_3l.cm", out / "realigned_3l.sto"), out / "cmbuild.log")
    summary.update(seed_cm=True, n_crops=len(crops), n_seed_sequences=len(realign) + 1,
                   seed_q10_bits=float(best.cm_bits.quantile(INCLUSION_QUANTILE)),
                   threshold_bits=inclusion_threshold(best.cm_bits))
    return summary


def build_final55(query_id: str, crops: pd.DataFrame, seed_dir: str | Path, out_dir: str | Path,
                  lineage: Mapping[str, str], infernal_bin: str | Path | None = None) -> dict[str, object]:
    """Search expansion crops with the seed CM, apply the inclusion rule and build final55.sto and final55.cm."""
    cmbuild = tool("cmbuild", infernal_bin)
    seed_dir, out = Path(seed_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    summary: dict[str, object] = {"query_id": query_id}
    if crops.empty:
        added = pd.DataFrame()
        accepted = pd.DataFrame(columns=["species", "lineage"])
    else:
        crops = crops.reset_index(drop=True).copy()
        crops["crop_id"] = [f"c{i:06d}" for i in range(1, len(crops) + 1)]
        write_fasta(zip(crops.crop_id, crops.crop_sequence), out / "crops.fa")
        shutil.copy2(seed_dir / "realigned_3l.cm", out / "seed20_3l.calibrated.cm")
        run_logged(cmcalibrate_command(tool("cmcalibrate", infernal_bin), out / "seed20_3l.calibrated.cm"),
                   out / "cmcalibrate.log")
        run_logged(cmsearch_command(tool("cmsearch", infernal_bin), out / "seed20_3l.calibrated.cm", out / "crops.fa",
                                    out / "cmsearch.tbl", out / "cmsearch.out"), out / "cmsearch.log")
        best = best_crop_hits(parse_tblout(out / "cmsearch.tbl"), crops)
        seed_bits = pd.read_parquet(seed_dir / "cmsearch_best_hits.parquet", columns=["cm_bits"]).cm_bits
        hits = apply_inclusion(best, seed_bits)
        hits["lineage"] = hits.species.map(lineage)
        hits.drop(columns="crop_sequence").to_parquet(out / "expansion_hits.parquet", index=False)
        accepted = hits.loc[hits.accepted]
        added = final_additions(hits, _stockholm_sequences(seed_dir / "realigned_3l.sto"))
        added.drop(columns="crop_sequence").to_parquet(out / "accepted_expansion_hits.parquet", index=False)
    if len(added):
        write_fasta(zip(added.final_sequence_id, added.cm_hit_sequence), out / "accepted_expansion.fa")
        run_logged(cmalign_command(tool("cmalign", infernal_bin), seed_dir / "realigned_3l.cm",
                                   out / "accepted_expansion.fa", out / "final55.sto",
                                   mapali=seed_dir / "realigned_3l.sto"), out / "cmalign.log")
    else:
        shutil.copy2(seed_dir / "realigned_3l.sto", out / "final55.sto")
    run_logged(cmbuild_command(cmbuild, out / "final55.cm", out / "final55.sto"), out / "cmbuild.log")
    summary.update(n_accepted_expansion_loci=len(accepted), n_accepted_expansion_species=int(accepted.species.nunique()),
                   n_added_unique_sequences=len(added),
                   n_final_alignment_sequences=len(_stockholm_sequences(out / "final55.sto")))
    for group in LINEAGES:
        summary[f"n_accepted_{group}_species"] = int(accepted.loc[accepted.lineage.eq(group), "species"].nunique())
    return summary


def _run_task(task: tuple) -> dict[str, object]:
    function, kwargs = task
    try:
        return {"status": "success", "error": "", **function(**kwargs)}
    except Exception as exc:  # keep failures in the manifest
        return {"query_id": kwargs.get("query_id"), "status": "failed", "error": repr(exc)}


def _run_parallel(tasks: list[tuple], workers: int) -> pd.DataFrame:
    """Run (function, kwargs) tasks in worker processes; one manifest row per task (status, error, summary)."""
    if not tasks:
        return _empty_frame({"query_id": "str", "status": "str", "error": "str"})
    with ProcessPoolExecutor(max_workers=workers) as executor:
        return pd.DataFrame(list(executor.map(_run_task, tasks)))


def failed_tasks(manifest: pd.DataFrame, what: str) -> str:
    """'' when every task succeeded, otherwise how many failed and the first error."""
    failed = manifest.loc[manifest.status.ne("success")]
    if failed.empty:
        return ""
    first = failed.iloc[0]
    return f"{what}: {len(failed):,} of {len(manifest):,} tasks failed; first: {first.query_id}: {first.error}"


def _check_manifest(manifest: pd.DataFrame, what: str, path: Path) -> None:
    """Report failed tasks on stderr; exit non-zero when no task succeeded."""
    message = failed_tasks(manifest, what)
    if not message:
        return
    message += f" (errors of all tasks in {path})"
    if manifest.status.ne("success").all():
        raise SystemExit(message)
    print(message, file=sys.stderr, flush=True)


def candidate_table(candidates: pd.DataFrame, counts: pd.DataFrame, final: pd.DataFrame) -> pd.DataFrame:
    """Candidates with a seed CM, their final55 species and loci and the high-copy flag (one row per query_id).

    `counts` comes from `locus_counts` of the seed-panel loci (indexed by query_id) and `final` holds the
    `build_final55` summaries. Species and loci of final55 are those of the seed loci plus accepted expansion loci.
    """
    table = candidates.merge(counts, left_on="query_id", right_index=True, how="left")
    columns = ["n_species_e01", "n_loci_e01", "n_species_hc", "n_loci_hc"]
    table[columns] = table[columns].fillna(0).astype(int)
    table = table.merge(final, on="query_id", how="inner", validate="one_to_one").sort_values("query_id")
    table["n_final55_species"] = table.n_species_hc + table.n_accepted_expansion_species.fillna(0).astype(int)
    table["n_final55_loci"] = table.n_loci_hc + table.n_accepted_expansion_loci.fillna(0).astype(int)
    table["high_copy"] = high_copy(table.n_loci_hc, table.n_final55_loci)
    return table.reset_index(drop=True)


# ---------------------------------------------------------------------------------------------------------------
# Structure-free homologues (input of covariation condition iii)


_HOMOLOGUE_ID = re.compile(r"[^A-Za-z0-9_.:-]+")


def _plain(sequence: str) -> str:
    return sequence.upper().replace("U", "T").replace("-", "").replace(".", "")


def homologue_records(core_id: str, seed_records: Sequence[tuple[str, str]],
                      expansion: pd.DataFrame | None = None) -> list[tuple[str, str]]:
    """The human core and its high-confidence homologues found without the structure (covariation condition iii).

    `seed_records` are the records of the seed `sequences.fa`: the human core (named `core_id`) and one sequence per
    distinct target of the high-confidence seed loci, cut at the Riboseek alignment. `expansion` holds the expansion
    crops (`crop_sequence`, `locus_offset_start`, `locus_offset_end`, `riboseek_evalue`, `riboseek_qcov`,
    `species`, `contig`, `locus_start`, `strand`, in crop order); loci whose Riboseek hit has E <= 1e-3 and query
    coverage >= 0.8 are cut at the merged locus interval in the orientation of the hit (not at the CM hit) and named
    `x{n:06d}_{species}_{contig}_{locus_start}_{strand}`. Sequences are uppercase DNA without gaps, names keep only
    letters, digits and `_.:-`, the human core comes first and exact duplicates (also of the human core) are removed.
    """
    records = [(name, sequence) for name, sequence in seed_records]
    if expansion is not None and len(expansion):
        hc = expansion.loc[expansion.riboseek_evalue.le(HC_EVALUE) & expansion.riboseek_qcov.ge(HC_QCOV)]
        for index, row in enumerate(hc.itertuples(index=False), start=1):
            sequence = str(row.crop_sequence)[int(row.locus_offset_start):int(row.locus_offset_end)]
            records.append((f"x{index:06d}_{row.species}_{row.contig}_{row.locus_start}_{row.strand}", sequence))
    human = [sequence for name, sequence in records if name == core_id]
    if len(human) != 1:
        raise ValueError(f"expected one human row named {core_id}")
    out, seen = [(core_id, _plain(human[0]))], {_plain(human[0])}
    for name, sequence in records:
        sequence = _plain(sequence)
        if name != core_id and sequence not in seen:
            out.append((_HOMOLOGUE_ID.sub("_", name), sequence))
            seen.add(sequence)
    return out


def write_homologues(core_id: str, seed_dir: str | Path, expand_dir: str | Path, output: str | Path) -> int:
    """Write `<core_id>.fa` for covariation condition (iii) from a candidate's `select seed` and `select expand` output.

    Reads `sequences.fa` of the seed directory and `crops.fa` with `expansion_hits.parquet` of the expansion
    directory; returns the number of records written.
    """
    seed_dir, expand_dir = Path(seed_dir), Path(expand_dir)
    expansion = None
    if (expand_dir / "expansion_hits.parquet").exists():
        expansion = pd.read_parquet(expand_dir / "expansion_hits.parquet")
        expansion["crop_sequence"] = expansion.crop_id.map(dict(iter_fasta(expand_dir / "crops.fa")))
    records = homologue_records(core_id, list(iter_fasta(seed_dir / "sequences.fa")), expansion)
    write_fasta(records, output)
    return len(records)


def write_all_homologues(query_ids: Iterable[str], seed_root: str | Path, expand_root: str | Path,
                         out_dir: str | Path) -> pd.DataFrame:
    """`<out_dir>/<core_id>.fa` for every candidate (`*_root/<core_id>/` are the per-candidate directories)."""
    rows = []
    for q in query_ids:
        n = write_homologues(q, Path(seed_root) / q, Path(expand_root) / q, Path(out_dir) / f"{q}.fa")
        rows.append({"query_id": q, "n_sequences": n})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------------
# Provisional families


def candidate_pairs(directed: pd.DataFrame) -> pd.DataFrame:
    """Undirected candidate pairs (unit_a < unit_b) from HMM hits (`source_query_id`, `target_query_id`)."""
    kept = directed.loc[directed.source_query_id.ne(directed.target_query_id)]
    ends = list(zip(kept.source_query_id.astype(str), kept.target_query_id.astype(str)))
    pairs = pd.DataFrame({"unit_a": [min(e) for e in ends], "unit_b": [max(e) for e in ends]})
    return pairs.drop_duplicates().sort_values(["unit_a", "unit_b"]).reset_index(drop=True)


def kmers(sequence: str, k: int = 4) -> set[str]:
    sequence = sequence.upper().replace("U", "T")
    if len(sequence) < k:
        return {sequence}
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 1.0


def sample_heldout(hits: pd.DataFrame, k: int = HELDOUT_PER_CANDIDATE) -> pd.DataFrame:
    """Up to k expansion sequences of one candidate for reciprocal scoring.

    One hit per species (the best by bits), then round-robin over the lineage groups, taking in each group the hit
    most different in 4-mers from those already chosen (then the higher bits, then the species name).
    """
    if hits.empty or k == 0:
        return hits.head(0).copy()
    candidates = hits.sort_values(["cm_bits", "species"], ascending=[False, True]).drop_duplicates("species").copy()
    candidates["_kmers"] = candidates.cm_hit_sequence.map(kmers)
    chosen, chosen_kmers = [], []
    by_group = {group: list(candidates.index[candidates.lineage.eq(group)]) for group in LINEAGES}
    while len(chosen) < min(k, len(candidates)):
        added = False
        for group in LINEAGES:
            remaining = [i for i in by_group[group] if i not in chosen]
            if not remaining:
                continue

            def key(i):
                diversity = (min(1.0 - jaccard(candidates.at[i, "_kmers"], x) for x in chosen_kmers)
                             if chosen_kmers else 1.0)
                return diversity, float(candidates.at[i, "cm_bits"]), str(candidates.at[i, "species"])

            best = max(remaining, key=key)
            chosen.append(best)
            chosen_kmers.append(candidates.at[best, "_kmers"])
            added = True
            if len(chosen) >= k:
                break
        if not added:
            break
    return candidates.loc[chosen].drop(columns="_kmers")


def direction_medians(scores: pd.DataFrame) -> pd.DataFrame:
    """b_{A->B}: median bits of each model on each target candidate's held-out sequences.

    `scores` has one row per (model_query_id, sample) with `target_query_id`, `sample_kind` ('heldout' for
    expansion sequences) and `cm_bits` (best hit of the model on that sequence).
    """
    heldout = scores.loc[scores.sample_kind.eq("heldout")]
    return heldout.groupby(["model_query_id", "target_query_id"]).cm_bits.agg(["median", "size"]).rename(
        columns={"median": "b", "size": "n_heldout"})


def relative_scores(pairs: pd.DataFrame, medians: pd.DataFrame, thresholds: Mapping[str, float]) -> pd.DataFrame:
    """r_{A->B} = (b_{A->B} - T_A) / T_A in both directions and their mean, for pairs scored in both directions."""
    out = pairs[["unit_a", "unit_b"]].copy()
    lookup = medians.b.to_dict()
    out["b_ab"] = [lookup.get((a, b), np.nan) for a, b in zip(out.unit_a, out.unit_b)]
    out["b_ba"] = [lookup.get((b, a), np.nan) for a, b in zip(out.unit_a, out.unit_b)]
    out = out.dropna(subset=["b_ab", "b_ba"]).copy()
    t_a, t_b = out.unit_a.map(thresholds), out.unit_b.map(thresholds)
    out["r_ab"] = (out.b_ab - t_a) / t_a
    out["r_ba"] = (out.b_ba - t_b) / t_b
    out["ratio_mean"] = (out.r_ab + out.r_ba) / 2
    return out.reset_index(drop=True)


def provisional_families(query_ids: Iterable[str], edges: pd.DataFrame,
                         threshold: float = FAMILY_THRESHOLD) -> pd.DataFrame:
    """Connected components of pairs with `ratio_mean` >= threshold (singletons included).

    Families are numbered C01, C02, ... by decreasing size, then by their smallest query_id; `rlmfam_id` is the
    published identifier of the family (RLMFam00001, RLMFam00002, ...; the paper refers to candidates by it).
    """
    parent = {q: q for q in query_ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    kept = edges.loc[edges.ratio_mean >= threshold]
    for a, b in zip(kept.unit_a, kept.unit_b):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[max(root_a, root_b)] = min(root_a, root_b)
    groups: dict[str, list[str]] = defaultdict(list)
    for q in parent:
        groups[find(q)].append(q)
    ordered = sorted(groups.values(), key=lambda members: (-len(members), min(members)))
    rows = [(q, f"C{rank:02d}", f"RLMFam{rank:05d}", len(members), rank)
            for rank, members in enumerate(ordered, start=1) for q in members]
    return pd.DataFrame(rows, columns=["query_id", "family_id", "rlmfam_id", "family_size", "family_rank"]).sort_values(
        ["family_rank", "query_id"]).reset_index(drop=True)


def _hmm_job(query_id: str, cm: str, targets: str, out_dir: str, infernal_bin) -> dict[str, object]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    run_logged(hmmsearch_command(tool("cmsearch", infernal_bin), cm, targets, out / "hmm.tbl"), out / "hmm.log")
    return {"query_id": query_id}


def _rescore_job(query_id: str, cm: str, fasta: str, out_dir: str, infernal_bin) -> dict[str, object]:
    out = Path(out_dir)
    run_logged(cmsearch_command(tool("cmsearch", infernal_bin), cm, fasta, out / "rescore.tbl"), out / "rescore.log")
    return {"query_id": query_id}


def score_families(candidates: pd.DataFrame, expand_dir: str | Path, out_dir: str | Path,
                   infernal_bin: str | Path | None = None, workers: int = 15) -> tuple[pd.DataFrame, pd.DataFrame]:
    """HMM pairs, held-out samples, reciprocal CM scores and provisional families; returns (families, edges).

    `candidates` is the table of `candidate_table` (query_id, sequence_u, threshold_bits); `expand_dir` holds the
    per-candidate output of `build_final55`. When the HMM search proposes no pair, every candidate is its own family
    and `edges` is empty.
    """
    if candidates.empty:
        raise ValueError("no candidates")
    expand_dir, out = Path(expand_dir), Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    candidates = candidates.sort_values("query_id")
    write_fasta(zip(candidates.query_id, candidates.sequence_u), out / "human_cores.fa", dna=False)
    cms = {q: str(expand_dir / q / "seed20_3l.calibrated.cm") for q in candidates.query_id}
    manifest = _run_parallel([(_hmm_job, {"query_id": q, "cm": cms[q], "targets": str(out / "human_cores.fa"),
                                          "out_dir": str(out / "pairs" / q), "infernal_bin": infernal_bin})
                              for q in candidates.query_id], workers)
    failure = failed_tasks(manifest, "HMM searches")
    if failure:
        raise RuntimeError(f"{failure} (logs under {out / 'pairs'})")
    directed = pd.concat([parse_tblout(out / "pairs" / q / "hmm.tbl").assign(source_query_id=q)
                          for q in candidates.query_id], ignore_index=True)
    pairs = candidate_pairs(directed.rename(columns={"target": "target_query_id"}))
    pairs.to_parquet(out / "candidate_pairs.parquet", index=False)
    if pairs.empty:
        return provisional_families(candidates.query_id, _empty_frame(EDGE_DTYPES)), _empty_frame(EDGE_DTYPES)

    units = sorted(set(pairs.unit_a) | set(pairs.unit_b))
    samples = []
    sequences = candidates.set_index("query_id").sequence_u
    for q in units:
        samples.append({"target_query_id": q, "sample_kind": "human", "species": "homo_sapiens",
                        "sequence": sequences[q].upper().replace("T", "U")})
        path = expand_dir / q / "accepted_expansion_hits.parquet"
        if path.exists():
            for hit in sample_heldout(pd.read_parquet(path)).itertuples(index=False):
                samples.append({"target_query_id": q, "sample_kind": "heldout", "species": hit.species,
                                "sequence": hit.cm_hit_sequence.upper().replace("T", "U")})
    samples = pd.DataFrame(samples)
    samples["sample_id"] = [f"s{i:08d}" for i in range(1, len(samples) + 1)]
    samples.to_parquet(out / "samples.parquet", index=False)

    neighbours: dict[str, set[str]] = defaultdict(set)
    for a, b in zip(pairs.unit_a, pairs.unit_b):
        neighbours[a].update((a, b))
        neighbours[b].update((a, b))
    tasks = []
    for q, targets in sorted(neighbours.items()):
        directory = out / "rescore" / q
        directory.mkdir(parents=True, exist_ok=True)
        subset = samples.loc[samples.target_query_id.isin(targets)]
        write_fasta(zip(subset.sample_id, subset.sequence), directory / "targets.fa", dna=False)
        tasks.append((_rescore_job, {"query_id": q, "cm": cms[q], "fasta": str(directory / "targets.fa"),
                                     "out_dir": str(directory), "infernal_bin": infernal_bin}))
    manifest = _run_parallel(tasks, workers)
    failure = failed_tasks(manifest, "reciprocal CM searches")
    if failure:
        raise RuntimeError(f"{failure} (logs under {out / 'rescore'})")
    lookup = samples.drop(columns="sequence")
    scores = []
    for q in neighbours:
        hits = parse_tblout(out / "rescore" / q / "rescore.tbl").rename(columns={"target": "sample_id"})
        hits = hits.sort_values("cm_bits", ascending=False).drop_duplicates("sample_id")
        scores.append(hits[["sample_id", "cm_bits"]].merge(lookup, on="sample_id").assign(model_query_id=q))
    scores = pd.concat(scores, ignore_index=True)
    scores.to_parquet(out / "reciprocal_scores.parquet", index=False)
    thresholds = candidates.set_index("query_id").threshold_bits.to_dict()
    edges = relative_scores(pairs, direction_medians(scores), thresholds)
    families = provisional_families(candidates.query_id, edges)
    return families, edges


# ---------------------------------------------------------------------------------------------------------------
# Command line: rlmhelix select <step>


def _cmd_cores(args) -> None:
    from .student import Student

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    windows = read_table(args.windows)
    if windows.empty:
        raise SystemExit(f"no windows in {args.windows}")
    top = top_windows(windows, args.fraction, args.p_column)
    write_table(top, out / "top_windows.parquet")
    loci = merge_loci(top, args.p_column, args.score_column)
    write_table(loci, out / "loci.parquet")
    print(f"windows={len(windows):,} top={len(top):,} loci={len(loci):,}", flush=True)
    student = Student(device=args.device, dtype=args.dtype)
    cores = extract_cores(loci, student, args.batch_size, args.n_jobs)
    write_table(cores, out / "cores.parquet")
    if cores.empty:
        raise SystemExit(f"no core: the peak windows of the {len(loci):,} loci decode no stem of >= {MIN_STEM} pairs")
    representatives, membership = collapse_cores(cores, loci, args.id_prefix)
    write_table(representatives, out / "representatives.parquet")
    write_table(membership, out / "membership.parquet")
    write_fasta(zip(representatives.query_id, representatives.canonical_sequence), out / "representatives_canonical.fa")
    write_fasta(zip(representatives.query_id, representatives.sequence_u), out / "representatives_source.fa")
    print(f"cores={len(cores):,} primary={len(membership):,} representatives={len(representatives):,}", flush=True)


def _cmd_chunk(args) -> None:
    for genome in args.genomes:
        species = args.species or Path(genome).name.split(".")[0].lower()
        print(chunk_genome(genome, species, args.output, args.chunk_size, args.overlap), flush=True)


def _species_list(args) -> list[str]:
    if args.species:
        return list(args.species)
    panel = load_panel(args.panel)
    return seed_species(panel) if args.stage == "seed" else expansion_species(panel)


def _riboseek_binary(path: str | None) -> str:
    """The Riboseek binary given with --riboseek (a path or a name on PATH), or `riboseek` from PATH."""
    if path is None:
        return require_tools(["riboseek"])["riboseek"]
    found = shutil.which(path)
    if found is None:
        raise SystemExit(f"riboseek binary {path} not found or not executable")
    return found


def _cmd_search(args) -> None:
    binary = _riboseek_binary(args.riboseek)
    species = _species_list(args)
    chunks = {s: Path(args.chunks) / f"{s}.10kb_o500.fa.gz" for s in species}
    missing = [str(p) for p in chunks.values() if not p.exists()]
    if missing:
        raise FileNotFoundError(f"missing chunk files: {missing[:3]}")
    try:
        status = run_riboseek(args.queries, chunks, args.output, binary, args.threads, force=args.force)
    except FileExistsError as exc:
        raise SystemExit(str(exc)) from None
    status.to_csv(Path(args.output) / "search_status.tsv", sep="\t", index=False)


def _cmd_routes(args) -> None:
    panel = load_panel(args.panel)
    representatives = read_table(args.representatives)
    frames = []
    for species in seed_species(panel):
        found = high_confidence_queries(Path(args.hits) / species / "hits.tsv")
        frames.append(pd.DataFrame({"query_id": sorted(found), "species": species}))
        print(f"{species}: {len(found):,} queries with a high-confidence hit", flush=True)
    counts = lineage_counts(pd.concat(frames, ignore_index=True), lineage_map(panel))
    routed = assign_routes(representatives, counts)
    selected = routed.loc[routed.selected].reset_index(drop=True)
    write_table(selected, args.output)
    write_fasta(zip(selected.query_id, selected.sequence_u), Path(args.output).with_suffix(".source.fa"))
    if selected.empty:
        raise SystemExit(f"no candidate of {len(representatives):,} meets a route")
    print(selected.route.value_counts().to_string(), flush=True)


def _seed_loci(hits_root: str | Path, species: list[str], query_ids: set[str]) -> pd.DataFrame:
    """Hit loci (E <= 0.1) of the queries in each species; species without loci add no rows."""
    frames = []
    for name in species:
        loci = hit_loci(Path(hits_root) / name / "hits.tsv", name, query_ids)
        frames.append(loci)
        print(f"{name}: {len(loci):,} loci", flush=True)
    found = [frame for frame in frames if len(frame)]
    return pd.concat(found, ignore_index=True) if found else frames[0]


def _cmd_seed(args) -> None:
    require_tools(INFERNAL_PROGRAMS, args.infernal_bin)
    panel = load_panel(args.panel)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    selected = read_table(args.selected)
    if selected.empty:
        raise SystemExit(f"no candidates in {args.selected}")
    loci = _seed_loci(args.hits, seed_species(panel), set(selected.query_id))
    write_table(loci.drop(columns=["qaln", "taln"]), out / "seed_loci.parquet")
    write_table(locus_counts(loci).reset_index(), out / "seed_locus_counts.parquet")
    hc = loci.loc[loci.any_high_confidence.astype(bool)]
    if hc.empty:
        raise SystemExit(f"no high-confidence locus (E <= {HC_EVALUE:g}, query coverage >= {HC_QCOV}) of the "
                         f"{len(selected):,} candidates in the seed panel ({args.hits}); no seed CM can be built")
    crops = [extract_crops(genome_path(args.genomes, s), part, "seed") for s, part in hc.groupby("species", sort=True)]
    crops = pd.concat(crops, ignore_index=True)
    hc_groups, crop_groups = dict(tuple(hc.groupby("query_id"))), dict(tuple(crops.groupby("query_id", sort=True)))
    tasks = [(build_seed, {"query_id": row.query_id, "query_sequence": row.sequence_u,
                           "structure": row.core_dot_bracket, "hc_loci": hc_groups.get(row.query_id, hc.head(0)),
                           "crops": crop_groups.get(row.query_id, crops.head(0)),
                           "out_dir": str(out / "candidates" / row.query_id), "infernal_bin": args.infernal_bin})
             for row in selected.itertuples(index=False)]
    manifest = _run_parallel(tasks, args.workers)
    write_table(manifest, out / "seed_manifest.parquet")
    _check_manifest(manifest, "select seed", out / "seed_manifest.parquet")
    print(f"seed CMs: {int(manifest.get('seed_cm', pd.Series(dtype=bool)).fillna(False).sum()):,} "
          f"of {len(manifest):,} candidates", flush=True)


def _cmd_expand(args) -> None:
    require_tools(INFERNAL_PROGRAMS, args.infernal_bin)
    panel = load_panel(args.panel)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    selected = read_table(args.selected)
    seed = read_table(Path(args.seed) / "seed_manifest.parquet")
    seed_cm = seed["seed_cm"].fillna(False).astype(bool) if "seed_cm" in seed else pd.Series(False, index=seed.index)
    built = set(seed.loc[seed.status.eq("success") & seed_cm, "query_id"])
    selected = selected.loc[selected.query_id.isin(built)]
    if selected.empty:
        raise SystemExit(f"no candidate with a seed CM in {Path(args.seed) / 'seed_manifest.parquet'}")
    loci = _seed_loci(args.hits, expansion_species(panel), set(selected.query_id))
    loci["model_length"] = loci.query_id.map(selected.set_index("query_id").core_length)
    crops = [extract_crops(genome_path(args.genomes, s), part, "expansion") for s, part in loci.groupby("species", sort=True)]
    crops = pd.concat(crops, ignore_index=True) if crops else pd.DataFrame(columns=["query_id"])
    crop_groups = dict(tuple(crops.groupby("query_id", sort=True)))
    lineage = lineage_map(panel)
    tasks = [(build_final55, {"query_id": q, "crops": crop_groups.get(q, crops.head(0)),
                              "seed_dir": str(Path(args.seed) / "candidates" / q),
                              "out_dir": str(out / "candidates" / q), "lineage": lineage,
                              "infernal_bin": args.infernal_bin}) for q in sorted(selected.query_id)]
    final = _run_parallel(tasks, args.workers)
    write_table(final, out / "final55_manifest.parquet")
    _check_manifest(final, "select expand", out / "final55_manifest.parquet")
    counts = read_table(Path(args.seed) / "seed_locus_counts.parquet").set_index("query_id")
    thresholds = seed.set_index("query_id")[["threshold_bits"]]
    table = candidate_table(selected.merge(thresholds, left_on="query_id", right_index=True),
                            counts, final.loc[final.status.eq("success")].drop(columns=["status", "error"]))
    table["msa_path"] = [str(out / "candidates" / q / "final55.sto") for q in table.query_id]
    table["cm_path"] = [str(out / "candidates" / q / "final55.cm") for q in table.query_id]
    write_table(table, out / "candidates.parquet")
    write_all_homologues(table.query_id, Path(args.seed) / "candidates", out / "candidates", out / "homologues")
    print(f"final55 candidates={len(table):,} median species={table.n_final55_species.median():g} "
          f"high-copy={int(table.high_copy.sum()):,}; homologues in {out / 'homologues'}", flush=True)


def _cmd_homologues(args) -> None:
    candidates = read_table(args.candidates)
    if candidates.empty:
        raise SystemExit(f"no candidates in {args.candidates}")
    counts = write_all_homologues(candidates.query_id, Path(args.seed) / "candidates",
                                  Path(args.expand) / "candidates", args.output)
    print(f"wrote {len(counts):,} FASTA files to {args.output} (median {counts.n_sequences.median():g} sequences)",
          flush=True)


def _cmd_families(args) -> None:
    require_tools(["cmsearch"], args.infernal_bin)
    candidates = read_table(args.candidates)
    if candidates.empty:
        raise SystemExit(f"no candidates in {args.candidates}")
    families, edges = score_families(candidates, Path(args.expand) / "candidates", args.output, args.infernal_bin,
                                     args.workers)
    write_table(edges, Path(args.output) / "edges.parquet")
    write_table(families, Path(args.output) / "families.tsv")
    if edges.empty:
        print("no candidate pair was scored (see candidate_pairs.parquet): every candidate is its own family",
              flush=True)
    print(f"families={families.family_id.nunique():,} (threshold {FAMILY_THRESHOLD})", flush=True)


PANEL_HELP = ("genome panel TSV (default: the 55-species panel shipped with the package, a copy of "
              "data/genome_panel_55.tsv)")


def register(subparsers) -> None:
    """Add `rlmhelix select <step>` (cores, chunk, search, routes, seed, expand, homologues, families) to a CLI."""
    parser = subparsers.add_parser("select", help="candidate selection from the scan output (see docs/pipeline_select.md)")
    steps = parser.add_subparsers(dest="step", required=True)

    cores = steps.add_parser("cores", help="top 1%% windows -> loci -> structural cores -> representatives (GPU)")
    cores.add_argument("windows", help="eligible scan windows (parquet/TSV) with sequence, coordinates, score, P")
    cores.add_argument("-o", "--output", required=True, help="output directory")
    cores.add_argument("--fraction", type=float, default=TOP_FRACTION)
    cores.add_argument("--p-column", help="calibrated P value column (default: p_value or p_min4_mixed)")
    cores.add_argument("--score-column", help="score column (default: score_min4 or score)")
    cores.add_argument("--id-prefix", default=ID_PREFIX, help="prefix of representative IDs")
    cores.add_argument("--device", help="torch device (default: cuda if available)")
    cores.add_argument("--dtype", choices=("bf16", "fp16", "fp32"), default="bf16")
    cores.add_argument("--batch-size", type=int, default=16)
    cores.add_argument("--n-jobs", type=int, default=12, help="CPU workers for RNAfold and the DP")
    cores.set_defaults(func=_cmd_cores)

    chunk = steps.add_parser("chunk", help="cut genomes into 10-kb chunks overlapping by 500 nt")
    chunk.add_argument("genomes", nargs="+")
    chunk.add_argument("--species", help="Ensembl slug (default: from the file name; one genome only)")
    chunk.add_argument("-o", "--output", required=True)
    chunk.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    chunk.add_argument("--overlap", type=int, default=CHUNK_OVERLAP)
    chunk.set_defaults(func=_cmd_chunk)

    search = steps.add_parser("search", help="Riboseek search of queries against chunked genomes")
    search.add_argument("queries", help="query FASTA")
    search.add_argument("--chunks", required=True, help="directory of <species>.10kb_o500.fa.gz")
    search.add_argument("--stage", choices=("seed", "expansion"), default="seed", help="species set of the panel")
    search.add_argument("--species", nargs="*", help="explicit species (overrides --stage)")
    search.add_argument("--panel", help=PANEL_HELP)
    search.add_argument("--riboseek", help="riboseek binary (default: from PATH)")
    search.add_argument("--threads", type=int, default=20)
    search.add_argument("--force", action="store_true",
                        help="rebuild the query database and search every species again (removes the done markers); "
                             "needed when the work directory holds the search of other queries")
    search.add_argument("-o", "--output", required=True, help="work directory (<species>/hits.tsv)")
    search.set_defaults(func=_cmd_search)

    routes = steps.add_parser("routes", help="high-confidence hits -> lineage counts -> four routes")
    routes.add_argument("representatives", help="representatives.parquet of 'select cores'")
    routes.add_argument("--hits", required=True, help="work directory of the seed-panel search")
    routes.add_argument("--panel", help=PANEL_HELP)
    routes.add_argument("-o", "--output", required=True, help="selected candidates (parquet/TSV)")
    routes.set_defaults(func=_cmd_routes)

    for name, text in (("seed", "seed alignments and CMs from the seed-panel search (Infernal)"),
                       ("expand", "expansion genomes -> inclusion rule -> final55 alignments and CMs (Infernal)")):
        step = steps.add_parser(name, help=text)
        step.add_argument("selected", help="selected candidates of 'select routes'")
        if name == "expand":
            step.add_argument("--seed", required=True, help="output directory of 'select seed'")
        step.add_argument("--hits", required=True, help=f"work directory of the source-strand {name}-panel search")
        step.add_argument("--genomes", required=True, help="directory of genome FASTA files (see genome_path)")
        step.add_argument("--panel", help=PANEL_HELP)
        step.add_argument("--infernal-bin", help="directory of the Infernal programs (default: from PATH)")
        step.add_argument("--workers", type=int, default=20)
        step.add_argument("-o", "--output", required=True)
        step.set_defaults(func=_cmd_seed if name == "seed" else _cmd_expand)

    homologues = steps.add_parser(
        "homologues", help="structure-free homologues per candidate (input of covariation condition iii)")
    homologues.add_argument("candidates", help="candidates.parquet of 'select expand'")
    homologues.add_argument("--seed", required=True, help="output directory of 'select seed'")
    homologues.add_argument("--expand", required=True, help="output directory of 'select expand'")
    homologues.add_argument("-o", "--output", required=True, help="directory for <core_id>.fa")
    homologues.set_defaults(func=_cmd_homologues)

    families = steps.add_parser("families", help="reciprocal CM scores -> provisional families (Infernal)")
    families.add_argument("candidates", help="candidates.parquet of 'select expand'")
    families.add_argument("--expand", required=True, help="output directory of 'select expand'")
    families.add_argument("--infernal-bin", help="directory of the Infernal programs (default: from PATH)")
    families.add_argument("--workers", type=int, default=15)
    families.add_argument("-o", "--output", required=True)
    families.set_defaults(func=_cmd_families)


def main(argv: list[str] | None = None) -> None:
    """Stand-alone entry point: `python -m rlmhelix.select <step> ...`."""
    import argparse

    parser = argparse.ArgumentParser(prog="rlmhelix")
    register(parser.add_subparsers(dest="command", required=True))
    args = parser.parse_args(["select", *(argv if argv is not None else sys.argv[1:])])
    args.func(args)


if __name__ == "__main__":
    main()
