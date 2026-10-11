"""Command-line interface.

    rlmhelix score sequences.fa -o scores.tsv                 # score, calibrated P value and structure
    rlmhelix score sequences.fa --shuffles 100 -o scores.tsv  # also compare each sequence with its own shuffles
    rlmhelix msa alignment.sto [...] -o msa.tsv               # RLMHelix-MSA per alignment
    rlmhelix scan prepare|score|fdr ...                       # genome scan (docs/pipeline_scan.md)
    rlmhelix select cores|chunk|search|routes|seed|expand|homologues|families ...
                                                              # candidate selection (docs/pipeline_select.md)
    rlmhelix catalogue candidates|icshape|covariation|build|concordance ...   # evidence and catalogue (catalogue/README.md)
"""
from __future__ import annotations

import argparse
import sys
from importlib import import_module
from pathlib import Path

import pandas as pd


PIPELINES = ("scan", "select", "catalogue")


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    records, name, chunks = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name, chunks = line[1:].split()[0] if line[1:].split() else "", []
        elif line.strip():
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def _score_one(name, sequence, apc):
    """One output row; n_pairs, n_stems and n_stems_scored mean what they mean in the `scan score` output."""
    from .score import rlmhelix_score, to_rna

    result = rlmhelix_score(apc, sequence)
    return {"id": name, "length": len(sequence), "score": result.score, "n_pairs": result.n_pairs,
            "n_stems": result.n_stems, "n_stems_scored": result.n_stems_scored, "sequence": to_rna(sequence),
            "structure": result.structure}


def score_command(args: argparse.Namespace) -> pd.DataFrame:
    from joblib import Parallel, delayed

    from .calibrate import Calibrator
    from .student import Student

    records = read_fasta(args.fasta)
    student = Student(device=args.device, dtype=args.dtype)
    rows = []
    for start in range(0, len(records), args.chunk):
        chunk = records[start:start + args.chunk]
        maps = student.predict([s for _, s in chunk], args.batch_size)
        rows.extend(Parallel(n_jobs=args.n_jobs)(delayed(_score_one)(n, s, m) for (n, s), m in zip(chunk, maps)))
    frame = pd.DataFrame(rows)
    frame.insert(3, "p_value", Calibrator.load().pvalues(frame.score, list(frame.sequence)))
    if args.shuffles:
        frame = pd.concat([frame, shuffle_columns(frame, student, args)], axis=1)
    return frame


def shuffle_columns(frame: pd.DataFrame, student, args: argparse.Namespace) -> pd.DataFrame:
    """Compare each sequence with its own dinucleotide-preserving shuffles (as in the paper's benchmarks)."""
    import numpy as np
    from joblib import Parallel, delayed

    from .shuffle import dinucleotide_shuffles

    rows = []
    for k, sequence in enumerate(frame.sequence):
        shuffles = dinucleotide_shuffles(sequence, args.shuffles, seed=args.seed + k)
        maps = student.predict(shuffles, args.batch_size)
        null = np.asarray([r["score"] for r in Parallel(n_jobs=args.n_jobs)(
            delayed(_score_one)("", s, m) for s, m in zip(shuffles, maps))])
        score = frame.score.iloc[k]
        sd = null.std(ddof=1)
        rows.append({"shuffle_mean": null.mean(), "shuffle_sd": sd,
                     "z": (score - null.mean()) / sd if sd > 0 else np.nan,
                     "empirical_p": (1 + int((null >= score).sum())) / (1 + len(null))})
    return pd.DataFrame(rows, index=frame.index)


def msa_command(args: argparse.Namespace) -> pd.DataFrame:
    from .msa import read_stockholm, score_alignment
    from .student import Student

    student = Student(device=args.device, dtype=args.dtype)
    rows = []
    for path in args.alignments:
        result = score_alignment(read_stockholm(path), student, args.batch_size)
        rows.append({"alignment": str(path), **result})
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    """The `rlmhelix` argument parser with every sub-command."""
    parser = argparse.ArgumentParser(
        prog="rlmhelix", description="RLMHelix scores of RNA sequences and alignments")
    sub = parser.add_subparsers(dest="command", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--output", help="output TSV (default: standard output)")
    common.add_argument("--device", help="torch device (default: cuda if available, else cpu)")
    common.add_argument("--dtype", choices=("bf16", "fp16", "fp32"), default="bf16", help="autocast type on GPU")
    common.add_argument("--batch-size", type=int, default=16, help="sequences per language-model batch")
    score = sub.add_parser("score", parents=[common], help="score single sequences (FASTA, 30-300 nt)")
    score.add_argument("fasta")
    score.add_argument("--n-jobs", type=int, default=4, help="CPU workers for RNAfold and the DP")
    score.add_argument("--chunk", type=int, default=1024, help="sequences held in memory at once")
    score.add_argument("--shuffles", type=int, default=0,
                       help="number of dinucleotide-preserving shuffles per sequence for an empirical P value and z")
    score.add_argument("--seed", type=int, default=0, help="random seed of the shuffles")
    msa = sub.add_parser("msa", parents=[common], help="score alignments (Stockholm) with RLMHelix-MSA")
    msa.add_argument("alignments", nargs="+")
    for name in PIPELINES:  # the pipelines of the paper add their own sub-commands, each setting `func`
        import_module(f"rlmhelix.{name}").register(sub)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command in PIPELINES:
        args.func(args)
        return
    frame = score_command(args) if args.command == "score" else msa_command(args)
    frame.to_csv(args.output or sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
