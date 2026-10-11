"""Command-line interface: every sub-command answers --help, and `rlmhelix score` uses the column names of the scan."""
import argparse

import numpy as np
import pytest

from rlmhelix import cli


def sub_commands(parser, prefix=()):
    """Every command path of the parser, e.g. ('select', 'homologues')."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            for name, child in action.choices.items():
                yield (*prefix, name)
                yield from sub_commands(child, (*prefix, name))


def test_every_command_has_help(capsys):
    commands = [(), *sub_commands(cli.build_parser())]
    assert {("score",), ("msa",), ("scan",), ("select", "homologues"), ("catalogue",)} <= set(commands)
    for command in commands:
        with pytest.raises(SystemExit) as stop:
            cli.main([*command, "--help"])
        assert stop.value.code == 0, command
        assert capsys.readouterr().out.startswith("usage: rlmhelix"), command


def test_score_columns_match_the_scan():
    pytest.importorskip("RNA")
    from rlmhelix.score import rlmhelix_score

    rng = np.random.default_rng(0)
    sequence = "".join(rng.choice(list("ACGU"), 80))
    apc = np.abs(rng.normal(size=(80, 80))).astype(np.float32)
    apc = (apc + apc.T) / 2
    row = cli._score_one("x", sequence, apc)
    result = rlmhelix_score(apc, sequence)
    assert list(row) == ["id", "length", "score", "n_pairs", "n_stems", "n_stems_scored", "sequence", "structure"]
    assert (row["n_pairs"], row["n_stems"], row["n_stems_scored"]) == (result.n_pairs, result.n_stems,
                                                                        result.n_stems_scored)
    assert row["n_stems"] >= row["n_stems_scored"] and row["n_pairs"] >= row["structure"].count("(")
