"""CLI parser tests for ``run_dimple.py``.

These stay at parser level so the command-line contract is covered without
running the full design pipeline.
"""

from run_dimple import build_parser


def test_dis_is_boolean_flag():
    args = build_parser().parse_args(["-geneFile", "gene.fa", "-dis"])
    assert args.dis is True


def test_dis_defaults_false():
    args = build_parser().parse_args(["-geneFile", "gene.fa"])
    assert args.dis is False
