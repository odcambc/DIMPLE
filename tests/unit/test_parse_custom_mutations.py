"""Unit tests for parse_custom_mutations (DIMPLE/utilities.py:52-72)."""

import pytest

from DIMPLE.utilities import parse_custom_mutations

_ALL = "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y"


@pytest.mark.parametrize(
    "mutation_text, expected",
    [
        # single position, single amino acid
        (["10:A"], {10: "A"}),
        # single position, All expands to the full 20-AA string
        (["10:All"], {10: _ALL}),
        # range, single amino acid — expands to one entry per position
        (["5-7:K"], {5: "K", 6: "K", 7: "K"}),
        # duplicate position: second entry is appended with comma
        (["10:A", "10:G"], {10: "A,G"}),
        # mixed: range+All followed by a single position
        (
            ["1-3:All", "5:M"],
            {1: _ALL, 2: _ALL, 3: _ALL, 5: "M"},
        ),
        # overlapping single-then-range: the range must accumulate onto the
        # single, not overwrite it (previously dropped the C at position 3).
        (
            ["3:C", "1-5:A"],
            {1: "A", 2: "A", 3: "C,A", 4: "A", 5: "A"},
        ),
        # overlapping range-then-single: symmetric — no residue lost either way.
        (
            ["1-5:A", "3:C"],
            {1: "A", 2: "A", 3: "A,C", 4: "A", 5: "A"},
        ),
        # duplicate residue across overlapping lines is de-duped, not doubled.
        (
            ["10:A", "8-12:A"],
            {8: "A", 9: "A", 10: "A", 11: "A", 12: "A"},
        ),
    ],
)
def test_parse_custom_mutations(mutation_text, expected):
    assert parse_custom_mutations(mutation_text) == expected


def test_parse_custom_mutations_skips_header_and_blank_lines():
    assert parse_custom_mutations(["Positions:Mutations\n", "\n", "1:All\n"]) == {1: _ALL}
