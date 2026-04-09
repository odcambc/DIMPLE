"""Unit tests for :mod:`DIMPLE.utilities`."""

from __future__ import annotations

import unittest
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from DIMPLE.utilities import findORF, parse_custom_mutations


class TestUtilities(unittest.TestCase):
    def test_parse_custom_mutations_merges_duplicate_positions(self) -> None:
        custom = parse_custom_mutations(["5:A", "5:G"])
        self.assertEqual(custom[5], "A,G")

    def test_parse_custom_mutations_expands_ranges(self) -> None:
        custom = parse_custom_mutations(["2-4:A", "7:W"])
        self.assertEqual(custom[2], "A")
        self.assertEqual(custom[3], "A")
        self.assertEqual(custom[4], "A")
        self.assertEqual(custom[7], "W")

    def test_parse_custom_mutations_all_keyword(self) -> None:
        all_aas = "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y"
        custom = parse_custom_mutations(["2-3:All", "9:All"])
        self.assertEqual(custom[2], all_aas)
        self.assertEqual(custom[3], all_aas)
        self.assertEqual(custom[9], all_aas)

    def test_find_orf_with_mocked_prompts(self) -> None:
        gene = SeqRecord(
            Seq("ATG" + "GCC" * 100 + "TAA"),
            id="orf_gene",
            name="orf_gene",
            description="orf_gene",
        )
        with patch("builtins.input", side_effect=["1", "y", "y"]):
            start, end = findORF(gene)
        self.assertEqual(start, 0)
        self.assertEqual(end, 303)


if __name__ == "__main__":
    unittest.main()
