"""Unit tests for findORF (DIMPLE/utilities.py)."""

from __future__ import annotations

import unittest
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from DIMPLE.utilities import findORF


class TestFindORF(unittest.TestCase):
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

    def test_find_orf_non_interactive_single_candidate(self) -> None:
        gene = SeqRecord(
            Seq("ATG" + "GCC" * 100 + "TAA"),
            id="orf_gene_single",
            name="orf_gene_single",
            description="orf_gene_single",
        )
        start, end = findORF(gene, non_interactive=True)
        self.assertEqual(start, 0)
        self.assertEqual(end, 303)

    def test_find_orf_non_interactive_requires_explicit_index_when_ambiguous(
        self,
    ) -> None:
        gene = SeqRecord(
            Seq("ATG" + "GCC" * 101 + "TAA" + "ATG" + "GCC" * 101 + "TAA"),
            id="orf_gene_multi",
            name="orf_gene_multi",
            description="orf_gene_multi",
        )
        with self.assertRaises(ValueError) as ctx:
            findORF(gene, non_interactive=True)
        self.assertIn("multiple orf candidates", str(ctx.exception).lower())


if __name__ == "__main__":
    unittest.main()
