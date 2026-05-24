"""Focused tests for post-QC behavior."""

from __future__ import annotations

import unittest
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from DIMPLE.DIMPLE import DIMPLE, post_qc, check_final_assembly
from DIMPLE.pool import Pool, DimpleRuntimeConfig


class TestQCPaths(unittest.TestCase):
    def setUp(self) -> None:
        self.config = DimpleRuntimeConfig()
        self.pool = Pool(self.config)
        self.gene = DIMPLE.__new__(DIMPLE)
        self.gene.geneid = "dummy_gene"
        self.gene.oligos = [
            SeqRecord(Seq("ATGCGTATGCGT"), id="dummy_fragment_000001", description="")
        ]
        self.gene.barPrimer = []
        self.gene.genePrimer = []
        self.gene.pool = self.pool
        self.pool.append(self.gene)

    def test_post_qc_calls_assembly_check_when_enzyme_is_set(self) -> None:
        self.config.enzyme = "BsaI"
        with patch("DIMPLE.qc.check_final_assembly") as mocked:
            post_qc(self.pool)
        mocked.assert_called_once_with(self.gene)

    def test_post_qc_skips_assembly_check_when_enzyme_is_none(self) -> None:
        self.config.enzyme = None
        with patch("DIMPLE.qc.check_final_assembly") as mocked:
            post_qc(self.pool)
        mocked.assert_not_called()

    def test_check_final_assembly_returns_none_when_enzyme_missing(self) -> None:
        self.config.enzyme = None
        result = check_final_assembly(self.gene)
        self.assertIsNone(result)

    def test_check_final_assembly_returns_none_for_unknown_enzyme(self) -> None:
        self.config.enzyme = "UnknownEnzyme"
        result = check_final_assembly(self.gene)
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
