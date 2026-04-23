"""Unit tests for addgene (DIMPLE/DIMPLE.py:44-67)."""

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import DIMPLE, addgene

# Minimal gene FASTA body: 200 bases, no restriction-site sequences.
_GENE_SEQ = "ATCG" * 50  # 200 bp, 50% GC, no BsmBI/BsaI sites


def _write_fasta(path, records):
    """Write a simple FASTA file; records is a list of (header, seq) tuples."""
    text = "".join(f">{header}\n{seq}\n" for header, seq in records)
    path.write_text(text)


def _minimal_state(dimple_state):
    """Configure just enough class state to create DIMPLE instances."""
    DIMPLE.avoid_sequence = []
    DIMPLE.stop_codon = True
    DIMPLE.synth_len = 230
    DIMPLE.maxfrag = 165
    DIMPLE.enzyme = None


class TestAddgeneStartEnd:
    def test_start_is_zero_indexed(self, tmp_path, dimple_state):
        """1-indexed start: in FASTA header is converted to 0-indexed."""
        _minimal_state(dimple_state)
        fa = tmp_path / "gene.fa"
        _write_fasta(fa, [("gene1 start:31 end:120", _GENE_SEQ)])
        OLS = addgene(str(fa))
        assert OLS[0].start == 30  # 31 - 1

    def test_end_is_unchanged(self, tmp_path, dimple_state):
        """end: in FASTA header is stored as-is (no offset)."""
        _minimal_state(dimple_state)
        fa = tmp_path / "gene.fa"
        _write_fasta(fa, [("gene1 start:31 end:120", _GENE_SEQ)])
        OLS = addgene(str(fa))
        assert OLS[0].end == 120

    def test_multi_record_returns_correct_count(self, tmp_path, dimple_state):
        """Multi-record FASTA produces one DIMPLE object per record."""
        _minimal_state(dimple_state)
        fa = tmp_path / "genes.fa"
        _write_fasta(
            fa,
            [
                ("geneA start:31 end:120", _GENE_SEQ),
                ("geneB start:31 end:120", _GENE_SEQ),
            ],
        )
        OLS = addgene(str(fa))
        assert len(OLS) == 2

    def test_geneid_matches_fasta_name(self, tmp_path, dimple_state):
        """DIMPLE.geneid is set from the FASTA record name field."""
        _minimal_state(dimple_state)
        fa = tmp_path / "gene.fa"
        _write_fasta(fa, [("Kir start:31 end:120", _GENE_SEQ)])
        OLS = addgene(str(fa))
        assert OLS[0].geneid == "Kir"
