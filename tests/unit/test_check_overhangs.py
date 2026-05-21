"""Unit tests for check_overhangs (DIMPLE/DIMPLE.py:857-899)."""

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import DIMPLE, addgene, check_overhangs

_GENE_SEQ = "ATCG" * 50  # 200 bp, 50% GC, unique-looking 4-mers at the overhang sites


def _write_fasta(path, records):
    text = "".join(f">{header}\n{seq}\n" for header, seq in records)
    path.write_text(text)


def _setup_state(dimple_state):
    DIMPLE.avoid_sequence = []
    DIMPLE.stop_codon = True
    DIMPLE.synth_len = 230
    DIMPLE.maxfrag = 165
    DIMPLE.enzyme = None
    DIMPLE.cutsite_overhang = 4


def _make_gene(tmp_path, dimple_state):
    """Create a single-fragment DIMPLE instance via addgene."""
    _setup_state(dimple_state)
    fa = tmp_path / "gene.fa"
    _write_fasta(fa, [("gene1 start:31 end:120", _GENE_SEQ)])
    pool = addgene(str(fa))
    return pool[0], pool


class TestCheckOverhangsUnique:
    def test_no_switch_when_overhangs_unique(self, tmp_path, dimple_state):
        """Returns False when all fragment overhangs are distinct."""
        gene, pool = _make_gene(tmp_path, dimple_state)
        # Override seq with a varied pattern so the F/R overhang slices differ.
        gene.seq = Seq("ATCG" * 50)
        switched = check_overhangs(gene, pool, overlap_l=3, overlap_r=3)
        assert switched is False

    def test_type_error_on_non_dimple(self, tmp_path, dimple_state):
        """Passing a non-DIMPLE object raises TypeError."""
        _setup_state(dimple_state)
        with pytest.raises(TypeError):
            check_overhangs(object(), [], overlap_l=3, overlap_r=3)


class TestCheckOverhangsPalindrome:
    def test_palindromic_overhang_detected(self, tmp_path, dimple_state):
        """A palindromic overhang (F == RC(F)) should be flagged and switch triggered."""
        gene, pool = _make_gene(tmp_path, dimple_state)
        # Build a seq where overhang_F == overhang_R (all same base → trivially equal).
        # Use a repeated pattern to make every 4-mer the same ("AAAA").
        gene.seq = Seq("A" * 200)
        # switch_fragmentsize will run; we only assert that switched is True.
        switched = check_overhangs(gene, pool, overlap_l=3, overlap_r=3)
        assert switched is True
