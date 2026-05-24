"""Unit tests for codon_usage (DIMPLE/utilities.py:75-100)."""

import pytest

from DIMPLE.utilities import codon_usage

# Standard genetic code: codon → amino acid (including stops as '*').
# Used to verify per-amino-acid frequency sums.
_CODON_TO_AA = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def _check_per_aa_sums(table, tol=0.02):
    """Assert that frequencies for each amino acid sum to ~1.0."""
    from collections import defaultdict

    sums = defaultdict(float)
    for codon, freq in table.items():
        sums[_CODON_TO_AA[codon]] += freq
    for aa, total in sums.items():
        assert (
            abs(total - 1.0) <= tol
        ), f"Amino acid {aa!r}: frequencies sum to {total:.4f}, expected ~1.0"


class TestEcoli:
    def test_returns_64_entries(self):
        table = codon_usage("ecoli")
        assert len(table) == 64

    def test_stop_codons_present(self):
        table = codon_usage("ecoli")
        for codon in ("TAA", "TAG", "TGA"):
            assert codon in table

    def test_spot_check_values(self):
        table = codon_usage("ecoli")
        assert table["TTT"] == pytest.approx(0.58)
        assert table["ATG"] == pytest.approx(1.0)

    def test_per_aa_frequency_sums(self):
        _check_per_aa_sums(codon_usage("ecoli"))


class TestHuman:
    def test_returns_64_entries(self):
        table = codon_usage("human")
        assert len(table) == 64

    def test_stop_codons_present(self):
        table = codon_usage("human")
        for codon in ("TAA", "TAG", "TGA"):
            assert codon in table

    def test_spot_check_values(self):
        table = codon_usage("human")
        assert table["TTT"] == pytest.approx(0.45)
        assert table["ATG"] == pytest.approx(1.0)

    def test_per_aa_frequency_sums(self):
        _check_per_aa_sums(codon_usage("human"))


def test_custom_dict_passthrough():
    custom = {"TTT": 1.0}
    assert codon_usage(custom) is custom
