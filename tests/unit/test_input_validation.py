"""Input-validation tests for ``addgene`` and ``DIMPLE.__init__``.

Most of these are not about interactivity -- they're about what happens when a
caller hands DIMPLE bad or ambiguous input. The one exception, in
``TestMissingORFCoords``, also verifies the ``non_interactive=True`` contract:
no ``input()`` prompt when the ORF cannot be uniquely identified.
"""

from __future__ import annotations

import pytest
from Bio.Seq import Seq

from DIMPLE.DIMPLE import addgene
from DIMPLE.pool import DimpleRuntimeConfig

_CLEAN_BODY = "ATCG" * 50

_BSMBI = Seq("CGTCTC")


def _write_fasta(path, records):
    path.write_text("".join(f">{header}\n{seq}\n" for header, seq in records))


def _minimal_config(**overrides) -> DimpleRuntimeConfig:
    defaults = dict(
        avoid_sequence=[],
        stop_codon=True,
        synth_len=230,
        maxfrag=165,
        enzyme=None,
        non_interactive=True,
    )
    defaults.update(overrides)
    return DimpleRuntimeConfig(**defaults)


class TestMissingORFCoords:
    """No ORF in FASTA header + ``non_interactive=True`` -> raise, don't prompt."""

    def test_no_orf_candidates_raises(self, tmp_path):
        fa = tmp_path / "no_orf.fa"
        _write_fasta(fa, [("no_orf", _CLEAN_BODY)])

        with pytest.raises(ValueError, match=r"(?i)no valid orf"):
            addgene(str(fa), _minimal_config())

    def test_multiple_orf_candidates_raises(self, tmp_path):
        ambiguous = "ATG" + "GCC" * 101 + "TAA" + "ATG" + "GCC" * 101 + "TAA"
        fa = tmp_path / "multi_orf.fa"
        _write_fasta(fa, [("multi_orf", ambiguous)])

        with pytest.raises(ValueError, match=r"(?i)multiple orf"):
            addgene(str(fa), _minimal_config())


class TestRestrictionSiteInORF:
    def test_bsmbi_site_in_gene_raises(self, tmp_path):
        with_site = "ATCG" * 20 + str(_BSMBI) + "ATCG" * 20
        fa = tmp_path / "with_cutsite.fa"
        _write_fasta(fa, [("with_cutsite start:1 end:" + str(len(with_site)), with_site)])

        with pytest.raises(ValueError, match=r"(?i)unwanted restriction"):
            addgene(str(fa), _minimal_config(avoid_sequence=[_BSMBI]))


class TestExplicitNonMultipleOfThreeORF:
    """Bad explicit coords trigger fallback to ``findORF``, not silent acceptance."""

    def test_bad_coords_not_silently_accepted(self, tmp_path):
        # Length 7, start:1 end:7 -> not mod 3. No clean ORF, so the fallback raises.
        body = "ATCGATC"
        fa = tmp_path / "bad_len.fa"
        _write_fasta(fa, [("bad_len start:1 end:7", body)])

        with pytest.raises(ValueError):
            addgene(str(fa), _minimal_config())

    def test_bad_coords_recovered_when_recoverable(self, tmp_path):
        # Header claims length 305 (not mod 3). The sequence has one clean ORF
        # of length 306, so fallback resolves to mod-3 coords.
        body = "ATG" + "GCC" * 100 + "TAA"
        fa = tmp_path / "recoverable.fa"
        _write_fasta(fa, [("recoverable start:1 end:305", body)])

        pool = addgene(str(fa), _minimal_config())
        gene = pool[0]
        assert (gene.end - gene.start) % 3 == 0


class TestInternalStopCodonInORF:
    """In-frame stop inside an explicit ORF: currently accepted, no validation.

    Pinned so any future tightening shows up as a deliberate change.
    """

    def test_internal_stop_currently_accepted(self, tmp_path):
        body = "ATG" + "GCC" * 20 + "TAA" + "GCC" * 20 + "TAA"
        fa = tmp_path / "internal_stop.fa"
        _write_fasta(fa, [("internal_stop start:1 end:" + str(len(body)), body)])

        pool = addgene(str(fa), _minimal_config())
        gene = pool[0]
        translated = str(Seq(body[gene.start : gene.end]).translate())
        assert "*" in translated[:-1], "fixture must contain an internal stop to be meaningful"
