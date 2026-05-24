"""Unit tests for addgene (the Pool factory)."""

from DIMPLE.DIMPLE import addgene
from DIMPLE.pool import DimpleRuntimeConfig

# Minimal gene FASTA body: 200 bases, no restriction-site sequences.
_GENE_SEQ = "ATCG" * 50  # 200 bp, 50% GC, no BsmBI/BsaI sites


def _write_fasta(path, records):
    """Write a simple FASTA file; records is a list of (header, seq) tuples."""
    text = "".join(f">{header}\n{seq}\n" for header, seq in records)
    path.write_text(text)


def _minimal_config() -> DimpleRuntimeConfig:
    """A runtime config with just enough set to create DIMPLE instances."""
    return DimpleRuntimeConfig(
        avoid_sequence=[],
        stop_codon=True,
        synth_len=230,
        maxfrag=165,
        enzyme=None,
    )


class TestAddgeneStartEnd:
    def test_start_is_zero_indexed(self, tmp_path):
        """1-indexed start: in FASTA header is converted to 0-indexed."""
        fa = tmp_path / "gene.fa"
        _write_fasta(fa, [("gene1 start:31 end:120", _GENE_SEQ)])
        pool = addgene(str(fa), _minimal_config())
        assert pool[0].start == 30  # 31 - 1

    def test_end_is_unchanged(self, tmp_path):
        """end: in FASTA header is stored as-is (no offset)."""
        fa = tmp_path / "gene.fa"
        _write_fasta(fa, [("gene1 start:31 end:120", _GENE_SEQ)])
        pool = addgene(str(fa), _minimal_config())
        assert pool[0].end == 120

    def test_multi_record_returns_correct_count(self, tmp_path):
        """Multi-record FASTA produces one DIMPLE object per record."""
        fa = tmp_path / "genes.fa"
        _write_fasta(
            fa,
            [
                ("geneA start:31 end:120", _GENE_SEQ),
                ("geneB start:31 end:120", _GENE_SEQ),
            ],
        )
        pool = addgene(str(fa), _minimal_config())
        assert len(pool) == 2

    def test_geneid_matches_fasta_name(self, tmp_path):
        """DIMPLE.geneid is set from the FASTA record name field."""
        fa = tmp_path / "gene.fa"
        _write_fasta(fa, [("Kir start:31 end:120", _GENE_SEQ)])
        pool = addgene(str(fa), _minimal_config())
        assert pool[0].geneid == "Kir"
