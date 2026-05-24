"""Shared pytest fixtures for the DIMPLE test suite.

Run configuration now lives on a :class:`~DIMPLE.pool.DimpleRuntimeConfig`
that each test builds explicitly, so the old ``dimple_state`` class-attribute
snapshot/restore fixture is no longer needed -- tests are isolated by
construction.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import pytest


# ---------------------------------------------------------------------------
# Human codon-usage table, exposed as a fixture so tests can opt in without
# redeclaring the table.
# ---------------------------------------------------------------------------
_HUMAN_USAGE: Dict[str, float] = {
    "TTT": 0.45, "TTC": 0.55, "TTA": 0.07, "TTG": 0.13,
    "TAT": 0.43, "TAC": 0.57, "TAA": 0.28, "TAG": 0.20,
    "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.41,
    "CAT": 0.41, "CAC": 0.59, "CAA": 0.25, "CAG": 0.75,
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16, "ATG": 1.00,
    "AAT": 0.46, "AAC": 0.54, "AAA": 0.42, "AAG": 0.58,
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    "GAT": 0.46, "GAC": 0.54, "GAA": 0.42, "GAG": 0.58,
    "TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06,
    "TGT": 0.45, "TGC": 0.55, "TGA": 0.52, "TGG": 1.00,
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21,
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    "AGT": 0.15, "AGC": 0.24, "AGA": 0.20, "AGG": 0.20,
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
}


@pytest.fixture
def dimple_human_usage() -> Dict[str, float]:
    """Return the human codon-usage table (a fresh copy per test).

    Tests pass this as ``DimpleRuntimeConfig(usage=...)``.
    """
    return dict(_HUMAN_USAGE)


@pytest.fixture
def kir_fa() -> Path:
    """Filesystem path to the canonical Kir FASTA test input."""
    new_path = Path(__file__).parent / "data" / "Kir.fa"
    old_path = Path(__file__).parent / "Kir.fa"
    return new_path if new_path.exists() else old_path


# ---------------------------------------------------------------------------
# --update-golden plumbing (exercised by the regression tests).
# ---------------------------------------------------------------------------
def pytest_addoption(parser: pytest.Parser) -> None:
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Overwrite golden files in tests/expected/ instead of asserting.",
    )


@pytest.fixture
def update_golden(request: pytest.FixtureRequest) -> bool:
    """True if pytest was invoked with --update-golden."""
    return bool(request.config.getoption("--update-golden"))
