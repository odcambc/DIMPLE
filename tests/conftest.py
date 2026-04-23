"""Shared pytest fixtures for the DIMPLE test suite.

Why this file exists
--------------------
The :class:`DIMPLE.DIMPLE.DIMPLE` class stores most of its configuration as
**class-level** attributes that tests mutate (e.g. ``DIMPLE.usage``,
``DIMPLE.cutsite``, ``DIMPLE.barcodeF``). Without isolation, tests pollute
each other in a single pytest session. The ``dimple_state`` fixture below
snapshots and restores those attributes around every test that requests it.

Caveat: the fixture protects against **rebinding** (``DIMPLE.foo = new``),
not against **in-place mutation** of shared mutable values
(``DIMPLE.foo.append(...)``). Tests that need that should deepcopy the
mutable value themselves before mutating.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict

import pytest

from DIMPLE.DIMPLE import DIMPLE


# Class attributes that integration/regression tests commonly override.
# Snapshot-and-restore these around each test using the `dimple_state` fixture.
# Note: `synth_len` is normally a @property descriptor; tests frequently replace
# it with a bare int. Reading via __dict__ below captures whichever is present.
_MANAGED_ATTRS: tuple[str, ...] = (
    "usage",
    "handle",
    "overlap",
    "synth_len",
    "maxfrag",
    "primerBuffer",
    "barcodeF",
    "barcodeR",
    "cutsite",
    "cutsite_buffer",
    "cutsite_overhang",
    "avoid_sequence",
    "dms",
    "stop_codon",
    "make_double",
    "maximize_nucleotide_change",
    "random_seed",
    "enzyme",
)

# Sentinel for "attribute was not set on the class before the test ran".
_MISSING = object()


def _snapshot_raw(cls: type, name: str) -> Any:
    """Return a deepcopy of the raw class-dict entry for ``name``, or _MISSING.

    Bypasses descriptor protocol so properties are preserved as-is. Falls back
    to a shallow copy or the bare reference if deepcopy fails (property
    objects, Seq objects with internal locks, etc.).
    """
    raw = cls.__dict__.get(name, _MISSING)
    if raw is _MISSING:
        return _MISSING
    try:
        return copy.deepcopy(raw)
    except Exception:
        try:
            return copy.copy(raw)
        except Exception:
            return raw


@pytest.fixture
def dimple_state():
    """Snapshot mutable ``DIMPLE`` class attributes; restore after the test."""
    snapshot: Dict[str, Any] = {
        name: _snapshot_raw(DIMPLE, name) for name in _MANAGED_ATTRS
    }

    yield DIMPLE

    for name, value in snapshot.items():
        if value is _MISSING:
            # Attribute wasn't set before the test; remove anything the test
            # added. Skip silently if the class doesn't allow deletion
            # (e.g. attribute is defined via a non-deletable descriptor).
            try:
                delattr(DIMPLE, name)
            except AttributeError:
                pass
        else:
            try:
                setattr(DIMPLE, name, value)
            except AttributeError:
                # Read-only descriptor; nothing sensible we can do.
                pass


# ---------------------------------------------------------------------------
# Human codon-usage table (moved verbatim from tests/test_dimple.py:19-84).
# Exposed as a fixture so unit tests can opt in without redeclaring the table.
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
def dimple_human_usage(dimple_state):
    """Assign the human codon-usage table to ``DIMPLE.usage`` for the test.

    Depends on ``dimple_state`` so the assignment is rolled back after the
    test runs. Returns the dict for convenience.
    """
    DIMPLE.usage = dict(_HUMAN_USAGE)  # copy to isolate test mutations
    return DIMPLE.usage


@pytest.fixture
def kir_fa() -> Path:
    """Filesystem path to the canonical Kir FASTA test input.

    Prefers the new ``tests/data/Kir.fa`` location (introduced in PR2) and
    falls back to the legacy ``tests/Kir.fa`` so PR1 remains green while the
    move is still pending.
    """
    new_path = Path(__file__).parent / "data" / "Kir.fa"
    old_path = Path(__file__).parent / "Kir.fa"
    return new_path if new_path.exists() else old_path


# ---------------------------------------------------------------------------
# --update-golden plumbing (exercised in PR2's regression tests).
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
