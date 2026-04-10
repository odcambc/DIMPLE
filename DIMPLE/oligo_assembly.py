"""Oligo assembly helpers for DIMPLE pipeline."""

from __future__ import annotations


def combine_fragments(tandem_oligos):
    """Compatibility wrapper for tandem oligo assembly."""
    from DIMPLE.DIMPLE import _legacy_combine_fragments

    return _legacy_combine_fragments(tandem_oligos)
