"""Primer design helpers for DIMPLE pipeline."""

from __future__ import annotations


def find_geneprimer(genefrag, start, end):
    """Compatibility wrapper for gene primer design."""
    from DIMPLE.DIMPLE import _legacy_find_geneprimer

    return _legacy_find_geneprimer(genefrag, start, end)


def find_fragment_primer(fragment, stop):
    """Compatibility wrapper for fragment primer design."""
    from DIMPLE.DIMPLE import _legacy_find_fragment_primer

    return _legacy_find_fragment_primer(fragment, stop)


def check_nonspecific(primer, fragment, point):
    """Compatibility wrapper for nonspecific-primer checks."""
    from DIMPLE.DIMPLE import _legacy_check_nonspecific

    return _legacy_check_nonspecific(primer, fragment, point)
