"""Output aggregation helpers for DIMPLE pipeline."""

from __future__ import annotations


def print_all(ols, folder="", config=None):
    """Compatibility wrapper for output aggregation."""
    from DIMPLE.DIMPLE import _legacy_print_all

    return _legacy_print_all(ols, folder=folder, config=config)
