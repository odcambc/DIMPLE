"""Fragment layout helpers for DIMPLE pipeline."""

from __future__ import annotations


def recalculate_num_fragments(gene):
    """Compatibility wrapper for fragment-count recalculation."""
    from DIMPLE.DIMPLE import _legacy_recalculate_num_fragments

    return _legacy_recalculate_num_fragments(gene)


def switch_fragmentsize(gene, detectedsite, ols):
    """Compatibility wrapper for fragment-size switching logic."""
    from DIMPLE.DIMPLE import _legacy_switch_fragmentsize

    return _legacy_switch_fragmentsize(gene, detectedsite, ols)


def check_overhangs(gene, ols, overlap_l, overlap_r):
    """Compatibility wrapper for overhang quality checks."""
    from DIMPLE.DIMPLE import _legacy_check_overhangs

    return _legacy_check_overhangs(gene, ols, overlap_l, overlap_r)
