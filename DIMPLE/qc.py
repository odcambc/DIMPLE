"""QC helpers for DIMPLE pipeline."""

from __future__ import annotations


def post_qc(ols, config=None):
    """Compatibility wrapper for post-run QC checks."""
    from DIMPLE.DIMPLE import _legacy_post_qc

    return _legacy_post_qc(ols, config=config)


def test_final_assembly(gene, config=None):
    """Compatibility wrapper for assembly simulation QC."""
    from DIMPLE.DIMPLE import _legacy_test_final_assembly

    return _legacy_test_final_assembly(gene, config=config)
