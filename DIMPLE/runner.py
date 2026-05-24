"""Shared pipeline driver for CLI and notebook entrypoints.

The :func:`build_runtime_config` helper takes the user-facing knobs that both
``run_dimple.py`` and ``DIMPLE.ipynb`` expose, threads them through the existing
``run_settings`` mutators, and returns a populated
:class:`~DIMPLE.pool.DimpleRuntimeConfig`. :func:`run_pipeline` then drives the
pipeline -- :func:`addgene` to build the :class:`~DIMPLE.pool.Pool`,
:func:`apply_instance_settings`, optional :func:`validate_insertions` /
:func:`align_genevariation`, :func:`generate_DMS_fragments`, :func:`post_qc`,
and :func:`print_all`.

Both entrypoints (CLI / notebook) collapse to one config call plus one
pipeline call, removing the ~70-line block that previously drifted between
them.
"""

from __future__ import annotations

import logging
import os
from typing import List, Optional, Sequence, Union

from DIMPLE.pool import DimpleRuntimeConfig, Pool, addgene
from DIMPLE.run_settings import (
    apply_barcode_start,
    apply_handle,
    apply_instance_settings,
    apply_random_seed,
    apply_restriction_settings,
    apply_runtime_policies,
    compute_overlaps_and_maxfrag,
    normalize_avoid_list,
    resolve_codon_usage,
    validate_insertions,
)

__all__ = ("build_runtime_config", "run_pipeline")


_logger = logging.getLogger(__name__)


def build_runtime_config(
    *,
    oligo_len: int,
    fragment_len: int,
    overlap: int,
    handle: str,
    restriction_sequence: str,
    avoid_sequence: Union[List[str], str],
    codon_usage: Union[str, dict],
    barcode_start: int,
    deletions: Union[bool, Sequence[int]],
    dms: bool,
    dis: bool,
    make_double: bool,
    stop_codon: bool,
    maximize_nucleotide_change: bool,
    non_interactive: bool,
    preferred_orf_index: Optional[int],
    link_policy: str,
    breaksite_change_policy: str,
    random_seed: Optional[int],
    logger: Optional[logging.Logger] = None,
) -> tuple[DimpleRuntimeConfig, int, int]:
    """Build a populated :class:`DimpleRuntimeConfig` from user-facing knobs.

    Mirrors the argparse fields of ``run_dimple.py`` and the ``#@param`` form
    fields of ``DIMPLE.ipynb``. Internally drives the existing ``run_settings``
    helpers, so behavior matches both entrypoints exactly.

    Args:
        oligo_len: Synthesized oligo length.
        fragment_len: Maximum fragment length; ``0`` means derive from oligo.
        overlap: Base pairs of overlap between adjacent fragments.
        handle: Domain-insertion genetic handle (nucleotides).
        restriction_sequence: Enzyme recognition string (e.g. ``CGTCTC(G)1/5``).
        avoid_sequence: Sequences to avoid (list or comma-separated string).
        codon_usage: Codon usage preset name, dict, or path.
        barcode_start: Skip the first N barcode primer pairs.
        deletions: Deletion lengths (list[int]) or ``False``; used here only to
            size the right-hand overlap so fragments survive symmetric
            deletions. The actual deletion list is passed to
            :func:`run_pipeline`.
        dms: Run a deep mutational scan.
        dis: Domain-insertion scan (stored on config for symmetry; consumed by
            :func:`generate_DMS_fragments` via :func:`run_pipeline`).
        make_double: Generate all double-mutant combinations within a fragment.
        stop_codon: Add stop codons to the scanning alphabet.
        maximize_nucleotide_change: Pick codons that maximize nt edits.
        non_interactive: Refuse interactive prompts (fail fast).
        preferred_orf_index: ORF to pick when several are found (1-based).
        link_policy: Policy for ``align_genevariation`` linking.
        breaksite_change_policy: Policy when breaksite endpoints change.
        random_seed: Seed for NumPy RNG; ``None`` is nondeterministic.
        logger: Optional logger forwarded to overlap / avoid helpers.

    Returns:
        ``(config, overlap_l, overlap_r)`` -- the populated config plus the
        left/right overlaps needed by :func:`generate_DMS_fragments`.
    """
    # `dis` is currently not stored on the config; it's plumbed through to
    # generate_DMS_fragments at pipeline time. Accept it here so callers have
    # a single place to declare mutation intent.
    del dis

    config = DimpleRuntimeConfig()
    apply_handle(handle, config=config)

    overlap_l, overlap_r = compute_overlaps_and_maxfrag(
        oligo_len,
        fragment_len,
        overlap,
        deletions if deletions else False,
        logger=logger,
        config=config,
    )

    apply_barcode_start(int(barcode_start), config=config)
    apply_restriction_settings(restriction_sequence, config=config)
    normalize_avoid_list(avoid_sequence, logger=logger, config=config)

    apply_runtime_policies(
        dms=dms,
        stop_codon=stop_codon,
        make_double=make_double,
        maximize_nucleotide_change=maximize_nucleotide_change,
        non_interactive=non_interactive,
        preferred_orf_index=preferred_orf_index,
        link_policy=link_policy,
        breaksite_change_policy=breaksite_change_policy,
        config=config,
    )

    apply_random_seed(random_seed, config=config)
    resolve_codon_usage(codon_usage, config=config)

    return config, overlap_l, overlap_r


def run_pipeline(
    target_file: str,
    work_dir: str,
    config: DimpleRuntimeConfig,
    overlap_l: int,
    overlap_r: int,
    *,
    include_synonymous: bool = False,
    custom_mutations=None,
    insertions: Union[bool, Sequence[str]] = False,
    deletions: Union[bool, Sequence[int]] = False,
    dis: bool = False,
    match_sequences: bool = False,
    aminoacids: Optional[Sequence[str]] = None,
) -> Pool:
    """Drive the DIMPLE pipeline end-to-end against *target_file*.

    Builds the :class:`Pool` from *target_file*, applies per-gene settings,
    optionally runs :func:`align_genevariation`, generates fragments, runs
    post-QC, and writes outputs to *work_dir*. Returns the populated pool.

    Raises:
        ValueError: if no mutation type is enabled (none of ``config.dms``,
            ``insertions``, ``deletions``, or ``dis``).
    """
    # Imported here to dodge a cycle: DIMPLE.DIMPLE re-exports addgene/Pool.
    from DIMPLE.DIMPLE import (
        align_genevariation,
        generate_DMS_fragments,
        post_qc,
        print_all,
    )

    if not any([config.dms, insertions, deletions, dis]):
        raise ValueError("Didn't select any mutations to generate")

    pool = addgene(os.path.join(work_dir, target_file).strip(), config)

    if aminoacids is not None:
        apply_instance_settings(pool, config=config, aminoacids=list(aminoacids))
    else:
        apply_instance_settings(pool, config=config)

    if insertions:
        validate_insertions(list(insertions), config)

    if match_sequences:
        align_genevariation(pool)

    _logger.info("Generating DMS fragments")
    generate_DMS_fragments(
        pool,
        overlap_l,
        overlap_r,
        include_synonymous,
        custom_mutations,
        config.dms,
        insertions,
        deletions,
        dis,
        work_dir,
    )

    post_qc(pool)
    print_all(pool, work_dir)

    return pool
