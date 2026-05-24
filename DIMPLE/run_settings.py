"""Shared CLI/GUI configuration for DIMPLE runs.

This module centralizes population of a :class:`DimpleRuntimeConfig` so the
command-line and GUI entrypoints stay aligned. Each helper takes the run's
config object and mutates it in place; the populated config is then handed to
:func:`addgene` to build the oligo :class:`~DIMPLE.pool.Pool`.

Random seed policy:
    * CLI: ``None`` unless ``-seed`` is passed (nondeterministic RNG per gene
      when ``None``).
    * GUI: defaults to ``1848`` for reproducibility unless overridden.
"""

from __future__ import annotations

import ast
import logging
import os
import re
import warnings
from typing import List, Optional, Union

from Bio.Seq import Seq

from DIMPLE.pool import PRIMER_BUFFER_BASE, DimpleRuntimeConfig
from DIMPLE.utilities import codon_usage

# Default GUI random seed (matches historical GUI behavior and tests).
DEFAULT_GUI_RANDOM_SEED: int = 1848


def configure_dimple_logging(log_file: str, level: int = logging.INFO) -> logging.Logger:
    """Configure the root logger for file output and return the ``run_settings`` logger.

    Args:
        log_file: Path to the log file (parent directory is created if missing).
        level: Logging level (default INFO).

    Returns:
        Logger for this module (callers may also use ``logging.getLogger(__name__)``).
    """
    log_dir = os.path.dirname(os.path.abspath(log_file))
    if log_dir and not os.path.exists(log_dir):
        os.makedirs(log_dir)
    logging.basicConfig(filename=log_file, level=level, force=True)
    return logging.getLogger(__name__)


def validate_handle(handle: str) -> None:
    """Raise ``ValueError`` if *handle* contains non-nucleotide characters."""
    invalid = [x for x in handle if x not in "ACGTacgt"]
    if invalid:
        raise ValueError("Genetic handle contains non-nucleic acid bases")


def apply_handle(handle: str, config: DimpleRuntimeConfig) -> None:
    """Validate *handle* and store it on *config*."""
    validate_handle(handle)
    config.handle = handle


_VALID_CODONS = frozenset(a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT")


def _validate_usage_table(table, source: str) -> None:
    """Raise ``ValueError`` if *table* is not a dict of all 64 codons in [0, 1]."""
    if not isinstance(table, dict):
        raise ValueError(
            f"Codon usage from {source} must be a dict of 64 codons; "
            f"got {type(table).__name__}."
        )
    missing = _VALID_CODONS - table.keys()
    extra = table.keys() - _VALID_CODONS
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"missing codons: {sorted(missing)}")
        if extra:
            parts.append(f"unrecognized keys: {sorted(extra)}")
        raise ValueError(f"Codon usage from {source} is malformed: {'; '.join(parts)}.")
    for codon, freq in table.items():
        if not isinstance(freq, (int, float)) or isinstance(freq, bool):
            raise ValueError(
                f"Codon usage from {source}: value for {codon!r} is not a number " f"({freq!r})."
            )
        if not 0 <= freq <= 1:
            raise ValueError(
                f"Codon usage from {source}: value for {codon!r} ({freq}) is " "outside [0, 1]."
            )


def resolve_codon_usage(usage_arg: Union[str, dict], config: DimpleRuntimeConfig) -> dict:
    """Resolve codon usage and assign it to ``config.usage``.

    * ``"ecoli"`` / ``"human"``: built-in tables via :func:`codon_usage`.
    * ``dict``: used directly (e.g. GUI custom table).
    * Path string (any other str): the file must exist; its contents are
      ``ast.literal_eval``'d to a dict.

    The resolved table is validated (64 codons present, all values numeric in
    ``[0, 1]``) before being assigned to ``config.usage``.

    Returns:
        The resolved codon usage table.
    """
    if isinstance(usage_arg, dict):
        config.usage = usage_arg
        source = "<custom dict>"
    elif usage_arg in ("ecoli", "human"):
        config.usage = codon_usage(usage_arg)
        source = f"built-in {usage_arg!r}"
    else:
        if not os.path.exists(usage_arg):
            raise ValueError(
                f"Codon usage {usage_arg!r}: not a built-in preset and no file "
                "exists at that path. Did you mean 'ecoli' or 'human'?"
            )
        with open(usage_arg, encoding="utf-8") as f:
            text = f.read()
        config.usage = ast.literal_eval(text.strip())
        source = f"file {usage_arg!r}"
    _validate_usage_table(config.usage, source)
    return config.usage


def apply_restriction_settings(restriction_sequence: str, config: DimpleRuntimeConfig) -> None:
    """Parse a restriction enzyme / site string into ``config``'s cutsite fields."""
    if re.match(r"[ACGT]+\([ACGT]\)\d+/\d+", restriction_sequence):
        tmp_cutsite = restriction_sequence.split("(")
        config.cutsite = Seq(tmp_cutsite[0])
        config.cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
        tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
        config.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])
        if (
            config.cutsite == Seq("GGTCTC")
            and config.cutsite_buffer == Seq("G")
            and config.cutsite_overhang == 4
        ):
            config.enzyme = "BsaI"
        elif (
            config.cutsite == Seq("CGTCTC")
            and config.cutsite_buffer == Seq("G")
            and config.cutsite_overhang == 4
        ):
            config.enzyme = "BsmBI"
        else:
            config.enzyme = None
    elif restriction_sequence.upper() in ("BSAI", "BSMBI"):
        if restriction_sequence.upper() == "BSAI":
            config.cutsite = Seq("GGTCTC")
            config.cutsite_buffer = Seq("G")
            config.cutsite_overhang = 4
            config.enzyme = "BsaI"
        else:
            config.cutsite = Seq("CGTCTC")
            config.cutsite_buffer = Seq("G")
            config.cutsite_overhang = 4
            config.enzyme = "BsmBI"
    else:
        raise ValueError(
            f"Restriction sequence {restriction_sequence!r} not recognized. " "Please check input."
        )


def compute_overlaps_and_maxfrag(
    oligo_len: int,
    fragment_len: int,
    overlap: int,
    deletions: Union[bool, List[int]],
    config: DimpleRuntimeConfig,
    logger: Optional[logging.Logger] = None,
) -> tuple[int, int]:
    """Set ``config``'s ``synth_len``, ``maxfrag``, and ``primer_buffer`` from overlap.

    Returns:
        ``(overlapL, overlapR)`` after deletion adjustment.
    """
    overlap_l = int(overlap)
    overlap_r = int(overlap)
    if deletions:
        overlap_r = max(int(x) for x in deletions) + overlap_r - 3

    config.synth_len = oligo_len
    if fragment_len != 0:
        config.maxfrag = fragment_len
        if logger:
            logger.info("Maximum fragment length: %s based on input", config.maxfrag)
    else:
        config.maxfrag = oligo_len - 64 - overlap_l - overlap_r
        if logger:
            logger.info(
                "Maximum fragment length: %s based on oligo length and overlap: 2 * %s",
                config.maxfrag,
                overlap,
            )

    config.primer_buffer = PRIMER_BUFFER_BASE + overlap_l
    return overlap_l, overlap_r


def apply_barcode_start(start: int, config: DimpleRuntimeConfig) -> None:
    """Drop the first *start* barcode primer pairs from ``config``."""
    idx = int(start)
    config.barcode_f = config.barcode_f[idx:]
    config.barcode_r = config.barcode_r[idx:]


def normalize_avoid_list(
    sequences: Union[List[str], str],
    config: DimpleRuntimeConfig,
    logger: Optional[logging.Logger] = None,
) -> List[Seq]:
    """Build ``config.avoid_sequence`` as ``Seq`` objects; strip whitespace (GUI comma lists)."""
    if isinstance(sequences, str):
        if sequences.strip() == "":
            raw: List[str] = []
        else:
            raw = [s.strip() for s in sequences.split(",")]
    else:
        raw = [str(s).strip() for s in sequences]

    config.avoid_sequence = [Seq(x) for x in raw]

    if config.cutsite not in config.avoid_sequence:
        config.avoid_sequence.append(config.cutsite)
        msg = (
            f"Restriction sequence {config.cutsite} was not included in the avoid list. "
            "Adding before continuing."
        )
        warnings.warn(msg)
        if logger:
            logger.warning(msg)

    return config.avoid_sequence


def apply_random_seed(seed: Optional[int], config: DimpleRuntimeConfig) -> None:
    """Set ``config.random_seed`` (``None`` means nondeterministic / NumPy default)."""
    config.random_seed = seed


def apply_runtime_policies(
    *,
    dms: bool,
    stop_codon: bool,
    make_double: bool,
    maximize_nucleotide_change: bool,
    non_interactive: bool,
    preferred_orf_index: Optional[int],
    link_policy: str,
    breaksite_change_policy: str,
    config: DimpleRuntimeConfig,
) -> None:
    """Apply runtime behavior policy flags to *config*."""
    config.dms = dms
    config.stop_codon = stop_codon
    config.make_double = make_double
    config.maximize_nucleotide_change = maximize_nucleotide_change
    config.non_interactive = non_interactive
    config.preferred_orf_index = preferred_orf_index
    config.link_policy = link_policy
    config.breaksite_change_policy = breaksite_change_policy


def validate_insertions(insertions: List[str], config: DimpleRuntimeConfig) -> None:
    """Validate insertion strings: ACGT only and must not contain the cutsite."""
    cutsite_str = str(config.cutsite)
    for insertion in insertions:
        if any(b not in "ACGTacgt" for b in insertion):
            raise ValueError(f"Insertions contain non-nucleic acid bases: {insertions!r}")
        if cutsite_str in insertion:
            raise ValueError(
                f"Insertions contain restriction sites ({cutsite_str}): {insertions!r}"
            )


def apply_instance_settings(
    instances: list,
    config: DimpleRuntimeConfig,
    aminoacids: Optional[List[str]] = None,
    doublefrag: Optional[int] = None,
    gene_primer_tm: Optional[tuple[int, int]] = None,
) -> None:
    """Apply per-gene and run-wide settings after :func:`addgene`.

    ``gene_primer_tm`` updates ``config.gene_primer_tm`` because primer logic
    reads run config, not per-instance values.

    Args:
        instances: The :class:`DIMPLE` gene instances (e.g. ``pool.genes``).
        config: The run configuration.
        aminoacids: Three-letter amino acid codes to scan (replaces ``gene.aminoacids``).
        doublefrag: ``0`` or ``1`` fragment-per-oligo layout.
        gene_primer_tm: Melting temperature bounds for gene primers.
    """
    if gene_primer_tm is not None:
        config.gene_primer_tm = gene_primer_tm
    for gene in instances:
        if aminoacids is not None:
            gene.aminoacids = [a.strip() for a in aminoacids]
        if doublefrag is not None:
            gene.doublefrag = int(doublefrag)
