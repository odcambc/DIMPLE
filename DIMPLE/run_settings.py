"""Shared CLI/GUI configuration for DIMPLE runs.

This module centralizes mutation of :class:`DIMPLE` class attributes so the
command-line and GUI entrypoints stay aligned.

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

from DIMPLE.DIMPLE import DIMPLE
from DIMPLE.utilities import codon_usage

# Base primer buffer before overlap extension (matches DIMPLE.primerBuffer default).
PRIMER_BUFFER_BASE: int = 30

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
    logging.basicConfig(filename=log_file, level=level)
    return logging.getLogger(__name__)


def validate_handle(handle: str) -> None:
    """Raise ``ValueError`` if *handle* contains non-nucleotide characters."""
    invalid = [x for x in handle if x not in "ACGTacgt"]
    if invalid:
        raise ValueError("Genetic handle contains non-nucleic acid bases")


def apply_handle(handle: str) -> None:
    """Set ``DIMPLE.handle`` after validation."""
    validate_handle(handle)
    DIMPLE.handle = handle


def resolve_codon_usage(usage_arg: Union[str, dict]) -> dict:
    """Resolve codon usage and assign ``DIMPLE.usage``.

    * ``\"ecoli\"`` / ``\"human\"``: built-in tables via :func:`codon_usage`.
    * Path string (any other str): file contents ``ast.literal_eval`` to a dict.
    * ``dict``: used directly (e.g. GUI custom table).

    Args:
        usage_arg: Named preset, path to a Python-literal dict file, or dict.

    Returns:
        The resolved codon usage table.
    """
    if isinstance(usage_arg, dict):
        DIMPLE.usage = usage_arg
        return DIMPLE.usage
    if usage_arg in ("ecoli", "human"):
        DIMPLE.usage = codon_usage(usage_arg)
        return DIMPLE.usage
    with open(usage_arg, encoding="utf-8") as f:
        text = f.read()
    DIMPLE.usage = ast.literal_eval(text.strip())
    return DIMPLE.usage


def apply_restriction_settings(restriction_sequence: str) -> None:
    """Parse restriction enzyme / site string and set cutsite-related class attrs."""
    if re.match(r"[ACGT]+\([ACGT]\)\d+/\d+", restriction_sequence):
        tmp_cutsite = restriction_sequence.split("(")
        DIMPLE.cutsite = Seq(tmp_cutsite[0])
        DIMPLE.cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
        tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
        DIMPLE.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])
        if (
            DIMPLE.cutsite == Seq("GGTCTC")
            and DIMPLE.cutsite_buffer == Seq("G")
            and DIMPLE.cutsite_overhang == 4
        ):
            DIMPLE.enzyme = "BsaI"
        elif (
            DIMPLE.cutsite == Seq("CGTCTC")
            and DIMPLE.cutsite_buffer == Seq("G")
            and DIMPLE.cutsite_overhang == 4
        ):
            DIMPLE.enzyme = "BsmBI"
        else:
            DIMPLE.enzyme = None
    elif restriction_sequence.upper() in ("BSAI", "BSMBI"):
        if restriction_sequence.upper() == "BSAI":
            DIMPLE.cutsite = Seq("GGTCTC")
            DIMPLE.cutsite_buffer = Seq("G")
            DIMPLE.cutsite_overhang = 4
            DIMPLE.enzyme = "BsaI"
        else:
            DIMPLE.cutsite = Seq("CGTCTC")
            DIMPLE.cutsite_buffer = Seq("G")
            DIMPLE.cutsite_overhang = 4
            DIMPLE.enzyme = "BsmBI"
    else:
        raise ValueError(
            f"Restriction sequence {restriction_sequence!r} not recognized. "
            "Please check input."
        )


def compute_overlaps_and_maxfrag(
    oligo_len: int,
    fragment_len: int,
    overlap: int,
    deletions: Union[bool, List[int]],
    logger: Optional[logging.Logger] = None,
) -> tuple[int, int]:
    """Set ``DIMPLE.synth_len``, ``DIMPLE.maxfrag``, and primer buffer from overlap.

    Resets ``DIMPLE.primerBuffer`` to ``PRIMER_BUFFER_BASE + overlapL`` so repeated
    configuration in the same process does not accumulate.

    Returns:
        ``(overlapL, overlapR)`` after deletion adjustment.
    """
    overlap_l = int(overlap)
    overlap_r = int(overlap)
    if deletions:
        overlap_r = max(int(x) for x in deletions) + overlap_r - 3

    DIMPLE.synth_len = oligo_len
    if fragment_len != 0:
        DIMPLE.maxfrag = fragment_len
        if logger:
            logger.info("Maximum fragment length: %s based on input", DIMPLE.maxfrag)
    else:
        DIMPLE.maxfrag = oligo_len - 64 - overlap_l - overlap_r
        if logger:
            logger.info(
                "Maximum fragment length: %s based on oligo length and overlap: 2 * %s",
                DIMPLE.maxfrag,
                overlap,
            )

    DIMPLE.primerBuffer = PRIMER_BUFFER_BASE + overlap_l
    return overlap_l, overlap_r


def apply_barcode_start(start: int) -> None:
    """Slice barcode primer lists from *start*."""
    idx = int(start)
    DIMPLE.barcodeF = DIMPLE.barcodeF[idx:]
    DIMPLE.barcodeR = DIMPLE.barcodeR[idx:]


def normalize_avoid_list(
    sequences: Union[List[str], str],
    logger: Optional[logging.Logger] = None,
) -> List[Seq]:
    """Build ``DIMPLE.avoid_sequence`` as ``Seq`` objects; strip whitespace (GUI comma lists)."""
    if isinstance(sequences, str):
        if sequences.strip() == "":
            raw: List[str] = []
        else:
            raw = [s.strip() for s in sequences.split(",")]
    else:
        raw = [str(s).strip() for s in sequences]

    DIMPLE.avoid_sequence = [Seq(x) for x in raw]

    if DIMPLE.cutsite not in DIMPLE.avoid_sequence:
        DIMPLE.avoid_sequence.append(DIMPLE.cutsite)
        msg = (
            f"Restriction sequence {DIMPLE.cutsite} was not included in the avoid list. "
            "Adding before continuing."
        )
        warnings.warn(msg)
        if logger:
            logger.warning(msg)

    return DIMPLE.avoid_sequence


def apply_random_seed(seed: Optional[int]) -> None:
    """Set ``DIMPLE.random_seed`` (``None`` means nondeterministic / NumPy default)."""
    DIMPLE.random_seed = seed


def validate_insertions(insertions: List[str]) -> None:
    """Validate insertion strings: ACGT only and must not contain the cutsite."""
    cutsite_str = str(DIMPLE.cutsite)
    for insertion in insertions:
        if any(b not in "ACGTacgt" for b in insertion):
            raise ValueError(
                f"Insertions contain non-nucleic acid bases: {insertions!r}"
            )
        if cutsite_str in insertion:
            raise ValueError(
                f"Insertions contain restriction sites ({cutsite_str}): {insertions!r}"
            )


def apply_instance_settings(
    instances: list,
    aminoacids: Optional[List[str]] = None,
    doublefrag: Optional[int] = None,
    gene_primer_tm: Optional[tuple[int, int]] = None,
) -> None:
    """Apply per-gene and class-level settings after :func:`addgene`.

    ``gene_primer_tm`` updates the **class** attribute ``DIMPLE.gene_primerTm`` because
    primer logic reads ``DIMPLE.gene_primerTm``, not per-instance values.

    Args:
        instances: List of :class:`DIMPLE` instances from ``addgene``.
        aminoacids: Three-letter amino acid codes to scan (replaces ``gene.aminoacids``).
        doublefrag: ``0`` or ``1`` fragment-per-oligo layout.
        gene_primer_tm: Melting temperature bounds for gene primers.
    """
    if gene_primer_tm is not None:
        DIMPLE.gene_primerTm = gene_primer_tm
    for gene in instances:
        if aminoacids is not None:
            gene.aminoacids = [a.strip() for a in aminoacids]
        if doublefrag is not None:
            gene.doublefrag = int(doublefrag)
