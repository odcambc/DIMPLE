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
from dataclasses import dataclass
from typing import List, Optional, Union

from Bio.Seq import Seq

from DIMPLE.DIMPLE import DIMPLE
from DIMPLE.utilities import codon_usage

# Base primer buffer before overlap extension (matches DIMPLE.primerBuffer default).
PRIMER_BUFFER_BASE: int = 30

# Default GUI random seed (matches historical GUI behavior and tests).
DEFAULT_GUI_RANDOM_SEED: int = 1848


@dataclass
class DimpleRuntimeConfig:
    """Explicit runtime settings used by CLI/GUI and core pipeline."""

    handle: str = ""
    usage: Optional[dict] = None
    cutsite: Optional[Seq] = None
    cutsite_buffer: Optional[Seq] = None
    cutsite_overhang: int = 0
    enzyme: Optional[str] = None
    synth_len: Optional[int] = None
    maxfrag: Optional[int] = None
    primer_buffer: int = PRIMER_BUFFER_BASE
    random_seed: Optional[int] = 0
    avoid_sequence: Optional[List[Seq]] = None
    barcode_f: Optional[list] = None
    barcode_r: Optional[list] = None
    dms: bool = False
    stop_codon: bool = False
    make_double: bool = False
    maximize_nucleotide_change: bool = False
    gene_primer_tm: tuple[int, int] = (58, 62)
    non_interactive: bool = False
    preferred_orf_index: Optional[int] = None
    link_policy: str = "prompt"
    breaksite_change_policy: str = "prompt"


_RUNTIME_CONFIG: Optional[DimpleRuntimeConfig] = None


def get_runtime_config() -> DimpleRuntimeConfig:
    """Return the active runtime config, creating one from current class state."""
    global _RUNTIME_CONFIG
    if _RUNTIME_CONFIG is None:
        _RUNTIME_CONFIG = DimpleRuntimeConfig(
            handle=getattr(DIMPLE, "handle", ""),
            usage=getattr(DIMPLE, "usage", None),
            cutsite=getattr(DIMPLE, "cutsite", None),
            cutsite_buffer=getattr(DIMPLE, "cutsite_buffer", None),
            cutsite_overhang=getattr(DIMPLE, "cutsite_overhang", 0),
            enzyme=getattr(DIMPLE, "enzyme", None),
            synth_len=getattr(DIMPLE, "synth_len", None),
            maxfrag=getattr(DIMPLE, "maxfrag", None),
            primer_buffer=getattr(DIMPLE, "primerBuffer", PRIMER_BUFFER_BASE),
            random_seed=getattr(DIMPLE, "random_seed", 0),
            avoid_sequence=list(getattr(DIMPLE, "avoid_sequence", [])),
            barcode_f=list(getattr(DIMPLE, "barcodeF", [])),
            barcode_r=list(getattr(DIMPLE, "barcodeR", [])),
            dms=getattr(DIMPLE, "dms", False),
            stop_codon=getattr(DIMPLE, "stop_codon", False),
            make_double=getattr(DIMPLE, "make_double", False),
            maximize_nucleotide_change=getattr(
                DIMPLE, "maximize_nucleotide_change", False
            ),
            gene_primer_tm=getattr(DIMPLE, "gene_primerTm", (58, 62)),
            non_interactive=getattr(DIMPLE, "non_interactive", False),
            preferred_orf_index=getattr(DIMPLE, "preferred_orf_index", None),
            link_policy=getattr(DIMPLE, "link_policy", "prompt"),
            breaksite_change_policy=getattr(
                DIMPLE, "breaksite_change_policy", "prompt"
            ),
        )
    return _RUNTIME_CONFIG


def reset_runtime_config() -> None:
    """Reset global runtime config cache (primarily for tests)."""
    global _RUNTIME_CONFIG
    _RUNTIME_CONFIG = None


def _sync_dimple_from_config(config: DimpleRuntimeConfig) -> None:
    """Dual-write helper during migration from class attrs to runtime config."""
    DIMPLE.handle = config.handle
    DIMPLE.usage = config.usage
    if config.cutsite is not None:
        DIMPLE.cutsite = config.cutsite
    if config.cutsite_buffer is not None:
        DIMPLE.cutsite_buffer = config.cutsite_buffer
    DIMPLE.cutsite_overhang = config.cutsite_overhang
    DIMPLE.enzyme = config.enzyme
    if config.synth_len is not None:
        DIMPLE.synth_len = config.synth_len
    if config.maxfrag is not None:
        DIMPLE.maxfrag = config.maxfrag
    DIMPLE.primerBuffer = config.primer_buffer
    DIMPLE.random_seed = config.random_seed
    if config.avoid_sequence is not None:
        DIMPLE.avoid_sequence = list(config.avoid_sequence)
    if config.barcode_f is not None:
        DIMPLE.barcodeF = list(config.barcode_f)
    if config.barcode_r is not None:
        DIMPLE.barcodeR = list(config.barcode_r)
    DIMPLE.dms = config.dms
    DIMPLE.stop_codon = config.stop_codon
    DIMPLE.make_double = config.make_double
    DIMPLE.maximize_nucleotide_change = config.maximize_nucleotide_change
    DIMPLE.gene_primerTm = config.gene_primer_tm
    DIMPLE.non_interactive = config.non_interactive
    DIMPLE.preferred_orf_index = config.preferred_orf_index
    DIMPLE.link_policy = config.link_policy
    DIMPLE.breaksite_change_policy = config.breaksite_change_policy


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


def apply_handle(handle: str, config: Optional[DimpleRuntimeConfig] = None) -> None:
    """Set ``DIMPLE.handle`` after validation."""
    validate_handle(handle)
    runtime = config or get_runtime_config()
    runtime.handle = handle
    _sync_dimple_from_config(runtime)


def resolve_codon_usage(
    usage_arg: Union[str, dict], config: Optional[DimpleRuntimeConfig] = None
) -> dict:
    """Resolve codon usage and assign ``DIMPLE.usage``.

    * ``\"ecoli\"`` / ``\"human\"``: built-in tables via :func:`codon_usage`.
    * Path string (any other str): file contents ``ast.literal_eval`` to a dict.
    * ``dict``: used directly (e.g. GUI custom table).

    Args:
        usage_arg: Named preset, path to a Python-literal dict file, or dict.

    Returns:
        The resolved codon usage table.
    """
    runtime = config or get_runtime_config()
    if isinstance(usage_arg, dict):
        runtime.usage = usage_arg
        _sync_dimple_from_config(runtime)
        return DIMPLE.usage
    if usage_arg in ("ecoli", "human"):
        runtime.usage = codon_usage(usage_arg)
        _sync_dimple_from_config(runtime)
        return DIMPLE.usage
    with open(usage_arg, encoding="utf-8") as f:
        text = f.read()
    runtime.usage = ast.literal_eval(text.strip())
    _sync_dimple_from_config(runtime)
    return DIMPLE.usage


def apply_restriction_settings(
    restriction_sequence: str, config: Optional[DimpleRuntimeConfig] = None
) -> None:
    """Parse restriction enzyme / site string and set cutsite-related class attrs."""
    runtime = config or get_runtime_config()
    if re.match(r"[ACGT]+\([ACGT]\)\d+/\d+", restriction_sequence):
        tmp_cutsite = restriction_sequence.split("(")
        runtime.cutsite = Seq(tmp_cutsite[0])
        runtime.cutsite_buffer = Seq(tmp_cutsite[1].split(")")[0])
        tmp_overhang = tmp_cutsite[1].split(")")[1].split("/")
        runtime.cutsite_overhang = int(tmp_overhang[1]) - int(tmp_overhang[0])
        if (
            runtime.cutsite == Seq("GGTCTC")
            and runtime.cutsite_buffer == Seq("G")
            and runtime.cutsite_overhang == 4
        ):
            runtime.enzyme = "BsaI"
        elif (
            runtime.cutsite == Seq("CGTCTC")
            and runtime.cutsite_buffer == Seq("G")
            and runtime.cutsite_overhang == 4
        ):
            runtime.enzyme = "BsmBI"
        else:
            runtime.enzyme = None
    elif restriction_sequence.upper() in ("BSAI", "BSMBI"):
        if restriction_sequence.upper() == "BSAI":
            runtime.cutsite = Seq("GGTCTC")
            runtime.cutsite_buffer = Seq("G")
            runtime.cutsite_overhang = 4
            runtime.enzyme = "BsaI"
        else:
            runtime.cutsite = Seq("CGTCTC")
            runtime.cutsite_buffer = Seq("G")
            runtime.cutsite_overhang = 4
            runtime.enzyme = "BsmBI"
    else:
        raise ValueError(
            f"Restriction sequence {restriction_sequence!r} not recognized. "
            "Please check input."
        )
    _sync_dimple_from_config(runtime)


def compute_overlaps_and_maxfrag(
    oligo_len: int,
    fragment_len: int,
    overlap: int,
    deletions: Union[bool, List[int]],
    logger: Optional[logging.Logger] = None,
    config: Optional[DimpleRuntimeConfig] = None,
) -> tuple[int, int]:
    """Set ``DIMPLE.synth_len``, ``DIMPLE.maxfrag``, and primer buffer from overlap.

    Resets ``DIMPLE.primerBuffer`` to ``PRIMER_BUFFER_BASE + overlapL`` so repeated
    configuration in the same process does not accumulate.

    Returns:
        ``(overlapL, overlapR)`` after deletion adjustment.
    """
    runtime = config or get_runtime_config()
    overlap_l = int(overlap)
    overlap_r = int(overlap)
    if deletions:
        overlap_r = max(int(x) for x in deletions) + overlap_r - 3

    runtime.synth_len = oligo_len
    if fragment_len != 0:
        runtime.maxfrag = fragment_len
        if logger:
            logger.info("Maximum fragment length: %s based on input", runtime.maxfrag)
    else:
        runtime.maxfrag = oligo_len - 64 - overlap_l - overlap_r
        if logger:
            logger.info(
                "Maximum fragment length: %s based on oligo length and overlap: 2 * %s",
                runtime.maxfrag,
                overlap,
            )

    runtime.primer_buffer = PRIMER_BUFFER_BASE + overlap_l
    _sync_dimple_from_config(runtime)
    return overlap_l, overlap_r


def apply_barcode_start(start: int, config: Optional[DimpleRuntimeConfig] = None) -> None:
    """Slice barcode primer lists from *start*."""
    runtime = config or get_runtime_config()
    idx = int(start)
    runtime.barcode_f = DIMPLE.barcodeF[idx:]
    runtime.barcode_r = DIMPLE.barcodeR[idx:]
    _sync_dimple_from_config(runtime)


def normalize_avoid_list(
    sequences: Union[List[str], str],
    logger: Optional[logging.Logger] = None,
    config: Optional[DimpleRuntimeConfig] = None,
) -> List[Seq]:
    """Build ``DIMPLE.avoid_sequence`` as ``Seq`` objects; strip whitespace (GUI comma lists)."""
    if isinstance(sequences, str):
        if sequences.strip() == "":
            raw: List[str] = []
        else:
            raw = [s.strip() for s in sequences.split(",")]
    else:
        raw = [str(s).strip() for s in sequences]

    runtime = config or get_runtime_config()
    runtime.avoid_sequence = [Seq(x) for x in raw]

    if runtime.cutsite not in runtime.avoid_sequence:
        runtime.avoid_sequence.append(runtime.cutsite)
        msg = (
            f"Restriction sequence {runtime.cutsite} was not included in the avoid list. "
            "Adding before continuing."
        )
        warnings.warn(msg)
        if logger:
            logger.warning(msg)

    _sync_dimple_from_config(runtime)
    return DIMPLE.avoid_sequence


def apply_random_seed(seed: Optional[int], config: Optional[DimpleRuntimeConfig] = None) -> None:
    """Set ``DIMPLE.random_seed`` (``None`` means nondeterministic / NumPy default)."""
    runtime = config or get_runtime_config()
    runtime.random_seed = seed
    _sync_dimple_from_config(runtime)


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
    config: Optional[DimpleRuntimeConfig] = None,
) -> None:
    """Apply runtime behavior policy flags using config dual-write."""
    runtime = config or get_runtime_config()
    runtime.dms = dms
    runtime.stop_codon = stop_codon
    runtime.make_double = make_double
    runtime.maximize_nucleotide_change = maximize_nucleotide_change
    runtime.non_interactive = non_interactive
    runtime.preferred_orf_index = preferred_orf_index
    runtime.link_policy = link_policy
    runtime.breaksite_change_policy = breaksite_change_policy
    _sync_dimple_from_config(runtime)


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
    config: Optional[DimpleRuntimeConfig] = None,
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
    runtime = config or get_runtime_config()
    if gene_primer_tm is not None:
        runtime.gene_primer_tm = gene_primer_tm
        _sync_dimple_from_config(runtime)
    for gene in instances:
        if aminoacids is not None:
            gene.aminoacids = [a.strip() for a in aminoacids]
        if doublefrag is not None:
            gene.doublefrag = int(doublefrag)
