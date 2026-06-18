"""Smoke test for the Tk GUI entrypoint (``run_dimple_gui.py``).

The GUI is never exercised by the rest of the suite, so it's the easiest
entrypoint to silently break when the shared pipeline API drifts. This test
covers the two failure modes that don't need a display:

1. The module imports cleanly -- catches stale ``from DIMPLE...`` imports the
   same way ``test_notebook_smoke.py`` catches stale notebook call sites.
2. ``run()`` maps GUI widget values into ``build_runtime_config`` against the
   current helper signature, producing a valid ``DimpleRuntimeConfig``.

No widgets are driven and no display is needed: importing the module doesn't
instantiate ``Application`` (the Tk root is created only under ``__main__``), and
``run()`` is fed a duck-typed fake ``app`` whose ``.get()`` accessors return the
GUI's own default values. ``addgene`` -- the first pipeline call after config
construction -- is stubbed to capture the config and halt, so the heavy pipeline
never runs.
"""

import types

import pytest
from Bio.Seq import Seq

import run_dimple_gui as gui
from DIMPLE.pool import DimpleRuntimeConfig
from DIMPLE.run_settings import DEFAULT_GUI_RANDOM_SEED


class _Var:
    """Stand-in for a Tk ``*Var`` / ``Entry``: ``.get()`` returns a fixed value."""

    def __init__(self, value):
        self._value = value

    def get(self):
        return self._value


class _Sink:
    """Stand-in for the GUI ``output_text`` widget; swallows ``.insert(...)``."""

    def insert(self, *args, **kwargs):
        pass


class _StopBeforePipeline(Exception):
    """Raised by the stubbed ``addgene`` to halt ``run()`` after config build."""


def _fake_app():
    """A duck-typed ``app`` carrying the GUI's default widget values (DMS run)."""
    return types.SimpleNamespace(
        # Mutation-type toggles -- DMS only, so the deletions/insertions
        # parsing branches stay quiet and validation passes.
        include_substitutions=_Var(1),
        delete=_Var(0),
        insert=_Var(0),
        dis=_Var(0),
        make_double=_Var(0),
        stop=_Var(0),
        max_mutations=_Var(0),
        synonymous=_Var(0),
        # Text/entry fields -- the real GUI defaults from Application.__init__.
        deletions=_Var("3,6"),
        insertions=_Var("GGC,GGCTCT,GGCTCTGGA"),
        fragmentLen=_Var("auto"),
        oligoLen=_Var("250"),
        overlap=_Var("4"),
        handle=_Var("AGCGGGAGACCGGGGTCTCTGAGC"),
        restriction_sequence=_Var("CGTCTC(G)1/5"),
        avoid_sequence=_Var("CGTCTC, GGTCTC"),
        barcode_start=_Var("0"),
        codon_usage="human",
        geneFile="/nonexistent/Kir.fa",
        wDir=None,
        output_text=_Sink(),
    )


def test_gui_module_imports():
    """Importing the GUI module resolves all shared-pipeline imports."""
    assert hasattr(gui, "run")
    assert hasattr(gui, "Application")


def test_gui_run_builds_valid_config(monkeypatch):
    """run() maps GUI defaults into a valid DimpleRuntimeConfig via the helper."""
    captured = {}

    def _stub_addgene(gene_file, config):
        captured["gene_file"] = gene_file
        captured["config"] = config
        raise _StopBeforePipeline

    # `app` is assigned only in the __main__ block, so it isn't a module
    # attribute at import time -- raising=False lets monkeypatch create it.
    monkeypatch.setattr(gui, "app", _fake_app(), raising=False)
    monkeypatch.setattr(gui, "addgene", _stub_addgene)
    # Avoid a real messagebox (no display) if an unexpected error path fires.
    monkeypatch.setattr(
        gui, "messagebox", types.SimpleNamespace(showerror=lambda *a, **k: None)
    )

    with pytest.raises(_StopBeforePipeline):
        gui.run()

    config = captured["config"]
    assert isinstance(config, DimpleRuntimeConfig)
    assert config.dms is True
    assert config.random_seed == DEFAULT_GUI_RANDOM_SEED
    # "CGTCTC(G)1/5" -> recognition site CGTCTC, 4-base overhang.
    assert config.cutsite == Seq("CGTCTC")
    assert config.cutsite_overhang == 4
