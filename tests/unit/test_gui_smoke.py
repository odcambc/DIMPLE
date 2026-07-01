"""Smoke test for the Tk GUI entrypoint (``run_dimple_gui.py``).

The GUI is never exercised by the rest of the suite, so it's the easiest
entrypoint to silently break when the shared pipeline API drifts. This module
covers:

1. The module imports cleanly -- catches stale ``from DIMPLE...`` imports the
   same way ``test_notebook_smoke.py`` catches stale notebook call sites.
2. ``run()`` maps GUI widget values into ``build_runtime_config`` against the
   current helper signature, producing a valid ``DimpleRuntimeConfig``
   (``addgene`` stubbed, so the heavy pipeline never runs).
3. ``run()`` drives the *whole* pipeline (addgene -> post_qc -> print_all) to
   completion on the real Kir fixture and emits DMS outputs -- guards
   shared-API drift the config-only test can't see (e.g. ``post_qc`` /
   ``print_all`` being handed a ``config=`` kwarg they don't accept). ``@slow``.
4. The real ``Application`` widget default has DMS selected (needs a Tk display;
   skipped where none is available).

Tests 1-3 need no display: importing the module doesn't instantiate
``Application`` (the Tk root is created only under ``__main__``), and ``run()``
is fed a duck-typed fake ``app`` whose ``.get()`` accessors return the GUI's own
default values.
"""

import types

import pytest
from Bio.Seq import Seq

import run_dimple_gui as gui
from DIMPLE.pool import DEFAULT_RANDOM_SEED, DimpleRuntimeConfig


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


def _fake_app(geneFile="/nonexistent/Kir.fa", wDir=None):
    """A duck-typed ``app`` carrying the GUI's default widget values (DMS run).

    Mirrors every ``.get()`` / attribute ``run()`` reads, including the ones
    past ``addgene`` (``substitutions``, Tm bounds, ``doublefrag``,
    ``avoid_breaksites`` ...), so the same fake can drive either a
    config-only run (with ``addgene`` stubbed) or the full pipeline.
    """
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
        # Fields read after addgene (apply_instance_settings + generate).
        substitutions=_Var(
            "Cys,Asp,Ser,Gln,Met,Asn,Pro,Lys,Thr,Phe,Ala,Gly,Ile,Leu,His,Arg,Trp,Val,Glu,Tyr"
        ),
        melting_temp_low=_Var("58"),
        melting_temp_high=_Var("62"),
        doublefrag=_Var(0),
        avoid_breaksites=_Var(0),
        matchSequences=_Var(0),
        custom_mutations={},
        avoid_others_list=_Var(""),
        geneFile=geneFile,
        wDir=wDir,
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
    monkeypatch.setattr(gui, "messagebox", types.SimpleNamespace(showerror=lambda *a, **k: None))

    with pytest.raises(_StopBeforePipeline):
        gui.run()

    config = captured["config"]
    assert isinstance(config, DimpleRuntimeConfig)
    assert config.dms is True
    assert config.random_seed == DEFAULT_RANDOM_SEED
    # "CGTCTC(G)1/5" -> recognition site CGTCTC, 4-base overhang.
    assert config.cutsite == Seq("CGTCTC")
    assert config.cutsite_overhang == 4


@pytest.mark.slow
def test_gui_run_completes_full_pipeline(monkeypatch, tmp_path, kir_fa):
    """run() drives the whole GUI pipeline to completion and emits DMS outputs.

    Unlike the config-only test, this does NOT stub addgene, so it reaches
    post_qc/print_all -- catching call-signature drift between the GUI and the
    shared pipeline (e.g. post_qc(pool, config=...) when post_qc takes only
    pool). run() re-raises after showing a messagebox, so any TypeError here
    fails the test.
    """
    wdir = str(tmp_path) + "/"
    monkeypatch.setattr(gui, "app", _fake_app(geneFile=str(kir_fa), wDir=wdir), raising=False)
    monkeypatch.setattr(gui, "messagebox", types.SimpleNamespace(showerror=lambda *a, **k: None))

    gui.run()  # must run clean through post_qc + print_all

    produced = {p.name for p in tmp_path.iterdir()}
    for name in ("Kir_DMS_Oligos.fasta", "Kir_designed_variants.csv", "All_Oligos.fasta"):
        assert name in produced, f"GUI run did not produce {name}; got {sorted(produced)}"


def test_gui_dms_selected_by_default():
    """The real DMS checkbox launches selected.

    Regression: a ``self.include_sub_check.deselect()`` at widget creation used
    to force ``include_substitutions`` back to 0, defeating the default-on
    behavior. The fake-app tests can't see this -- it only shows on the real
    ``Application``, so this instantiates it (skipped where no display exists).
    """
    tk = pytest.importorskip("tkinter")
    try:
        root = tk.Tk()
    except tk.TclError as exc:
        pytest.skip(f"no Tk display available: {exc}")
    try:
        real_app = gui.Application(master=root)
        assert (
            real_app.include_substitutions.get() == 1
        ), "DMS checkbox should launch selected (default-on)"
    finally:
        root.destroy()
