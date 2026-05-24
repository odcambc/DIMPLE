"""Smoke test for the Colab notebook (``DIMPLE.ipynb``).

Parses the notebook as JSON, locates the pipeline cell by header marker, and
execs it against the Kir FASTA fixture in a temp working directory. Doesn't
assert on output byte-equivalence — just that the cell runs to completion
against the current API and produces the expected DMS outputs.

Catches stale API call sites the next time the notebook drifts from the public
surface. Uses stdlib json (the .ipynb format is stable JSON) so this test has
no dev-dep dependency beyond pytest itself.
"""

import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = REPO_ROOT / "DIMPLE.ipynb"
PIPELINE_CELL_MARKER = "3. Configure parameters and run DIMPLE"


def _find_pipeline_cell(nb_dict):
    for cell in nb_dict["cells"]:
        if cell["cell_type"] != "code":
            continue
        src = cell["source"]
        if isinstance(src, list):
            src = "".join(src)
        if PIPELINE_CELL_MARKER in src:
            return src
    raise AssertionError(f"No code cell containing {PIPELINE_CELL_MARKER!r} in {NOTEBOOK}")


@pytest.mark.slow
def test_notebook_pipeline_cell_runs(tmp_path, kir_fa, monkeypatch):
    nb = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    source = _find_pipeline_cell(nb)

    work = tmp_path / "workspace"
    work.mkdir()
    (work / "Kir.fa").write_bytes(kir_fa.read_bytes())

    monkeypatch.chdir(REPO_ROOT)

    namespace = {"target_file": "Kir.fa", "directory": str(work) + "/"}
    exec(compile(source, str(NOTEBOOK), "exec"), namespace)

    assert (work / "Kir_DMS_Oligos.fasta").exists()
    assert (work / "Kir_designed_variants.csv").exists()
    assert (work / "Kir_mutations.csv").exists()
    assert (work / "All_Oligos.fasta").exists()
    assert (work / "All_Primers.fasta").exists()
