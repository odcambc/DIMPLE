"""Smoke test for the command-line entrypoint (``run_dimple.py``).

Invokes the published CLI as a subprocess against the Kir FASTA fixture and
asserts the expected DMS output files appear. Doesn't check byte-equivalence --
just that the real ``python run_dimple.py`` path runs to completion against the
current API and writes its outputs.

Complements ``test_cli.py`` (parser-level) and ``test_notebook_smoke.py``: this
is the only test that exercises ``run_dimple.py``'s actual ``main()`` -- argparse
wiring, ``build_runtime_config`` + ``run_pipeline``, and the cwd-relative log
path. It runs out-of-process so NumPy's global RNG and the ``logs/`` directory
stay isolated from the rest of the suite.

Subprocess cwd is the temp work dir, so outputs and ``logs/`` land there rather
than polluting the repo. ``run_dimple.py`` is invoked by absolute path: Python
puts the script's directory on ``sys.path[0]``, so ``import DIMPLE`` resolves
against the repo root regardless of cwd.
"""

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_DIMPLE = REPO_ROOT / "run_dimple.py"

# Files a single-gene Kir DMS run is expected to produce.
_EXPECTED_OUTPUTS = [
    "Kir_DMS_Oligos.fasta",
    "Kir_DMS_Gene_Primers.fasta",
    "Kir_DMS_Oligo_Primers.fasta",
    "Kir_mutations.csv",
    "Kir_designed_variants.csv",
    "All_Oligos.fasta",
    "All_Primers.fasta",
]


@pytest.mark.slow
def test_cli_dms_run_writes_outputs(tmp_path, kir_fa):
    """``run_dimple.py -DMS`` on Kir produces the expected output files."""
    work = tmp_path / "workspace"
    work.mkdir()
    (work / "Kir.fa").write_bytes(kir_fa.read_bytes())

    result = subprocess.run(
        [
            sys.executable,
            str(RUN_DIMPLE),
            "-geneFile",
            "Kir.fa",
            "-DMS",
            "--non_interactive",
        ],
        cwd=work,
        capture_output=True,
        text=True,
        timeout=300,
    )

    assert result.returncode == 0, (
        f"run_dimple.py exited {result.returncode}\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )

    missing = [name for name in _EXPECTED_OUTPUTS if not (work / name).exists()]
    assert not missing, f"missing expected outputs: {missing}"
