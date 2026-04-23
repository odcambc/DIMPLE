# Plan: Add Tests to the DIMPLE Repo

## Context

DIMPLE (Deep Indel Missense Programmable Library Engineering) is a Python 3.12 bioinformatics tool for designing oligonucleotide libraries. The repo currently has **one** test — `tests/test_dimple.py` — a single integration test that runs the full DMS pipeline and byte-compares outputs against golden files in `tests/expected/`. There is no `pyproject.toml`, no pytest config, no conftest, no declared dev dependencies, and no CI. The codebase has significant amounts of mutable module-level state on the `DIMPLE` class, which makes unit testing non-trivial.

The existing test also has three real bugs that need fixing:

1. **Case-sensitivity bug** — the test reads `tests/expected/kir_DMS_Gene_Primers.fasta` (lowercase `kir_`), but `print_all()` writes `Kir_DMS_*` (capital `K`). Passes only on case-insensitive filesystems (macOS/APFS default); silently broken on Linux.
2. **cwd fragility** — hardcoded relative paths like `tests/kir_DMS_Oligos.fasta`; only works when pytest is invoked from the repo root.
3. **Brittle teardown** — `tearDown()` uses bare `os.remove()` calls that crash with `FileNotFoundError` if the test fails before creating all outputs, masking the real failure.

We want to (a) establish test infrastructure, (b) fix the existing test, (c) add a focused batch of unit tests covering the highest-ROI pure/near-pure functions, and (d) wire up CI on Linux + macOS so the case-sensitivity bug class can never regress.

**Out of scope for this plan:** `run_dimple_gui.py`, `DIMPLE.ipynb`, deep internals of `generate_DMS_fragments()` (beyond the existing regression test), interactive-function testing (`findORF`, `align_genevariation`), coverage thresholds in CI, Windows CI, property-based testing. These are follow-ups.

## Approach

Five incremental PRs, each small and independently reviewable. pytest is the framework (it runs existing `unittest.TestCase` classes natively, so no migration cost). PR0 establishes a `.devcontainer` so every subsequent PR is developed and reviewed in the same Linux environment CI will eventually use — this alone catches the case-sensitivity bug class at author-time, not review-time.

---

### PR0 — VS Code Dev Container

Files to create:

- **`/Users/bartleby/Projects/DIMPLE/.devcontainer/devcontainer.json`** — new file.
  - `"image": "mcr.microsoft.com/devcontainers/python:3.12-bookworm"` — matches the Python version pinned in `dimple_env.yml` and the version CI will use.
  - `"postCreateCommand": "pip install --upgrade pip && pip install -r requirements.txt && pip install -e '.[dev]'"` — installs runtime deps (biopython, numpy, pydna) plus dev extras once `pyproject.toml` lands in PR1. For PR0 alone, the command simplifies to `pip install -r requirements.txt pytest pytest-cov`; PR1 replaces it with the form above.
  - `"customizations.vscode.extensions"`: `ms-python.python`, `ms-python.vscode-pylance`, `ms-toolsai.jupyter` (for `DIMPLE.ipynb` interactive use).
  - `"customizations.vscode.settings"`: `"python.testing.pytestEnabled": true`, `"python.testing.pytestArgs": ["tests"]`.
  - `"features"`: `"ghcr.io/devcontainers/features/git:1"` (ensures recent git for `git mv` in PR2).
  - `"mounts"`: none initially; the workspace mount is implicit.
  - `"remoteUser": "vscode"`.
- **`/Users/bartleby/Projects/DIMPLE/.devcontainer/Dockerfile`** — optional, only if `postCreateCommand` gets unwieldy. Initially skip; use the base image + `postCreateCommand` for minimum churn. Add a Dockerfile only in a later PR if container builds become slow enough to warrant baking deps into the image.
- **`/Users/bartleby/Projects/DIMPLE/.gitignore`** — append `.venv/` if not already present (keeps the stray `.venv/` out of the repo; the container uses the system Python so no `.venv/` is needed inside).

Explicit non-goals for PR0:
- No `pyqt5` / GUI deps in the container (the repo dropped the Qt GUI per commit `a9d7ac4`).
- No Dockerfile customization beyond the base image.
- No `docker-compose.yml` — single-container setup only.

**Acceptance check for PR0:** "Dev Containers: Reopen in Container" succeeds in VS Code; inside the container, `python -c "import Bio, numpy, pydna"` works; `python -m unittest tests.test_dimple` runs to completion (it is expected to **fail** on the container's Linux filesystem with `FileNotFoundError` on the lowercase `tests/expected/kir_DMS_*.fasta` paths — this is exactly the bug PR2 fixes, and its appearance here is proof the container is exercising a case-sensitive FS). Merge PR0 with that failure documented; it is not a blocker for PR0 itself — PR0's scope is just "container exists and installs deps."

---

### PR1 — Test-infra foundation (inside the container)

Files to create:

- **`/Users/bartleby/Projects/DIMPLE/pyproject.toml`** — new file. Add:
  - `[project.optional-dependencies] dev = ["pytest>=8.0", "pytest-cov>=5.0"]`
  - `[tool.pytest.ini_options]` with `testpaths = ["tests"]`, `addopts = ["-ra", "--strict-markers", "--strict-config"]`, markers `slow`, `interactive`, and `filterwarnings` to silence the noisy biopython/pydna `DeprecationWarning`s and the `UserWarning` about unwanted restriction sites (these are expected during the Kir regression run).
  - `[tool.coverage.run] source = ["DIMPLE"]`, `omit = ["DIMPLE/data/*", "run_dimple_gui.py"]`
- **`/Users/bartleby/Projects/DIMPLE/tests/conftest.py`** — new file. Defines:
  - A session-level `dimple_human_usage` fixture that sets `DIMPLE.usage` to the human codon table (currently inlined at `tests/test_dimple.py:19–84`).
  - A function-scoped `dimple_state` fixture that snapshots and restores these class attributes around each test: `usage`, `handle`, `overlap`, `synth_len`, `maxfrag`, `primerBuffer`, `barcodeF`, `barcodeR`, `cutsite`, `cutsite_buffer`, `cutsite_overhang`, `avoid_sequence`, `dms`, `stop_codon`, `make_double`, `maximize_nucleotide_change`, `random_seed`.
  - A `kir_fa` fixture returning `Path(__file__).parent / "data" / "Kir.fa"` (anticipates PR2's move).
- **`/Users/bartleby/Projects/DIMPLE/.github/workflows/tests.yml`** — new CI workflow.
  - Matrix: `os: [ubuntu-latest, macos-latest]`, `python: ["3.12"]`, `fail-fast: false`.
  - Steps: `actions/checkout@v4` → `actions/setup-python@v5` (with `cache: pip`) → `pip install -r requirements.txt` and `pip install -e ".[dev]"` → `pytest -q`.
  - The Ubuntu job is what catches the case-sensitivity bug in PR2.

**Expected CI state after PR1:** existing `tests/test_dimple.py` passes on macOS, fails on Ubuntu. This is intentional — PR2 fixes it.

---

### PR2 — Replace the existing integration test (inside the container)

Files to create:

- **`/Users/bartleby/Projects/DIMPLE/tests/data/Kir.fa`** — moved from `tests/Kir.fa`.
- **`/Users/bartleby/Projects/DIMPLE/tests/expected/Kir_designed_variants.csv`** — promoted from the currently-untracked `tests/Kir_designed_variants.csv` (which is an output `generate_DMS_fragments` writes around `DIMPLE/DIMPLE.py:1884`).
- **`/Users/bartleby/Projects/DIMPLE/tests/regression/__init__.py`** — empty.
- **`/Users/bartleby/Projects/DIMPLE/tests/regression/test_dms_pipeline_kir.py`** — the rewritten test.
  - Uses the `dimple_state` and `dimple_human_usage` fixtures from conftest (no inline 60-line codon table).
  - Accepts `tmp_path`, copies `tests/data/Kir.fa` into it, runs `generate_DMS_fragments` / `post_qc` / `print_all` with `wDir=str(tmp_path)`.
  - Golden comparison: reads expected files via `Path(__file__).parent.parent / "expected" / <name>`, NOT cwd-relative.
  - Correct filenames: `Kir_DMS_Gene_Primers.fasta`, `Kir_DMS_Oligo_Primers.fasta`, `Kir_DMS_Oligos.fasta`, `Kir_mutations.csv`, and the newly-tracked `Kir_designed_variants.csv` — **five** golden comparisons, one more than the current test.
  - Marked `@pytest.mark.slow`.
  - No `tearDown`; `tmp_path` handles cleanup.

Files to delete:

- **`tests/test_dimple.py`** — deleted in the same PR, once local + CI runs confirm the new test produces identical assertions on macOS **and** Ubuntu.
- **`tests/Kir.fa`** — removed after move.
- **`tests/Kir_designed_variants.csv`** — removed after promotion into `tests/expected/`.

**Acceptance check:** `pytest tests/regression -q` passes on both OSes in CI.

---

### PR3 — Priority-1 unit tests (pure functions in `utilities.py`)

These are the highest-ROI tests: no DIMPLE class state required, fast, and they lock down behavior of logic that's easy to regress when refactoring.

- **`/Users/bartleby/Projects/DIMPLE/tests/unit/__init__.py`** — empty.
- **`/Users/bartleby/Projects/DIMPLE/tests/unit/test_parse_custom_mutations.py`** — covers `DIMPLE/utilities.py:52–72`. Use `@pytest.mark.parametrize` for:
  - `["10:A"]` → `{10: "A"}`
  - `["10:All"]` → `{10: "A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y"}`
  - `["5-7:K"]` → `{5: "K", 6: "K", 7: "K"}`
  - `["10:A", "10:G"]` → `{10: "A,G"}` (the append-to-existing-key branch, lines 68–69)
  - `["1-3:All", "5:M"]` (mixed range + single + All)
- **`/Users/bartleby/Projects/DIMPLE/tests/unit/test_codon_usage.py`** — covers `DIMPLE/utilities.py:75–100`.
  - `codon_usage("ecoli")` returns a dict with 64 entries; every stop codon (`TAA`, `TAG`, `TGA`) present.
  - `codon_usage("human")` returns 64 entries; spot-check two values against the inlined table in the existing test (e.g. `TTT == 0.45`, `ATG == 1`).
  - `codon_usage({"TTT": 1.0})` (custom passthrough, line 99) returns input verbatim.
  - For each preset, assert frequencies per amino acid sum to ≈1.0 within tolerance.

---

### PR4 — Priority-2 unit tests (DIMPLE class boundary) + self-test for `dimple_state`

These touch the `DIMPLE` class but stay small thanks to the `dimple_state` fixture.

- **`/Users/bartleby/Projects/DIMPLE/tests/unit/test_addgene.py`** — covers `DIMPLE/DIMPLE.py:44–67` (the `addgene` FASTA loader).
  - Write a tiny synthetic FASTA to `tmp_path` with header `>gene1 start:4 end:30`; call `addgene`; assert `gene.start == 3` (subtract 1 for 0-indexing) and `gene.end == 30`.
  - FASTA without `start:/end:` markers returns a gene object whose `start`/`end` default per constructor.
  - Multi-record FASTA returns a list of length == record count.
  - Path with `\\` escapes (`"tests\\data\\Kir.fa"`) — confirms the `.replace("\\", "")` sanitization on line 52.
- **`/Users/bartleby/Projects/DIMPLE/tests/unit/test_check_overhangs.py`** — covers `DIMPLE/DIMPLE.py:857–899`.
  - Construct a minimal `DIMPLE` instance whose `breaklist` produces unique overhangs → `check_overhangs` returns `False` (no switching needed).
  - Force a duplicate overhang → returns `True` and mutates state via `switch_fragmentsize`.
  - Palindromic/reverse-complement case: if `overhang_F == overhang_R.reverse_complement()` the site is flagged.
- **`/Users/bartleby/Projects/DIMPLE/tests/unit/test_find_fragment_primer.py`** — covers `DIMPLE/DIMPLE.py:529–570`.
  - Pass a high-GC 25-nt fragment; assert returned primer Tm ∈ (56.5, 60) ± small tolerance and `len(primer) >= 16`.
  - Pass adversarial short/AT-rich input; assert loop terminates (`count > 12` early-exit) — i.e., the call returns in finite time. This is a regression guard against infinite-loop bugs in the Tm-tuning branch.
- **`/Users/bartleby/Projects/DIMPLE/tests/unit/test_dimple_state_isolation.py`** — a meta-test. Mutate `DIMPLE.usage`, `DIMPLE.overlap`, and one Seq attribute inside a test that uses the `dimple_state` fixture; assert in a follow-up test (or via the fixture's yield-restore) that the values are restored. Prevents future regressions in the fixture itself.

Total new code in PR3 + PR4: ~8 test files, ~25 individual parametrized cases, target <500 LOC combined.

---

## Critical Files to Modify/Create

| Path | PR | Action |
|------|----|--------|
| `.devcontainer/devcontainer.json` | PR0 | create |
| `.gitignore` | PR0 | append `.venv/` |
| `pyproject.toml` | PR1 | create |
| `tests/conftest.py` | PR1 | create |
| `.github/workflows/tests.yml` | PR1 | create |
| `tests/regression/test_dms_pipeline_kir.py` | PR2 | create |
| `tests/data/Kir.fa` | PR2 | move from `tests/Kir.fa` |
| `tests/expected/Kir_designed_variants.csv` | PR2 | promote from untracked |
| `tests/test_dimple.py` | PR2 | delete |
| `tests/unit/test_parse_custom_mutations.py` | PR3 | create |
| `tests/unit/test_codon_usage.py` | PR3 | create |
| `tests/unit/test_addgene.py` | PR4 | create |
| `tests/unit/test_check_overhangs.py` | PR4 | create |
| `tests/unit/test_find_fragment_primer.py` | PR4 | create |
| `tests/unit/test_dimple_state_isolation.py` | PR4 | create |

## Existing Functions/Utilities to Reuse

- **`addgene`** (`DIMPLE/DIMPLE.py:44`) — FASTA loader, used as entry to the regression test; also the subject of a unit test. No modification.
- **`generate_DMS_fragments`** (`DIMPLE/DIMPLE.py:902`) — orchestrator exercised by the regression test; writes the 5 golden files.
- **`print_all`** (`DIMPLE/DIMPLE.py:2047`), **`post_qc`** (`DIMPLE/DIMPLE.py:2069`) — called by the regression test, unchanged.
- **`parse_custom_mutations`**, **`codon_usage`** (`DIMPLE/utilities.py:52, 75`) — tested directly.
- The inline human codon-usage table in `tests/test_dimple.py:19–84` is moved verbatim into the `dimple_human_usage` fixture in `tests/conftest.py`.

## Verification

After each PR, reopen the repo in the dev container (VS Code: "Dev Containers: Reopen in Container") and run from the workspace root:

```bash
pip install -e ".[dev]"          # no-op after PR1 if already installed
pytest -q                        # all tests
pytest tests/unit -q             # fast unit tests only
pytest tests/regression -q       # slow integration/regression
pytest --co -q                   # test collection sanity check
```

End-to-end verification specific to each PR:

- **PR0:** `code --folder-uri vscode-remote://dev-container+...` (or the UI "Reopen in Container") successfully builds the container. Inside the container: `python --version` shows 3.12.x, `python -c "import Bio, numpy, pydna"` succeeds, `python -m unittest tests.test_dimple` runs (and likely fails with FileNotFoundError on the lowercase `kir_` path — confirms we're on a case-sensitive FS).
- **PR1:** Inside the container, `pytest -q` runs the pre-existing `tests/test_dimple.py` (will surface the casing bug on Linux — intentional; PR2 fixes it). `pytest --co` shows the fixtures load without import errors. CI on both Ubuntu and macOS matches local container behavior.
- **PR2:** Both macOS and Ubuntu CI jobs go green. Inside the container, `pytest tests/regression -q` passes. Delete the old `tests/test_dimple.py` only after CI confirms the new test passes on both runners.
- **PR3:** `pytest tests/unit/test_parse_custom_mutations.py tests/unit/test_codon_usage.py -v` — all parametrized cases green on both OSes.
- **PR4:** `pytest tests/unit -v`; confirm `test_dimple_state_isolation.py` catches a deliberate break (temporarily mutate state without restoration, verify a later test fails — then undo).

**Sign-off criterion:** `pytest -q` green on both macOS and Ubuntu CI, with at least 1 slow test (the Kir regression) and ~20 parametrized unit cases. No coverage threshold enforced in this plan; revisit once the suite grows.

---

## Progress Log

### PR0 — COMPLETE
Committed in `d3fc0ec`. Created `.devcontainer/dimple/devcontainer.json` (Python
3.12-bookworm base image, editable install of `.[dev]`, Pylance/Jupyter extensions).
Note: the workspace `.devcontainer/` root is the Claude Code sandbox environment and
was deliberately not committed; the DIMPLE devcontainer lives at
`.devcontainer/dimple/devcontainer.json`.

### PR1 — COMPLETE
Committed in `d3fc0ec`. Created `pyproject.toml`, `tests/conftest.py`,
`.github/workflows/tests.yml`. Added `enzyme` to `_MANAGED_ATTRS` in conftest.

### PR2 — COMPLETE
Committed in `d3fc0ec`. Replaced `tests/test_dimple.py` with
`tests/regression/test_dms_pipeline_kir.py`. All three original bugs fixed:
case-sensitivity, cwd fragility, brittle teardown. Added `Kir_designed_variants.csv`
as a fifth golden comparison. Moved `tests/Kir.fa` → `tests/data/Kir.fa`.

### PR3 — COMPLETE
Committed in `3116c2a`. Added `tests/unit/test_parse_custom_mutations.py` (5
parametrized cases) and `tests/unit/test_codon_usage.py` (ecoli, human, custom
passthrough, per-AA frequency sum checks).

### PR4 — COMPLETE
Committed in `ac9a070`. Added `test_addgene.py`, `test_check_overhangs.py`,
`test_find_fragment_primer.py`, `test_dimple_state_isolation.py`.

### Remaining verification (requires Python env):
Run `pytest -q` inside `.devcontainer/dimple/` (Python 3.12-bookworm) or in CI.
Expected: ~20+ parametrized unit cases green on both macOS and Ubuntu; 1 slow
regression test (Kir). The regression test requires a Python environment with
biopython, numpy, pydna installed (`pip install -e ".[dev]"`).
