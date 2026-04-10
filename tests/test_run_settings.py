"""Unit tests for :mod:`DIMPLE.run_settings`."""

from __future__ import annotations

import logging
import os
import tempfile
import unittest
import warnings

from Bio import SeqIO
from Bio.Seq import Seq

from DIMPLE.DIMPLE import DIMPLE
from DIMPLE.run_settings import (
    DimpleRuntimeConfig,
    PRIMER_BUFFER_BASE,
    apply_barcode_start,
    apply_instance_settings,
    apply_random_seed,
    apply_restriction_settings,
    apply_runtime_policies,
    compute_overlaps_and_maxfrag,
    configure_dimple_logging,
    get_runtime_config,
    normalize_avoid_list,
    reset_runtime_config,
    resolve_codon_usage,
    validate_handle,
    validate_insertions,
)
from DIMPLE.utilities import codon_usage


def _reload_barcodes() -> None:
    """Restore barcode lists from packaged data (tests may slice them)."""
    data_dir = os.path.join(
        os.path.dirname(__file__), "..", "DIMPLE", "data"
    )
    fwd = os.path.join(data_dir, "forward_finalprimers.fasta")
    rev = os.path.join(data_dir, "reverse_finalprimers.fasta")
    with open(fwd, encoding="utf-8") as handle:
        DIMPLE.barcodeF = list(SeqIO.parse(handle, "fasta"))
    with open(rev, encoding="utf-8") as handle:
        DIMPLE.barcodeR = list(SeqIO.parse(handle, "fasta"))


class _FakeGene:
    """Minimal stand-in for :class:`DIMPLE` instances in ``apply_instance_settings``."""

    def __init__(self) -> None:
        self.aminoacids = ["Cys", "Asp"]
        self.doublefrag = 0


class TestRunSettings(unittest.TestCase):
    def setUp(self) -> None:
        reset_runtime_config()
        _reload_barcodes()
        DIMPLE.primerBuffer = PRIMER_BUFFER_BASE
        DIMPLE.enzyme = None
        DIMPLE.random_seed = 0
        DIMPLE.gene_primerTm = (58, 62)
        DIMPLE.non_interactive = False
        DIMPLE.breaksite_change_policy = "prompt"
        DIMPLE.dms = False

    def test_validate_handle_accepts_acgt(self) -> None:
        validate_handle("ACGTacgt")

    def test_validate_handle_rejects_non_nucleotide(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            validate_handle("ACGTN")
        self.assertIn("non-nucleic", str(ctx.exception).lower())

    def test_resolve_codon_usage_ecoli_matches_utilities(self) -> None:
        resolve_codon_usage("ecoli")
        self.assertEqual(DIMPLE.usage, codon_usage("ecoli"))
        self.assertEqual(DIMPLE.usage["TTT"], 0.58)

    def test_resolve_codon_usage_human_matches_utilities(self) -> None:
        resolve_codon_usage("human")
        self.assertEqual(DIMPLE.usage, codon_usage("human"))
        self.assertEqual(DIMPLE.usage["TTT"], 0.45)

    def test_resolve_codon_usage_dict(self) -> None:
        custom = {"TTT": 0.5, "TTC": 0.5}
        resolve_codon_usage(custom)
        self.assertIs(DIMPLE.usage, custom)

    def test_resolve_codon_usage_from_file(self) -> None:
        table = {"TTT": 0.11, "TTC": 0.89}
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".txt",
            delete=False,
            encoding="utf-8",
        ) as f:
            f.write(str(table))
            path = f.name
        try:
            resolve_codon_usage(path)
            self.assertEqual(DIMPLE.usage["TTT"], 0.11)
            self.assertEqual(DIMPLE.usage["TTC"], 0.89)
        finally:
            os.unlink(path)

    def test_apply_restriction_settings_pattern_bsmbi(self) -> None:
        apply_restriction_settings("CGTCTC(G)1/5")
        self.assertEqual(DIMPLE.cutsite, Seq("CGTCTC"))
        self.assertEqual(DIMPLE.cutsite_buffer, Seq("G"))
        self.assertEqual(DIMPLE.cutsite_overhang, 4)
        self.assertEqual(DIMPLE.enzyme, "BsmBI")

    def test_apply_restriction_settings_keyword_bsai(self) -> None:
        apply_restriction_settings("BsaI")
        self.assertEqual(DIMPLE.enzyme, "BsaI")
        self.assertEqual(DIMPLE.cutsite, Seq("GGTCTC"))

    def test_apply_restriction_settings_invalid(self) -> None:
        with self.assertRaises(ValueError):
            apply_restriction_settings("NOT_AN_ENZYME")

    def test_compute_overlaps_and_maxfrag(self) -> None:
        overlap_l, overlap_r = compute_overlaps_and_maxfrag(
            oligo_len=230,
            fragment_len=0,
            overlap=4,
            deletions=False,
        )
        self.assertEqual(overlap_l, 4)
        self.assertEqual(overlap_r, 4)
        self.assertEqual(DIMPLE.synth_len, 230)
        self.assertEqual(DIMPLE.maxfrag, 230 - 64 - 4 - 4)
        self.assertEqual(DIMPLE.primerBuffer, PRIMER_BUFFER_BASE + 4)

    def test_compute_overlaps_and_maxfrag_deletions_adjust_overlap_r(self) -> None:
        overlap_l, overlap_r = compute_overlaps_and_maxfrag(
            oligo_len=230,
            fragment_len=100,
            overlap=4,
            deletions=[9, 3, 6],
        )
        self.assertEqual(overlap_l, 4)
        self.assertEqual(overlap_r, 9 + 4 - 3)
        self.assertEqual(DIMPLE.maxfrag, 100)

    def test_apply_barcode_start_slices(self) -> None:
        n0 = len(DIMPLE.barcodeF)
        apply_barcode_start(2)
        self.assertEqual(len(DIMPLE.barcodeF), n0 - 2)
        self.assertEqual(len(DIMPLE.barcodeR), n0 - 2)

    def test_normalize_avoid_list_strips_whitespace_in_comma_string(self) -> None:
        apply_restriction_settings("CGTCTC(G)1/5")
        with warnings.catch_warnings(record=True) as wrec:
            warnings.simplefilter("always")
            seqs = normalize_avoid_list("CGTCTC, GGTCTC")
        self.assertEqual(len(seqs), 2)
        self.assertEqual(seqs[0], Seq("CGTCTC"))
        self.assertEqual(seqs[1], Seq("GGTCTC"))

    def test_normalize_avoid_list_cli_list_appends_cutsite(self) -> None:
        apply_restriction_settings("CGTCTC(G)1/5")
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            normalize_avoid_list(["GGTCTC"])
        self.assertIn(Seq("CGTCTC"), DIMPLE.avoid_sequence)

    def test_apply_random_seed(self) -> None:
        apply_random_seed(42)
        self.assertEqual(DIMPLE.random_seed, 42)
        apply_random_seed(None)
        self.assertIsNone(DIMPLE.random_seed)

    def test_get_runtime_config_reflects_current_class_state(self) -> None:
        DIMPLE.non_interactive = True
        DIMPLE.link_policy = "never"
        cfg = get_runtime_config()
        self.assertIsInstance(cfg, DimpleRuntimeConfig)
        self.assertEqual(cfg.non_interactive, DIMPLE.non_interactive)
        self.assertEqual(cfg.link_policy, DIMPLE.link_policy)

    def test_apply_runtime_policies_dual_write(self) -> None:
        cfg = get_runtime_config()
        apply_runtime_policies(
            dms=True,
            stop_codon=True,
            make_double=False,
            maximize_nucleotide_change=True,
            non_interactive=True,
            preferred_orf_index=2,
            link_policy="never",
            breaksite_change_policy="warn",
            config=cfg,
        )
        self.assertTrue(DIMPLE.dms)
        self.assertTrue(DIMPLE.non_interactive)
        self.assertEqual(DIMPLE.preferred_orf_index, 2)
        self.assertEqual(DIMPLE.link_policy, "never")
        self.assertEqual(DIMPLE.breaksite_change_policy, "warn")

    def test_validate_insertions_ok(self) -> None:
        apply_restriction_settings("CGTCTC(G)1/5")
        validate_insertions(["GAC", "GGGCCC"])

    def test_validate_insertions_bad_bases(self) -> None:
        apply_restriction_settings("CGTCTC(G)1/5")
        with self.assertRaises(ValueError) as ctx:
            validate_insertions(["GAC", "GAX"])
        self.assertIn("non-nucleic", str(ctx.exception).lower())

    def test_validate_insertions_contains_cutsite(self) -> None:
        apply_restriction_settings("CGTCTC(G)1/5")
        with self.assertRaises(ValueError) as ctx:
            validate_insertions(["CGTCTCGGG"])
        self.assertIn("restriction", str(ctx.exception).lower())

    def test_apply_instance_settings(self) -> None:
        g = _FakeGene()
        apply_instance_settings(
            [g],
            aminoacids=["Ala", " Gly "],
            doublefrag=1,
            gene_primer_tm=(55, 65),
        )
        self.assertEqual(g.aminoacids, ["Ala", "Gly"])
        self.assertEqual(g.doublefrag, 1)
        self.assertEqual(DIMPLE.gene_primerTm, (55, 65))

    def test_apply_instance_settings_noops_when_none(self) -> None:
        g = _FakeGene()
        prior_aa = list(g.aminoacids)
        apply_instance_settings([g])
        self.assertEqual(g.aminoacids, prior_aa)

    def test_configure_dimple_logging_creates_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            log_path = os.path.join(tmp, "sub", "dimple.log")
            configure_dimple_logging(log_path)
            self.assertTrue(os.path.isfile(log_path))
            logger = logging.getLogger("tests.configure_dimple_logging")
            logger.info("test")
            self.assertGreater(os.path.getsize(log_path), 0)

    def test_breaksites_non_interactive_error_policy_raises(self) -> None:
        gene = DIMPLE.__new__(DIMPLE)
        gene.geneid = "dummy"
        gene._DIMPLE__breaksites = [30, 60, 90]
        DIMPLE.non_interactive = True
        DIMPLE.breaksite_change_policy = "error"
        with self.assertRaises(ValueError):
            gene.breaksites = [33, 60, 90]

    def test_breaksites_non_interactive_warn_policy_allows(self) -> None:
        gene = DIMPLE.__new__(DIMPLE)
        gene.geneid = "dummy"
        gene._DIMPLE__breaksites = [30, 60, 90]
        DIMPLE.non_interactive = True
        DIMPLE.breaksite_change_policy = "warn"
        gene.breaksites = [33, 60, 90]
        self.assertEqual(gene.breaksites, [33, 60, 90])


if __name__ == "__main__":
    unittest.main()
