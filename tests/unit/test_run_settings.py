"""Unit tests for :mod:`DIMPLE.run_settings`.

Each ``run_settings`` helper populates a :class:`DimpleRuntimeConfig`; the
tests build a fresh config, apply a helper, and assert on the config object.
"""

from __future__ import annotations

import logging
import os
import tempfile
import unittest
import warnings

from Bio.Seq import Seq

from DIMPLE.DIMPLE import DIMPLE
from DIMPLE.pool import DimpleRuntimeConfig, Pool, PRIMER_BUFFER_BASE
from DIMPLE.run_settings import (
    apply_barcode_start,
    apply_instance_settings,
    apply_random_seed,
    apply_restriction_settings,
    apply_runtime_policies,
    compute_overlaps_and_maxfrag,
    configure_dimple_logging,
    normalize_avoid_list,
    resolve_codon_usage,
    validate_handle,
    validate_insertions,
)
from DIMPLE.utilities import codon_usage


class _FakeGene:
    """Minimal stand-in for :class:`DIMPLE` instances in ``apply_instance_settings``."""

    def __init__(self) -> None:
        self.aminoacids = ["Cys", "Asp"]
        self.doublefrag = 0


class TestRunSettings(unittest.TestCase):
    def test_validate_handle_accepts_acgt(self) -> None:
        validate_handle("ACGTacgt")

    def test_validate_handle_rejects_non_nucleotide(self) -> None:
        with self.assertRaises(ValueError) as ctx:
            validate_handle("ACGTN")
        self.assertIn("non-nucleic", str(ctx.exception).lower())

    def test_resolve_codon_usage_ecoli_matches_utilities(self) -> None:
        cfg = DimpleRuntimeConfig()
        resolve_codon_usage("ecoli", cfg)
        self.assertEqual(cfg.usage, codon_usage("ecoli"))
        self.assertEqual(cfg.usage["TTT"], 0.58)

    def test_resolve_codon_usage_human_matches_utilities(self) -> None:
        cfg = DimpleRuntimeConfig()
        resolve_codon_usage("human", cfg)
        self.assertEqual(cfg.usage, codon_usage("human"))
        self.assertEqual(cfg.usage["TTT"], 0.45)

    def test_resolve_codon_usage_dict(self) -> None:
        custom = {"TTT": 0.5, "TTC": 0.5}
        cfg = DimpleRuntimeConfig()
        resolve_codon_usage(custom, cfg)
        self.assertIs(cfg.usage, custom)

    def test_resolve_codon_usage_from_file(self) -> None:
        table = {"TTT": 0.11, "TTC": 0.89}
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False, encoding="utf-8"
        ) as f:
            f.write(str(table))
            path = f.name
        try:
            cfg = DimpleRuntimeConfig()
            resolve_codon_usage(path, cfg)
            self.assertEqual(cfg.usage["TTT"], 0.11)
            self.assertEqual(cfg.usage["TTC"], 0.89)
        finally:
            os.unlink(path)

    def test_apply_restriction_settings_pattern_bsmbi(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("CGTCTC(G)1/5", cfg)
        self.assertEqual(cfg.cutsite, Seq("CGTCTC"))
        self.assertEqual(cfg.cutsite_buffer, Seq("G"))
        self.assertEqual(cfg.cutsite_overhang, 4)
        self.assertEqual(cfg.enzyme, "BsmBI")

    def test_apply_restriction_settings_keyword_bsai(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("BsaI", cfg)
        self.assertEqual(cfg.enzyme, "BsaI")
        self.assertEqual(cfg.cutsite, Seq("GGTCTC"))

    def test_apply_restriction_settings_invalid(self) -> None:
        with self.assertRaises(ValueError):
            apply_restriction_settings("NOT_AN_ENZYME", DimpleRuntimeConfig())

    def test_compute_overlaps_and_maxfrag(self) -> None:
        cfg = DimpleRuntimeConfig()
        overlap_l, overlap_r = compute_overlaps_and_maxfrag(
            oligo_len=230, fragment_len=0, overlap=4, deletions=False, config=cfg
        )
        self.assertEqual(overlap_l, 4)
        self.assertEqual(overlap_r, 4)
        self.assertEqual(cfg.synth_len, 230)
        self.assertEqual(cfg.maxfrag, 230 - 64 - 4 - 4)
        self.assertEqual(cfg.primer_buffer, PRIMER_BUFFER_BASE + 4)

    def test_compute_overlaps_and_maxfrag_deletions_adjust_overlap_r(self) -> None:
        cfg = DimpleRuntimeConfig()
        overlap_l, overlap_r = compute_overlaps_and_maxfrag(
            oligo_len=230, fragment_len=100, overlap=4, deletions=[9, 3, 6], config=cfg
        )
        self.assertEqual(overlap_l, 4)
        self.assertEqual(overlap_r, 9 + 4 - 3)
        self.assertEqual(cfg.maxfrag, 100)

    def test_apply_barcode_start_slices(self) -> None:
        cfg = DimpleRuntimeConfig()
        n0 = len(cfg.barcode_f)
        apply_barcode_start(2, cfg)
        self.assertEqual(len(cfg.barcode_f), n0 - 2)
        self.assertEqual(len(cfg.barcode_r), n0 - 2)

    def test_normalize_avoid_list_strips_whitespace_in_comma_string(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("CGTCTC(G)1/5", cfg)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            seqs = normalize_avoid_list("CGTCTC, GGTCTC", cfg)
        self.assertEqual(len(seqs), 2)
        self.assertEqual(seqs[0], Seq("CGTCTC"))
        self.assertEqual(seqs[1], Seq("GGTCTC"))

    def test_normalize_avoid_list_cli_list_appends_cutsite(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("CGTCTC(G)1/5", cfg)
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            normalize_avoid_list(["GGTCTC"], cfg)
        self.assertIn(Seq("CGTCTC"), cfg.avoid_sequence)

    def test_apply_random_seed(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_random_seed(42, cfg)
        self.assertEqual(cfg.random_seed, 42)
        apply_random_seed(None, cfg)
        self.assertIsNone(cfg.random_seed)

    def test_apply_runtime_policies(self) -> None:
        cfg = DimpleRuntimeConfig()
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
        self.assertTrue(cfg.dms)
        self.assertTrue(cfg.non_interactive)
        self.assertEqual(cfg.preferred_orf_index, 2)
        self.assertEqual(cfg.link_policy, "never")
        self.assertEqual(cfg.breaksite_change_policy, "warn")

    def test_validate_insertions_ok(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("CGTCTC(G)1/5", cfg)
        validate_insertions(["GAC", "GGGCCC"], cfg)

    def test_validate_insertions_bad_bases(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("CGTCTC(G)1/5", cfg)
        with self.assertRaises(ValueError) as ctx:
            validate_insertions(["GAC", "GAX"], cfg)
        self.assertIn("non-nucleic", str(ctx.exception).lower())

    def test_validate_insertions_contains_cutsite(self) -> None:
        cfg = DimpleRuntimeConfig()
        apply_restriction_settings("CGTCTC(G)1/5", cfg)
        with self.assertRaises(ValueError) as ctx:
            validate_insertions(["CGTCTCGGG"], cfg)
        self.assertIn("restriction", str(ctx.exception).lower())

    def test_apply_instance_settings(self) -> None:
        cfg = DimpleRuntimeConfig()
        g = _FakeGene()
        apply_instance_settings(
            [g], cfg, aminoacids=["Ala", " Gly "], doublefrag=1, gene_primer_tm=(55, 65)
        )
        self.assertEqual(g.aminoacids, ["Ala", "Gly"])
        self.assertEqual(g.doublefrag, 1)
        self.assertEqual(cfg.gene_primer_tm, (55, 65))

    def test_apply_instance_settings_noops_when_none(self) -> None:
        cfg = DimpleRuntimeConfig()
        g = _FakeGene()
        prior_aa = list(g.aminoacids)
        apply_instance_settings([g], cfg)
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
        config = DimpleRuntimeConfig(
            non_interactive=True, breaksite_change_policy="error"
        )
        gene = DIMPLE.__new__(DIMPLE)
        gene.geneid = "dummy"
        gene.pool = Pool(config)
        gene._DIMPLE__breaksites = [30, 60, 90]
        with self.assertRaises(ValueError):
            gene.breaksites = [33, 60, 90]

    def test_breaksites_non_interactive_warn_policy_allows(self) -> None:
        config = DimpleRuntimeConfig(
            non_interactive=True, breaksite_change_policy="warn"
        )
        gene = DIMPLE.__new__(DIMPLE)
        gene.geneid = "dummy"
        gene.pool = Pool(config)
        gene._DIMPLE__breaksites = [30, 60, 90]
        gene.breaksites = [33, 60, 90]
        self.assertEqual(gene.breaksites, [33, 60, 90])


if __name__ == "__main__":
    unittest.main()
