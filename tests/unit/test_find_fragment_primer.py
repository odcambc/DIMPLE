"""Unit tests for find_fragment_primer (DIMPLE/DIMPLE.py:529-570).

These tests verify:
  1. A normal-length fragment returns a primer within the configured Tm window.
  2. An adversarially short/AT-rich input terminates in finite time rather than
     looping indefinitely — this is the main regression guard.
"""

from Bio.Seq import Seq

from DIMPLE.DIMPLE import DIMPLE, find_fragment_primer

# 25-nt sequence with ~50% GC — Tm_NN should land close to the primerTm window.
_BALANCED_25MER = Seq("ATCGATCGATCGATCGATCGATCGA")

# 10-nt AT-rich sequence below the 16-nt minimum — exercises the early-exit path.
_SHORT_ATRICH = Seq("ATATATATAT")


class TestFindFragmentPrimerNormal:
    def test_returns_primer_and_float_tm(self):
        """Return type is (Seq, float)."""
        primer, tm = find_fragment_primer(_BALANCED_25MER, stop=25)
        assert isinstance(tm, float)
        assert len(primer) >= 1

    def test_primer_within_minimum_length(self):
        """Returned primer is at least 16 nt when stop allows it."""
        primer, tm = find_fragment_primer(_BALANCED_25MER, stop=25)
        # If the function found a valid Tm primer, len >= 16; if it fell back to
        # the count>12 exit, len == stop (25 >= 16 still).  Either way >= 16.
        assert len(primer) >= 16

    def test_tm_in_configured_window(self):
        """For a balanced-GC 25-mer the Tm should lie within primerTm bounds.

        Tolerance is ±2 °C to account for rounding across the two Tm tables
        (DNA_NN2 and DNA_NN4) and any borderline sequences.
        """
        lo, hi = DIMPLE.primerTm
        primer, tm = find_fragment_primer(_BALANCED_25MER, stop=25)
        assert lo - 2 <= tm <= hi + 2, f"Tm {tm:.2f} °C is well outside primerTm window {lo}–{hi}"


class TestFindFragmentPrimerAdversarial:
    def test_short_input_terminates(self):
        """An AT-rich fragment shorter than the minimum primer length terminates."""
        primer, tm = find_fragment_primer(_SHORT_ATRICH, stop=len(_SHORT_ATRICH))
        # The function must return — no infinite loop.
        assert primer is not None

    def test_stop_less_than_minimum_terminates(self):
        """stop < 16 triggers the count>12 / end>stop exit; must not hang."""
        primer, tm = find_fragment_primer(_BALANCED_25MER, stop=10)
        assert len(primer) <= 10
