"""Meta-tests: verify that the dimple_state fixture properly restores DIMPLE class state.

These tests are deliberately order-sensitive (the second test in each pair must run
after the first) and serve as a regression guard on the fixture itself.  If dimple_state
fails to restore an attribute, the second test in the pair will catch it.

DIMPLE.random_seed is used as the canary because it has a known class-level default (0),
making the expected post-restore value deterministic.
"""

import pytest

from DIMPLE.DIMPLE import DIMPLE

_SENTINEL = 0xDEAD_BEEF  # Unlikely to be the real class default for any attribute.


class TestDimpleStateRestoresScalars:
    def test_mutation_is_visible_inside_test(self, dimple_state):
        """State mutations made inside a test are visible within that test."""
        assert DIMPLE.random_seed == 0  # class-level default
        DIMPLE.random_seed = _SENTINEL
        assert DIMPLE.random_seed == _SENTINEL

    def test_state_is_restored_after_previous_test(self, dimple_state):
        """After the previous test exits, dimple_state restores random_seed to 0."""
        assert DIMPLE.random_seed != _SENTINEL, (
            "dimple_state fixture did not restore DIMPLE.random_seed — "
            "state leaked from a previous test"
        )
        assert DIMPLE.random_seed == 0


class TestDimpleStateRestoresNewAttributes:
    def test_new_attribute_added_inside_test(self, dimple_state):
        """An attribute that did not exist before the test can be set and seen."""
        DIMPLE.overlap = 42
        assert DIMPLE.overlap == 42

    def test_new_attribute_removed_after_test(self, dimple_state):
        """An attribute added by the previous test is removed after restoration."""
        # If dimple_state failed to delete the attribute it added, this would
        # see DIMPLE.overlap == 42 rather than raising AttributeError.
        # We accept either outcome: attribute absent OR value != 42.
        try:
            val = DIMPLE.overlap
            assert val != 42, (
                "dimple_state did not remove DIMPLE.overlap that was added by the previous test"
            )
        except AttributeError:
            pass  # Attribute was correctly removed — this is the expected path.
