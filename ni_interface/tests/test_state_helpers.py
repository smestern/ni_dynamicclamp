"""Phase 2: hardware-free state-helper tests.

These tests build a tiny Brian2 cpp_standalone network that never calls
``step_clamp`` — so DAQmx is never opened — but does call the
test-helper accessors registered in ``ni_test_helpers``. The accessors
operate on plain process globals; linking ``interface.cpp`` does not
touch the NI card at static-init time, so this test runs anywhere
``libnidaqmx.so`` is loadable. If the library itself is missing the
test is skipped.
"""
import os
import sys

import numpy as np
import pytest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS_DIR)

from ni_test_helpers import (  # noqa: E402,F401
    ni_steps_taken,
    ni_total_debt,
    ni_last_read_t,
    ni_reset_state,
)


pytestmark = pytest.mark.skipif(
    not os.path.exists(os.path.join(os.path.dirname(THIS_DIR), "libnidaqmx.so")),
    reason="libnidaqmx.so not present; cannot link interface.cpp",
)


def test_reset_state_clears_counters(brian2_standalone):
    """After ``ni_reset_state()`` the three observable counters read zero."""
    from brian2 import NeuronGroup, ms, second, run

    G = NeuronGroup(
        1,
        """
        did_reset : 1
        steps : 1
        debt : second
        last_t : second
        """,
        method="euler",
    )
    # Run once at t=0 to reset, then sample on subsequent steps.
    G.run_regularly("did_reset = ni_reset_state(t)", dt=1 * ms, when="start")
    G.run_regularly(
        """
        steps = ni_steps_taken(t)
        debt = ni_total_debt(t)
        last_t = ni_last_read_t(t)
        """,
        dt=1 * ms,
        when="end",
    )

    run(2 * ms)

    # Read back recorded state. ``steps`` and ``debt`` are counters that
    # must be zero after a clean reset when ``step_clamp`` was never
    # called. ``last_t`` is intentionally re-anchored to wall-clock
    # ``ClockGetTime()`` by ``reset_state``, so it should be a large
    # positive monotonic-clock reading, not zero.
    assert int(G.steps[0]) == 0
    assert float(G.debt[0] / second) == 0.0
    last_t_s = float(G.last_t[0] / second)
    assert last_t_s > 0.0, "last_read_t was not re-anchored by reset_state"
