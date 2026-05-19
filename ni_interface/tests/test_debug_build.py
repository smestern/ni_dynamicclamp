"""Phase 4: exercise the ``-DDEBUG`` build of interface.cpp.

When compiled with ``-DDEBUG``, ``step_clamp`` records the measured
real time of every 1000th step into ``read_times``, and ``clean_up``
dumps the array to a file at the path in the ``NI_READ_TIMES_PATH``
environment variable (default ``./read_times.txt``).

The ``brian2_standalone_debug`` fixture adds ``-DDEBUG`` to Brian2's
``extra_compile_args_gcc`` so the standalone build picks it up
automatically.
"""
import os
import sys

import pytest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NI_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, NI_DIR)

from ni_test_helpers import (  # noqa: E402,F401
    ni_steps_taken,
    ni_last_read_t,
    ni_total_debt,
)
from ni_brian2 import init_neuron_device, step_clamp  # noqa: E402,F401

pytestmark = pytest.mark.hardware


def test_read_times_dump(brian2_standalone_debug, tmp_path, monkeypatch):
    """A long enough -DDEBUG run writes a parseable read_times dump."""
    from brian2 import NeuronGroup, ms, pA, defaultclock, run

    dump_path = tmp_path / "read_times.txt"
    monkeypatch.setenv("NI_READ_TIMES_PATH", str(dump_path))

    dt = 0.1 * ms
    defaultclock.dt = dt
    init_neuron_device(brian2_standalone_debug, dt=dt, runtime=1.0)

    G = NeuronGroup(
        1,
        """
        v : volt
        I_in : amp
        """,
        method="euler",
    )
    G.I_in = 0 * pA
    G.run_regularly("v = step_clamp(t, I_in)", dt=dt, when="start")

    # Need > 1000 steps for read_times to record at least one sample
    # (it samples every 1000 steps; clean_up writes ``steps_taken/1000``
    # lines). Use 3000 to leave headroom against the off-by-one Brian2
    # scheduling.
    N = 3000
    run(N * dt)

    assert dump_path.exists(), f"{dump_path} was not written"
    lines = dump_path.read_text().strip().splitlines()
    assert len(lines) >= 2, f"expected ≥2 samples, got {len(lines)}: {lines}"
    # Each line is a `long double` formatted with %Lf — must parse as float.
    for line in lines:
        v = float(line)
        # Values are wall-clock elapsed times per sampled step; bounded
        # by sane real-world limits.
        assert 0.0 <= v < 10.0, f"implausible read_time {v}"
