"""Phase 3: hardware-gated clamp-timing tests.

Each test builds a tiny Brian2 cpp_standalone network that drives
``step_clamp`` via ``run_regularly`` (using ``ni_brian2.init_neuron_device``
to open the NI card) and captures helper accessors via
``StateMonitor`` -recorded NeuronGroup variables.

Skipped automatically when the NI Dev1 is not reachable — see
``conftest.py::pytest_collection_modifyitems``.
"""
import os
import sys

import numpy as np
import pytest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NI_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, NI_DIR)

# These must be importable at module scope so they are present in the
# calling frame when Brian2's ``run_regularly`` resolves identifiers.
from ni_test_helpers import (  # noqa: E402,F401
    ni_steps_taken,
    ni_last_read_t,
    ni_total_debt,
)
from ni_brian2 import init_neuron_device, step_clamp  # noqa: E402,F401

pytestmark = pytest.mark.hardware


def _make_clamped_neuron(dt, capture_helpers=True):
    """Build a 1-neuron NeuronGroup whose ``v`` is driven by step_clamp.

    Returns (G, monitor) where the monitor records all sampled variables
    each timestep. Caller is responsible for ``run(...)``.
    """
    from brian2 import NeuronGroup, StateMonitor, ms, mV, pA, second

    model = """
        v : volt
        I_in : amp
        steps_seen : 1
        last_t_seen : second
        debt_seen : second
    """
    G = NeuronGroup(1, model=model, method="euler")
    G.I_in = 0 * pA
    # step_clamp must run before the helper sampling so we observe
    # post-step state.
    G.run_regularly("v = step_clamp(t, I_in)", dt=dt, when="start")
    if capture_helpers:
        G.run_regularly(
            """
            steps_seen = ni_steps_taken(t)
            last_t_seen = ni_last_read_t(t)
            debt_seen = ni_total_debt(t)
            """,
            dt=dt,
            when="end",
        )
    mon = StateMonitor(G, ["steps_seen", "last_t_seen", "debt_seen", "v"], record=True)
    return G, mon


def test_step_count_matches_calls(brian2_standalone):
    """After N timesteps, ``get_steps_taken()`` reports ~N.

    Brian2 may schedule the ``when='start'`` slot one fewer time than
    ``N = run_duration / dt`` depending on how the network's clock
    aligns with the run-window boundary, so we tolerate ±2.
    """
    from brian2 import ms, defaultclock, run

    dt = 0.1 * ms
    defaultclock.dt = dt
    init_neuron_device(brian2_standalone, dt=dt, runtime=0.2)

    G, mon = _make_clamped_neuron(dt)
    N = 100
    run(N * dt)

    final_steps = int(mon.steps_seen[0][-1])
    # Monotone, positive, close to N.
    assert final_steps > 0
    assert abs(final_steps - N) <= 2, (
        f"expected ~{N} step_clamp invocations, got {final_steps}"
    )
    # Monitor recorded a sample per timestep — confirm count too.
    assert mon.steps_seen.shape[1] >= N - 2


def test_last_read_t_anchored_to_target(brian2_standalone):
    """``LAST_READ_T`` should advance by exactly dt each step.

    Validates the target-time anchoring (``LAST_READ_T += step_time_net``)
    documented in interface.cpp. We assert the per-step delta matches dt
    within 50 µs across the run.
    """
    from brian2 import ms, defaultclock, second, run

    dt = 0.1 * ms
    defaultclock.dt = dt
    init_neuron_device(brian2_standalone, dt=dt, runtime=0.2)

    G, mon = _make_clamped_neuron(dt)
    N = 500
    run(N * dt)

    last_t = np.asarray(mon.last_t_seen[0] / second)  # in seconds
    # Drop the first sample: it is captured before any step_clamp call
    # has executed, so it still reflects reset_state()'s anchor.
    last_t = last_t[1:]
    deltas = np.diff(last_t)
    dt_s = float(dt / second)
    drift = np.abs(deltas - dt_s)
    assert drift.max() < 50e-6, (
        f"per-step drift {drift.max()*1e6:.2f} µs exceeds 50 µs tolerance"
    )


def test_debt_stays_bounded_for_short_run(brian2_standalone, on_rt_kernel):
    """Total debt should remain small over a short, in-rate run.

    On PREEMPT_RT this should be ≤ a few ms over 1000 dt=0.1ms steps.
    On non-RT kernels we relax to ≤ 50ms — the goal is just to catch
    catastrophic regressions where debt explodes.
    """
    from brian2 import ms, defaultclock, second, run

    dt = 0.1 * ms
    defaultclock.dt = dt
    init_neuron_device(brian2_standalone, dt=dt, runtime=0.2)

    G, mon = _make_clamped_neuron(dt)
    N = 1000
    run(N * dt)

    debt = float(mon.debt_seen[0][-1] / second)
    # On RT kernels with sub-µs precision, debt should stay well under the
    # total run wall-time. On stock kernels, DAQmx overhead alone can take
    # several ms per step, so we only assert debt grows at most linearly
    # with N and never exceeds ~10× the total network time — anything
    # beyond that signals a regression (e.g. unbounded busy-wait).
    network_time_s = float((N * dt) / second)
    limit = 5e-3 if on_rt_kernel else max(50e-3, 10 * network_time_s)
    assert 0.0 <= debt < limit, (
        f"debt {debt*1e3:.2f} ms exceeds {limit*1e3:.1f} ms limit "
        f"(rt_kernel={on_rt_kernel}, N*dt={network_time_s*1e3:.1f} ms)"
    )


def test_jitter_histogram(brian2_standalone, on_rt_kernel):
    """Inter-step wall-time jitter measured via last_read_t deltas.

    Note: ``last_read_t`` is target-anchored, not measured. To get a
    *real* jitter estimate we use ``total_debt`` deltas instead, which
    track ``step_time_real - step_time_net`` whenever the loop falls
    behind. p99 should be modest on a real-time kernel; on stock
    kernels we only assert it is finite (informational).
    """
    from brian2 import ms, defaultclock, second, run

    dt = 0.1 * ms
    defaultclock.dt = dt
    init_neuron_device(brian2_standalone, dt=dt, runtime=2.0)

    G, mon = _make_clamped_neuron(dt)
    N = 5000
    run(N * dt)

    debt = np.asarray(mon.debt_seen[0] / second)
    debt_deltas = np.diff(debt)  # per-step "overrun" in seconds
    # Always-positive: when on schedule, delta is ~0; only positive
    # contributions accumulate.
    debt_deltas = np.clip(debt_deltas, 0.0, None)

    p50 = float(np.percentile(debt_deltas, 50))
    p99 = float(np.percentile(debt_deltas, 99))
    print(f"[jitter] p50={p50*1e6:.2f}µs p99={p99*1e6:.2f}µs rt={on_rt_kernel}")

    if on_rt_kernel:
        assert p99 < 200e-6, f"p99 jitter {p99*1e6:.2f}µs exceeds 200µs on RT"
    else:
        # Informational on non-RT kernels — only fail on grossly broken
        # timing (>50 ms per step).
        assert p99 < 50e-3, f"p99 jitter {p99*1e3:.2f}ms is implausibly large"
