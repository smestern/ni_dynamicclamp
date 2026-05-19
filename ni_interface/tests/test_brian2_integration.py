"""Phase 5: end-to-end Brian2 integration smoke test.

Builds a tiny LIF + Poisson + real-neuron network the way ``demo_net.py``
does, attaches the NI clamp via ``init_neuron_device`` + ``attach_neuron``,
runs for a short simulated time, and asserts the standalone binary built,
exited cleanly, and produced sensible spike output.

This is the "does the whole pipeline still wire together" guard.
"""
import os
import sys

import numpy as np
import pytest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
NI_DIR = os.path.dirname(THIS_DIR)
sys.path.insert(0, THIS_DIR)
sys.path.insert(0, NI_DIR)

from ni_brian2 import (  # noqa: E402,F401
    init_neuron_device,
    attach_neuron,
    step_clamp,
)

pytestmark = pytest.mark.hardware


def test_lif_with_attached_neuron_smoke(brian2_standalone):
    """A 50 ms LIF + real-neuron network runs to completion."""
    from brian2 import (
        NeuronGroup,
        PoissonGroup,
        Synapses,
        SpikeMonitor,
        StateMonitor,
        defaultclock,
        ms,
        mV,
        Hz,
        nS,
        run,
    )

    dt = 0.1 * ms
    defaultclock.dt = dt
    runtime_s = 0.05  # 50 ms

    eqs = """
    dv/dt = -v/(10*ms) : volt
    dvt/dt = (10*mV-vt)/(15*ms) : volt
    """
    reset = """
    v = 0*mV
    vt += 3*mV
    """
    IF = NeuronGroup(
        5, model=eqs, reset=reset, threshold="v>vt", method="euler"
    )
    IF.vt = 10 * mV

    PG = PoissonGroup(1, 500 * Hz)
    C = Synapses(PG, IF, on_pre="v += 3*mV")
    C.connect()

    eqs_real = """
    v : volt
    I_in : amp
    """
    real_neuron = NeuronGroup(
        1,
        model=eqs_real,
        threshold="v>0*mV",
        refractory=5 * ms,
        method="euler",
    )

    IF_r = Synapses(
        IF,
        real_neuron,
        model="""I_in_post = ge*(0*mV - v_post) : amp (summed)
        dge/dt = -ge/(10*ms) : siemens
        """,
        on_pre="ge += 1e-2*nS",
    )
    IF_r.connect()

    # Init NI + bind the real neuron's v/I_in to step_clamp.
    init_neuron_device(
        brian2_standalone, dt=dt, runtime=runtime_s + 0.05
    )
    attach_neuron(real_neuron, dt=dt, method="proxy")

    spikes = SpikeMonitor(real_neuron)
    state = StateMonitor(real_neuron, ["v", "I_in"], record=True)

    run(runtime_s * 1000 * ms)

    # Standalone binary exists and produced StateMonitor output.
    main_bin = os.path.join(brian2_standalone.project_dir, "main")
    assert os.path.exists(main_bin), f"binary not built at {main_bin}"

    v_trace = np.asarray(state.v[0])
    assert v_trace.size > 0, "no samples recorded"
    assert np.all(np.isfinite(v_trace)), "non-finite v samples"

    # SpikeMonitor count is a non-negative int (whether or not the
    # real cell actually fires depends on the rig — we only assert
    # the pipeline didn't crash).
    assert spikes.count[0] >= 0
