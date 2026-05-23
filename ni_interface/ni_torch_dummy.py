"""No-DAQ stub for ``ni_torch``. Mirrors the public API so unit tests
and CI machines without the NI card can import and run shape-level
checks. Uses a trivial leaky-integrator model in pure numpy to produce
a plausible voltage trace; do **not** rely on it for physiology.

Use::

    from ni_interface import ni_torch_dummy as ni_torch
    ni_torch.init_ni(dt_ms=0.1)
    V, spk = ni_torch.run_step_loop(I_array, dt_ms=0.1)
"""
import numpy as np


class NiTorchError(RuntimeError):
    pass


_state = {
    "initialised": False,
    "dt_ms": None,
    "scale_in": None,
    "scale_out": None,
    "ai_chan": None,
    "ao_chan": None,
    "v": -0.070,           # membrane voltage in volts (matches `data` default)
    "proxy_spike": False,
    "vthresh": 0.0,
    "vreset": 0.0,
    "last_spike": 0,
}

# Simple passive-cell time constants (seconds).
_TAU = 20e-3
_R = 100e6  # ohms, matches 100 MΩ headstage assumption
_EL = -0.070


def init_ni(dt_ms, scalefactor_in=0.1, scalefactor_out=1 / 0.5,
            runtime=1.0, ai_chan="Dev1/ai0", ao_chan="Dev1/ao0"):
    if _state["initialised"]:
        raise NiTorchError("init_ni() already called; call clean_up() first.")
    _state.update(
        initialised=True,
        dt_ms=float(dt_ms),
        scale_in=float(scalefactor_in),
        scale_out=float(scalefactor_out),
        ai_chan=ai_chan,
        ao_chan=ao_chan,
        v=_EL,
        last_spike=0,
    )


def turn_on_proxy_spike(vthresh_mV=-30.0, vreset_mV=-70.0):
    _state["proxy_spike"] = True
    _state["vthresh"] = float(vthresh_mV) / 1000.0
    _state["vreset"] = float(vreset_mV) / 1000.0


def step_clamp(t_s, I_pA):
    if not _state["initialised"]:
        raise NiTorchError("init_ni() must be called first.")
    dt_s = _state["dt_ms"] / 1000.0
    v = _state["v"]
    I_amps = I_pA * 1e-12
    dv = (-(v - _EL) + I_amps * _R) / _TAU * dt_s
    v = v + dv
    _state["v"] = v
    # Note: returns V in volts pre-scaling, matching `data * SF_IN` in
    # interface.cpp where data is in volts. The Brian2 path multiplies
    # by SF_IN; here we return v * SF_IN for parity.
    return v * _state["scale_in"]


def run_step_loop(I, dt_ms, t0=0.0, return_spikes=False):
    if not _state["initialised"]:
        raise NiTorchError("init_ni() must be called before run_step_loop().")
    I = np.ascontiguousarray(I, dtype=np.float64)
    if I.ndim != 1:
        raise ValueError(f"I must be 1-D, got shape {I.shape}")

    V = np.empty_like(I)
    spk = np.empty_like(I) if return_spikes else None

    dt_s = float(dt_ms) / 1000.0
    v = _state["v"]
    last_spike = _state["last_spike"]
    vthresh = _state["vthresh"]
    vreset = _state["vreset"]
    sf_in = _state["scale_in"]
    proxy = _state["proxy_spike"]

    for i in range(I.shape[0]):
        I_amps = I[i] * 1e-12
        v = v + (-(v - _EL) + I_amps * _R) / _TAU * dt_s
        if proxy:
            if v > vthresh:
                if last_spike == 0:
                    last_spike = 1
                else:
                    v = vreset
            else:
                last_spike = 0
        V[i] = v * sf_in
        if spk is not None:
            spk[i] = float(last_spike)

    _state["v"] = v
    _state["last_spike"] = last_spike
    return V, spk


def clean_up():
    if not _state["initialised"]:
        return 0.0
    _state["initialised"] = False
    return 0.0


def is_initialised():
    return _state["initialised"]


__all__ = [
    "NiTorchError",
    "init_ni",
    "turn_on_proxy_spike",
    "step_clamp",
    "run_step_loop",
    "clean_up",
    "is_initialised",
]
