"""ctypes wrapper around ``libni.so`` for the PyTorch dynamic-clamp path.

Unlike the Brian2 integration (which has ``cpp_standalone`` codegen
compile ``interface.cpp`` inline with the generated network — see
``ni_brian2.py``), PyTorch needs a **pre-built** shared library that we
hot-load at import time. Build the library once with::

    cd ni_interface && bash compile_cpp

Then import this module:

    from ni_interface import ni_torch
    ni_torch.init_ni(dt_ms=0.1)
    V, spk = ni_torch.run_step_loop(I_array)
    ni_torch.clean_up()

This module intentionally has **no torch import** at module top so it
stays usable from non-torch code paths and from unit tests. The
torch-facing ``autograd.Function`` / ``nn.Module`` wrapper lives in
``snn-nidaq/ni_torch_layer.py``.
"""
import ctypes
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_LIB_PATH = os.path.join(_HERE, "libni.so")


class NiTorchError(RuntimeError):
    """Raised for libni.so load failures or DAQmx errors surfaced from C."""


def _load_lib():
    if not os.path.exists(_LIB_PATH):
        raise NiTorchError(
            f"libni.so not found at {_LIB_PATH}. Build it with "
            f"`cd ni_interface && bash compile_cpp` (or "
            f"`python setuptools_generic.py build_ext --inplace`)."
        )
    try:
        lib = ctypes.CDLL(_LIB_PATH)
    except OSError as exc:
        raise NiTorchError(
            f"Failed to load {_LIB_PATH}: {exc}. Ensure libnidaqmx.so is "
            f"discoverable (a copy lives in {_HERE})."
        ) from exc

    # The libni.so exports two parallel symbol sets:
    #   * C++-mangled names (init_ni, step_clamp, ...) consumed by
    #     Brian2's cpp_standalone codegen via interface.h.
    #   * C-linkage ``ni_*`` forwarding shims (interface_torch_export.cpp)
    #     bound below for the PyTorch / ctypes hot-load path.

    # ---- ni_init_ni(dt_ms, scalein, scaleout, runtime_s, aI, aO) -> int ----
    lib.ni_init_ni.argtypes = [
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_char_p, ctypes.c_char_p,
    ]
    lib.ni_init_ni.restype = ctypes.c_int

    # ---- ni_turn_on_proxy_spike(vthresh_mV, vreset_mV) ----
    # NB: the C signature uses ``long double`` which is 80-bit on x86_64
    # under SysV. ctypes has no portable long-double type, but on x86_64
    # GCC long double passes via x87 registers, so ctypes.c_longdouble
    # works at the ABI level. If this ever becomes a portability issue,
    # we can add a thin ``double``-typed wrapper to interface.cpp.
    lib.ni_turn_on_proxy_spike.argtypes = [
        ctypes.c_longdouble, ctypes.c_longdouble,
    ]
    lib.ni_turn_on_proxy_spike.restype = None

    # ---- ni_step_clamp(t_s, I_pA) -> mV ----
    lib.ni_step_clamp.argtypes = [ctypes.c_double, ctypes.c_double]
    lib.ni_step_clamp.restype = ctypes.c_double

    # ---- ni_run_step_loop(I, V_out, spk_out_or_NULL, n, dt_ms, t0_s) -> int
    _dbl_p = ctypes.POINTER(ctypes.c_double)
    lib.ni_run_step_loop.argtypes = [
        _dbl_p, _dbl_p, _dbl_p,
        ctypes.c_int, ctypes.c_double, ctypes.c_double,
    ]
    lib.ni_run_step_loop.restype = ctypes.c_int

    # ---- ni_clean_up() -> double ----
    lib.ni_clean_up.argtypes = []
    lib.ni_clean_up.restype = ctypes.c_double

    return lib


_lib = None
_initialised = False


def _lib_or_load():
    global _lib
    if _lib is None:
        _lib = _load_lib()
    return _lib


def init_ni(dt_ms, scalefactor_in=0.1, scalefactor_out=1 / 0.5,
            runtime=1.0, ai_chan="Dev1/ai0", ao_chan="Dev1/ao0"):
    """Initialise the NI-DAQ. Mirrors ``ni_brian2.init_neuron_device`` args.

    Parameters
    ----------
    dt_ms : float
        Per-step time advance in milliseconds (matches Brian2's
        ``defaultclock.dt / ms``).
    scalefactor_in, scalefactor_out : float
        Amplifier gain conversions. Defaults assume a 100 MΩ headstage.
    runtime : float
        Expected total runtime in seconds. Only used to preallocate the
        DEBUG-build ``read_times`` buffer.
    ai_chan, ao_chan : str
        DAQmx channel strings. Empty string is treated as "keep default"
        by the C side (see ``init_ni`` in ``interface.cpp``).
    """
    global _initialised
    if _initialised:
        raise NiTorchError(
            "init_ni() already called; libni.so holds NI task handles in "
            "file-scope globals, so only one clamp may be open per process. "
            "Call clean_up() before re-initialising."
        )
    lib = _lib_or_load()
    rc = lib.ni_init_ni(
        float(dt_ms), float(scalefactor_in), float(scalefactor_out),
        float(runtime),
        ai_chan.encode("ascii") if ai_chan else b"",
        ao_chan.encode("ascii") if ao_chan else b"",
    )
    if rc != 0:
        raise NiTorchError(f"init_ni returned {rc} — see stderr from libni.so")
    _initialised = True


def turn_on_proxy_spike(vthresh_mV=-30.0, vreset_mV=-70.0):
    """Enable proxy-spike mode (see interface.cpp). vthresh/vreset in mV."""
    lib = _lib_or_load()
    lib.ni_turn_on_proxy_spike(
        ctypes.c_longdouble(vthresh_mV), ctypes.c_longdouble(vreset_mV),
    )


def step_clamp(t_s, I_pA):
    """Single-step entry point. Returns membrane voltage in V * SF_IN."""
    return _lib_or_load().ni_step_clamp(float(t_s), float(I_pA))


def run_step_loop(I, dt_ms, t0=0.0, return_spikes=False):
    """Drive the dynamic clamp for one batch.

    Parameters
    ----------
    I : np.ndarray
        1-D float64 input currents in pA. Will be made C-contiguous if
        not already.
    dt_ms : float
        Per-step advance in milliseconds (must match init_ni's value).
    t0 : float
        Starting network time in seconds. Informational only — pacing is
        driven by libni's internal ``LAST_READ_T``, not by ``t``.
    return_spikes : bool
        If True, return a second array with the proxy-spike flag per
        step (0/1). Requires ``turn_on_proxy_spike`` was called.

    Returns
    -------
    V : np.ndarray
        Membrane voltage trace, same shape/dtype as ``I``.
    spk : np.ndarray or None
        Proxy-spike trace if ``return_spikes=True``, else ``None``.
    """
    if not _initialised:
        raise NiTorchError("init_ni() must be called before run_step_loop().")

    I = np.ascontiguousarray(I, dtype=np.float64)
    if I.ndim != 1:
        raise ValueError(f"I must be 1-D, got shape {I.shape}")
    n = I.shape[0]

    V = np.empty(n, dtype=np.float64)
    spk = np.empty(n, dtype=np.float64) if return_spikes else None

    _dbl_p = ctypes.POINTER(ctypes.c_double)
    I_p = I.ctypes.data_as(_dbl_p)
    V_p = V.ctypes.data_as(_dbl_p)
    spk_p = spk.ctypes.data_as(_dbl_p) if spk is not None else _dbl_p()

    rc = _lib_or_load().ni_run_step_loop(
        I_p, V_p, spk_p, ctypes.c_int(n),
        ctypes.c_double(dt_ms), ctypes.c_double(t0),
    )
    if rc != 0:
        raise NiTorchError(f"run_step_loop returned {rc}")
    return V, spk


def clean_up():
    """Stop tasks and reset the device. Safe to call multiple times."""
    global _initialised
    if not _initialised:
        return 0.0
    rc = _lib_or_load().ni_clean_up()
    _initialised = False
    return rc


def is_initialised():
    return _initialised


__all__ = [
    "NiTorchError",
    "init_ni",
    "turn_on_proxy_spike",
    "step_clamp",
    "run_step_loop",
    "clean_up",
    "is_initialised",
]
