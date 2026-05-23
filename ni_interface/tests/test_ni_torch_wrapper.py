"""ctypes wrapper smoke tests for ``ni_torch.py`` (the PyTorch hot-load
path against ``libni.so``).

These tests build ``libni.so`` on demand if missing, then exercise the
real C entry points. Anything that requires the NI card is gated on
``@pytest.mark.hardware`` (auto-skipped on CI by ``conftest._probe_ni_hardware``).
"""
import os
import subprocess
import sys

import numpy as np
import pytest

_NI_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
_LIBNI = os.path.join(_NI_DIR, "libni.so")
if _NI_DIR not in sys.path:
    sys.path.insert(0, _NI_DIR)


@pytest.fixture(scope="module")
def libni_built():
    """Ensure libni.so exists; rebuild if missing."""
    if not os.path.exists(_LIBNI):
        subprocess.check_call(
            ["bash", os.path.join(_NI_DIR, "compile_cpp")],
            cwd=_NI_DIR,
        )
    assert os.path.exists(_LIBNI), "compile_cpp did not produce libni.so"
    return _LIBNI


def test_libni_exports_required_symbols(libni_built):
    """Hardware-free: verify the shared library exports the C ABI we
    bind to in ``ni_torch.py``."""
    out = subprocess.check_output(["nm", "-D", libni_built], text=True)
    required = {
        "ni_init_ni",
        "ni_clean_up",
        "ni_step_clamp",
        "ni_run_step_loop",
        "ni_turn_on_proxy_spike",
    }
    found = {line.split()[-1] for line in out.splitlines() if " T " in line}
    missing = required - found
    assert not missing, f"libni.so missing symbols: {missing}"


def test_ni_torch_module_loads(libni_built):
    """Hardware-free: ctypes binding succeeds for all entry points."""
    import importlib
    import ni_torch
    importlib.reload(ni_torch)  # ensure clean state if other tests touched it
    lib = ni_torch._load_lib()
    assert lib.ni_init_ni.restype is not None
    assert lib.ni_run_step_loop.restype is not None
    assert not ni_torch.is_initialised()


@pytest.mark.hardware
def test_run_step_loop_smoke(libni_built):
    """Hardware-gated: init -> 1000-step ramp -> clean_up returns finite V."""
    import importlib
    import ni_torch
    importlib.reload(ni_torch)

    ni_torch.init_ni(dt_ms=0.1, runtime=1.0)
    try:
        I = np.linspace(-50.0, 50.0, 1000).astype(np.float64)
        V, spk = ni_torch.run_step_loop(I, dt_ms=0.1, return_spikes=False)
        assert V.shape == I.shape
        assert np.all(np.isfinite(V))
        assert spk is None
    finally:
        ni_torch.clean_up()


@pytest.mark.hardware
def test_run_step_loop_with_proxy_spikes(libni_built):
    """Hardware-gated: proxy_spike path returns a (0/1) spike trace."""
    import importlib
    import ni_torch
    importlib.reload(ni_torch)

    ni_torch.init_ni(dt_ms=0.1, runtime=1.0)
    try:
        ni_torch.turn_on_proxy_spike(vthresh_mV=-30.0, vreset_mV=-70.0)
        I = np.zeros(500, dtype=np.float64)
        V, spk = ni_torch.run_step_loop(I, dt_ms=0.1, return_spikes=True)
        assert spk is not None and spk.shape == I.shape
        # All values must be 0 or 1.
        assert set(np.unique(spk)).issubset({0.0, 1.0})
    finally:
        ni_torch.clean_up()


def test_init_ni_singleton_guard(libni_built, monkeypatch):
    """Hardware-free: double-init raises rather than corrupting state.

    We monkeypatch the ctypes ``init_ni`` to a no-op success so this
    works without an NI card.
    """
    import importlib
    import ctypes
    import ni_torch
    importlib.reload(ni_torch)

    lib = ni_torch._lib_or_load()
    monkeypatch.setattr(
        lib, "ni_init_ni", lambda *a, **kw: 0,
    )
    monkeypatch.setattr(lib, "ni_clean_up", lambda *a, **kw: 0.0)

    ni_torch.init_ni(dt_ms=0.1)
    try:
        with pytest.raises(ni_torch.NiTorchError):
            ni_torch.init_ni(dt_ms=0.1)
    finally:
        # reset internal flag for subsequent tests
        ni_torch._initialised = False
