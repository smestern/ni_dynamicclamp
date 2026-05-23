"""Shape / contract tests for ``NiDAQClampLayer``.

All tests force ``use_dummy=True`` so they run anywhere without an NI
card. The dummy backend mirrors the ctypes API of ``ni_torch`` so any
shape regression caught here applies to the real-DAQ build too.
"""
import os
import sys
import warnings

import numpy as np
import pytest
import torch

_HERE = os.path.dirname(os.path.abspath(__file__))
_SNN = os.path.abspath(os.path.join(_HERE, os.pardir))
if _SNN not in sys.path:
    sys.path.insert(0, _SNN)

from ni_torch_layer import NiDAQClampLayer  # noqa: E402


@pytest.fixture
def fresh_layer():
    """Yield a fresh dummy-backed layer and guarantee teardown.

    ``proxy_spike=True`` so ``last_spikes`` is populated for the shape
    contract tests.
    """
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    layer = NiDAQClampLayer(
        dt_ms=0.1, grad_mode="detach", use_dummy=True, proxy_spike=True,
    )
    try:
        yield layer
    finally:
        layer.close()


def test_1d_input_shape_contract(fresh_layer):
    I = torch.linspace(-20.0, 20.0, 100)
    V = fresh_layer(I)
    assert V.shape == (100,)
    assert torch.isfinite(V).all()
    assert fresh_layer.last_spikes.shape == (100,)


def test_2d_input_shape_contract(fresh_layer):
    I = torch.zeros(50, 4)
    V = fresh_layer(I)
    assert V.shape == (50, 4)
    assert fresh_layer.last_spikes.shape == (50, 4)


def test_detach_grad_mode_blocks_grad(fresh_layer):
    I = torch.full((32,), 10.0, requires_grad=True)
    V = fresh_layer(I)
    V.sum().backward()
    assert I.grad is None


def test_surrogate_grad_mode_passes_through():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    layer = NiDAQClampLayer(
        dt_ms=0.1, grad_mode="surrogate", use_dummy=True,
    )
    try:
        I = torch.full((32,), 10.0, requires_grad=True)
        V = layer(I)
        V.sum().backward()
        assert I.grad is not None
        # straight-through estimator: grad of sum wrt I_i is 1.
        assert torch.allclose(I.grad, torch.ones_like(I))
    finally:
        layer.close()


def test_singleton_enforcement():
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    a = NiDAQClampLayer(dt_ms=0.1, use_dummy=True)
    try:
        with pytest.raises(RuntimeError, match="already active"):
            NiDAQClampLayer(dt_ms=0.1, use_dummy=True)
    finally:
        a.close()


def test_time_cursor_advances(fresh_layer):
    t0 = fresh_layer._t_cursor
    _ = fresh_layer(torch.zeros(100))
    expected = t0 + 100 * (0.1 / 1000.0)
    assert fresh_layer._t_cursor == pytest.approx(expected, rel=1e-9)


def test_reset_time_zeros_cursor(fresh_layer):
    _ = fresh_layer(torch.zeros(10))
    fresh_layer.reset_time(0.0)
    assert fresh_layer._t_cursor == 0.0


def test_backend_name_is_dummy(fresh_layer):
    assert "dummy" in fresh_layer.backend_name.lower()
