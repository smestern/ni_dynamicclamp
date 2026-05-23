"""PyTorch wrapper around the precompiled ``libni.so`` dynamic-clamp backend.

This module bridges PyTorch tensors to the C ``run_step_loop`` entry
point in ``interface.cpp``. The DAQ I/O itself is non-differentiable; we
expose two gradient modes:

* ``"detach"`` (default) — backward returns no gradient for the input
  current. Matches the original ``snn-nidaq/nidaq_SNN.py`` pattern of
  computing the loss against a detached real-cell trace.
* ``"surrogate"`` — backward passes the upstream gradient straight back
  to the input (identity / pass-through). This treats the DAQ as an
  identity I→V map for gradient purposes, which lets gradients flow to
  upstream layers while keeping the forward pass faithful to the real
  cell.

Usage::

    from ni_torch_layer import NiDAQClampLayer
    clamp = NiDAQClampLayer(dt_ms=0.1, grad_mode="surrogate",
                            proxy_spike=True)
    V = clamp(I_batch)          # I_batch shape (T,) or (T, B)

The C side keeps task handles in file-scope globals, so **only one
``NiDAQClampLayer`` may be live per process**. Constructing a second
one before the first is destroyed raises ``RuntimeError``.
"""
from __future__ import annotations

import atexit
import os
import sys
from typing import Optional

import numpy as np
import torch
import torch.nn as nn


# Resolve the ni_interface package directory. Layout: this file lives in
# snn-nidaq/ at the repo root; the interface lives in ../ni_interface/.
_HERE = os.path.dirname(os.path.abspath(__file__))
_NI_DIR = os.path.abspath(os.path.join(_HERE, os.pardir, "ni_interface"))
if _NI_DIR not in sys.path:
    sys.path.insert(0, _NI_DIR)


def _select_backend(use_dummy: bool):
    """Return the (module, name) tuple for the chosen backend."""
    if use_dummy:
        import ni_torch_dummy as ni  # type: ignore
        return ni, "ni_torch_dummy"
    import ni_torch as ni  # type: ignore
    return ni, "ni_torch"


# Process-wide singleton bookkeeping. The C side has no notion of
# multiple clamps, so we enforce single-instance here.
_active_layer: Optional["NiDAQClampLayer"] = None


def _atexit_cleanup():
    global _active_layer
    if _active_layer is not None:
        try:
            _active_layer.close()
        except Exception:
            pass


atexit.register(_atexit_cleanup)


class _NiDAQClampFn(torch.autograd.Function):
    """autograd.Function that calls ``run_step_loop`` once per batch elem.

    Forward
    -------
    Input: ``I`` shape ``(T,)`` or ``(T, B)``, float tensor, CPU or GPU.
    Output: ``V`` same shape on the same device.

    Backward
    --------
    ``detach`` mode -> returns ``None`` for ``I``'s gradient.
    ``surrogate`` mode -> returns ``grad_V`` unchanged.
    """

    @staticmethod
    def forward(ctx, I: torch.Tensor, layer: "NiDAQClampLayer"):
        ctx.grad_mode = layer.grad_mode
        ctx.save_for_backward(I)

        if I.dim() == 1:
            I_2d = I.unsqueeze(1)        # (T, 1)
            squeeze_out = True
        elif I.dim() == 2:
            I_2d = I
            squeeze_out = False
        else:
            raise ValueError(
                f"NiDAQClampLayer expects 1-D or 2-D input (T,) or (T, B); "
                f"got shape {tuple(I.shape)}"
            )

        T, B = I_2d.shape
        I_cpu = I_2d.detach().to(torch.float64).cpu().contiguous().numpy()
        V_out = np.empty((T, B), dtype=np.float64)
        spk_out = np.empty((T, B), dtype=np.float64) if layer.proxy_spike else None

        for b in range(B):
            V_b, spk_b = layer._ni.run_step_loop(
                I_cpu[:, b],
                dt_ms=layer.dt_ms,
                t0=layer._t_cursor,
                return_spikes=layer.proxy_spike,
            )
            V_out[:, b] = V_b
            if spk_b is not None:
                spk_out[:, b] = spk_b
            # Advance the network-time cursor so successive batches keep
            # a monotonic ``t``. Pacing is driven by libni's LAST_READ_T,
            # so this is informational, but it keeps proxy-spike DEBUG
            # printouts readable.
            layer._t_cursor += T * (layer.dt_ms / 1000.0)

        V = torch.from_numpy(V_out).to(I.dtype).to(I.device)
        spk_tensor = (
            torch.from_numpy(spk_out).to(I.dtype).to(I.device)
            if spk_out is not None else None
        )
        if squeeze_out:
            V = V.squeeze(1)
            if spk_tensor is not None:
                spk_tensor = spk_tensor.squeeze(1)
        layer._last_spikes = spk_tensor

        return V

    @staticmethod
    def backward(ctx, grad_V: torch.Tensor):
        if ctx.grad_mode == "surrogate":
            # Identity pass-through: dV/dI = 1.
            return grad_V, None
        # "detach" mode: no gradient w.r.t. I.
        return None, None


class NiDAQClampLayer(nn.Module):
    """``nn.Module`` wrapping the dynamic-clamp DAQ as a forward pass.

    Parameters
    ----------
    dt_ms : float
        Per-step time advance in milliseconds; **must** match the
        upstream network's integration step.
    scalefactor_in, scalefactor_out : float
        Amplifier-gain conversions (see ``init_ni`` docstring).
    runtime : float
        Expected total wall time (seconds); only used for the DEBUG
        build's ``read_times`` preallocation.
    ai_chan, ao_chan : str
        DAQmx channel strings.
    proxy_spike : bool
        Whether to enable libni's proxy-spike workaround for noisy
        real-cell voltages. When True, the per-step spike flag is
        captured into ``self.last_spikes`` after each forward pass.
    vthresh_mV, vreset_mV : float
        Proxy-spike thresholds (used only when ``proxy_spike=True``).
    grad_mode : {"detach", "surrogate"}
        Backward behaviour. See module docstring.
    use_dummy : bool
        If True, use ``ni_torch_dummy`` (pure-numpy stub) instead of
        loading ``libni.so``. Useful for CI/dev machines without an
        NI card. Auto-detect: if libni.so is missing, fall back to
        dummy and emit a warning.
    """

    _VALID_GRAD_MODES = ("detach", "surrogate")

    def __init__(
        self,
        dt_ms: float,
        *,
        scalefactor_in: float = 0.1,
        scalefactor_out: float = 1 / 0.5,
        runtime: float = 1.0,
        ai_chan: str = "Dev1/ai0",
        ao_chan: str = "Dev1/ao0",
        proxy_spike: bool = False,
        vthresh_mV: float = -30.0,
        vreset_mV: float = -70.0,
        grad_mode: str = "detach",
        use_dummy: bool = False,
    ):
        super().__init__()
        if grad_mode not in self._VALID_GRAD_MODES:
            raise ValueError(
                f"grad_mode must be one of {self._VALID_GRAD_MODES}, "
                f"got {grad_mode!r}"
            )
        global _active_layer
        if _active_layer is not None:
            raise RuntimeError(
                "Another NiDAQClampLayer is already active in this process. "
                "libni.so holds NI task handles in file-scope globals; "
                "call .close() on the existing layer first."
            )

        # Pick backend. Auto-fall back to dummy if libni.so is missing.
        backend_dummy = use_dummy
        if not backend_dummy:
            if not os.path.exists(os.path.join(_NI_DIR, "libni.so")):
                import warnings
                warnings.warn(
                    f"libni.so not found in {_NI_DIR}; falling back to "
                    f"ni_torch_dummy. Build with "
                    f"`cd ni_interface && bash compile_cpp` for real DAQ.",
                    RuntimeWarning,
                )
                backend_dummy = True

        self._ni, self._backend_name = _select_backend(backend_dummy)

        self.dt_ms = float(dt_ms)
        self.proxy_spike = bool(proxy_spike)
        self.grad_mode = grad_mode

        self._ni.init_ni(
            dt_ms=self.dt_ms,
            scalefactor_in=scalefactor_in,
            scalefactor_out=scalefactor_out,
            runtime=runtime,
            ai_chan=ai_chan,
            ao_chan=ao_chan,
        )
        if self.proxy_spike:
            self._ni.turn_on_proxy_spike(vthresh_mV, vreset_mV)

        self._t_cursor: float = 0.0
        self._last_spikes: Optional[torch.Tensor] = None
        self._closed = False
        _active_layer = self

    # ---- properties ----

    @property
    def last_spikes(self) -> Optional[torch.Tensor]:
        """Proxy-spike flag tensor from the most recent forward pass.

        Shape matches the forward output; ``None`` if ``proxy_spike=False``.
        """
        return self._last_spikes

    @property
    def backend_name(self) -> str:
        return self._backend_name

    # ---- core ----

    def forward(self, I: torch.Tensor) -> torch.Tensor:
        if self._closed:
            raise RuntimeError("NiDAQClampLayer has been closed.")
        return _NiDAQClampFn.apply(I, self)

    def reset_time(self, t0: float = 0.0):
        """Reset the network-time cursor passed to ``run_step_loop``."""
        self._t_cursor = float(t0)

    def close(self):
        global _active_layer
        if self._closed:
            return
        try:
            self._ni.clean_up()
        finally:
            self._closed = True
            if _active_layer is self:
                _active_layer = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass


__all__ = ["NiDAQClampLayer"]
