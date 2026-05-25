"""Hardware-free PyTorch demo for the dynamic-clamp layer.

Runs a toy training loop that drives the ``NiDAQClampLayer`` with the
``ni_torch_dummy`` backend (pure-numpy LIF stub — no NI card required).
Useful for:

* Smoke-testing the autograd plumbing without ``libni.so``.
* Sanity-checking shapes / batching before plugging in real hardware.
* CI on machines without DAQmx.

Run::

    cd snn-nidaq
    python demo_torch_dummy.py

Switch to the real backend by setting ``use_dummy=False`` (and making
sure ``libni.so`` is built — see ``ni_interface/compile_cpp``).
"""
from __future__ import annotations

import torch
import torch.nn as nn

from ni_torch_layer import NiDAQClampLayer


def main():
    torch.manual_seed(0)

    # --- network-time grid ---------------------------------------------------
    dt_ms = 0.1
    T = 500           # 50 ms per trial
    B = 4             # batch size

    # --- a trivial upstream module: linear projection -> input current -------
    # Maps a learned weight vector to a (T, B) current waveform. In a real
    # SNN this would be your spiking layers feeding synaptic current.
    upstream = nn.Linear(B, B, bias=False)

    # --- the dynamic-clamp layer (dummy backend) -----------------------------
    clamp = NiDAQClampLayer(
        dt_ms=dt_ms,
        proxy_spike=True,
        grad_mode="surrogate",   # let gradients flow back through I
        use_dummy=True,          # <-- no DAQ hardware needed
    )

    opt = torch.optim.Adam(upstream.parameters(), lr=1e-2)

    # A made-up target trace: ramp from -70 mV to -50 mV.
    V_target = torch.linspace(-0.070, -0.050, T).unsqueeze(1).expand(T, B)

    # Constant drive shaped (T, B); upstream learns to scale it.
    base = torch.ones(T, B) * 1e-10  # 100 pA-ish in amps

    for step in range(20):
        opt.zero_grad()
        I = upstream(base)                  # (T, B)
        V = clamp(I)                        # (T, B) — runs through dummy LIF
        loss = ((V - V_target) ** 2).mean()
        loss.backward()
        opt.step()
        if step % 5 == 0 or step == 19:
            print(
                f"step {step:3d}  loss={loss.item():.4e}  "
                f"V mean={V.mean().item():+.4f} V  "
                f"spikes={int(clamp.last_spikes.sum().item()) if clamp.last_spikes is not None else 0}"
            )

    clamp.close()
    print("done.")


if __name__ == "__main__":
    main()
