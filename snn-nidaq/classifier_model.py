"""Shared model + data utilities for the DAQ-in-the-loop MNIST classifier.

Architecture (PoC):

    image (28x28)
        -> rate-encode to (T, B, 784) spikes
        -> Linear(784, num_hidden)
        -> snn.Synaptic LIF (membrane state mem1)
        -> Linear(num_hidden, 1) * scale_pA      <-- 1-channel bottleneck
        -> NiDAQClampLayer (real or dummy)        <-- the cell
        -> (mean, std, max, min of V, spike count) -> Linear(5, 32) -> ReLU
        -> Linear(32, num_classes)

The DAQ layer is the **only** non-differentiable step; we rely on
``grad_mode="surrogate"`` in :class:`NiDAQClampLayer` to pass the
upstream gradient straight through, so the encoder + bottleneck train
end-to-end against a cross-entropy loss on the classifier head.

The bottleneck is intentionally 1-channel because the rig has one cell;
batches are serialized through the same DAQ inside ``_NiDAQClampFn``.
"""
from __future__ import annotations

import os
from typing import Optional, Sequence, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import snntorch as snn
from snntorch import spikegen


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def make_mnist_subset(
    digits: Sequence[int] = (0, 1, 7),
    n_per_class_train: int = 400,
    n_per_class_test: int = 100,
    batch_size: int = 16,
    root: str = "../data/mnist",
    seed: int = 0,
) -> Tuple[torch.utils.data.DataLoader, torch.utils.data.DataLoader, Tuple[int, ...]]:
    """Return train/test DataLoaders restricted to ``digits``.

    Labels are remapped to ``0..len(digits)-1`` so the classifier head
    can be a small ``Linear(_, len(digits))``.
    """
    from torchvision import datasets, transforms

    tfm = transforms.Compose([
        transforms.ToTensor(),                # (1, 28, 28) in [0, 1]
        transforms.Normalize((0.0,), (1.0,)), # keep [0, 1] for rate encoding
    ])

    train = datasets.MNIST(root=root, train=True, download=True, transform=tfm)
    test = datasets.MNIST(root=root, train=False, download=True, transform=tfm)

    digits = tuple(digits)
    remap = {d: i for i, d in enumerate(digits)}

    def _subset(ds, n_per_class):
        rng = np.random.default_rng(seed)
        targets = ds.targets.numpy() if hasattr(ds.targets, "numpy") else np.asarray(ds.targets)
        keep_idx = []
        for d in digits:
            idx = np.where(targets == d)[0]
            rng.shuffle(idx)
            keep_idx.extend(idx[:n_per_class].tolist())
        rng.shuffle(keep_idx)

        class _Remapped(torch.utils.data.Dataset):
            def __len__(self_inner):
                return len(keep_idx)

            def __getitem__(self_inner, i):
                x, y = ds[keep_idx[i]]
                return x.view(-1), remap[int(y)]   # flatten to (784,)

        return _Remapped()

    train_loader = torch.utils.data.DataLoader(
        _subset(train, n_per_class_train),
        batch_size=batch_size, shuffle=True, drop_last=False,
    )
    test_loader = torch.utils.data.DataLoader(
        _subset(test, n_per_class_test),
        batch_size=batch_size, shuffle=False, drop_last=False,
    )
    return train_loader, test_loader, digits


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class SNN_DAQ_Classifier(nn.Module):
    """Encoder SNN -> 1-channel DAQ bottleneck -> MLP readout.

    Parameters
    ----------
    daq : NiDAQClampLayer
        Pre-built DAQ layer (real or dummy). The classifier does NOT own
        its lifecycle; the caller is responsible for ``daq.close()``.
        ``daq.grad_mode`` should be ``"surrogate"`` for end-to-end
        training; ``"detach"`` works too but only trains the head.
    num_inputs : int
        Flattened input dim (784 for MNIST).
    num_hidden : int
        Width of the SNN hidden layer.
    num_classes : int
        Output dim of the classifier head.
    num_steps : int
        Time steps fed through the encoder + DAQ per sample. At
        ``daq.dt_ms=0.1`` and ``num_steps=200`` this is 20 ms of real
        cell time per sample (rig wall-clock).
    init_scale_pA : float
        Initial value of the learnable scalar that converts the
        bottleneck projection into picoamps before it hits the DAQ.
    max_abs_pA : float
        Hard physiological clamp (pA) applied via tanh to the
        bottleneck current before it is sent to the DAQ. Prevents
        runaway currents (e.g. >>1 nA) that would kill a real cell.
    beta, alpha : float
        Initial values for ``snn.Synaptic`` membrane / synaptic decays
        (both learnable, matching the convention in
        ``nidaq_SNN.py``).
    """

    def __init__(
        self,
        daq,
        *,
        num_inputs: int = 28 * 28,
        num_hidden: int = 64,
        num_classes: int = 3,
        num_steps: int = 50,
        init_scale_pA: float = 50.0,
        max_abs_pA: float = 500.0,
        beta: float = 0.9,
        alpha: float = 0.8,
    ):
        super().__init__()
        self.daq = daq
        self.num_inputs = num_inputs
        self.num_hidden = num_hidden
        self.num_classes = num_classes
        self.num_steps = num_steps

        # Encoder ------------------------------------------------------
        self.fc1 = nn.Linear(num_inputs, num_hidden)
        self.beta_1 = nn.Parameter(torch.tensor(float(beta)))
        self.alpha_1 = nn.Parameter(torch.tensor(float(alpha)))
        self.lif1 = snn.Synaptic(
            beta=self.beta_1,
            alpha=self.alpha_1,
            init_hidden=False,
            reset_mechanism="zero",
            learn_alpha=True,
            learn_beta=True,
        )

        # Bottleneck ---------------------------------------------------
        # One scalar current per timestep per batch element. Bias on so
        # the cell can sit at a non-zero baseline drive.
        self.fc_bn = nn.Linear(num_hidden, 1, bias=True)
        # Learnable amplitude in pA. tanh-squashed below so it stays
        # in a physiologically sane range and doesn't blow the
        # amplifier headroom.
        # Physiological saturation limit (pA). Fixed, not learnable —
        # this is a safety rail for the real cell, not a free
        # parameter. tanh-squashed in _run_encoder.
        self.register_buffer(
            "max_abs_pA", torch.tensor(float(max_abs_pA))
        )
        self.scale_pA = nn.Parameter(torch.tensor(float(init_scale_pA)))

        # Head ---------------------------------------------------------
        # 5 hand-picked summary stats of the V trace / spikes per
        # batch element. Keeps the head tiny so the gradient path
        # through the DAQ is short and well-conditioned. LayerNorm
        # is critical here: V trace magnitudes depend on scale_in
        # (real ~0.07 V, dummy ~7e-3 V), and the spike-count feature
        # lives on a totally different scale; without normalisation
        # the head ignores all but the dominant feature.
        self.feat_norm = nn.LayerNorm(5)
        self.head = nn.Sequential(
            nn.Linear(5, 32),
            nn.ReLU(),
            nn.Linear(32, num_classes),
        )

    # ---------- helpers ----------

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Rate-encode a flattened image batch ``(B, 784)`` to spikes
        ``(T, B, 784)``."""
        return spikegen.rate(x, num_steps=self.num_steps)

    def _run_encoder(self, x: torch.Tensor) -> torch.Tensor:
        """Return the bottleneck current ``(T, B)`` in picoamps."""
        spk_in = self.encode(x)                        # (T, B, 784)
        T, B, _ = spk_in.shape
        syn, mem = self.lif1.init_synaptic()
        cur_b = []
        for t in range(T):
            cur1 = self.fc1(spk_in[t])
            _spk1, syn, mem = self.lif1(cur1, syn, mem)
            # Bottleneck: project membrane potential (continuous) to a
            # single scalar current. Using mem (not spk1) keeps the
            # signal differentiable in the encoder.
            scalar = self.fc_bn(mem).squeeze(-1)        # (B,)
            cur_b.append(scalar)
        cur = torch.stack(cur_b, dim=0)                 # (T, B)
        # Squash into ±max_abs_pA so the real cell never sees more
        # current than its amplifier / biology can handle. The
        # learnable scale_pA still controls the *slope* near zero
        # (gain), but the asymptote is hard-bounded.
        return self.max_abs_pA * torch.tanh(
            cur * self.scale_pA / self.max_abs_pA
        )

    def _summarise(self, V: torch.Tensor, spikes: Optional[torch.Tensor]) -> torch.Tensor:
        """V ``(T, B)`` -> features ``(B, 5)``."""
        feats = [
            V.mean(dim=0),
            V.std(dim=0),
            V.max(dim=0).values,
            V.min(dim=0).values,
        ]
        if spikes is not None:
            feats.append(spikes.sum(dim=0))
        else:
            feats.append(torch.zeros_like(feats[0]))
        return torch.stack(feats, dim=1)                # (B, 5)

    # ---------- forward ----------

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.dim() == 3:                                # (B, 28, 28)
            x = x.view(x.shape[0], -1)
        I = self._run_encoder(x)                        # (T, B) pA
        # Reset the DAQ's network-time cursor for each forward so
        # successive batches don't bake in unbounded t. Real-time
        # pacing is driven by LAST_READ_T in the C side; this is just
        # bookkeeping for debug prints.
        self.daq.reset_time(0.0)
        V = self.daq(I)                                 # (T, B) volts*sf_in
        spikes = self.daq.last_spikes
        feats = self._summarise(V, spikes)
        feats = self.feat_norm(feats)
        return self.head(feats)


__all__ = ["make_mnist_subset", "SNN_DAQ_Classifier"]
