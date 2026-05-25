"""Hardware-free demo: train the DAQ-in-the-loop MNIST classifier on
the dummy backend.

Trains :class:`SNN_DAQ_Classifier` end-to-end (encoder + bottleneck +
post-DAQ head) on a 3-digit MNIST subset, using
``NiDAQClampLayer(use_dummy=True)`` so no NI card is required. The
DAQ contributes a passive-cell V trace whose summary stats feed the
classifier head; the surrogate gradient pass-through in
``_NiDAQClampFn.backward`` lets the encoder learn.

Run::

    cd snn-nidaq
    python demo_classifier_dummy.py

This is the offline counterpart to ``demo_classifier_real.py``; both
share the model definition in :mod:`classifier_model`.
"""
from __future__ import annotations

import time

import torch
import torch.nn as nn

from ni_torch_layer import NiDAQClampLayer
from classifier_model import SNN_DAQ_Classifier, make_mnist_subset


def evaluate(model: SNN_DAQ_Classifier, loader) -> float:
    model.eval()
    correct = total = 0
    with torch.no_grad():
        for x, y in loader:
            logits = model(x)
            pred = logits.argmax(dim=1)
            correct += (pred == y).sum().item()
            total += y.numel()
    model.train()
    return correct / max(total, 1)


def main():
    torch.manual_seed(0)

    # ---- config ----
    digits = (0, 1, 7)
    epochs = 15
    batch_size = 16
    num_steps = 5000          # 500 ms of simulated time at dt=0.1 ms
    lr = 3e-3
    init_scale_pA = 300.0
    # Proxy-spike thresholds chosen so the passive dummy LIF can
    # actually cross threshold given the above current scale; on the
    # real rig you'd use the standard -30 / -70 mV pair (the defaults).
    vthresh_mV = -55.0
    vreset_mV = -70.0

    # ---- data ----
    train_loader, test_loader, digits = make_mnist_subset(
        digits=digits,
        n_per_class_train=200,
        n_per_class_test=50,
        batch_size=batch_size,
    )
    print(f"[data] train={len(train_loader.dataset)}  test={len(test_loader.dataset)}  "
          f"classes={digits}")

    # ---- DAQ ----
    clamp = NiDAQClampLayer(
        dt_ms=0.1,
        proxy_spike=True,
        vthresh_mV=vthresh_mV,
        vreset_mV=vreset_mV,
        grad_mode="surrogate",
        use_dummy=True,
    )

    try:
        # ---- model ----
        model = SNN_DAQ_Classifier(
            daq=clamp,
            num_classes=len(digits),
            num_steps=num_steps,
            num_hidden=64,
            init_scale_pA=init_scale_pA,
        )
        opt = torch.optim.Adam(model.parameters(), lr=lr)
        loss_fn = nn.CrossEntropyLoss()

        # ---- train ----
        for epoch in range(epochs):
            t0 = time.time()
            running_loss = 0.0
            n_batches = 0
            for x, y in train_loader:
                opt.zero_grad()
                logits = model(x)
                loss = loss_fn(logits, y)
                loss.backward()
                opt.step()
                running_loss += loss.item()
                n_batches += 1
            train_acc = evaluate(model, train_loader)
            test_acc = evaluate(model, test_loader)
            print(
                f"epoch {epoch+1:2d}/{epochs}  "
                f"loss={running_loss / max(n_batches, 1):.4f}  "
                f"train_acc={train_acc:.3f}  test_acc={test_acc:.3f}  "
                f"({time.time() - t0:.1f}s)"
            )
    finally:
        clamp.close()

    print("done.")


if __name__ == "__main__":
    main()
