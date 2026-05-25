"""Diagnostic plot: show the bottleneck current `I(t)` and the
neuron's V response (with spike markers) for a few sample images.

Useful for sanity-checking the rig wiring or the dummy backend
before/after training. Example::

    cd snn-nidaq
    python debug_plot.py                       # untrained dummy
    python debug_plot.py --init-from runs/<ts>_classifier.pt
    python debug_plot.py --backend real        # uses libni.so

For each requested class, picks one sample from the test set, runs it
through the encoder + DAQ, and plots:

  row 1 : input image
  row 2 : bottleneck current I(t)  (pA)
  row 3 : DAQ-returned V(t) (volts * scale_in) + spike markers if any
"""
from __future__ import annotations

import argparse
import os

import torch

from ni_torch_layer import NiDAQClampLayer
from classifier_model import SNN_DAQ_Classifier, make_mnist_subset


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--backend", choices=("dummy", "real"), default="dummy")
    p.add_argument("--digits", type=int, nargs="+", default=[0, 1, 7])
    p.add_argument("--num-steps", type=int, default=5000,
                   help="Time steps per sample (dt=0.1 ms). 5000 = 500 ms.")
    p.add_argument("--num-hidden", type=int, default=64)
    p.add_argument("--init-scale-pA", type=float, default=300.0)
    p.add_argument("--init-from", type=str, default=None,
                   help="Optional checkpoint to load before plotting.")
    p.add_argument("--vthresh-mV", type=float, default=None,
                   help="Proxy-spike threshold. Defaults: -55 (dummy), -30 (real).")
    p.add_argument("--vreset-mV", type=float, default=-70.0)
    p.add_argument("--out", type=str, default=None,
                   help="If set, save the figure here instead of showing.")
    return p.parse_args()


def main():
    args = parse_args()
    use_dummy = args.backend == "dummy"
    if args.vthresh_mV is None:
        args.vthresh_mV = -55.0 if use_dummy else -30.0

    import matplotlib
    if args.out is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    torch.manual_seed(0)

    here = os.path.dirname(os.path.abspath(__file__))
    _, test_loader, digits = make_mnist_subset(
        digits=tuple(args.digits),
        n_per_class_train=1,    # unused; we only sample from test
        n_per_class_test=8,
        batch_size=1,
        root=os.path.join(here, os.pardir, "data", "mnist"),
    )

    # Grab one sample per class from the test set.
    seen = {}
    for x, y in test_loader:
        c = int(y.item())
        if c not in seen:
            seen[c] = x[0]                  # (784,)
        if len(seen) == len(digits):
            break
    sample_xs = torch.stack([seen[c] for c in sorted(seen)], dim=0)   # (C, 784)
    sample_ys = sorted(seen)

    dt_ms = 0.1
    clamp = NiDAQClampLayer(
        dt_ms=dt_ms,
        proxy_spike=True,
        vthresh_mV=args.vthresh_mV,
        vreset_mV=args.vreset_mV,
        grad_mode="surrogate",
        use_dummy=use_dummy,
        runtime=max(2.0, args.num_steps * dt_ms / 1000.0 * len(sample_xs) * 1.5),
    )
    print(f"[daq] backend={clamp.backend_name}  vthresh={args.vthresh_mV} mV")

    try:
        model = SNN_DAQ_Classifier(
            daq=clamp,
            num_classes=len(digits),
            num_steps=args.num_steps,
            num_hidden=args.num_hidden,
            init_scale_pA=args.init_scale_pA,
        )
        if args.init_from:
            sd = torch.load(args.init_from, map_location="cpu")
            missing, unexpected = model.load_state_dict(sd, strict=False)
            print(f"[init] {args.init_from}  "
                  f"missing={len(missing)}  unexpected={len(unexpected)}")
        model.eval()

        with torch.no_grad():
            # Hand-roll the forward so we can grab I separately from V.
            I = model._run_encoder(sample_xs)             # (T, C) pA
            clamp.reset_time(0.0)
            V = clamp(I)                                  # (T, C)
            spikes = clamp.last_spikes                    # (T, C) or None
            feats = model._summarise(V, spikes)
            feats_n = model.feat_norm(feats)
            logits = model.head(feats_n)
            pred = logits.argmax(dim=1)

    finally:
        clamp.close()

    T, C = I.shape
    t_ms = torch.arange(T) * dt_ms

    fig, axes = plt.subplots(3, C, figsize=(3.4 * C, 7.0), squeeze=False)
    for c in range(C):
        # row 0: image
        ax = axes[0, c]
        ax.imshow(sample_xs[c].view(28, 28).numpy(), cmap="gray")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(f"digit={digits[sample_ys[c]]}  pred={digits[int(pred[c])]}")

        # row 1: input current
        ax = axes[1, c]
        ax.plot(t_ms.numpy(), I[:, c].numpy(), lw=0.6)
        ax.set_ylabel("I (pA)" if c == 0 else "")
        ax.axhline(0, color="k", lw=0.4, alpha=0.4)
        ax.set_xlim(0, t_ms[-1].item())

        # row 2: V response
        ax = axes[2, c]
        ax.plot(t_ms.numpy(), V[:, c].numpy(), lw=0.6, color="C1")
        ax.set_xlabel("t (ms)")
        ax.set_ylabel("V (volts * sf_in)" if c == 0 else "")
        ax.set_xlim(0, t_ms[-1].item())
        if spikes is not None:
            # snntorch-style: just mark rising edges of the spike flag.
            spk = spikes[:, c].numpy()
            edges = ((spk[1:] - spk[:-1]) > 0).nonzero()[0] + 1
            n_spk = int(spk.sum().item()) if hasattr(spk, "sum") else 0
            ax.eventplot(t_ms[edges].numpy(),
                         lineoffsets=V[:, c].max().item(),
                         linelengths=(V[:, c].max() - V[:, c].min()).item() * 0.15,
                         colors="r")
            ax.text(0.98, 0.05,
                    f"spikes={int(spk.sum())} (rises={len(edges)})",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8, color="r")

    fig.suptitle(f"Bottleneck I and DAQ V response  "
                 f"[{clamp.backend_name}, num_steps={T}, dt={dt_ms} ms]")
    fig.tight_layout()

    if args.out:
        fig.savefig(args.out, dpi=120)
        print(f"[save] {args.out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
