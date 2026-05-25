"""On-rig demo: train the DAQ-in-the-loop MNIST classifier end-to-end
against a **real cell** via NI-DAQ.

Run from the rig host (NI card present, ``libni.so`` built, real-time
kernel, etc.)::

    cd snn-nidaq
    python demo_classifier_real.py --epochs 1 --batch-size 1 --num-steps 5000

This is the rig counterpart of ``demo_classifier_dummy.py``: identical
model (see :mod:`classifier_model`), just with ``use_dummy=False`` on
the :class:`NiDAQClampLayer`. Each sample is a real-time recording, so
defaults are aggressively small to keep an epoch under a minute of
wall-clock.

Wall-clock budget (rough): ``epochs * train_size * num_steps * dt_ms``
of real-time DAQ I/O, plus per-batch Python overhead. With the
defaults (epochs=1, batch=1, num_steps=5000, n_per_class=20, 3 classes)
that's ~30 s of DAQ time.

CLI::

    --dummy            Use the dummy backend (for a sanity dry run).
    --epochs N         Training epochs.
    --batch-size B     Mini-batch size (serialised through the cell).
    --num-steps T      Time steps per sample at dt=0.1 ms.
    --n-per-class N    Samples per class in the training set.
    --lr LR            Adam learning rate.
    --init-from PATH   Load a state_dict (e.g. from the dummy demo) to
                       warm-start training and shorten rig time.
    --no-train         Skip training; only evaluate on the test set.
    --out-dir DIR      Where to save weights + metrics CSV.

Saves on completion (or Ctrl-C):
    runs/<timestamp>_classifier.pt    -- model state_dict
    runs/<timestamp>_metrics.csv      -- per-epoch loss / acc
"""
from __future__ import annotations

import argparse
import csv
import datetime
import os
import sys
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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dummy", action="store_true",
                   help="Use ni_torch_dummy instead of libni.so (dry run).")
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=1)
    p.add_argument("--num-steps", type=int, default=5000,
                   help="Time steps per sample (dt=0.1 ms). 5000 = 500 ms.")
    p.add_argument("--n-per-class", type=int, default=20)
    p.add_argument("--n-per-class-test", type=int, default=10)
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--num-hidden", type=int, default=64)
    p.add_argument("--init-scale-pA", type=float, default=300.0)
    p.add_argument("--digits", type=int, nargs="+", default=[0, 1, 7])
    p.add_argument("--init-from", type=str, default=None,
                   help="Load a state_dict (e.g. snn-nidaq/runs/<ts>_classifier.pt).")
    p.add_argument("--no-train", action="store_true",
                   help="Skip training; evaluate only (useful with --init-from).")
    p.add_argument("--out-dir", type=str, default="runs")
    p.add_argument("--ai-chan", type=str, default="Dev1/ai0")
    p.add_argument("--ao-chan", type=str, default="Dev1/ao0")
    # Use the real rig's defaults when --dummy is NOT set; the dummy
    # backend's passive cell needs a lower threshold to actually spike
    # (mirrors demo_classifier_dummy.py).
    p.add_argument("--vthresh-mV", type=float, default=-30.0)
    p.add_argument("--vreset-mV", type=float, default=-70.0)
    return p.parse_args()


def save_checkpoint(out_dir: str, ts: str, model: nn.Module,
                    metrics: list[dict]) -> tuple[str, str]:
    os.makedirs(out_dir, exist_ok=True)
    pt_path = os.path.join(out_dir, f"{ts}_classifier.pt")
    csv_path = os.path.join(out_dir, f"{ts}_metrics.csv")
    torch.save(model.state_dict(), pt_path)
    if metrics:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(metrics[0].keys()))
            writer.writeheader()
            writer.writerows(metrics)
    return pt_path, csv_path


def main():
    args = parse_args()
    torch.manual_seed(0)

    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = args.out_dir if os.path.isabs(args.out_dir) \
        else os.path.join(here, args.out_dir)

    # ---- data ----
    train_loader, test_loader, digits = make_mnist_subset(
        digits=tuple(args.digits),
        n_per_class_train=args.n_per_class,
        n_per_class_test=args.n_per_class_test,
        batch_size=args.batch_size,
        root=os.path.join(here, os.pardir, "data", "mnist"),
    )
    print(f"[data] train={len(train_loader.dataset)}  "
          f"test={len(test_loader.dataset)}  classes={digits}")

    # Wall-clock budget estimate (DAQ time only, ignoring Python overhead).
    daq_secs = (args.epochs
                * len(train_loader.dataset)
                * args.num_steps * 0.0001)
    runtime_budget = max(daq_secs * 1.5, 2.0)
    print(f"[budget] estimated DAQ time: {daq_secs:.1f} s "
          f"(runtime arg: {runtime_budget:.1f} s)")

    # ---- DAQ ----
    clamp = NiDAQClampLayer(
        dt_ms=0.1,
        proxy_spike=True,
        vthresh_mV=args.vthresh_mV,
        vreset_mV=args.vreset_mV,
        grad_mode="surrogate",
        use_dummy=args.dummy,
        runtime=runtime_budget,
        ai_chan=args.ai_chan,
        ao_chan=args.ao_chan,
    )
    print(f"[daq] backend={clamp.backend_name}")

    metrics: list[dict] = []
    try:
        # ---- model ----
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
            print(f"[init] loaded {args.init_from}  "
                  f"missing={len(missing)} unexpected={len(unexpected)}")

        opt = torch.optim.Adam(model.parameters(), lr=args.lr)
        loss_fn = nn.CrossEntropyLoss()

        if args.no_train:
            test_acc = evaluate(model, test_loader)
            print(f"[eval-only] test_acc={test_acc:.3f}")
            metrics.append(dict(epoch=0, loss=float("nan"),
                                train_acc=float("nan"), test_acc=test_acc,
                                wall_s=0.0))
        else:
            for epoch in range(args.epochs):
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
                wall = time.time() - t0
                avg_loss = running_loss / max(n_batches, 1)
                print(
                    f"epoch {epoch+1:2d}/{args.epochs}  "
                    f"loss={avg_loss:.4f}  "
                    f"train_acc={train_acc:.3f}  test_acc={test_acc:.3f}  "
                    f"({wall:.1f}s)"
                )
                metrics.append(dict(epoch=epoch + 1, loss=avg_loss,
                                    train_acc=train_acc, test_acc=test_acc,
                                    wall_s=wall))
    except KeyboardInterrupt:
        print("\n[interrupt] saving partial checkpoint...", file=sys.stderr)
    finally:
        # Always release the DAQ task handles, even on Ctrl-C / exception.
        try:
            pt, csvp = save_checkpoint(out_dir, ts, model, metrics)
            print(f"[save] {pt}\n[save] {csvp}")
        except NameError:
            pass    # model never got built
        clamp.close()

    print("done.")


if __name__ == "__main__":
    main()
