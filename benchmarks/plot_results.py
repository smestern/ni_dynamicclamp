"""Plot benchmark sweep results.

Reads a CSV produced by ``sweep.py`` and renders two stacked subplots:

* mean step time (with std error band and p99 marker) vs. axis
* real-time factor + total_debt (twin axis) vs. axis

Optional ``--read-times-dir`` overlays a separate figure with one
histogram per benchmark point from saved ``*_read_times.txt`` files
(see ``bench_runner.py --save-read-times``).

Usage::

    python benchmarks/plot_results.py --csv benchmarks/results/size_*.csv \\
        --x n --out benchmarks/results/size.png
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np


def _load_rows(csv_path: Path) -> list[dict]:
    with csv_path.open() as fh:
        return list(csv.DictReader(fh))


def _f(row: dict, key: str) -> float:
    val = row.get(key, "")
    try:
        return float(val)
    except (TypeError, ValueError):
        return float("nan")


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", type=Path, required=True)
    p.add_argument("--x", choices=["n", "dt_ms"], required=True,
                   help="Which CSV column drives the x-axis.")
    p.add_argument("--out", type=Path, required=True,
                   help="Output PNG path.")
    p.add_argument("--read-times-dir", type=Path, default=None,
                   help="If set, also render a histogram figure to "
                        "<out>_hist.png from saved read_times files.")
    p.add_argument("--logx", action="store_true")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = _load_rows(args.csv)
    if not rows:
        print(f"ERROR: no rows in {args.csv}", file=sys.stderr)
        return 1

    rows.sort(key=lambda r: _f(r, args.x))
    xs = np.array([_f(r, args.x) for r in rows])
    mean = np.array([_f(r, "mean_step_s") for r in rows]) * 1e6  # us
    std = np.array([_f(r, "std_step_s") for r in rows]) * 1e6
    p99 = np.array([_f(r, "p99_step_s") for r in rows]) * 1e6
    rt = np.array([_f(r, "rt_factor") for r in rows])
    debt = np.array([_f(r, "total_debt_s") for r in rows]) * 1e3  # ms
    dt_ms = np.array([_f(r, "dt_ms") for r in rows])

    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    ax_top.fill_between(xs, mean - std, mean + std, alpha=0.25, label="±1σ")
    ax_top.plot(xs, mean, "o-", label="mean step (µs)")
    ax_top.plot(xs, p99, "s--", label="p99 step (µs)", color="tab:red")
    # Show the per-point dt budget as a horizontal reference if it
    # varies across the sweep.
    if args.x == "dt_ms":
        ax_top.plot(xs, xs * 1e3, ":", color="gray", label="dt budget")
    elif np.allclose(dt_ms, dt_ms[0]):
        ax_top.axhline(dt_ms[0] * 1e3, color="gray", linestyle=":",
                       label=f"dt budget ({dt_ms[0]} ms)")
    ax_top.set_ylabel("Per-step time (µs)")
    ax_top.set_title(f"Dynamic-clamp benchmark: {args.csv.name}")
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(loc="best", fontsize=8)

    ax_bot.plot(xs, rt, "o-", color="tab:blue", label="real-time factor")
    ax_bot.axhline(1.0, color="gray", linestyle=":", linewidth=1)
    ax_bot.set_ylabel("Real-time factor (wall / sim)", color="tab:blue")
    ax_bot.tick_params(axis="y", labelcolor="tab:blue")
    ax_bot.set_xlabel({"n": "Network size N", "dt_ms": "dt (ms)"}[args.x])
    ax_bot.grid(True, alpha=0.3)

    ax_debt = ax_bot.twinx()
    ax_debt.plot(xs, debt, "s--", color="tab:red", label="total_debt (ms)")
    ax_debt.set_ylabel("Total debt (ms)", color="tab:red")
    ax_debt.tick_params(axis="y", labelcolor="tab:red")

    if args.logx:
        ax_top.set_xscale("log")
        ax_bot.set_xscale("log")

    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"[plot] wrote {args.out}")

    if args.read_times_dir is not None:
        _plot_histograms(rows, args, plt)

    return 0


def _plot_histograms(rows, args, plt) -> None:
    import matplotlib.pyplot as _plt
    files = []
    for r in rows:
        path = r.get("read_times_path", "")
        if path:
            p = Path(path)
            if p.exists():
                files.append((r, p))
    if not files:
        return
    n = len(files)
    fig, axes = _plt.subplots(n, 1, figsize=(8, 2.0 * n + 1), sharex=True)
    if n == 1:
        axes = [axes]
    for ax, (r, p) in zip(axes, files):
        try:
            data = np.loadtxt(p) * 1e6  # us
        except Exception:
            continue
        ax.hist(data, bins=80, color="tab:blue", alpha=0.8)
        ax.set_ylabel("count")
        label = f"N={r.get('n')} dt={r.get('dt_ms')}ms (mean={data.mean():.1f}µs)"
        ax.set_title(label, fontsize=9)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("Sampled per-step interval (µs)")
    fig.tight_layout()
    out = args.out.with_name(args.out.stem + "_hist.png")
    fig.savefig(out, dpi=120)
    print(f"[plot] wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
