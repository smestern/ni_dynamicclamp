"""Sweep driver for the dynamic-clamp benchmark.

Fans out one ``bench_runner.py`` subprocess per benchmark point. Two
sweep modes:

* ``--grid size`` — vary network size (``--sizes``) at fixed ``--dt-ms``.
* ``--grid rate`` — vary clock dt (``--dts-ms``) at fixed ``--n``.
* ``--grid both`` — run both sweeps back-to-back into separate CSVs.

Each subprocess gets a fresh build dir under ``--scratch`` and appends
its row to the shared ``--out-csv``. Failed points are logged but do
not abort the sweep.

Usage::

    python benchmarks/sweep.py --grid size --runtime-s 1.5
    python benchmarks/sweep.py --grid rate --quick --runtime-s 2
"""
from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
RUNNER = THIS_DIR / "bench_runner.py"
RESULTS_DIR = THIS_DIR / "results"


DEFAULT_SIZES = [100, 250, 500, 1000, 2000, 5000]
DEFAULT_DTS_MS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0]
QUICK_SIZES = [100, 500, 1000]
QUICK_DTS_MS = [0.1, 0.2, 0.5]


def _csv_int_list(s: str) -> list[int]:
    return [int(x.strip()) for x in s.split(",") if x.strip()]


def _csv_float_list(s: str) -> list[float]:
    return [float(x.strip()) for x in s.split(",") if x.strip()]


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grid", choices=["size", "rate", "both"], default="size")
    p.add_argument("--runtime-s", type=float, default=1.5)
    p.add_argument("--p-conn", type=float, default=0.02)
    p.add_argument("--sizes", type=_csv_int_list, default=None,
                   help="Comma-separated neuron counts (size grid).")
    p.add_argument("--dts-ms", type=_csv_float_list, default=None,
                   help="Comma-separated dt values in ms (rate grid).")
    p.add_argument("--fixed-dt-ms", type=float, default=0.1,
                   help="dt used for the size grid.")
    p.add_argument("--fixed-n", type=int, default=1000,
                   help="Network size used for the rate grid.")
    p.add_argument("--quick", action="store_true",
                   help="Use the small QUICK_* defaults for fast iteration.")
    p.add_argument("--out-csv", type=Path, default=None,
                   help="Output CSV (default: results/<grid>_<timestamp>.csv).")
    p.add_argument("--scratch", type=Path,
                   default=Path(tempfile.gettempdir()) / "ni_bench",
                   help="Parent dir for per-point Brian2 build dirs.")
    p.add_argument("--keep-build-dirs", action="store_true",
                   help="Don't delete build dirs after each point completes.")
    p.add_argument("--ai-chan", default="Dev1/ai0")
    p.add_argument("--ao-chan", default="Dev1/ao0")
    return p.parse_args(argv)


def _run_point(label: str, n: int, dt_ms: float, args, out_csv: Path,
               point_idx: int) -> bool:
    build_dir = args.scratch / f"{label}_{point_idx:03d}_n{n}_dt{dt_ms}"
    save_rt = RESULTS_DIR / f"{out_csv.stem}_n{n}_dt{dt_ms}_read_times.txt"

    cmd = [
        sys.executable, str(RUNNER),
        "--n", str(n),
        "--dt-ms", str(dt_ms),
        "--runtime-s", str(args.runtime_s),
        "--p-conn", str(args.p_conn),
        "--out-csv", str(out_csv),
        "--build-dir", str(build_dir),
        "--save-read-times", str(save_rt),
        "--ai-chan", args.ai_chan,
        "--ao-chan", args.ao_chan,
    ]
    print(f"\n[sweep] [{label} {point_idx+1}] N={n} dt={dt_ms}ms -> {build_dir}",
          flush=True)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd)
    elapsed = time.perf_counter() - t0

    if not args.keep_build_dirs and build_dir.exists():
        shutil.rmtree(build_dir, ignore_errors=True)

    if proc.returncode != 0:
        print(f"[sweep] FAILED point N={n} dt={dt_ms}ms (exit {proc.returncode}, "
              f"{elapsed:.1f}s)", file=sys.stderr)
        return False
    print(f"[sweep] OK N={n} dt={dt_ms}ms in {elapsed:.1f}s", flush=True)
    return True


def _print_summary(out_csv: Path) -> None:
    if not out_csv.exists():
        print(f"[sweep] no CSV produced at {out_csv}", file=sys.stderr)
        return
    with out_csv.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        print(f"[sweep] CSV empty: {out_csv}")
        return
    print(f"\n[sweep] Summary ({out_csv}):")
    cols = ["n", "dt_ms", "wall_s", "rt_factor", "steps_taken",
            "mean_step_s", "p99_step_s", "total_debt_s"]
    print("  " + "  ".join(f"{c:>13}" for c in cols))
    for r in rows:
        print("  " + "  ".join(f"{r.get(c, ''):>13}" for c in cols))


def _resolve_out_csv(args, label: str) -> Path:
    if args.out_csv is not None and args.grid != "both":
        return args.out_csv
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    return RESULTS_DIR / f"{label}_{stamp}.csv"


def run_size_grid(args) -> Path:
    sizes = args.sizes or (QUICK_SIZES if args.quick else DEFAULT_SIZES)
    out_csv = _resolve_out_csv(args, "size")
    if out_csv.exists():
        out_csv.unlink()  # fresh CSV per sweep
    args.scratch.mkdir(parents=True, exist_ok=True)
    print(f"[sweep] size grid: sizes={sizes} dt={args.fixed_dt_ms}ms "
          f"runtime={args.runtime_s}s -> {out_csv}")
    for i, n in enumerate(sizes):
        _run_point("size", n, args.fixed_dt_ms, args, out_csv, i)
    _print_summary(out_csv)
    return out_csv


def run_rate_grid(args) -> Path:
    dts = args.dts_ms or (QUICK_DTS_MS if args.quick else DEFAULT_DTS_MS)
    out_csv = _resolve_out_csv(args, "rate")
    if out_csv.exists():
        out_csv.unlink()
    args.scratch.mkdir(parents=True, exist_ok=True)
    print(f"[sweep] rate grid: dts_ms={dts} N={args.fixed_n} "
          f"runtime={args.runtime_s}s -> {out_csv}")
    for i, dt in enumerate(dts):
        _run_point("rate", args.fixed_n, dt, args, out_csv, i)
    _print_summary(out_csv)
    return out_csv


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.grid in ("size", "both"):
        run_size_grid(args)
    if args.grid in ("rate", "both"):
        run_rate_grid(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
