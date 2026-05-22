"""Single-point benchmark runner for the dynamic-clamp pipeline.

Builds a parameterized LIF E/I network (Davis-2021-style), wires neuron 0
to the in-vitro cell via ``attach_neuron(... method='proxy')``, runs it
through Brian2 ``cpp_standalone`` against the real NI-DAQ, and writes a
single CSV row of timing metrics.

Designed to be invoked as a subprocess by ``benchmarks/sweep.py`` so each
benchmark point gets a clean Brian2 device and a fresh build dir.

Usage::

    python -m benchmarks.bench_runner \
        --n 1000 --dt-ms 0.1 --runtime-s 1.5 \
        --out-csv benchmarks/results/run.csv \
        --build-dir /tmp/bench_build_n1000

Hardware-only: exits non-zero with a clear message if Dev1 is not
reachable. Always builds with ``-DDEBUG`` so the C++ backend dumps
``read_times.txt`` for jitter analysis.
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import os
import shutil
import sys
import time
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent
NI_DIR = REPO_ROOT / "ni_interface"
TESTS_DIR = NI_DIR / "tests"

# Make ni_brian2 + the test-helper Brian2 wrappers importable.
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(NI_DIR))
sys.path.insert(0, str(TESTS_DIR))
from ni_brian2 import * # noqa: F401

CSV_COLUMNS = [
    "n",
    "dt_ms",
    "runtime_s",
    "p_conn",
    "wall_s",
    "sim_wall_s",
    "wall_outer_s",
    "rt_factor",
    "steps_taken",
    "total_debt_s",
    "mean_step_s",
    "std_step_s",
    "p50_step_s",
    "p99_step_s",
    "max_step_s",
    "read_times_path",
]


def probe_ni_hardware() -> bool:
    """Return True iff ``DAQmxResetDevice("Dev1")`` succeeds."""
    libpath = NI_DIR / "libnidaqmx.so"
    try:
        lib = ctypes.CDLL(str(libpath))
    except OSError:
        return False
    try:
        lib.DAQmxResetDevice.argtypes = [ctypes.c_char_p]
        lib.DAQmxResetDevice.restype = ctypes.c_int
        return lib.DAQmxResetDevice(b"Dev1") == 0
    except Exception:
        return False


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, required=True, help="Total neuron count.")
    p.add_argument("--dt-ms", type=float, required=True,
                   help="defaultclock.dt in milliseconds.")
    p.add_argument("--runtime-s", type=float, default=1.5,
                   help="Simulated runtime in seconds.")
    p.add_argument("--p-conn", type=float, default=0.02,
                   help="Random connection probability (E->all and I->all).")
    p.add_argument("--out-csv", type=Path, required=True,
                   help="Append-mode CSV path (header written if missing).")
    p.add_argument("--build-dir", type=Path, required=True,
                   help="Brian2 cpp_standalone output dir (will be wiped).")
    p.add_argument("--ai-chan", default="Dev1/ai0")
    p.add_argument("--ao-chan", default="Dev1/ao0")
    p.add_argument("--save-read-times", type=Path, default=None,
                   help="If set, copy read_times.txt to this path post-run.")
    p.add_argument("--seed", type=int, default=43)
    p.add_argument("--skip-hardware-probe", action="store_true",
                   help="(Debug only) skip the Dev1 probe.")
    return p.parse_args(argv)


def build_network(args, defaultclock, NeuronGroup, Synapses, Equations,
                  np_, units, attach_neuron):
    """Construct the parameterized E/I LIF network. Returns the population P."""
    ms = units["ms"]
    mV = units["mV"]
    pA = units["pA"]
    nS = units["nS"]
    pF = units["pF"]
    nF = units["nF"]
    second = units["second"]

    N = args.n
    if N < 2:
        raise SystemExit(f"--n must be >= 2 (got {N})")

    # Davis-2021-style parameters, simplified (no AdEx adaptation,
    # plain conductance-based exponential synapses).
    C = 200 * pF
    gL = 10 * nS
    EL = -65 * mV
    Ee = 0 * mV
    Ei = -80 * mV
    vth_default = -50 * mV
    vreset = -70 * mV
    taue = 5 * ms
    taui = 10 * ms
    we = 5 * nS
    wi = 25 * nS
    refrac = 5 * ms

    eqs = Equations(
        """
        dv/dt = (gL*(EL - v) + ge*(Ee - v) + gi*(Ei - v)) / C * neur_dyn : volt (unless refractory)
        d_I = clip(ge*(Ee - v) + gi*(Ei - v), -550*pA, 550*pA) : amp
        dge/dt = -ge / taue : siemens
        dgi/dt = -gi / taui : siemens
        vth : volt
        neur_dyn : 1
        C : farad
        gL : siemens
        EL : volt
        Ee : volt
        Ei : volt
        taue : second
        taui : second
        we : siemens
        wi : siemens
        vreset : volt
        refrac : second
        """
    )

    P = NeuronGroup(N, eqs, threshold="v > vth", reset="v = vreset",
                    refractory=refrac, method="euler")
    P.v = EL
    P.vth = vth_default
    P.neur_dyn = 1
    P.C = C
    P.gL = gL
    P.EL = EL
    P.Ee = Ee
    P.Ei = Ei
    P.taue = taue
    P.taui = taui
    P.we = we
    P.wi = wi
    P.vreset = vreset
    P.refrac = refrac

    Ne = max(1, int(N * 0.8))
    Pe = P[:Ne]
    Pi = P[Ne:]

    Ce = Synapses(Pe, P, on_pre="ge += we", delay=0.3 * ms)
    Ce.connect(p=args.p_conn)
    Ci = Synapses(Pi, P, on_pre="gi += wi", delay=0.3 * ms)
    Ci.connect(p=args.p_conn)

    # Wire neuron 0 to the in-vitro cell. Disable its intrinsic dynamics
    # and threshold so step_clamp owns its membrane voltage entirely.
    neuron, _ = attach_neuron(
        P, idx=0, i_mem_var="d_I",
        dt=defaultclock.dt, when="start", method="proxy",
    )
    neuron.neur_dyn = 0
    neuron.vth = 999 * mV

    return P, Ce, Ci, neuron


def main(argv=None) -> int:
    args = parse_args(argv)

    if not args.skip_hardware_probe and not probe_ni_hardware():
        print("ERROR: NI Dev1 not reachable; cannot run hardware benchmark.",
              file=sys.stderr)
        return 2

    # Wipe build dir for a clean cpp_standalone build.
    if args.build_dir.exists():
        shutil.rmtree(args.build_dir)
    args.build_dir.mkdir(parents=True, exist_ok=True)

    read_times_path = args.build_dir / "read_times.txt"
    os.environ["NI_READ_TIMES_PATH"] = str(read_times_path)

    # Heavy imports happen here so probing/CLI errors are fast.
    from brian2 import (
        NeuronGroup, Synapses, StateMonitor, Equations, defaultclock,
        ms, mV, pA, nS, pF, nF, second, Hz, set_device, prefs, run, seed,
    )
    import numpy as np

   
    # Imported for their @implementation side effect (Brian2 needs the
    # functions registered at module-import time).
    from ni_test_helpers import (  # noqa: F401
        ni_steps_taken, ni_total_debt, ni_total_rate, ni_last_read_t, ni_run_time,
    )

    # Force -DDEBUG so interface.cpp populates read_times.txt.
    base_args = list(getattr(prefs.codegen.cpp, "extra_compile_args_gcc", []) or [])
    if "-DDEBUG" not in base_args:
        base_args.append("-DDEBUG")
    prefs.codegen.cpp.extra_compile_args_gcc = base_args

    seed(args.seed)
    np.random.seed(args.seed)

    set_device("cpp_standalone", build_on_run=False, directory=str(args.build_dir))
    from brian2 import device  # rebind after set_device

    defaultclock.dt = args.dt_ms * ms

    units = dict(ms=ms, mV=mV, pA=pA, nS=nS, pF=pF, nF=nF, second=second, Hz=Hz)
    P, Ce, Ci, neuron = build_network(
        args, defaultclock, NeuronGroup, Synapses, Equations,
        np, units, attach_neuron,
    )

    # Capture the test-helper accessors at every dt so we can read the
    # final values after the run. One-element NeuronGroup keeps the
    # StateMonitor cheap.
    probe_model = """
        steps_seen : 1
        debt_seen : second
        rate_seen : second
        run_time_seen : second
    """
    probe = NeuronGroup(1, probe_model)
    probe.run_regularly(
        """
        steps_seen = ni_steps_taken(t)
        debt_seen = ni_total_debt(t)
        rate_seen = ni_total_rate(t)
        run_time_seen = ni_run_time(t)
        """,
        dt=defaultclock.dt,
        when="end",
    )
    mon = StateMonitor(probe, ["steps_seen", "debt_seen", "rate_seen", "run_time_seen"], record=True)

    device = init_neuron_device(
        device, dt=defaultclock.dt,
        runtime=args.runtime_s,
        aI=args.ai_chan, aO=args.ao_chan,
    )

    print(f"[bench_runner] N={args.n} dt={args.dt_ms}ms runtime={args.runtime_s}s "
          f"p_conn={args.p_conn} build_dir={args.build_dir}", flush=True)

    # Build (compile but don't run) so we can time only the standalone
    # binary execution as wall clock.
    run(args.runtime_s * second)
    device.build(directory=str(args.build_dir), compile=True, run=False)

    main_bin = args.build_dir / "main"
    if not main_bin.exists():
        print(f"ERROR: standalone binary not found at {main_bin}", file=sys.stderr)
        return 3

    # Use ``device.run(...)`` (rather than a raw subprocess) so Brian2
    # marks the run as completed and reloads StateMonitor results from
    # the standalone results/ dir. Otherwise accessing ``mon.steps_seen``
    # raises "Cannot retrieve the values of state variables in standalone
    # code before the simulation has been run."
    t0 = time.perf_counter()
    try:
        device.run(
            directory=str(args.build_dir),
            with_output=False,
            run_args=[],
        )
    except Exception as e:
        print(f"ERROR: standalone run failed: {e}", file=sys.stderr)
        return 3
    wall_outer = time.perf_counter() - t0

    # Brian2 publishes two finer-grained timings written by the standalone
    # binary itself; prefer them over our Python-side perf_counter so
    # build/compile time and Python subprocess overhead are excluded.
    #   * ``timers['run_binary']`` — wall time of the binary subprocess
    #     (includes init_ni + clean_up + Network::run loop).
    #   * ``_last_run_time``      — time spent inside Network::run loop.
    binary_wall = float(getattr(device, "timers", {}).get("run_binary", wall_outer))
    sim_wall = float(getattr(device, "_last_run_time", binary_wall) or binary_wall)

    try:
        steps_arr = np.asarray(mon.steps_seen[0])
        debt_arr = np.asarray(mon.debt_seen[0])
        rate_arr = np.asarray(mon.rate_seen[0])
        run_time_arr = np.asarray(mon.run_time_seen[0])
    except Exception as e:
        print(f"ERROR: failed to read StateMonitor results: {e}", file=sys.stderr)
        return 4

    if steps_arr.size == 0:
        print("ERROR: StateMonitor produced no samples", file=sys.stderr)
        return 5

    final_steps = float(steps_arr[-1])
    final_debt_s = float(debt_arr[-1])
    final_rate_s = float(rate_arr[-1])
    final_run_time_s = float(run_time_arr[-1])
    # Authoritative wall: measured inside step_clamp via CLOCK_MONOTONIC_RAW
    # (LAST_READ_T - full_run_time). Excludes init_ni + clean_up + subprocess
    # overhead. Falls back to Brian2's binary timer if for some reason the
    # accessor returned 0 (e.g. on a zero-step run).
    wall = final_run_time_s if final_run_time_s > 0 else binary_wall
    mean_step_s = (final_rate_s / final_steps) if final_steps > 0 else float("nan")

    # Jitter stats from the DEBUG-mode 1-in-1000 sample buffer.
    rt_samples = np.array([], dtype=float)
    if read_times_path.exists():
        try:
            rt_samples = np.loadtxt(read_times_path)
            if rt_samples.ndim == 0:
                rt_samples = rt_samples.reshape(1)
        except Exception as e:
            print(f"WARN: failed to parse {read_times_path}: {e}", file=sys.stderr)

    if rt_samples.size > 0:
        std_step_s = float(np.std(rt_samples))
        p50_step_s = float(np.percentile(rt_samples, 50))
        p99_step_s = float(np.percentile(rt_samples, 99))
        max_step_s = float(np.max(rt_samples))
    else:
        std_step_s = p50_step_s = p99_step_s = max_step_s = float("nan")

    saved_rt_path = ""
    if args.save_read_times is not None and read_times_path.exists():
        args.save_read_times.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(read_times_path, args.save_read_times)
        saved_rt_path = str(args.save_read_times)

    rt_factor = wall / args.runtime_s if args.runtime_s > 0 else float("nan")

    row = {
        "n": args.n,
        "dt_ms": args.dt_ms,
        "runtime_s": args.runtime_s,
        "p_conn": args.p_conn,
        "wall_s": f"{wall:.6f}",
        "sim_wall_s": f"{sim_wall:.6f}",
        "wall_outer_s": f"{wall_outer:.6f}",
        "rt_factor": f"{rt_factor:.6f}",
        "steps_taken": int(final_steps),
        "total_debt_s": f"{final_debt_s:.9f}",
        "mean_step_s": f"{mean_step_s:.9f}",
        "std_step_s": f"{std_step_s:.9f}",
        "p50_step_s": f"{p50_step_s:.9f}",
        "p99_step_s": f"{p99_step_s:.9f}",
        "max_step_s": f"{max_step_s:.9f}",
        "read_times_path": saved_rt_path,
    }

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    write_header = not args.out_csv.exists() or args.out_csv.stat().st_size == 0
    with args.out_csv.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=CSV_COLUMNS)
        if write_header:
            w.writeheader()
        w.writerow(row)

    print(f"[bench_runner] DONE wall={wall:.3f}s rt_factor={rt_factor:.3f} "
          f"steps={int(final_steps)} debt={final_debt_s*1e3:.2f}ms "
          f"mean_step={mean_step_s*1e6:.2f}us p99={p99_step_s*1e6:.2f}us",
          flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
