## Dynamic Clamp scripts for National Instruments Devices

In this repo are C/C++ Demo scripts for creating (dynamic clamp)[https://www.sciencedirect.com/science/article/pii/S0166223604000554]-based bio-hybrid networks between in vitro neurons and an arbitrary neural network. Also included is a python module for easy integration with Brian2 / Pytorch.


### Hardware
Currently, this hardware has been tested with   
- BNC-2110  
- PCIE-6321  
- AMD Ryzen™ 7 5700 (Running Ubuntu LTS with real-time kernel)  
- 64GB memory  

### Achievable Speeds
Currently clamping rates up to 20 kHz have been achieved. I am hoping to push up to 50kHz and would appreciate any tips / pull requests on performance.

## Usage: Brian2

The main focus of this project is integration with the neural / ode modeling software [brian2](https://brian2.readthedocs.io/en/stable/). This project enables the "dynamic clamping" of pretty much any Brian2 implemented spiking model. 
With 2 caveats:
- The dynamic clamp current is configured for ["current-clamp" mode](https://www.moleculardevices.com/applications/patch-clamp-electrophysiology/what-current-clamp-method): E.g. the network outputs a current, the current is applied to the real neuron, the voltage is fed back into the network. Voltage clamp is possible but will require some modifications that I have yet to apply.
- The network must be able to run in [brian2's standalone mode](https://brian2.readthedocs.io/en/2.9.0/user/computation.html#standalone-code-generation). This may require some modifications of your script, but the standalone version of brian2 supports most if not all functions of the pythonic / cpython version.

The package hooks itself into brian2 standalone code generation, and cross compiles, enabling ultra-fast interface.

## Benchmarking

A small benchmark suite under `benchmarks/` measures how clamp performance scales with **network size** and **clock rate** (`defaultclock.dt`). Each benchmark point compiles a parameterized LIF E/I network through the standard Brian2 + NI-DAQ pipeline and reports mean step time, jitter (std / p99), real-time factor, and accumulated `total_debt` from `interface.cpp`. Hardware (Dev1) must be reachable.

```bash
# Sweep network size at fixed dt (0.1 ms by default).
python benchmarks/sweep.py --grid size --runtime-s 1.5

# Sweep dt at fixed N (1000 by default).
python benchmarks/sweep.py --grid rate --runtime-s 1.5

# Quick smoke-test grids.
python benchmarks/sweep.py --grid both --quick --runtime-s 1.0

# Plot a CSV produced by sweep.py.
python benchmarks/plot_results.py \
    --csv benchmarks/results/size_<timestamp>.csv \
    --x n --out benchmarks/results/size.png
```

Each sweep writes `benchmarks/results/<grid>_<timestamp>.csv` plus per-point `*_read_times.txt` files (sampled every 1000 steps in DEBUG builds) that `plot_results.py --read-times-dir ...` can render as jitter histograms. Benchmarks are *not* recordings — they always build with `-DDEBUG`; do not use the same build for actual experiments.
