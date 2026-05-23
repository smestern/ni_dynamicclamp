// =========================================================================
// interface_torch_export.cpp
//
// C-linkage forwarding shims for the PyTorch / ctypes hot-load path
// (`ni_torch.py` → `libni.so`). This TU is compiled ONLY into
// `libni.so`; it is NOT consumed by Brian2's `cpp_standalone` codegen,
// which only sees `interface.cpp`.
//
// Why a separate TU? Wrapping `interface.cpp`'s public API in
// `extern "C"` was observed to perturb Brian2's real-time runs. Keeping
// the C++ linkage in `interface.cpp` and routing PyTorch through `ni_*`
// shims here gives both consumers what they need without cross-talk.
//
// All exported symbols are prefixed `ni_` to avoid ODR collision with
// the C++-mangled originals that also live in `libni.so`.
// =========================================================================
#include "interface.h"

// `last_spike` is a file-scope `int` defined in interface.cpp with
// default (external) C++ linkage; forward-declare it here.
extern int last_spike;

extern "C"
{

int ni_init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout,
               float64 runtime, char *aI, char *aO)
{
    return init_ni(net_clock_dt, scalein, scaleout, runtime, aI, aO);
}

double ni_step_clamp(double t, double I)
{
    return step_clamp(t, I);
}

double ni_clean_up(void)
{
    return clean_up();
}

void ni_turn_on_proxy_spike(long double vthresh, long double vreset)
{
    turn_on_proxy_spike(vthresh, vreset);
}

// Batch entry point. Delegates each step to step_clamp(), so all
// real-time pacing, overrun-debt accounting, and proxy-spike logic
// are reused unchanged.
//
// * step_clamp's pacing is driven by LAST_READ_T (set inside
//   step_clamp); the `t` argument is informational only.
// * No allocations, no printf on the hot path. spk_out is optional.
// * Per-step real-time pacing is preserved inside this loop, so the
//   Python caller pays one FFI hop per batch instead of one per sample.
// * step_clamp's read_sample/write_sample wrappers log DAQmx errors in
//   DEBUG builds but do not propagate them through a return code;
//   surfacing per-step DAQmx errors here would require a hot-path
//   branch that we intentionally avoid. Callers needing strict error
//   reporting should use the single-step step_clamp() entry point.
int ni_run_step_loop(double *I_in, double *V_out, double *spk_out,
                     int n_steps, double dt_ms, double t0)
{
    const double dt_s = dt_ms / 1000.0;
    for (int i = 0; i < n_steps; i++)
    {
        double t = t0 + (double)i * dt_s;
        V_out[i] = step_clamp(t, I_in[i]);
        if (spk_out != NULL)
        {
            spk_out[i] = (double)last_spike;
        }
    }
    return 0;
}

} // extern "C"
