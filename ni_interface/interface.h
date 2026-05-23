#include <stdio.h>
extern "C" {
#include <NIDAQmx.h>
}

#ifdef __cplusplus
extern "C" float64 data;
#else
extern float64 data;
#endif


int nidaqrec(void);
int read_sample();
int write_sample(float64 val);

// ---- Public API (C++ linkage) ----
// Default C++ linkage. Brian2's cpp_standalone codegen includes this
// header and links against `interface.cpp` directly, so it sees the
// same mangled names it always has. The PyTorch / ctypes path consumes
// these via `ni_*` C-linkage forwarding shims defined in
// interface_torch_export.cpp.
int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout, float64 runtime, char *aI, char *aO);
void turn_on_proxy_spike(long double vthresh, long double vreset);
double clean_up();
double step_clamp(double t, double I);

int set_thread_priority_max();
