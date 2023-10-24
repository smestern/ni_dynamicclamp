
#include <NIDAQmx.h>
#include <stdio.h>

float64 data[1];
TaskHandle taskHandle;
TaskHandle taskHandleWrite;
int nidaqrec();
void read_sample();
int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout);
double clean_up();
double step_clamp(double t, double I);
int set_thread_priority_max();
int run_step_loop(double *I, double *out, int size);
