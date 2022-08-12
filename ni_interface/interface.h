#include <stdio.h>
#include <NIDAQmx.h>

int nidaqrec(void);
int init_ni(float64 net_clock_dt);
double clean_up();
double step_clamp(double t, double I);