#include <stdio.h>
extern "C" {
#include <NIDAQmx.h>
}

#ifdef __cplusplus
extern "C" float64 data;
#else
extern float64 data;
#endif


extern "C" {
int nidaqrec(void);
void read_sample();

}
int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout);
double clean_up();
double step_clamp(double t, double I);
int set_thread_priority_max();
