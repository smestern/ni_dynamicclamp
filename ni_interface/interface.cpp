// ---- Includes ----
#include <stdio.h>
#include <time.h>
#include <sched.h>
#include <pthread.h>
#include <chrono>
#include <thread>

extern "C" {
// NI-DAQmx C API is C-style; wrap in extern "C" to avoid name mangling.
#include <NIDAQmx.h>
#include <cstdlib>
#include <cstring>
}

// ---- Macros ----
// needs -lrt (real-time lib) for clock_gettime

#define DAQmxErrChk(functionCall)                \
    if (DAQmxFailed(error = (functionCall)))     \
        goto Error;                              \
    else

#define tscmp(a, b, CMP) \
    (((a)->tv_sec == (b)->tv_sec) ? ((a)->tv_nsec CMP (b)->tv_nsec) : ((a)->tv_sec CMP (b)->tv_sec))

#define tsadd(a, b, result)                                  \
    do {                                                     \
        (result)->tv_sec  = (a)->tv_sec  + (b)->tv_sec;      \
        (result)->tv_nsec = (a)->tv_nsec + (b)->tv_nsec;     \
        if ((result)->tv_nsec >= 1000000000) {               \
            ++(result)->tv_sec;                              \
            (result)->tv_nsec -= 1000000000;                 \
        }                                                    \
    } while (0)

#define tssub(a, b, result)                                  \
    do {                                                     \
        (result)->tv_sec  = (a)->tv_sec  - (b)->tv_sec;      \
        (result)->tv_nsec = (a)->tv_nsec - (b)->tv_nsec;     \
        if ((result)->tv_nsec < 0) {                         \
            --(result)->tv_sec;                              \
            (result)->tv_nsec += 1000000000;                 \
        }                                                    \
    } while (0)

// ---- Forward declarations ----
long double ClockGetTime();

// ---- NI globals ----
// `data` is exported via interface.h, so it needs C linkage.
extern "C" {
float64 data = -0.070; // default value to prevent brian2 from freaking out
}

char errBuff[2048] = {'\0'};   // buffer for DAQmx error messages
char aIChan[256] = "Dev1/ai0"; // default analog input channel
char aOChan[256] = "Dev1/ao0"; // default analog output channel

int SAMPLE_RATE = 100000;      // in Hz
int32 error = 0;
TaskHandle taskHandle = 0;
TaskHandle taskHandleWrite = 0;
static int totalRead = 0;
int32 read_ni = 0;
float64 point;
float64 SF_IN;                 // scale factor for input
float64 SF_OUT;                // scale factor for output

// ---- Timing state ----
// TOLERANCE: well above clock_gettime ns resolution; controls busy-wait exit precision.
const long double TOLERANCE = 5e-7;         // in seconds (~500ns)
int LAST_READ = 0;                          // last read sample count
long double LAST_READ_T = ClockGetTime();   // last read time
long double new_read_T = 0;                 // new read time
long double now = ClockGetTime();           // current time
long double full_run_time = ClockGetTime(); // full run time
long double step_time_real;                 // time step in real time, in seconds
long double step_time_net;                  // time step in neural network time, in seconds
long double LAST_NET_T = 0;                 // last network time
long double total_debt = 0;                 // total debt in seconds
long int steps_taken = 0;                   // total number of steps taken
long double total_rate = 0;                 // total rate in seconds

// ---- Proxy-spike state ----
// Variables for the proxy spike mechanism: a hack to allow brian2 to run without
// freaking out about voltages that are too high (which would otherwise look like
// many spurious spikes).
int last_spike = 0;        // last spike time
long double vthresh = 0.0; // voltage threshold
long double vreset = 0.0;  // reset voltage
bool proxy_spike = false;  // proxy spike flag

// ---- Debug state ----
long double *read_times; // array for storing read times in debug mode

// ---- Time helpers ----

// Returns wall-clock time in seconds with full nanosecond precision.
// Uses CLOCK_MONOTONIC_RAW to avoid NTP adjtime slewing.
long double ClockGetTime()
{
    struct timespec ts; // local to avoid thread-safety issues
    clock_gettime(CLOCK_MONOTONIC_RAW, &ts);
    return (long double)ts.tv_sec + (long double)ts.tv_nsec / 1e9L;
}

int busySleep(uint32_t nanoseconds)
{
    struct timespec now;
    struct timespec then;
    struct timespec start;
    struct timespec sleep;
    if (nanoseconds > 999999999)
    {
        return 1;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    now = start;
    sleep.tv_sec = 0;
    sleep.tv_nsec = nanoseconds;
    tsadd(&start, &sleep, &then);
    while (tscmp(&now, &then, <))
    {
        clock_gettime(CLOCK_MONOTONIC_RAW, &now);
    }
    return 0;
}

// ---- NI low-level wrappers ----

extern "C" {

int nidaqrec(void)
{
    // DAQmx analog voltage channel and timing parameters

    // input task
    DAQmxErrChk(DAQmxCreateTask("", &taskHandle));
    // output task
    DAQmxErrChk(DAQmxCreateTask("", &taskHandleWrite));

    // Analog input channel
    DAQmxErrChk(DAQmxCreateAIVoltageChan(taskHandle, aIChan, "", DAQmx_Val_RSE, -1.0, 1.0, DAQmx_Val_Volts, NULL));
    // Analog output channel
    DAQmxErrChk(DAQmxCreateAOVoltageChan(taskHandleWrite, aOChan, "", -2.0, 2.0, DAQmx_Val_Volts, NULL));

    // DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"",SAMPLE_RATE,DAQmx_Val_Rising,DAQmx_Val_ContSamps,1000));

    // Ensure we only read and write the number of samples we expect
    DAQmxSetSampTimingType(taskHandleWrite, DAQmx_Val_OnDemand);
    DAQmxSetSampTimingType(taskHandle, DAQmx_Val_OnDemand);
    DAQmxSetReadOverWrite(taskHandle, DAQmx_Val_OverwriteUnreadSamps);

    // DAQmx Start Code
    DAQmxErrChk(DAQmxStartTask(taskHandle));
    DAQmxErrChk(DAQmxStartTask(taskHandleWrite));

    return 0; // success

Error:
    if (DAQmxFailed(error))
        DAQmxGetExtendedErrorInfo(errBuff, 2048);
    if (taskHandle != 0)
    {
        DAQmxStopTask(taskHandle);
        DAQmxClearTask(taskHandle);
        taskHandle = 0;
    }
    if (taskHandleWrite != 0)
    {
        DAQmxStopTask(taskHandleWrite);
        DAQmxClearTask(taskHandleWrite);
        taskHandleWrite = 0;
    }
    if (DAQmxFailed(error))
        printf("DAQmx Error: %s\n", errBuff);
    return -1; // failure
}

int read_sample()
{
    int32 err = DAQmxReadAnalogF64(taskHandle, 1, 1.0e-6, DAQmx_Val_GroupByScanNumber, &data, 1, &read_ni, NULL);
    if (DAQmxFailed(err))
    {
#ifdef DEBUG
        printf("DAQmx read error: %d\n", (int)err);
#endif
        return -1;
    }
    return 0;
}

int write_sample(float64 val)
{
    val = val * SF_OUT;
    int32 err = DAQmxWriteAnalogF64(taskHandleWrite, 1, 1, 1.0e-6, DAQmx_Val_GroupByScanNumber, &val, NULL, NULL);
    if (DAQmxFailed(err))
    {
#ifdef DEBUG
        printf("DAQmx write error: %d\n", (int)err);
#endif
        return -1;
    }
    return 0;
}

void clean_up_ni()
{
    // clear the task and reset the device
    DAQmxStopTask(taskHandle);
    DAQmxClearTask(taskHandle);
    DAQmxStopTask(taskHandleWrite);
    DAQmxClearTask(taskHandleWrite);
    DAQmxResetDevice("Dev1");
}

} // extern "C"

// ---- RT scheduling ----

int set_thread_priority_max()
{
    struct sched_param param;
    int ret;

    // Pin this thread to a dedicated core (core 4) so the busy-wait
    // doesn't compete with NI driver and kernel IRQ threads on other cores.
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(4, &cpuset); // pin to core 4 (0-indexed); adjust for your system
    ret = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (ret != 0)
    {
        printf("Warning: failed to set CPU affinity (error %d).\n", ret);
    }
    else
    {
        printf("Thread pinned to CPU core 4\n");
    }

    // Use mid-range FIFO priority (~50) so NI driver IRQ threads
    // (often priority 50-90) can still preempt during I/O waits.
    param.sched_priority = 50;
    ret = pthread_setschedparam(pthread_self(), SCHED_FIFO, &param);
    if (ret != 0)
    {
        printf("Warning: failed to set RT priority (error %d). Run as root or with CAP_SYS_NICE.\n", ret);
        return -1;
    }
    printf("Thread set to SCHED_FIFO with priority %d\n", param.sched_priority);
    return 0;
}

// ---- Public API ----

int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout, float64 runtime, char *aI, char *aO)
{
    // attempt RT scheduling (requires root or CAP_SYS_NICE)
    set_thread_priority_max();

    // set the sample rate to the network clock rate
    SAMPLE_RATE = 1 / (net_clock_dt / 1000);

    // set the scale factors
    SF_IN = scalein;
    SF_OUT = scaleout;

    // Only override defaults when a non-empty channel string is provided.
    if (aI != NULL && aI[0] != '\0')
    {
        strncpy(aIChan, aI, sizeof(aIChan) - 1);
        aIChan[sizeof(aIChan) - 1] = '\0';
    }
    if (aO != NULL && aO[0] != '\0')
    {
        strncpy(aOChan, aO, sizeof(aOChan) - 1);
        aOChan[sizeof(aOChan) - 1] = '\0';
    }

    printf("Using channels AI='%s' AO='%s'\n", aIChan, aOChan);

    // initialize the NI card
    if (nidaqrec() != 0)
    {
        printf("Error: NI card initialization failed\n");
        return -1;
    }
    printf("NI card initialized\n with scale factors: %f and %f\n", SF_IN, SF_OUT);

#ifdef DEBUG
    // initialize an empty array for storing read times (sampled every 1000 steps)
    int total_steps = (int)(runtime + 6) * SAMPLE_RATE;
    total_steps = total_steps / 1000;
    printf("Total steps: %d\n for a runtime of %d\n", total_steps, (int)runtime);
    read_times = (long double *)malloc((total_steps) * sizeof(long double));
    printf("Read times: %p\n", read_times);
    if (read_times == NULL)
    {
        printf("Error allocating memory for read times\n");
        return -1;
    }
#endif

    return 0;
}

void turn_on_proxy_spike(long double _vthresh, long double _vreset)
{
    // call this function prior to running the simulation to set the proxy spike values
    proxy_spike = true;
    vthresh = _vthresh / 1000; // convert to volts
    vreset = _vreset / 1000;
    printf("Proxy spike turned on with vthresh: %Lf and vreset: %Lf\n", vthresh, vreset);
}

double step_clamp(double t, double I)
{
    // t in seconds, I in pA

    step_time_net = (t - LAST_NET_T); // step in neural network time, in seconds
    if (step_time_net <= 0.0)
    {
        // if for some reason the network time is negative or zero, do nothing
        // and return the last value
        LAST_NET_T = t;
        full_run_time = ClockGetTime(); // reset the full run time
        LAST_READ_T = ClockGetTime();   // reset the last read time
        return data;
    }
    else
    {
        // write the sample to the NI card; first check how much time has
        // passed since last read. If the network time is ahead of real time,
        // wait; otherwise proceed.
        step_time_real = ClockGetTime() - LAST_READ_T;
        if (step_time_net < step_time_real)
        {
            // simulation fell behind real time — read but skip write
            total_debt += (step_time_real - step_time_net);
            // read the sample from the NI card (previous output remains)
            read_sample();
        }
        else
        {
            // busy-wait until real time catches up to network time
            while ((step_time_net - step_time_real) > TOLERANCE)
            {
                step_time_real = (ClockGetTime() - LAST_READ_T);
            }

            // inverted reading due to loop placement
            // brian2 integrates forward -> send current -> apply current ->
            // get signal -> volt used for next step
            write_sample(I * 1e9); // write the current to the NI card
            read_sample();
        }
    }

    LAST_NET_T = t;
    // Anchor to target time, not measured time: prevents overshoot from
    // accumulating across steps. Each deadline is exactly dt after the previous.
    LAST_READ_T += step_time_net;
    steps_taken++;
    // use target step size for rate reporting (debt tracked separately)
    total_rate += step_time_net;

#ifdef DEBUG
    // only store the read times every 1000 steps
    if (steps_taken % 1000 == 0)
    {
        read_times[steps_taken / 1000 - 1] = step_time_real;
    }
#endif

    if (proxy_spike)
    {
        // to trick brian2 we need to let the neuron go over the threshold for
        // one step; use last_spike to track threshold crossings
        if (data * SF_IN > vthresh)
        {
#ifdef DEBUG
            printf("Proxy spike at time: %Lf\n", t);
#endif
            if (last_spike == 0)
            {
                last_spike = 1;
            }
            else
            {
                data = vreset / SF_IN;
            }
        }
        else
        {
            last_spike = 0;
        }
    }

    return data * SF_IN;
}

double clean_up()
{
    long double rate = (total_rate / (long double)steps_taken) * 1000.0; // in ms

    printf("Average clamp rate: %Lf ms\n", rate);
    printf("Run time: %Lf with total delay debt of: %Lf\n", (LAST_READ_T - full_run_time), total_debt);

#ifdef DEBUG
    // dump the read times to a file
    FILE *fp;
    fp = fopen("/home/smestern/Dropbox/PVN_MODELLING_WORK/CADEX_MODEL/output/read_times.txt", "w");
    if (fp == NULL)
    {
        printf("Error opening file for writing\n");
        return -1;
    }
    for (int i = 0; i < steps_taken / 1000; i++)
    {
        fprintf(fp, "%Lf\n", read_times[i]);
    }
    fclose(fp);
    printf("Read times written to file\n");
    free(read_times);
#endif

    clean_up_ni();
    return 0.0;
}
