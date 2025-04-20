#include <stdio.h>
#include <chrono>
#include <thread>
#include <pthread.h>


extern "C" { 
#include <NIDAQmx.h>
#include <cstdlib>
#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else;

#ifdef HHDEBUG
#include <hh_tester.h>
#endif


}




extern "C" {float64 data = -0.070;}

#include <time.h>
#include <sys/timeb.h>
// needs -lrt (real-time lib)
// 1970-01-01 epoch UTC time, 1 mcs resolution (divide by 1M to get time_t)
struct timespec ts;
long double ClockGetTime()
{
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (long double)(ts.tv_sec * 1000000LL + ts.tv_nsec / 1000LL) / 1e6; // in seconds
}
# define tscmp(a, b, CMP)                             \
  (((a)->tv_sec == (b)->tv_sec) ?                         \
   ((a)->tv_nsec CMP (b)->tv_nsec) :                          \
   ((a)->tv_sec CMP (b)->tv_sec))
# define tsadd(a, b, result)                              \
  do {                                        \
    (result)->tv_sec = (a)->tv_sec + (b)->tv_sec;                 \
    (result)->tv_nsec = (a)->tv_nsec + (b)->tv_nsec;                  \
    if ((result)->tv_nsec >= 1000000000)                          \
      {                                       \
    ++(result)->tv_sec;                           \
    (result)->tv_nsec -= 1000000000;                          \
      }                                       \
  } while (0)
# define tssub(a, b, result)                              \
  do {                                        \
    (result)->tv_sec = (a)->tv_sec - (b)->tv_sec;                 \
    (result)->tv_nsec = (a)->tv_nsec - (b)->tv_nsec;                  \
    if ((result)->tv_nsec < 0) {                          \
      --(result)->tv_sec;                             \
      (result)->tv_nsec += 1000000000;                        \
    }                                         \
  } while (0)

int busySleep( uint32_t nanoseconds )
{
    struct timespec now;
    struct timespec then;
    struct timespec start;
    struct timespec sleep;
    if ( nanoseconds > 999999999 )
    {
        return 1;
    }
    clock_gettime( CLOCK_MONOTONIC_RAW, &start);
    now = start;
    sleep.tv_sec = 0;
    sleep.tv_nsec = nanoseconds;
    tsadd( &start, &sleep, &then );
    while ( tscmp( &now, &then, < )  )
    {
        clock_gettime( CLOCK_MONOTONIC_RAW, &now);
    }
    return 0;
}




int SAMPLE_RATE = 100000; //in Hz
int LAST_READ = 0; //last 
const long double TOLERANCE = 1e-12; //in seconds, the tolerance for the time difference between the network time and the code time
int32       error=0;
TaskHandle  taskHandle=0;
TaskHandle taskHandleWrite=0;
static int  totalRead=0;
int32       read_ni=0;
float64 point;
float64 SF_IN; //scale factor for input
float64 SF_OUT; //scale factor for output
long double LAST_READ_T = ClockGetTime(); //last read time
long double new_read_T = 0; //new read time
long double now = ClockGetTime(); //current time
long double full_run_time = ClockGetTime(); //full run time
long double step_time_real; //time steps in real time, in seconds
long double step_time_net; //time steps in neural network time, in seconds
long double LAST_NET_T = 0; //last network time
long double total_debt = 0; //total debt in seconds
long int steps_taken = 0; //total number of steps taken
long double total_rate = 0; //total rate in seconds
long double *read_times; //array for storing read times in debug mode
#ifdef HHDEBUG
HHModel model; //HH model

void _write_sample_HH(float64 val){
        //write the sample to the HH model
        hh_model_step(&model, val, (1.0/SAMPLE_RATE)*1000.0); //step the model forward in time
        data = model.Vm; //get the voltage from the model
}
void _read_sample_HH(){
        //read the sample from the HH model
        data = model.Vm; //get the voltage from the model
}
void _init_HH(float64 net_clock_dt, float64 scalein, float64 scaleout, float64 runtime){
        //set the sample rate to the network clock rate
        SAMPLE_RATE = 1/(net_clock_dt/1000);
        //set the scale factors
        SF_IN = scalein;
        SF_OUT = scaleout;
        //initialize the HH model
        hh_model_init(&model, 0.0); //initialize the model with a voltage offset of 0.0
}

#endif

extern "C" { //NI interface functions, in C style


int nidaqrec(void)
{
    
        char        errBuff[2048]={'\0'};

        // DAQmx analog voltage channel and timing parameters

        //input task
        DAQmxErrChk (DAQmxCreateTask("", &taskHandle));
        //output task
        DAQmxErrChk (DAQmxCreateTask("", &taskHandleWrite));
        
        //Analog input channel
        DAQmxErrChk(DAQmxCreateAIVoltageChan(taskHandle, "Dev1/ai0", "", DAQmx_Val_RSE, -1.0, 1.0, DAQmx_Val_Volts, NULL));
        //Analog output channel
        DAQmxErrChk(DAQmxCreateAOVoltageChan(taskHandleWrite, "Dev1/ao0", "", -2.0, 2.0, DAQmx_Val_Volts, NULL));

        //DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"",SAMPLE_RATE,DAQmx_Val_Rising,DAQmx_Val_ContSamps,1000));

        //Ensure we only read and write the number of samples we expect
        DAQmxSetSampTimingType(taskHandleWrite, DAQmx_Val_OnDemand);
        
        DAQmxSetSampTimingType(taskHandle,DAQmx_Val_OnDemand);
        DAQmxSetReadOverWrite(taskHandle, DAQmx_Val_OverwriteUnreadSamps);
        /*********************************************/
        // DAQmx Start Code
        /*********************************************/
        DAQmxErrChk (DAQmxStartTask(taskHandle));
        DAQmxErrChk (DAQmxStartTask(taskHandleWrite));



   
        // DAQmx Read Code

        //DAQmxErrChk(DAQmxReadAnalogF64(taskHandle, -1, -1, DAQmx_Val_GroupByChannel, data, 1000, &read, NULL));

        // Stop and clear task

        Error:
        if( DAQmxFailed(error) )
        DAQmxGetExtendedErrorInfo(errBuff,2048);
        if( taskHandle!=0 )  {
        //DAQmxStopTask(taskHandle);
        //DAQmxClearTask(taskHandle);
        }
        if( DAQmxFailed(error) )
        printf("DAQmx Error: %s\n",errBuff);
        return 50;
}


void read_sample(){
        DAQmxReadAnalogF64(taskHandle,1,1.0e-6,DAQmx_Val_GroupByScanNumber,&data,1,&read_ni,NULL);
}

void write_sample(float64 val){
        val = val*SF_OUT;
        DAQmxWriteAnalogF64(taskHandleWrite,1, 1, 1.0e-6, DAQmx_Val_GroupByScanNumber,&val,NULL,NULL);
}

void clean_up_ni(){
        //clear the task and reset the device
        DAQmxStopTask(taskHandle);
        DAQmxClearTask(taskHandle);
        DAQmxStopTask(taskHandleWrite);
        DAQmxClearTask(taskHandleWrite);
        DAQmxResetDevice("Dev1");
}



}

double clean_up(){
        long double rate = (total_rate/(long double)steps_taken) * 1000.0; //in ms
        

        printf("Average clamp rate: %Lf ms\n", rate);
        printf("Run time: %Lf with total delay debt of: %Lf\n", (LAST_READ_T - full_run_time), total_debt);

        #ifdef DEBUG
        //dump the read times to a file
        FILE *fp;
        fp = fopen("/home/smestern/Dropbox/PVN_MODELLING_WORK/CADEX_MODEL/output/read_times.txt", "w");
        if (fp == NULL) {
                printf("Error opening file for writing\n");
                return -1;
        }
        for (int i = 0; i < steps_taken/1000; i++) {
                fprintf(fp, "%Lf\n", read_times[i]);
        }
        fclose(fp);
        printf("Read times written to file\n");
        free(read_times);
        #endif

        clean_up_ni();  //clean up NI
        return 0.0;
}

int set_thread_priority_max(){
        int policy;
        struct sched_param param;
        pthread_attr_t attr;
        int ret;

        ret = pthread_attr_setschedpolicy(&attr, SCHED_FIFO);
        printf("ret: %d\n", ret);
        pthread_getschedparam(pthread_self(), &policy, &param);
        param.sched_priority = sched_get_priority_max(SCHED_FIFO); //set to RT priority
        ret = pthread_setschedparam(pthread_self(), policy, &param);
        return 0;
}



int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout, float64 runtime){
        //set the sample rate to the network clock rate
        set_thread_priority_max();
        SAMPLE_RATE = 1/(net_clock_dt/1000);
        //set the scale factors
        SF_IN = scalein;
        SF_OUT = scaleout;
        //initialize the NI card
        nidaqrec();
        //intialize an empty array for storing read times
        #ifdef DEBUG
        int total_steps = (int)(runtime+6)*SAMPLE_RATE; //total number of steps to take
        //actually only sample every 1000 steps, so divide by 1000
        total_steps = total_steps / 1000; // Corrected division
        printf("Total steps: %d\n for a runtime of %d\n", total_steps, (int)runtime);
        //allocate memory for the read times
        read_times = (long double*)malloc((total_steps)*sizeof(long double));
        printf("Read times: %p\n", read_times);
        if (read_times == NULL) {
                printf("Error allocating memory for read times\n");
                return -1;
        }
        #endif

        #ifdef HHDEBUG
        //initialize the HH model
        _init_HH(net_clock_dt, scalein, scaleout, runtime);

        //override the read and write functions to use the HH model
        write_sample = []_write_sample_HH;
        read_sample = _read_sample_HH;
        #else
        #endif


        return 0;
}

double step_clamp(double t, double I) {
        // t in seconds , I in pA

        step_time_net = (t - LAST_NET_T); //time steps in neural network time, in seconds
        if (step_time_net <= 0.0){
                //if for some reason the network time is negative, or zero, do nothing and return the last value
                LAST_NET_T = t;
                full_run_time = ClockGetTime();; //reset the full run time
                LAST_READ_T = ClockGetTime(); //reset the last read time
                return data;

        } else {
                

                //write the sample to the NI card
                //check how much time has passed since last read
                step_time_real = ClockGetTime() - LAST_READ_T; //time steps since last call of read, also in seconds
                 //check if the network time is ahead of the code time, if so, wait, otherwise proceed
                if (step_time_net < step_time_real) { //if neural network time is ahead of code time, wait, otherwise proceed
                        //printf("Code running slower than real time with a delay of: %lf\n with a network step of : %lf\n and a real time of: %lf", 1000*(step_time_real-step_time_net), 1000*(step_time_net), 1000*(step_time_net));
                        //dont write just read and return
                        total_debt += (step_time_real - step_time_net);
                } else {
                        
                        
                        write_sample(I*1e9); //write the current to the NI card
                        //force wait to slow the network down to match code time
                        //these funcions seem to be inaccurate, so we will use a busy wait instead, 
                        //stealing idea from https://github.com/CompEphys-team/stdpc/blob/17fa31b760e5f9210e22c84b2a702ce182971be4/src/drivers/Clock.cpp#L68
                        while ((step_time_net - step_time_real)>TOLERANCE){ //busy wait until the network time is behind the code time
                                step_time_real = (ClockGetTime() - LAST_READ_T);
                        }
                        //printf("Network running faster than real time with a step diff of: %Lf\n with a network step of : %Lf\n and a real time of: %Lf", 1000*(step_time_net - step_time_real), 1000*(step_time_real), 1000*(step_time_net));
                        
                }
                
                //read the sample from the NI card
                read_sample();
        }
       
       
        LAST_NET_T = t;
        LAST_READ_T = step_time_real + LAST_READ_T;
        steps_taken++;
        total_rate+=step_time_real;
        #ifdef DEBUG
        //printf("Step time: %Lf\n", step_time_real);
        //printf("Network time: %Lf\n", step_time_net);
        //printf("steps taken: %ld\n", steps_taken);
        //only store the read times every 1000 steps
        if (steps_taken % 1000 == 0) {
                read_times[steps_taken/1000-1] = step_time_real;
        }
        
        #endif
        return data*SF_IN;
}


// int run_step_loop(double *I, double *out){

//         for (int i = 0; i < sizeof(I)/sizeof(double); i++){
//                 out[i] = step_clamp(i*0.0001, I[i]);
//         }
//         return 0;
// }
