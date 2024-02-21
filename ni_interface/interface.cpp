#include <stdio.h>
#include <chrono>
#include <thread>
#include <pthread.h>


extern "C" { 
#include <NIDAQmx.h>

#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else;
}


extern "C" {float64 data = -0.070;}

#include <time.h>
#include <sys/timeb.h>
// needs -lrt (real-time lib)
// 1970-01-01 epoch UTC time, 1 mcs resolution (divide by 1M to get time_t)
struct timespec ts;
long double ClockGetTime()
{
    clock_gettime(CLOCK_REALTIME, &ts);
    return (long double)(ts.tv_sec * 1000000LL + ts.tv_nsec / 1000LL) / 1e6; // in seconds
}

int SAMPLE_RATE = 100000; //in Hz
int LAST_READ = 0; //last 
const long double TOLERANCE = 1e-9; //in seconds, the tolerance for the time difference between the network time and the code time
int32       error=0;
TaskHandle  taskHandle=0;
TaskHandle taskHandleWrite=0;
static int  totalRead=0;
int32       read=0;
float64 point;
float64 SF_IN;
float64 SF_OUT;
long double LAST_READ_T = ClockGetTime();
long double now = ClockGetTime();
long double full_run_time = ClockGetTime();
long double step_time_real;
long double step_time_net;
long double LAST_NET_T = 0;



extern "C" {


int nidaqrec(void)
{
    
        char        errBuff[2048]={'\0'};

        // DAQmx analog voltage channel and timing parameters

        //input task
        DAQmxErrChk (DAQmxCreateTask("", &taskHandle));
        //output task
        DAQmxErrChk (DAQmxCreateTask("", &taskHandleWrite));
        
        //Analog input channel
        DAQmxErrChk(DAQmxCreateAIVoltageChan(taskHandle, "Dev2/ai0", "", DAQmx_Val_RSE, -1.0, 1.0, DAQmx_Val_Volts, NULL));
        //Analog output channel
        DAQmxErrChk(DAQmxCreateAOVoltageChan(taskHandleWrite, "Dev2/ao0", "", -1.0, 1.0, DAQmx_Val_Volts, NULL));

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
        DAQmxReadAnalogF64(taskHandle,1,1.0e-6,DAQmx_Val_GroupByScanNumber,&data,1,&read,NULL);
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
        DAQmxResetDevice("Dev2");
}



}

double clean_up(){
        printf("%Lf/n", (LAST_READ_T - full_run_time));
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



int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout){
        //set the sample rate to the network clock rate
        set_thread_priority_max();
        SAMPLE_RATE = 1/(net_clock_dt/1000);
        //set the scale factors
        SF_IN = scalein;
        SF_OUT = scaleout;
        //initialize the NI card
        nidaqrec();
       
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
                //read the sample from the NI card
                read_sample();
                //write the sample to the NI card
                //check how much time has passed since last read
                step_time_real = ClockGetTime() - LAST_READ_T; //time steps since last call of read, also in seconds
                
                //time steps since last call of read
        }
        
        //check if the network time is ahead of the code time, if so, wait, otherwise proceed
        if (step_time_net < step_time_real) { //if neural network time is ahead of code time, wait, otherwise proceed
                //printf("Code running slower than real time with a delay of: %lf\n with a network step of : %lf\n and a real time of: %lf", 1000*(step_time_real-step_time_net), 1000*(step_time_net), 1000*(step_time_net));
                //dont write just read and return
        } else {
                
                //printf("Network running faster than real time with a step diff of: %Lf\n with a network step of : %Lf\n and a real time of: %Lf", 1000*(step_time_net - step_time_real), 1000*(step_time_real), 1000*(step_time_net));
                write_sample(I*1e9); //write the current to the NI card
                //force wait to slow the network down to match code time
                //these funcions seem to be inaccurate, so we will use a busy wait instead, 
                //stealing idea from https://github.com/CompEphys-team/stdpc/blob/17fa31b760e5f9210e22c84b2a702ce182971be4/src/drivers/Clock.cpp#L68
                while ((step_time_net - step_time_real)>TOLERANCE){ //busy wait until the network time is behind the code time
                        step_time_real = (ClockGetTime() - LAST_READ_T);
                }
                //printf("Network running faster than real time with a step diff of: %Lf\n with a network step of : %Lf\n and a real time of: %Lf", 1000*(step_time_net - step_time_real), 1000*(step_time_real), 1000*(step_time_net));
                
        }
        LAST_NET_T = t;
        LAST_READ_T = ClockGetTime();
        
        return data*SF_IN;
}


// int run_step_loop(double *I, double *out){

//         for (int i = 0; i < sizeof(I)/sizeof(double); i++){
//                 out[i] = step_clamp(i*0.0001, I[i]);
//         }
//         return 0;
// }
