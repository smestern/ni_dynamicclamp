#include <stdio.h>
#include <NIDAQmx.h>
#include <time.h>
#include "interface_c.h"
#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else;

float64 data[1];
int LAST_READ = 0;
int32 error=0;
TaskHandle taskHandle=0;
TaskHandle taskHandleWrite=0;
int totalRead=0;
int32 read[1];
float64 point;
float64 SF_IN;
float64 SF_OUT;
double LAST_READ_T = 0;
double LAST_NET_T = 0;
double now = 0;
double step_time_real;
double step_time_net;
double full_run_time = 0;
double SAMPLE_RATE = 100000;

int nidaqrec(void)
{
    
        //char  errBuff[2048]={'\0'};

        // DAQmx analog voltage channel and timing parameters

        //input task
        DAQmxCreateTask("", &taskHandle);
        //output task
        DAQmxCreateTask("", &taskHandleWrite);
        
        //Analog input channel
        DAQmxCreateAIVoltageChan(taskHandle, "Dev2/ai0", "", DAQmx_Val_RSE, -1.0, 1.0, DAQmx_Val_Volts, NULL);
        //Analog output channel
        DAQmxCreateAOVoltageChan(taskHandleWrite, "Dev2/ao0", "", -1.0, 1.0, DAQmx_Val_Volts, NULL);

        //DAQmxErrChk (DAQmxCfgSampClkTiming(taskHandle,"",SAMPLE_RATE,DAQmx_Val_Rising,DAQmx_Val_ContSamps,1000));

        //Ensure we only read and write the number of samples we expect
        DAQmxSetSampTimingType(taskHandleWrite, DAQmx_Val_OnDemand);
        
        DAQmxSetSampTimingType(taskHandle,DAQmx_Val_OnDemand);
        DAQmxSetReadOverWrite(taskHandle, DAQmx_Val_OverwriteUnreadSamps);
        /*********************************************/
        // DAQmx Start Code
        /*********************************************/
        DAQmxStartTask(taskHandle);
        DAQmxStartTask(taskHandleWrite);



   
        // DAQmx Read Code

        //DAQmxErrChk(DAQmxReadAnalogF64(taskHandle, -1, -1, DAQmx_Val_GroupByChannel, data, 1000, &read, NULL));

        // Stop and clear task

        Error:
        if( DAQmxFailed(error) )
        //DAQmxGetExtendedErrorInfo(errBuff,2048);
        if( taskHandle!=0 )  {
        //DAQmxStopTask(taskHandle);
        //DAQmxClearTask(taskHandle);
        }
        if( DAQmxFailed(error) )
        //printf("DAQmx Error: %s\n",errBuff);
        return 50;
}


void read_sample(){
        int32 newread[1];
        int32 read_code;
        //print the task handle
        
        //read the s
        read_code = DAQmxReadAnalogF64(taskHandle,1,1,DAQmx_Val_GroupByChannel,data,1,newread,NULL);
        
} 

void write_sample(float64 val){
        float64 newval[1];
        int32 write_code[1];
        newval[0] = val*SF_OUT;
        DAQmxWriteAnalogF64(taskHandleWrite,1, 1, 1, DAQmx_Val_GroupByChannel,newval,write_code,NULL);
}

void clean_up_ni(){
        //clear the task and reset the device
        DAQmxStopTask(taskHandle);
        DAQmxClearTask(taskHandle);
        DAQmxStopTask(taskHandleWrite);
        DAQmxClearTask(taskHandleWrite);
        DAQmxResetDevice("Dev1");
}




double clean_up(){
        clean_up_ni();  //clean up NI

        return 0.0;
}



int init_ni(float64 net_clock_dt, float64 scalein, float64 scaleout){
        //set the sample rate to the network clock rate
        SAMPLE_RATE = 1/(net_clock_dt/1000);
        nidaqrec();
        full_run_time = 0;
        SF_IN = scalein;
        SF_OUT = scaleout;
        printf("NI initialized");
        return 0;
}

double step_clamp(double t, double I) {
        // t in seconds , I in pA
         
        step_time_net = (t - LAST_NET_T); //time steps in neural network time
        
        if (step_time_net <= 0.0){
                //if for some reason the network time is negative, or zero, do nothing and return the last value
                printf("step_time_net: %lf\n", step_time_net);
                printf("t: %lf\n", t);
                printf("LAST_NET_T: %lf\n", LAST_NET_T);
                LAST_NET_T = t;
                return data[0]*SF_IN;

        } else {
                printf("step_time_net: %lf\n", step_time_net);
                //read the sample from the NI card
                read_sample(); 
                //write the sample to the NI card
                write_sample(I*1e9);
                //check how much time has passed since last read
                step_time_real = (double)time(NULL) - LAST_READ_T;
                //time steps since last call of read
        }
        
        
        
        if (step_time_net < step_time_real) { //if neural network time is ahead of code time, wait, otherwise proceed
                printf("Code running slower than real time with a delay of: %lf\n", -1*(step_time_net - step_time_real));
        } else {
                //force wait to slow the network down to match code time
                //sleep for the difference in time in miliseconds
                sleep((step_time_net - step_time_real));
                printf("Network running faster than real time with a step diff of: %lf\n", (step_time_net - step_time_real));

        }
        LAST_NET_T = t;
        LAST_READ_T = (double)time(NULL);
        //printf("data[0]: %lf\n", data[0]);
        return data[0]*SF_IN;
    }

int run_step_loop(double *I, double *out, int size){
        //loop over I and out and call step_clamp
        double dt = 1/SAMPLE_RATE;
       
        for (int i = 0; i < size; i++){
                out[i] = step_clamp(i*dt, I[i]);
                printf("out[%d]: %lf\n", i, out[i]);
        }
        return 0;

}