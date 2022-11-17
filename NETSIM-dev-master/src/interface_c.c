#include <stdio.h>



#include <NIDAQmx.h>

#define DAQmxErrChk(functionCall) if( DAQmxFailed(error=(functionCall)) ) goto Error; else;




float64 data;

int SAMPLE_RATE = 100000; //in Hz
int LAST_READ = 0;
int32       error=0;
TaskHandle  taskHandle=0;
TaskHandle taskHandleWrite=0;
static int  totalRead=0;
int32       read=0;
float64 point;
float64 SF_IN;
float64 SF_OUT;
double LAST_READ_T = 0;
double LAST_NET_T = 0;
double now = 0;
double step_time_real;
double step_time_net;
double full_run_time = 0;




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
        DAQmxErrChk(DAQmxCreateAOVoltageChan(taskHandleWrite, "Dev1/ao0", "", -1.0, 1.0, DAQmx_Val_Volts, NULL));

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
        DAQmxReadAnalogF64(taskHandle,1,1.0e-4,DAQmx_Val_GroupByScanNumber,&data,1,&read,NULL);
}

void write_sample(float64 val){
        val = val*SF_OUT;
        DAQmxWriteAnalogF64(taskHandleWrite,1, 1, 1.0e-4, DAQmx_Val_GroupByScanNumber,&val,NULL,NULL);
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
                LAST_NET_T = t;
                read_sample();
                return data*SF_IN;

        } else {
                //read the sample from the NI card
                read_sample();
                //write the sample to the NI card
                write_sample(I*1e9);
                //check how much time has passed since last read
                step_time_real = 0;//0std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - LAST_READ_T).count();
                //time steps since last call of read
        }
        
        
        
        if (step_time_net < step_time_real) { //if neural network time is ahead of code time, wait, otherwise proceed
                printf("Code running slower than real time with a delay of: %lf\n", -1*(step_time_net - step_time_real));
        } else {
                //force wait to slow the network down to match code time
                //std::this_thread::sleep_for(std::chrono::duration<double>((step_time_net - step_time_real))); //sleep for the difference in time in miliseconds
                //printf("Network running faster than real time with a step diff of: %lf\n", (step_time_net - step_time_real));

        }
        LAST_NET_T = t;
        LAST_READ_T = 0;//std::chrono::high_resolution_clock::now();
        
        return data*SF_IN;
    }

