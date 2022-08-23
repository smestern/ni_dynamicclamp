from brian2 import *
from brian2.equations.equations import SUBEXPRESSION, PARAMETER, DIFFERENTIAL_EQUATION
import os
import numpy
from brian2.utils.topsort import topsort
from collections import namedtuple

#Create the equation object to bypass the default equations check

#before this we want to attach our real neuron
current_dir = os.path.abspath(os.path.dirname(__file__))

#Here we link to our source and header files. Not sure exactly what fixed this but its currently working
@implementation('cpp', '''//''',  sources=[os.path.join(current_dir,
                                      'interface.cpp'), os.path.join(current_dir,'libnidaqmx.so')],
                headers=['"interface.h"', '"NIDAQmx.h"'],
                include_dirs=[current_dir], )
@check_units(t=second, I=pA, result=mV)
def step_clamp(t, I):
    return -999*mV/second


def init_neuron_device(device, dt=defaultclock.dt, scalefactor_in=0.1, scalefactor_out=1/0.5):
    #one of these works cant figure out what fixed it
    prefs.codegen.cpp.include_dirs = [current_dir]
    prefs.codegen.cpp.library_dirs = [current_dir]
    prefs.codegen.cpp.headers = ['"interface.h"', '"NIDAQmx.h"']

    #here we call our function to intialize the NIDAQ and start the recording
    device.insert_code('after_start', f'init_ni({defaultclock.dt/ms}, {scalefactor_in}, {scalefactor_out});') #we want to make sure that the sampling rate is the same as the defaultclock.dt

    device.insert_code('before_end', 'clean_up();') #clean up the NIDAQ, log the time
    return device

def attach_neuron(neurongroup, idx=0, v_mem_var='v', i_mem_var='I_in'):
    '''
    Helper function to subsitute neuron in brian2 neuron group with an invitro neuron. 
    '''
    #slice the neuron group to only include the neuron we want to attach
   
    neuron = neurongroup[idx]
    
    #add a run regularly statement to the neuron group
    neuron.run_regularly(f'{v_mem_var} = step_clamp(t, {i_mem_var})', dt=defaultclock.dt)
    #update the resetter of the full group to exclude the neuron we are attaching
    #replace the remaining idxs with a gapjunction like synapse


    return neuron, neurongroup