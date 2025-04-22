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
                                      'interface_adex.cpp'), os.path.join(current_dir,'libnidaqmx.so')],
                headers=['"interface.h"', '"NIDAQmx.h"'],
                include_dirs=[current_dir])
@check_units(t=second, I=pA, result=mV)
def step_clamp(t, I):
    raise NotImplementedError('step_clamp should not be called directly, this function is replaced by the C++ code')


def init_neuron_device(device, dt=defaultclock.dt, scalefactor_in=0.1, scalefactor_out=1/0.5, runtime=1.0):
    prefs.codegen.cpp.include_dirs = [current_dir]
    prefs.codegen.cpp.library_dirs = [current_dir]
    #dummy script for compatibility with the old code, don't worry about this
    return device

def attach_neuron(neurongroup, idx=0, v_mem_var='v', i_mem_var='I_in', dt=None, when='before_thresholds'):
    '''
    Helper function to subsitute neuron in brian2 neuron group with an invitro neuron. 
    '''
    #slice the neuron group to only include the neuron we want to attach
   
    neuron = neurongroup[idx]
    
    #add a run regularly statement to the neuron group
    neuron.run_regularly(f'{v_mem_var} = step_clamp(t, {i_mem_var})', dt=dt, when=when)
    #update the resetter of the full group to exclude the neuron we are attaching
    #replace the remaining idxs with a gapjunction like synapse


    return neuron, neurongroup