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
                include_dirs=[current_dir])
@check_units(t=second, I=pA, result=mV)
def step_clamp(t, I):
    raise NotImplementedError('step_clamp should not be called directly, this function is replaced by the C++ code')


def init_neuron_device(device, dt=defaultclock.dt, scalefactor_in=0.1, scalefactor_out=1/0.5, runtime=1.0, proxy_spike=False):
    #one of these works cant figure out what fixed it
    prefs.codegen.cpp.include_dirs = [current_dir]
    prefs.codegen.cpp.library_dirs = [current_dir]
    prefs.codegen.cpp.headers = ['"interface.h"', '"NIDAQmx.h"']

    #here we call our function to intialize the NIDAQ and start the recording
    device.insert_code('after_start', f'init_ni({dt/ms}, {scalefactor_in}, {scalefactor_out}, {runtime});') #we want to make sure that the sampling rate is the same as the defaultclock.dt
    if proxy_spike:
        device.insert_code('after_start', f'turn_on_proxy_spike(-30., -70.);') #
    device.insert_code('before_end', 'clean_up();') #clean up the NIDAQ, log the time
    return device

def attach_neuron_brute_force(neurongroup, eqs, idx=0, v_mem_var='v', i_mem_var='I_in', dt=None, when='before_thresholds'):
    '''
    Helper function to subsitute neuron in brian2 neuron group with an invitro neuron. 
    '''
    #slice the neuron group to only include the neuron we want to attach
    #if eqs is a equation object we need to bruteforce sub
    
    neuron = NeuronGroup(1, model=eqs, threshold='v>Vcut', refractory=f'v>(Vcut - 3*mV)', method='euler')
    
    #add a run regularly statement to the neuron group
    neuron.run_regularly(f'{v_mem_var} = step_clamp(t, {i_mem_var})', dt=dt, when=when)
    #update the resetter of the full group to exclude the neuron we are attaching
    

    return neuron, neurongroup

def attach_neuron_proxy(neurongroup, idx=0, v_mem_var='v', i_mem_var='I_in', dt=None, when='before_thresholds'):
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

def attach_neuron(neurongroup, eqs=None, idx=0, v_mem_var='v', i_mem_var='I_in', dt=None, when='before_thresholds', method='proxy'):
    '''
    Helper function to subsitute neuron in brian2 neuron group with an invitro neuron. 
    '''
    #brute force method
    if method == 'brute_force':
        return attach_neuron_brute_force(neurongroup, eqs, idx, v_mem_var, i_mem_var, dt, when)
    if method == 'proxy':
        return attach_neuron_proxy(neurongroup, idx, v_mem_var, i_mem_var, dt, when)
    else :
        raise ValueError('Method must be either brute_force or substitution')
    
    return neuron, neurongroup