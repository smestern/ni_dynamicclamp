from brian2 import *
from brian2.equations.equations import SUBEXPRESSION, PARAMETER, DIFFERENTIAL_EQUATION
import os
import numpy
from brian2.utils.topsort import topsort
from collections import namedtuple

#Create the equation object to bypass the default equations check

#before this we want to attach our real neuron
current_dir = os.path.abspath(os.path.dirname(__file__))

#Here we link to our source and header files. 

@implementation('cpp', '''//''',  sources=[os.path.join(current_dir,
                                      'interface.cpp'), os.path.join(current_dir,'libnidaqmx.so')], 
                headers=['"interface.h"', '"NIDAQmx.h"'], 
                include_dirs=[current_dir])
@check_units(t=second, I=pA, result=mV)
def step_clamp(t, I):
    raise NotImplementedError('step_clamp should not be called directly, this function is replaced by the C++ code')

# Not sure exactly what fixed this but its currently working, otherwise it struggles to find the header files 
# even though they are in the same directory
# IDK if I can even legally distribute this code with the header files included
# but im just a poor grad student pls don't sue NI 


def init_neuron_device(device, dt=defaultclock.dt, scalefactor_in=0.1, scalefactor_out=1/0.5,  runtime=1.0, aI=None, aO=None, proxy_spike=False):
    """Helper function to initialize the NIDAQ and start the recording, this should be called before running the simulation
    Args:
        device: the brian2 device to insert the code into
        dt: the time step of the simulation, this should match the defaultclock.dt
        scalefactor_in: the scale factor for the input current, this is used to convert from pA to V, defaults to 0.1 (1V = 10pA), this should be set to 0.1 for a 100 MOhm electrode
        scalefactor_out: the scale factor for the output voltage, this is used to convert from V to pA, defaults to 1/0.5 (1pA = 0.5V) for a 100 MOhm electrode
        aI: the analog input channel, if None it will use the default "Dev1/ai0"
        aO: the analog output channel, if None it will use the default "Dev1/ao0"
        runtime: the total runtime of the simulation in seconds, this is used to preallocate memory for storing read times in debug mode
        proxy_spike: whether to turn on the proxy spike mechanism, this is a hack to allow brian2 to run without freaking out about voltages that are too high.
    returns:
        device: the brian2 device with the NIDAQ code inserted
    """
    #one of these works cant figure out what fixed it
    prefs.codegen.cpp.include_dirs = [current_dir]
    prefs.codegen.cpp.library_dirs = [current_dir]
    prefs.codegen.cpp.headers = ['"interface.h"', '"NIDAQmx.h"'] #again link the header files, not sure if this is necessary but it seems to be

    #here we call our function to intialize the NIDAQ and start the recording
    init_ni_str = f'init_ni({dt/ms}, {scalefactor_in}, {scalefactor_out}, {runtime}, "{aI if aI is not None else ""}", "{aO if aO is not None else ""}");'
    device.insert_code('after_start', init_ni_str) #we want to make sure that the sampling rate is the same as the defaultclock.dt
    if proxy_spike:
        device.insert_code('after_start', f'turn_on_proxy_spike(-30., -70.);') #
    device.insert_code('before_end', 'clean_up();') #clean up the NIDAQ, log the time
    return device

def attach_neuron_brute_force(neurongroup, eqs, idx=0, v_mem_var='v', i_mem_var='I_in', dt=None, when='before_thresholds'):
    '''
    Helper function to subsitute neuron in brian2 neuron group with an invitro neuron. 
    '''
    #Now we create a new neuron group with the same equations as the original group, but with only one neuron, and we replace the equations with our own equations that include the step_clamp function,
    #sometimes this helps with speed in brian2. But can be a bit hacky and may not work in all cases, especially if the original equations are complex or have a lot of state variables.
    
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