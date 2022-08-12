from brian2 import *
import os
import time
defaultclock.dt = 0.05*ms
set_device('cpp_standalone', build_on_run=True)
eqs = '''
dv/dt = -v/(10*ms) : volt
dvt/dt = (10*mV-vt)/(15*ms) : volt
'''

reset = '''
v = 0*mV
vt += 3*mV
'''

IF = NeuronGroup(1000, model=eqs, reset=reset, threshold='v>vt',
                 method='euler')
IF.vt = 10*mV
PG = PoissonGroup(1, 500 * Hz)

C = Synapses(PG, IF, on_pre='v += 3*mV')
C.connect()

eqs_real_neuron = '''
c = step_clamp(t, 0*pA) : volt'''

real_neuron = NeuronGroup(1, model=eqs_real_neuron, threshold='c>0*mV', method='euler')

Mv = StateMonitor(real_neuron, ['c'], record=True)
#Mvt = StateMonitor(IF, 'vt', record=True)
# Record the value of v when the threshold is crossed
#M_crossings = SpikeMonitor(IF, variables='v')
#run(2*second, report='text')



#before this we want to attach our real neuron
current_dir = os.path.abspath(os.path.dirname(__file__))
current_dir = "/home/smestern/Dropbox/RTXI/ni_interface"
#Here we link to our source and header files. Not sure exactly what fixed this but its currently working
@implementation('cpp', '''//''',  sources=[os.path.join(current_dir,
                                      'interface.cpp'), '/home/smestern/Dropbox/RTXI/ni_interface/libnidaqmx.so'],
                headers=['"interface.h"', '"NIDAQmx.h"'],
                include_dirs=[current_dir], )
@check_units(t=second, I=pA, result=mV)
def step_clamp(t, I):
    return -999*mV



#one of these works cant figure out what fixed it
prefs.codegen.cpp.include_dirs = [current_dir]
prefs.codegen.cpp.library_dirs = [current_dir]
prefs.codegen.cpp.headers = ['"interface.h"', '"NIDAQmx.h"']

#here we call our function to intialize the NIDAQ and start the recording
device.insert_code('after_start', f'init_ni({defaultclock.dt/ms});') #we want to make sure that the sampling rate is the same as the defaultclock.dt

device.insert_code('before_end', 'clean_up();') #clean up the NIDAQ, log the time

time_start = time.time()
run(20*second, report='text')
#print(time.time() - time_start)
plot(Mv.t/ms, Mv.c[0], label='v')
show()
