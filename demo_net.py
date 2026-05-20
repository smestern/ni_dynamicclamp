from brian2 import *
from ni_interface.ni_brian2 import *
import os
import time

DYN_CLAMP = True
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

IF = NeuronGroup(10, model=eqs, reset=reset, threshold='v>vt',
                 method='euler')
IF.vt = 10*mV
PG = PoissonGroup(1, 500 * Hz)

C = Synapses(PG, IF, on_pre='v += 3*mV')
C.connect()

eqs_real_neuron = '''
c = step_clamp(t, I_in) : volt
I_in : amp
'''

real_neuron = NeuronGroup(1, model=eqs_real_neuron, threshold='c>0*mV', refractory=5*ms, method='euler')

r_IF = Synapses(real_neuron, IF, on_pre='v += 3*mV')
r_IF.connect()

IF_r = Synapses(IF, real_neuron, 
model='''I_in_post = ge*(0*mV - c_post) : amp (summed)
dge/dt = -ge/(10*ms) : siemens
''', on_pre='ge +=  1e-2*nS')
IF_r.connect()






Mv = StateMonitor(real_neuron, ['c', 'I_in'], record=True)
Mv2 = StateMonitor(IF, 'v', record=True)
# Record the value of v when the threshold is crossed
M_crossings = SpikeMonitor(real_neuron, variables='c')
#run(2*second, report='text')



#before this we want to attach our real neuron
#here we call our helper to initialize the NIDAQ and start the recording.
#`step_clamp` and the C++ source/header wiring are provided by ni_interface.ni_brian2.
if DYN_CLAMP:
    device = init_neuron_device(device=device, dt=defaultclock.dt)

time_start = time.time()
run(20*second, report='text')
#print(time.time() - time_start)
plot(Mv.t/ms, Mv.c[0]/mV, label='v')
twinx()
plot(Mv2.t/ms, Mv.I_in[0]/pA, label='v2', c='r', alpha=0.1)
scatter(M_crossings.t/ms, np.full(len(M_crossings.t), 0), c='r', marker='x')
show()
