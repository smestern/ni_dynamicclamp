from brian2 import *
from ni_interface.ni_brian2 import *
import os
import time
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use('agg')
seed(43)
plt.ioff()  
DYN_CLAMP = True
defaultclock.dt = 0.1*ms
set_device('cpp_standalone', build_on_run=False)

start_scope()
start_time = time.time()
seed(43)

#ramp G to show ringing
g_values = np.linspace(0, 100, 20)

g_values = [[0, g] for g in g_values]
g_values = np.hstack(g_values)

g_timed = TimedArray(np.hstack([g_values, g_values])*nS, dt=500*ms)

# Membrane Equation + y*(Ee-v)
eqs = Equations('''
dv/dt = (v_rest - v) / (10*ms) : volt
d_I = (g*(Ee-v)) : amp
g = g_timed(t) : siemens
Ee : volt
v_rest : volt''')

neurons = NeuronGroup(1, eqs)
neurons.v_rest = -70*mV

neurons.Ee = 0*mV

neurons_dyn, neurons = attach_neuron(neurons, 0, i_mem_var="d_I", dt=defaultclock.dt, when='start')



Mv = StateMonitor(neurons_dyn, ['v', 'd_I', 'g'], record=True)

device = init_neuron_device(device=device, dt=defaultclock.dt, runtime=10)

run(len(g_values)*500*ms)
neurons.Ee = -80*mV
#same but with Ee = -80 mV
run(len(g_values)*500*ms)

device.build()

print("Simulation time:", time.time() - start_time)
plt.figure()
plt.plot(Mv.t/ms, Mv.v[0]/mV, color='blue', label='v')
plt.twinx()
plt.plot(Mv.t/ms, Mv.d_I[0]/nA, color='red', label='d_I')
plt.plot(Mv.t/ms, Mv.g[0]/nS, color='green', label='g')
plt.legend(loc='best')
plt.savefig("membrane_potential.png")
plt.show()