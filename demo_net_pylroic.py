
# # Example 1 (Pyloric network of the crustacean stomatogastric ganglion) – biophysically detailed model
# 
# This example is equivalent to the [simplified model](example_1_pyloric_network.ipynb), but uses a more biophysically detailed neuron model. For a detailed explanation of the modelling approach and the conversion of equations into Brian syntax, see the comments in that file.
# 
# For details about the neuron model, see:
# 
# Golowasch, J., Casey, M., Abbott, L. F., & Marder, E. (1999).  
# Network Stability from Activity-Dependent Regulation of Neuronal Conductances.  
# Neural Computation, 11(5), 1079-1096.  
# https://doi.org/10.1162/089976699300016359

# %%
from brian2.only import *
from ni_interface.ni_brian2 import *
import numpy as np

set_device('cpp_standalone', build_on_run=False)
defaultclock.dt = 0.01*ms

init_time = 0.5*second
observe_time = 3*second
adapt_time = 24 * second

### Class-independent constants
# Reversal potentials
E_L = -68*mV
E_Na = 20*mV
E_K = -80*mV
E_Ca = 120*mV
E_proc = -10*mV

# Capacitance
C_s = 0.2*nF
C_a = 0.02*nF

# maximal conductances
g_E = 10*nS
g_La = 7.5*nS
g_Na = 300*nS
g_Kd = 4*uS
G_Ca = 0.2*uS
G_K = 16*uS

# time constants (independent of V)
tau_h_Ca = 150*ms
tau_m_A = 0.1*ms
tau_h_A = 50*ms
tau_m_proc = 6*ms
tau_m_Na = 0.025*ms
tau_z = 5*second

# Synapses
s_fast = 0.2/mV
V_fast = -50*mV
s_slow = 1/mV
V_slow = -55*mV
E_syn = -75*mV
k_1 = 1/ms

### Neuronal equations
eqs = '''
# somatic compartment
dV_s/dt = (-I_syn - I_L - I_Ca - I_K - I_A - I_proc - g_E*(V_s - V_a))/C_s * neur_dyn: volt
I_total = clip(I_syn, -150*pA, 150*pA) : amp
I_L = g_Ls*(V_s - E_L) : amp
I_K = g_K*m_K**4*(V_s - E_K) : amp
I_A = g_A*m_A**3*h_A*(V_s - E_K) : amp
I_proc = g_proc*m_proc*(V_s - E_proc) : amp

I_syn = I_fast + I_slow: amp
I_fast : amp
I_slow : amp

I_Ca = g_Ca*m_Ca**3*h_Ca*(V_s - E_Ca)            : amp
dm_Ca/dt = (m_Ca_inf - m_Ca)/tau_m_Ca            : 1
m_Ca_inf = 1/(1 + exp(0.205/mV*(-61.2*mV - V_s))): 1
tau_m_Ca = 30*ms -5*ms/(1 + exp(0.2/mV*(-65*mV - V_s))) : second
dh_Ca/dt = (h_Ca_inf - h_Ca)/tau_h_Ca            : 1
h_Ca_inf = 1/(1 + exp(-0.15/mV*(-75*mV - V_s)))  : 1

dm_K/dt = (m_K_inf - m_K)/tau_m_K                : 1
m_K_inf = 1/(1 + exp(0.1/mV*(-35*mV - V_s)))     : 1
tau_m_K = 2*ms + 55*ms/(1 + exp(-0.125/mV*(-54*mV - V_s))) : second

dm_A/dt = (m_A_inf - m_A)/tau_m_A                : 1
m_A_inf = 1/(1 + exp(0.2/mV*(-60*mV - V_s)))     : 1
dh_A/dt = (h_A_inf - h_A)/tau_h_A                : 1
h_A_inf = 1/(1 + exp(-0.18/mV*(-68*mV - V_s)))   : 1

dm_proc/dt = (m_proc_inf - m_proc)/tau_m_proc    : 1
m_proc_inf = 1/(1 + exp(0.2/mV*(-55*mV - V_s)))  : 1

# axonal compartment
dV_a/dt = (-g_La*(V_a - E_L) - g_Na*m_Na**3*h_Na*(V_a - E_Na)
           -g_Kd*m_Kd**4*(V_a - E_K) - g_E*(V_a - V_s))/C_a : volt

dm_Na/dt = (m_Na_inf - m_Na)/tau_m_Na            : 1
m_Na_inf = 1/(1 + exp(0.1/mV*(-42.5*mV - V_a)))  : 1
dh_Na/dt = (h_Na_inf - h_Na)/tau_h_Na            : 1
h_Na_inf = 1/(1 + exp(-0.13/mV*(-50*mV - V_a)))  : 1
tau_h_Na = 10*ms/(1 + exp(0.12/mV*(-77*mV - V_a))) : second

dm_Kd/dt = (m_Kd_inf - m_Kd)/tau_m_Kd            : 1
m_Kd_inf = 1/(1 + exp(0.2/mV*(-41*mV - V_a)))    : 1
tau_m_Kd = 12.2*ms + 10.5*ms/(1 + exp(-0.05/mV*(58*mV - V_a))) : second

# class-specific fixed maximal conductances
g_Ls   : siemens (constant)
g_A    : siemens (constant)
g_proc : siemens (constant)

# Adaptive conductances
g_Ca = G_Ca/2*(1 + tanh(z)) : siemens
g_K = G_K/2*(1 - tanh(z))   : siemens
I_diff = (I_target + I_Ca) : amp
dz/dt = tanh(I_diff/nA)/tau_z : 1
I_target : amp (constant)

# neuron class
label : integer (constant)
# neuron dynamics
neur_dyn : 1 (constant)
'''
ABPD, LP, PY = 0, 1, 2

circuit = NeuronGroup(3, eqs, method='euler',
                      threshold='m_Na > 0.5', refractory='m_Na > 0.5')

# class-specific constants
circuit.label = [ABPD, LP, PY]
circuit.I_target = [0.4, 0.3, 0.5]*nA
circuit.g_Ls = [30, 25, 15]*nS
circuit.g_A = [450, 100, 250]*nS
circuit.g_proc = [6, 8, 0]*nS

# Initial conditions
circuit.V_s = E_L
circuit.V_a = E_L
circuit.m_Ca = 'm_Ca_inf'
circuit.h_Ca = 'h_Ca_inf'
circuit.m_K = 'm_K_inf'
circuit.m_A = 'm_A_inf'
circuit.h_A = 'h_A_inf'
circuit.m_proc = 'm_proc_inf'
circuit.m_Na = 'm_Na_inf'
circuit.h_Na = 'h_Na_inf'
circuit.m_Kd = 'm_Kd_inf'
circuit.neur_dyn = 1

#circuit[1].run_regularly(f'V_s = -70*mV')
dyn_clamp_neuron,_ = attach_neuron(circuit, 1, 'V_s', 'I_total', dt=0.1*ms)
dyn_clamp_neuron.neur_dyn = 0
# %% [markdown]
# The definition of the synapses and the simulation protocol are identical to the simplified model:

# %%
eqs_fast = '''
g_fast : siemens (constant)
I_fast_post = g_fast*(V_s_post - E_syn)/(1+exp(s_fast*(V_fast-V_s_pre))) : amp (summed)
'''
fast_synapses = Synapses(circuit, circuit, model=eqs_fast)
fast_synapses.connect('label_pre != label_post and not (label_pre == PY and label_post == ABPD)')
fast_synapses.g_fast['label_pre == ABPD and label_post == LP'] = 0.015*uS
fast_synapses.g_fast['label_pre == ABPD and label_post == PY'] = 0.005*uS
fast_synapses.g_fast['label_pre == LP and label_post == ABPD'] = 0.01*uS
fast_synapses.g_fast['label_pre == LP and label_post == PY']   = 0.02*uS
fast_synapses.g_fast['label_pre == PY and label_post == LP']   = 0.005*uS

eqs_slow = '''
k_2 : 1/second (constant)
g_slow : siemens (constant)
I_slow_post = g_slow*m_slow*(V_s_post-E_syn) : amp (summed)
dm_slow/dt = k_1*(1-m_slow)/(1+exp(s_slow*(V_slow-V_s_pre))) - k_2*m_slow : 1 (clock-driven)
'''
slow_synapses = Synapses(circuit, circuit, model=eqs_slow, method='exact')
slow_synapses.connect('label_pre == ABPD and label_post != ABPD')
slow_synapses.g_slow['label_post == LP'] = 0.025*uS
slow_synapses.k_2['label_post == LP']    = 0.03/ms
slow_synapses.g_slow['label_post == PY'] = 0.015*uS
slow_synapses.k_2['label_post == PY']    = 0.008/ms


M = StateMonitor(circuit, 'V_s', record=True)
spikes = SpikeMonitor(circuit)

M.active = False
run(init_time, report='text')
M.active = True
run(observe_time, report='text')
M.active = False
run(adapt_time, report='text')
M.active = True
run(observe_time, report='text')
device = init_neuron_device(device, dt=0.1*ms, scalefactor_out=2.5)
device.build(directory='example_1_complex')

# %% [markdown]
# We use the same analysis/plotting procedure as in the simplified model:

# %%
spike_trains = spikes.spike_trains()

# Plot

import numpy as np

from plotly import tools
import plotly.graph_objs as go

from brian2.units import second, mV

def do_pyloric_net_plot(spike_trains, times, membrane_potential, varname,
                        init_time, observe_time, adapt_time):
    fig = tools.make_subplots(rows=7, cols=2, shared_xaxes=True, shared_yaxes=True,
                              start_cell='bottom-left', subplot_titles=['initial', 'adapted',
                                                                        '', '', '', '', '', ''],
                             specs=[[{}, {}],
                                    [{'rowspan': 2}, {'rowspan': 2}],
                                    [None, None],
                                    [{'rowspan': 2}, {'rowspan': 2}],
                                    [None, None],
                                    [{'rowspan': 2}, {'rowspan': 2}],
                                    [None, None]],
                             print_grid=False)
    
    traces = []
    before_adaptation = (times>=init_time) & (times < (init_time + observe_time))
    after_adapt_time = init_time + observe_time + adapt_time
    after_adaptation = (times>=after_adapt_time)
    for idx, (label, color) in enumerate(zip(['AB/PD', 'LP', 'PY'],
                                            ['#1f77b4', '#ff7f03', '#2ca02c'])):
        
        trace = go.Scattergl(x=(times[before_adaptation] - init_time) / second,
                           y=membrane_potential[idx][before_adaptation] / mV,
                           marker={'color': color},
                           showlegend=False, name=label)
        fig.append_trace(trace, 2+idx*2, 1)
        if spike_trains is not None:
            spike_times = spike_trains[idx]
            spike_times = spike_times[(spike_times >= init_time) & (spike_times < (init_time + observe_time))]
            spike_trace = go.Scattergl(x=(spike_times - init_time) / second,
                                     y=np.ones(len(spike_times))*(3-idx),
                                     marker={'symbol': 'line-ns', "line": {"width": 2, 'color': color}},
                                     mode='markers', showlegend=False, name=label)
            fig.append_trace(spike_trace, 1, 1)
        trace = go.Scattergl(x=(times[after_adaptation] - after_adapt_time) / second,
                           y=membrane_potential[idx][after_adaptation] / mV,
                           marker={'color': color},
                           name=label)
        fig.append_trace(trace, 2+idx*2, 2)
        if spike_trains is not None:
            spike_times = spike_trains[idx]
            spike_times = spike_times[spike_times >= after_adapt_time]
            spike_trace = go.Scattergl(x=(spike_times - after_adapt_time) / second,
                             y=np.ones(len(spike_times))*(3-idx),
                             marker={'symbol': 'line-ns', "line": {"width": 2, 'color': color}},
                             mode='markers', showlegend=False, name=label)
            fig.append_trace(spike_trace, 1, 2)
    fig['layout'].update(xaxis1={'range': (0, observe_time/second),
                                 'title': 'time (in s)',
                                 'zeroline': False},
                         xaxis2={'range': (0, observe_time/second),
                                 'title': 'time (in s)',
                                 'zeroline': False},
                         yaxis3={'title': 'v (in mV)'},
                         yaxis1={'showline': False,
                                 'showticklabels': False})
    fig.show()

do_pyloric_net_plot(spike_trains, M.t, M.V_s, 'V_s',
                    init_time, observe_time, adapt_time);

print('f')


