#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 14:45:48 2018

@author: tdesbordes
"""

from brian2 import *
from scipy.io import savemat, loadmat
from ni_interface.ni_brian2 import *
seed(4323)
defaultclock.dt = 0.01*ms
set_device('cpp_standalone', build_on_run=True)
N=10000
K=100

T=10.0*second
dt=0.05*ms

we = 3.3276023360e-9*siemens
wi = 8e-8*siemens
taum = 20e-3*second
vreset = -70e-3*volt
vth = -50e-3*volt
taur = 5e-3*second
Ee = 0*volt
Ei = -80e-3*volt
taue = 5e-3*second
taui = 5e-3*second
EL = -65e-3*volt
Cm = 200e-12*farad

defaultclock.dt = dt

gL = Cm/taum

vm_mean = -0.0680901*volt
vm_sigma = 0.0057296*volt
ge_mean = 1.14244e-08*siemens
ge_sigma = 4.30869e-09*siemens
gi_mean = 7.95233e-08*siemens
gi_sigma = 5.80938e-08*siemens

# Membrane Equation
eqs = Equations('''
dv/dt = ( gL*(EL-v) + ge*(Ee-v) + gi*(Ei-v) ) * (1./Cm) : volt (unless refractory)
d_I = clip(( ge*(Ee-v) + gi*(Ei-v)), -250*pA, 250*pA) : amp
dge/dt = -ge*(1./taue) : siemens
dgi/dt = -gi*(1./taui) : siemens
vth : volt
''')

# Set up the network
P = NeuronGroup(N, model=eqs, threshold='v>vth', reset='v=vreset', refractory=taur, method='euler')
Ne = int(N*8/10)
Ni = N-Ne
Pe = P[:Ne]
Pi = P[Ne:]
neuron = P[:1]
P.vth = vth
neuron.run_regularly(f'v = step_clamp(t, d_I)', dt=defaultclock.dt)
neuron.vth = -0*mV
taudel = 300e-6*second

Ce = Synapses( Pe, P, on_pre='ge+=we', delay=taudel )
Ci = Synapses( Pi, P, on_pre='gi+=wi', delay=taudel )

# data = loadmat('data/ij.mat')
# exc_i = int32(data['exc_i']).flatten()
# exc_j = int32(data['exc_j']).flatten()
# inh_i = int32(data['inh_i']).flatten()
# inh_j = int32(data['inh_j']).flatten()

Ce.connect( p=0.01 )
Ci.connect( p=0.01 )

# Initialization = gaussian sampling with mean and sigma given by param file
P.v = np.random.normal( vm_mean, vm_sigma, len(P) ) * volt
P.ge = np.random.normal( ge_mean, ge_sigma, len(P) ) * siemens
P.gi = np.random.normal( gi_mean, gi_sigma, len(P) ) * siemens

# Record the number of spikes
M = SpikeMonitor(P)
Mv = StateMonitor(P, 'v', record=[0])
Mr = PopulationRateMonitor(P)
device = init_neuron_device(device=device, scalefactor_out=2.5)
# Run simulation
run( T, report='stdout' )

si = M.i[:] + 1
st = M.t/ms
st = st/1000

mdic = {"si":si, "st":st, "N":N, "T":T, "Ne":double(Ne), "dt":dt/second}
savemat("brian2.mat",mdic)
