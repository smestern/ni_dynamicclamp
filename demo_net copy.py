from brian2 import *
from ni_interface.ni_brian2 import *
import os
import time
seed(4323)
defaultclock.dt = 0.1*ms
set_device('cpp_standalone', build_on_run=True)
# state 1 (network-bursting): a = 400 * nS; we = 0.1 * nS; wi = 100 * nS
# state 2 (network-reduced-inhibition): a = 400 * nS; we = 0.1 * nS; wi = 10 * nS
# state 3 (network-tonic): a = 4 * nS; we = 10 * nS; wi = 100 * nS
import time
start_scope()
start_time = time.time()
seed(4323)
# network parameters
N = 1e3
p_ei = 0.02
p_ie = 0.04 # e->i & i->e connection probability

# neuron parameters
C = 21.64577686 * pF
taum = 26.45276208*ms
gL = C/taum #30 * nS


DeltaT = 12.0032413*mV
pr = 1
EL = -67.91702674 * mV
VT = -55.19868804 * mV
VR = -71.40471178 * mV
LTS = {'tauw':(92.07978488 + np.random.normal(loc=10, scale=5, size=500))*ms, 'a': np.random.normal(loc=0.090873665, scale=0.0001, size=500)*nS, 
       'b': np.abs(np.random.normal(loc=38.65665726, scale=3, size=500))*pA } # low-threshold spiking(?)
FS = {'tauw': 13.5686673*ms, 'a': 0.*nS, 'b': 3*pA }
#{'tauw': 144*ms, 'a': 0.1*nS, 'b': 0*nA } # fast spiking(?)
#{'N': 1, '_run_time': 2, 'Cm': 20.21738406, 'EL': -67.44385, 'VT': -45.55570857, 'VR': -53.08817532,
#'taum': 20.54444831, 'tauw': 13.5686673, 'a': 0.08555308306564564, 'b': 0.17, 'DeltaT': 8.543155448, 'dt': 0.1}
# synapse parameters
Ee = 0 * mvolt
Ei = -80 * mvolt
taue = 18.24179238* msecond
tauCRH = 250.7198441 * msecond
taun = 1 * msecond
taui = 11.96379858 * msecond
wCRH = 0.005620078 * nS
we = 5.371492379 * nS
wn = .000001 * nS
wi = 12.54548072 * nS
taup = 4 *second
taubr = 8 * second



#sin_ = np.sin(0.005 * np.arange(0, 2500)) * 10 * pA
#sin_ = (chirp(np.arange(0, 10, 0.001),f0=0,t1=10,f1=3) * 1)* pA
#input_curr = TimedArray(sin_, dt=1*ms)
# 
#-pr*(1./taup) 
# Membrane Equation + y*(Ee-v)
eqs = Equations('''
dv/dt = ( gL*(EL-v) + gL*DeltaT*exp((v - VT)/DeltaT) + ge*(Ee-v) + gi*(Ei-v) - w) * (1./C) : volt (unless refractory)
d_I = clip(( ge*(Ee-v) + gi*(Ei-v)), -250*pA, 250*pA) : amp 
dw/dt = ( a*(v - EL) - w ) / tauw : amp 
dgi/dt = -gi*(1./taui) : siemens
dy/dt = -y*(1./taue) : siemens
dbd/dt = -bd*(1./taubr): amp
dge/dt = (y-ge)/taue : siemens
dpr/dt = (1./taup) * pr * (1 - (pr/1)): 1
p_w = clip(pr, 0.01,1): 1
br = b + bd : amp
tauw : second
taue : second
taui : second
a : siemens
I : amp
b: amp
Vcut: volt
VT : volt
refrac : second
input_1 : 1
input_2 : 1
input_3 : 1
VR : volt
v_gap : volt''')


# build network
P = NeuronGroup( N, model=eqs, threshold='v>Vcut', reset='v=VR; w+=br', refractory='refrac', method='euler')
P.pr = 1
P.VT = VT
P.VR = VR
P.Vcut = (VT + 5 * DeltaT )
P.taui = taui
#print(P.Vcut)
CRH = P[:int(floor(0.5*N))]; GABA = P[int(floor(0.5*N)):]
CRH.tauw = LTS['tauw']; CRH.a = LTS['a']; CRH.b = LTS['b']
GABA.tauw = FS['tauw']; GABA.a = FS['a']; GABA.b = FS['b']
GABA.taue = tauCRH; CRH.taue = taue
GABA.refrac = 1*ms
CRH.refrac = (1 + np.abs(np.random.normal(loc=1.5, scale=1, size=500))) *ms

neuron = P[:1]
neuron.run_regularly(f'v = step_clamp(t, d_I*int(t/second>3))', dt=defaultclock.dt)
neuron.Vcut = 0*mV


add_neuron = P[1:30]
add_neuron.run_regularly(f'v = v_gap + 5*mV*randn()', dt=defaultclock.dt)
gap_syn = Synapses(neuron, add_neuron, '''v_gap_post = v_pre : volt (summed)''', method='euler')
gap_syn.connect()
#neuron.r
#CRH[:1].Rrsetter.up
#presynaptic indices
GABA_TO_0 = np.array([  5,  22, 161, 193, 222, 249, 278, 301, 400, 413, 439, 456, 458,
       462])




# connect
EI = Synapses( CRH, GABA, on_pre='y=clip((y + wCRH), 0*nS, 30*nS)' ) #*(rand()<p_w)
EI.connect( True, p=p_ei )
IE = Synapses( GABA, CRH, on_pre='gi+=wi*(rand()<p_w)' ) #*(rand()<p_w)
IE.connect( True, p=p_ie )
P.input_1 = 1
P.input_2 = 0
P.input_3 = 0
# poisson input
input = []
sub = []
input=PoissonInput(CRH, 'ge', 1, 50.53785123*Hz, weight='we*input_1')
input2 = PoissonInput( CRH, 'ge', 1, 50*Hz, weight='we*input_2')
input3 = PoissonInput( CRH, 'ge', 1, 75*Hz, weight='we*input_3')
#input2 = PoissonInput( CRH, 'gi', 1, 1*Hz, weight=wi )
# init
P.v = ((randn(len(P)) * (2* mV))) + EL
P.ge = (randn(len(P)) * 2 + 5) * nS
P.gi = ((randn(len(P)) * 2 + 5) + 0) * nS


# run simulation

print("=== Net Sim Start ===")

Mv = StateMonitor(CRH, ['v', 'd_I'], record=[0])
#record the neurons 
Mv2 = StateMonitor(GABA, 'v', record=GABA_TO_0)
# Record the value of v when the threshold is crossed
M_crossings = SpikeMonitor(CRH, 'v', record=[0])
#run(2*second, report='text')

device = init_neuron_device(device=device, dt=defaultclock.dt, scalefactor_out=2)



time_start = time.time()
run(600*second, report='text')
#print(time.time() - time_start)
subplot(211)
plot(Mv.t/ms, Mv.v[0]/mV, label='v', c='k', alpha=0.5)
spike_dict = M_crossings.spike_trains()
scatter(spike_dict[0]/ms, np.full(len(spike_dict[0]), 10), c='r', marker='x')
twinx()
plot(Mv.t/ms, Mv.d_I[0]/pA, label='v2', c='r', alpha=0.5)
subplot(212)
#plot the gaba responders
for i,_ in enumerate(GABA_TO_0):
    plot(Mv2.t/ms, Mv2.v[i]/mV, label='v', c='k', alpha=0.5)
show()
print('p')
