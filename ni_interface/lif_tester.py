from brian2 import *

def binarize_spikes(spike_times, dt=1e-4, duration=1):
    #spike_times is a list of spike times in seconds
    #dt is the timestep in seconds
    #duration is the duration of the signal in seconds
    #return the binarized spikes
    num_steps = int(duration/dt)
    spikes = np.zeros(num_steps)
    spike_indices = (np.array(spike_times)/dt).astype(int)
    spikes[spike_indices] = 1
    return spikes


def lif_model(I, dt=1e-4):
    #simple LIFF neuron model
    #I is the input current in pA
    #dt is the timestep in seconds
    #return the voltage in mV
    I_curr = TimedArray(I*pA, dt=dt*second)
    eqs = '''
    dv/dt = (gL*(EL-v)+I_curr(t))/C : volt (unless refractory)
    tau : second
    EL : volt
    gL : siemens
    C: farad
    '''
    G = NeuronGroup(1, eqs, threshold='v>-50*mV', reset='v=-60*mV', refractory=1*ms, method='exact')
    G.C = 50*pF
    G.EL = -70*mV
    G.gL = 1/(0.2*Gohm)
    G.v = -70*mV
    M = StateMonitor(G, 'v', record=True)
    spk = SpikeMonitor(G)
    run(len(I)*(dt)*second)
    return M.v[0]/mV, binarize_spikes(spk.t/second, dt=dt, duration=len(I)*dt)