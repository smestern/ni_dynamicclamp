import ot
import torch
import snntorch as snn
from snntorch._neurons.neurons import _SpikeTorchConv, _SpikeTensor

#make a new neuron model with biologically plausible parameters
class sLIFin(snn.SpikingNeuron):
    #here we want parameters for the LIF neuron, G, C, EL, tau
    #we also want the parameters for the synapse, alpha and beta
    #we also want the threshold
    def __init__(self, R=0.1, EL=-70, tau=0.1, ge=3, taue=0.005, threshold=-40, dt=1e-5, init_hidden=True, **kwargs):
        super(snn.SpikingNeuron, self).__init__(
                
            )
        #in gigohms convert to ohms
        self.R = R * 1e9
        self.EL = EL * 1e-3
        #in seconds leave as is
        self.tau = tau 
        #in nanosiemens convert to siemens
        self.ge = ge * 1e-9
        #in seconds leave as is
        self.taue = taue
        #in mV convert to V
        self.threshold = threshold * 1e-3
        #in seconds leave as is
        self.dt = dt
        self.init_hidden = init_hidden
      
        for param in [self.R, self.EL, self.tau, self.ge, self.taue, self.threshold]:
            if not isinstance(param, torch.Tensor):
                param = torch.as_tensor(param)
            
            self._param_register_buffer(param, True)

        if self.init_hidden:
            self.syn, self.mem = self.init_sLIF()

    def init_sLIF(self):
        spk = _SpikeTensor(init_flag=False)
        syn = _SpikeTensor(init_flag=False)
        mem = _SpikeTensor(init_flag=False)
        return syn, mem

    def forward(self, input_, syn=False, mem=False):

        if hasattr(syn, "init_flag") or hasattr(
            mem, "init_flag"
        ):  # only triggered on first-pass
            syn, mem = _SpikeTorchConv(syn, mem, input_=input_)
        elif mem is False and hasattr(
            self.mem, "init_flag"
        ):  # init_hidden case
            self.syn, self.mem = _SpikeTorchConv(
                self.syn, self.mem, input_=input_
            )

        if not self.init_hidden:
            self.reset = self.mem_reset(mem)
            syn, mem = self._build_state_function(input_, syn, mem)

            if self.state_quant:
                syn = self.state_quant(syn)
                mem = self.state_quant(mem)

            if self.inhibition:
                spk = self.fire_inhibition(mem.size(0), mem)
            else:
                spk = self.fire(mem)

            return spk, syn, mem

        # intended for truncated-BPTT where instance variables are
        # hidden states
        if self.init_hidden:
            self._synaptic_forward_cases(mem, syn)
            self.reset = self.mem_reset(self.mem)
            self.syn, self.mem = self._build_state_function_hidden(input_)

            if self.state_quant:
                self.syn = self.state_quant(self.syn)
                self.mem = self.state_quant(self.mem)

            if self.inhibition:
                self.spk = self.fire_inhibition(self.mem.size(0), self.mem)
            else:
                self.spk = self.fire(self.mem)

            if self.output:
                return self.spk, self.syn, self.mem
            else:
                return self.spk



    def _base_state_function(self, input_, syn, mem):
        base_fn_syn = self.alpha.clamp(0, 1) * syn + input_
        base_fn_mem = self.beta.clamp(0, 1) * mem + base_fn_syn
        return base_fn_syn, base_fn_mem

    def _base_state_reset_zero(self, input_, syn, mem):
        base_fn_syn = self.alpha.clamp(0, 1) * syn + input_
        base_fn_mem = self.beta.clamp(0, 1) * mem + base_fn_syn
        return 0, base_fn_mem

    def _build_state_function(self, input_, syn, mem):
        if self.reset_mechanism_val == 0:  # reset by subtraction
            state_fn = tuple(
                map(
                    lambda x, y: x - y,
                    self._base_state_function(input_, syn, mem),
                    (0, self.reset * self.threshold),
                )
            )
        elif self.reset_mechanism_val == 1:  # reset to zero
            state_fn = tuple(
                map(
                    lambda x, y: x - self.reset * y,
                    self._base_state_function(input_, syn, mem),
                    self._base_state_reset_zero(input_, syn, mem),
                )
            )
        elif self.reset_mechanism_val == 2:  # no reset, pure integration
            state_fn = self._base_state_function(input_, syn, mem)
        return state_fn

    def _base_state_function_hidden(self, input_):
        base_fn_syn = self.alpha.clamp(0, 1) * self.syn + input_
        base_fn_mem = self.beta.clamp(0, 1) * self.mem + base_fn_syn
        return base_fn_syn, base_fn_mem

    def _base_state_reset_zero_hidden(self, input_):
        base_fn_syn = self.alpha.clamp(0, 1) * self.syn + input_
        base_fn_mem = self.beta.clamp(0, 1) * self.mem + base_fn_syn
        return 0, base_fn_mem

    def _build_state_function_hidden(self, input_):
        if self.reset_mechanism_val == 0:  # reset by subtraction
            state_fn = tuple(
                map(
                    lambda x, y: x - y,
                    self._base_state_function_hidden(input_),
                    (0, self.reset * self.threshold),
                )
            )
        elif self.reset_mechanism_val == 1:  # reset to zero
            state_fn = tuple(
                map(
                    lambda x, y: x - self.reset * y,
                    self._base_state_function_hidden(input_),
                    self._base_state_reset_zero_hidden(input_),
                )
            )
        elif self.reset_mechanism_val == 2:  # no reset, pure integration
            state_fn = self._base_state_function_hidden(input_)
        return state_fn

    def _param_register_buffer(self, alpha, learn_alpha):
        if not isinstance(alpha, torch.Tensor):
            alpha = torch.as_tensor(alpha)
        if learn_alpha:
            self.alpha = torch.nn.Parameter(alpha)
        else:
            self.register_buffer("alpha", alpha)

    def _synaptic_forward_cases(self, mem, syn):
        if mem is not False or syn is not False:
            raise TypeError(
                "When `init_hidden=True`, Synaptic expects 1 input argument."
            )
        


class EMD():
    def __init__(self, dt=1e-5, duration=1, bin_size=-1, batch_axis=-1):
        self.dt = dt
        self.duration = duration
        self.bin_size = bin_size
        self.batch_axis = batch_axis
    def __call__(self, x, y):
        if self.bin_size > 0:
            x = torch.stack([self.bin(x[:, i]) for i in range(x.shape[1])], dim=1)
            y = torch.stack([self.bin(y[:, i]) for i in range(y.shape[1])], dim=1)
            _coords = torch.arange(start=0.0,end=self.duration, step=self.bin_size, device=x.device)
        else:
            _coords = torch.arange(x.shape[0], device=x.device)
        losses = []
        for i in range(x.shape[1]):
            losses.append(ot.wasserstein_1d(_coords, _coords, x[:, i], y[:, i],  p=2))
        error = torch.mean(torch.abs(torch.stack(losses, dim=0)))
        return  error

    def bin(self, x):
        #bin the spikes, in a window
        #x is the spikes, already binned AT dt resolution
        #dt is the timestep in seconds
        #duration is the duration of the signal in seconds
        #return the binned spikes
        average_window = int(self.bin_size/self.dt)
        num_steps = int(self.duration/self.bin_size)
        x = x.reshape(num_steps, average_window)
        x = torch.mean(x, dim=1)
        return x