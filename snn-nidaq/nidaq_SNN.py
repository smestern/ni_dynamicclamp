from snntorch._neurons.neurons import *
import nidaqmx
#define our custom LIF CLASS representing the real neuron
class niDAQ(LIF):
    def __init__(
        self,
        beta,
        threshold=1.0,
        spike_grad=None,
        init_hidden=False,
        inhibition=False,
        learn_beta=False,
        learn_threshold=False,
        reset_mechanism="subtract",
        state_quant=False,
        output=False,
    ):
        super(niDAQ, self).__init__(
            beta,
            threshold,
            spike_grad,
            init_hidden,
            inhibition,
            learn_beta,
            learn_threshold,
            reset_mechanism,
            state_quant,
            output,
        )

        if self.init_hidden:
            self.mem = self.init_leaky()
            self.state_fn = self._build_state_function_hidden
        else:
            self.state_fn = self._build_state_function

        self.task_write = nidaqmx.Task()
        self.task_write.ao_channels.add_ao_voltage_chan("Dev2/ao0")
        self.task_write.start()
        self.task_read = nidaqmx.Task()
        self.task_read.ai_channels.add_ai_voltage_chan("Dev2/ai0")
        self.task_read.start()

    def forward(self, input_, mem=False):

        if hasattr(mem, "init_flag"):  # only triggered on first-pass
            mem = _SpikeTorchConv(mem, input_=input_)
        elif mem is False and hasattr(self.mem, "init_flag"):  # init_hidden case
            self.mem = _SpikeTorchConv(self.mem, input_=input_)

        # TO-DO: alternatively, we could do torch.exp(-1 / self.beta.clamp_min(0)),
        # giving actual time constants instead of values in [0, 1] as initial beta
        # beta = self.beta.clamp(0, 1)

        if not self.init_hidden:
            self.reset = self.mem_reset(mem)
            mem = self.state_fn(input_, mem)

            if self.state_quant:
                mem = self.state_quant(mem)

            if self.inhibition:
                spk = self.fire_inhibition(mem.size(0), mem)  # batch_size
            else:
                spk = self.fire(mem)

            return spk, mem

        # intended for truncated-BPTT where instance variables are hidden states
        if self.init_hidden:
            self._niDAQ_forward_cases(mem)
            self.reset = self.mem_reset(self.mem)
            self.mem = self.state_fn(input_)

            if self.state_quant:
                self.mem = self.state_quant(self.mem)

            if self.inhibition:
                self.spk = self.fire_inhibition(self.mem.size(0), self.mem)
            else:
                self.spk = self.fire(self.mem)

            if self.output:  # read-out layer returns output+states
                return self.spk, self.mem
            else:  # hidden layer e.g., in nn.Sequential, only returns output
                return self.spk


    def internal_nidaq(self, input_):
        self.task_write.write(input_.detach().numpy()[:,0])
        return self.task_read.read()

    def _base_state_function(self, input_, mem):
        out = self.internal_nidaq(input_)
        base_fn = self.beta.clamp(0, 1) * mem + input_
        return base_fn

    def _build_state_function(self, input_, mem):
        if self.reset_mechanism_val == 0:  # reset by subtraction
            state_fn = (
                self._base_state_function(input_, mem - self.reset * self.threshold)
            )
        elif self.reset_mechanism_val == 1:  # reset to zero
            state_fn = self._base_state_function(
                input_, mem
            ) - self.reset * self._base_state_function(input_, mem)
        elif self.reset_mechanism_val == 2:  # no reset, pure integration
            state_fn = self._base_state_function(input_, mem)
        return state_fn

    def _base_state_function_hidden(self, input_):
        base_fn = self.beta.clamp(0, 1) * self.mem + input_
        return base_fn

    def _build_state_function_hidden(self, input_):
        if self.reset_mechanism_val == 0:  # reset by subtraction
            state_fn = (
                self._base_state_function_hidden(input_) - self.reset * self.threshold
            )
        elif self.reset_mechanism_val == 1:  # reset to zero
            state_fn = self._base_state_function_hidden(
                input_
            ) - self.reset * self._base_state_function_hidden(input_)
        elif self.reset_mechanism_val == 2:  # no reset, pure integration
            state_fn = self._base_state_function_hidden(input_)
        return state_fn

    def _niDAQ_forward_cases(self, mem):
        if mem is not False:
            raise TypeError("When `init_hidden=True`, niDAQ expects 1 input argument.")

    @classmethod
    def detach_hidden(cls):
        """Returns the hidden states, detached from the current graph.
        Intended for use in truncated backpropagation through time where hidden state variables are instance variables."""

        for layer in range(len(cls.instances)):
            if isinstance(cls.instances[layer], niDAQ):
                cls.instances[layer].mem.detach_()

    @classmethod
    def reset_hidden(cls):
        """Used to clear hidden state variables to zero.
        Intended for use where hidden state variables are instance variables.
        Assumes hidden states have a batch dimension already."""
        for layer in range(len(cls.instances)):
            if isinstance(cls.instances[layer], niDAQ):
                cls.instances[layer].mem = _SpikeTensor(init_flag=False)