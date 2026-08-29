import ot
import math
import torch
import snntorch as snn

#make a new neuron model with biologically plausible parameters
class sLIFin(snn.SpikingNeuron):
    """Trainable synaptic-LIF proxy for a real neuron's *spiking*.

    This model is optimised to predict *when* a real cell spikes, not to
    reproduce its absolute membrane trace. Dynamics are therefore kept in a
    normalised, O(1) range (``mem`` starts at 0) so a downstream squashing
    non-linearity (e.g. ``tanh``) and the surrogate spike gradient both stay
    in their sensitive regions during training.

    snntorch 1.0 removed ``_SpikeTensor`` / ``_SpikeTorchConv`` and the
    ``init_flag`` hidden-state protocol. State is now held in lazily-sized
    buffers (``syn`` / ``mem``) that are resized to match the input on the
    first forward pass, mirroring :class:`snntorch.Synaptic`.

    The physical time constants seed the learnable per-step decay factors
    ``alpha = exp(-dt/taue)`` (synaptic) and ``beta = exp(-dt/tau)``
    (membrane); ``threshold`` is learnable so the proxy can shift its spike
    boundary to match the target cell. ``R`` / ``EL`` / ``ge`` are retained
    as reference attributes but deliberately kept out of the (normalised)
    membrane update.

    Alongside the spike prediction the neuron emits ``psc`` — a
    post-synaptic-current trace filtered from its *output* spikes
    (``psc = alpha*psc + spk``, reusing the synaptic decay). This is the
    STEP 2.5 -> STEP 3 signal that drives the next in-silico layer. Pass the
    real neuron's measured spikes via ``real_spk`` to run the trace
    straight-through: the forward ``psc`` is built from the real spikes
    while gradients still flow back through the proxy's predicted spikes.
    """

    def __init__(self, R=0.1, EL=-70, tau=0.1, ge=3, taue=0.005,
                 threshold=1.0, dt=1e-5, init_hidden=False,
                 reset_mechanism="subtract", learn_decay=True,
                 learn_threshold=True, **kwargs):
        super().__init__(
            threshold=threshold,
            init_hidden=init_hidden,
            reset_mechanism=reset_mechanism,
            learn_threshold=learn_threshold,
            **kwargs,
        )
        # kept for reference / bookkeeping, not used in the normalised update
        self.R = R
        self.EL = EL
        self.ge = ge
        self.tau = tau
        self.taue = taue
        self.dt = dt

        # per-step decay factors derived from the time constants
        alpha = torch.as_tensor(math.exp(-dt / taue), dtype=torch.float)
        beta = torch.as_tensor(math.exp(-dt / tau), dtype=torch.float)
        if learn_decay:
            self.alpha = torch.nn.Parameter(alpha)
            self.beta = torch.nn.Parameter(beta)
        else:
            self.register_buffer("alpha", alpha)
            self.register_buffer("beta", beta)

        self._init_mem()

    def _init_mem(self):
        self.register_buffer("syn", torch.zeros(0), False)
        self.register_buffer("mem", torch.zeros(0), False)
        self.register_buffer("psc", torch.zeros(0), False)

    def reset_mem(self):
        self.syn = torch.zeros_like(self.syn, device=self.syn.device)
        self.mem = torch.zeros_like(self.mem, device=self.mem.device)
        self.psc = torch.zeros_like(self.psc, device=self.psc.device)
        return self.syn, self.mem, self.psc

    def init_synaptic(self):
        """Backwards-compatible alias for :meth:`reset_mem`."""
        return self.reset_mem()

    def forward(self, input_, syn=None, mem=None, psc=None, real_spk=None):
        if syn is not None:
            self.syn = syn
        if mem is not None:
            self.mem = mem
        if psc is not None:
            self.psc = psc

        if self.init_hidden and (
            syn is not None or mem is not None or psc is not None
        ):
            raise TypeError(
                "`mem`, `syn` or `psc` should not be passed as an argument "
                "while `init_hidden=True`"
            )

        if self.syn.shape != input_.shape:
            self.syn = torch.zeros_like(input_, device=self.syn.device)
        if self.mem.shape != input_.shape:
            self.mem = torch.zeros_like(input_, device=self.mem.device)
        if self.psc.shape != input_.shape:
            self.psc = torch.zeros_like(input_, device=self.psc.device)

        self.reset = self.mem_reset(self.mem)
        self.syn, self.mem = self._build_state_function(input_)

        if self.state_quant:
            self.syn = self.state_quant(self.syn)
            self.mem = self.state_quant(self.mem)

        if self.inhibition:
            spk = self.fire_inhibition(self.mem.size(0), self.mem)
        else:
            spk = self.fire(self.mem)

        # STEP 2.5 -> STEP 3: post-synaptic trace from the OUTPUT spikes.
        # With `real_spk` supplied, run straight-through -- the forward trace
        # follows the real cell's measured spikes while gradients keep
        # flowing through the proxy's predicted spikes.
        if real_spk is not None:
            spk = spk + (real_spk - spk).detach()
        self.psc = self.alpha.clamp(0, 1) * self.psc + spk

        if self.output:
            return spk, self.syn, self.mem, self.psc
        elif self.init_hidden:
            return spk
        else:
            return spk, self.syn, self.mem, self.psc

    def _base_state_function(self, input_):
        base_fn_syn = self.alpha.clamp(0, 1) * self.syn + input_
        base_fn_mem = self.beta.clamp(0, 1) * self.mem + base_fn_syn
        return base_fn_syn, base_fn_mem

    def _build_state_function(self, input_):
        syn, mem = self._base_state_function(input_)
        if self.reset_mechanism_val == 0:  # reset by subtraction
            mem = mem - self.reset * self.threshold
        elif self.reset_mechanism_val == 1:  # reset to zero
            mem = mem - self.reset * mem
        return syn, mem



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


class VanRossum:
    """Van Rossum spike-train distance (differentiable, spike-only).

    Each spike train is convolved with a causal exponential kernel
    ``h(t) = exp(-t/tau)`` and the distance is the L2 difference of the
    filtered traces, ``D^2 = (dt/tau) * sum_t (f_x - f_y)^2``. Unlike
    :class:`EMD` it needs no binning and stays fully differentiable, so it
    is the natural objective when fitting the proxy purely on spike timing.
    """

    def __init__(self, tau=5e-3, dt=1e-4):
        self.tau = tau
        self.dt = dt

    def __call__(self, x, y):
        # x, y: (T, B) spike trains at dt resolution.
        decay = math.exp(-self.dt / self.tau)
        fx = self._filter(x, decay)
        fy = self._filter(y, decay)
        d2 = ((fx - fy) ** 2).sum(dim=0) * (self.dt / self.tau)
        return d2.mean()

    def _filter(self, s, decay):
        out = []
        trace = torch.zeros_like(s[0])
        for t in range(s.shape[0]):
            trace = decay * trace + s[t]
            out.append(trace)
        return torch.stack(out, dim=0)


def fit_proxy(proxy, real_neuron, input_gen, *, n_iter=600, lr=5e-2,
              spike_loss="vanrossum", mem_weight=1.0, tau=5e-3, dt=1e-4,
              duration=0.25, bin_size=2e-3, device=None, use_scheduler=True,
              callback=None, verbose=True):
    """Fit a differentiable spiking ``proxy`` to a black-box ``real_neuron``.

    Parameters
    ----------
    proxy : torch.nn.Module
        ``proxy(cur)`` -> ``(spk, mem)`` each shaped ``(T, B)``, differentiable
        in ``proxy.parameters()``.
    real_neuron : callable
        ``real_neuron(cur)`` -> ``(spk, mem)`` shaped ``(T, B)``. Treated as the
        (detached) target; ``cur`` is the ``(B, T)`` current batch.
    input_gen : callable or torch.Tensor
        ``input_gen()`` -> current batch ``(B, T)``, or a fixed tensor reused
        every iteration.
    spike_loss : {"vanrossum", "emd"} or callable
        Spike-train objective. ``"vanrossum"`` is the spike-only default.
    mem_weight : float
        Weight on the membrane MSE term; set to ``0`` for a spike-only fit.
    callback : callable(it, loss, tensors) or None
        Optional hook (e.g. live plotting); ``tensors`` holds the current
        ``cur`` / ``proxy_spk`` / ``proxy_mem`` / ``real_spk`` / ``real_mem``.

    Returns
    -------
    list[float]
        Per-iteration loss history.
    """
    if device is None:
        device = next(proxy.parameters()).device

    if spike_loss == "vanrossum":
        spk_loss_fn = VanRossum(tau=tau, dt=dt)
    elif spike_loss == "emd":
        spk_loss_fn = EMD(dt=dt, duration=duration, bin_size=bin_size)
    elif callable(spike_loss):
        spk_loss_fn = spike_loss
    else:
        raise ValueError(
            "spike_loss must be 'vanrossum', 'emd', or a callable; got "
            f"{spike_loss!r}"
        )
    mem_loss_fn = torch.nn.MSELoss()

    opt = torch.optim.Adam(proxy.parameters(), lr=lr)
    sched = (
        torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=n_iter, eta_min=1e-9)
        if use_scheduler else None
    )

    hist = []
    for it in range(n_iter):
        cur = input_gen() if callable(input_gen) else input_gen
        cur = cur.to(device)

        proxy_spk, proxy_mem = proxy(cur)
        with torch.no_grad():
            real_spk, real_mem = real_neuron(cur)
        real_spk = real_spk.to(device)
        real_mem = real_mem.to(device)

        loss = torch.mean(spk_loss_fn(proxy_spk, real_spk))
        if mem_weight:
            loss = loss + mem_weight * mem_loss_fn(proxy_mem, real_mem)

        opt.zero_grad()
        loss.backward()
        opt.step()
        if sched is not None:
            sched.step()

        hist.append(float(loss.detach().cpu()))
        if callback is not None:
            callback(it, hist[-1], dict(
                cur=cur, proxy_spk=proxy_spk, proxy_mem=proxy_mem,
                real_spk=real_spk, real_mem=real_mem,
            ))
        if verbose:
            print(f"[fit_proxy] iter {it:4d}  loss={hist[-1]:.4f}")

    return hist