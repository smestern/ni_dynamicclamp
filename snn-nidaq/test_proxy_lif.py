import torch, torch.nn as nn
#import snntorch as snn
import sys
import os
import argparse
import faulthandler

faulthandler.enable() #to debug seg faults and timeouts
batch_size = 4
data_path='./data/mnist'
sys.path.append('/home/smestern/Dropbox/ni_dynamicclamp/ni_interface/')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
#import ni_generic as ni
from lif_tester import lif_model
import utils_snn as utils

import matplotlib.pyplot as plt
import numpy as np 

# Optional real-DAQ path: same proxy_net but swap the offline lif_model
# call for a NiDAQClampLayer forward. Use --use-daq on the command line.
parser = argparse.ArgumentParser()
parser.add_argument("--use-daq", action="store_true",
                    help="Drive the real NI card via NiDAQClampLayer instead "
                         "of the offline lif_model proxy.")
parser.add_argument("--use-dummy", action="store_true",
                    help="With --use-daq, force the ni_torch_dummy backend "
                         "(no hardware required).")
_args, _ = parser.parse_known_args()

_daq_layer = None
if _args.use_daq:
    from ni_torch_layer import NiDAQClampLayer
    _daq_layer = NiDAQClampLayer(
        dt_ms=0.1, grad_mode="detach", use_dummy=_args.use_dummy,
    )
pA = 1e-12
device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")


#make a really generic dense network of 3 layers

dt = 0.0001
time_steps = int(0.25 * (1/dt))
num_steps = time_steps
beta = 0.8 #membrane potential decay rate in seconds 
alpha = 0.05 #synaptic decay rate in seconds
# Network Architecture
class proxy_net(nn.Module):
    def __init__(self) -> None:
        super().__init__()
        self.scale_fator = nn.Parameter(torch.tensor(1.0))
        self.lif1 = utils.sLIFin(0.1, 30, -50, 0.0002, 0.1, 0.1, -20)
        self.tanh = nn.Tanh()

    def forward(self, x):
        """Run the proxy over a current batch ``(B, T)``.

        Returns ``(spikes, mem)`` each shaped ``(T, B)``.
        """
        syn1, mem1, psc1 = self.lif1.reset_mem()
        spikes_rec = []
        mem_rec = []
        for step in range(x.shape[1]):
            cur1 = x[:, step]
            spk1, syn1, mem1, psc1 = self.lif1(cur1 * self.scale_fator, syn1, mem1, psc1)
            mem1 = self.tanh(mem1)
            spikes_rec.append(spk1)
            mem_rec.append(mem1)
        return torch.stack(spikes_rec, dim=0), torch.stack(mem_rec, dim=0)


def real_neuron(cur):
    """Target neuron driven by current batch ``cur`` ``(B, T)``.

    Uses the offline brian2 ``lif_model`` or, with ``--use-daq``, the live
    rig via ``NiDAQClampLayer``. Returns ``(spikes, mem)`` each ``(T, B)``
    with ``mem`` normalised to ~[0, 1] over the [-70, -40] mV span.
    """
    out_mem = []
    spk_ = []
    for b in range(cur.shape[0]):
        series = cur[b]
        if _daq_layer is None:
            # Offline proxy via brian2 LIF model.
            resp, spikes = lif_model(series.detach().cpu().numpy() * 1000, dt=dt)
        else:
            # Real DAQ via NiDAQClampLayer. Convert the (units-free)
            # network signal to pA the same way lif_model does internally
            # (input * 1000) and feed it through the clamp. Spikes are
            # taken from libni's proxy-spike flag if enabled, else zeros.
            I_pA = (series.detach().cpu() * 1000.0).to(torch.float32)
            with torch.no_grad():
                v_volts = _daq_layer(I_pA).cpu().numpy()
            resp = v_volts * 1000.0  # V -> mV (parity with lif_model)
            spikes = (
                _daq_layer.last_spikes.cpu().numpy()
                if _daq_layer.last_spikes is not None
                else np.zeros_like(resp)
            )
        out_mem.append(torch.as_tensor(resp, dtype=torch.float32, device=device))
        spk_.append(torch.as_tensor(spikes, dtype=torch.float32, device=device))

    real_mem = torch.stack(out_mem, dim=1)          # (T, B)
    real_mem = (real_mem - -70) / (-40 - -70)
    real_spk = torch.stack(spk_, dim=1)             # (T, B)
    return real_spk, real_mem


#daq = ni.init_ni(0.1, 0.1, 1/0.05)

# Define Network
network = proxy_net().to(device)


def input_gen():
    inp = torch.randn(batch_size, num_steps, device=device)
    inp[:, 200:-200] += 0.25   # step-current epoch in the middle
    return inp


# Live-plotting callback (keeps fit_proxy itself headless).
loss_hist = []


def _plot_cb(it, loss, t):
    loss_hist.append(loss)
    plt.figure(num=1)
    plt.clf()
    plt.plot(t["proxy_mem"].detach().cpu().numpy()[:, 0], label="proxy")
    plt.plot(t["real_mem"].detach().cpu().numpy()[:, 0], label="real")
    plt.legend()
    plt.pause(0.01)
    plt.figure(num=3)
    plt.clf()
    plt.plot(t["proxy_spk"].detach().cpu().numpy()[:, 0], label="proxy")
    plt.plot(t["real_spk"].detach().cpu().numpy()[:, 0] + 1, label="real")
    plt.legend()
    plt.pause(0.01)
    plt.figure(num=2)
    plt.clf()
    plt.plot(loss_hist)
    plt.pause(0.01)


# Fit the proxy to the real neuron: spike-only van Rossum distance plus a
# membrane-MSE regulariser (set mem_weight=0 for a pure spike-timing fit).
utils.fit_proxy(
    network, real_neuron, input_gen,
    n_iter=600, lr=5e-2, spike_loss="vanrossum", mem_weight=1.0,
    tau=5e-3, dt=dt, device=device, callback=_plot_cb,
)

#clean up the DAQ
#ni.clean_up()
#plot the loss hist
plt.figure()
plt.plot(loss_hist)
plt.show()