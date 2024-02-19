import torch, torch.nn as nn
import snntorch as snn
import sys
import faulthandler

faulthandler.enable() #to debug seg faults and timeouts
batch_size = 4
data_path='./data/mnist'
sys.path.append('/home/smestern/Dropbox/RTXI/ni_interface')
#import ni_generic as ni
from lif_tester import lif_model
import utils_snn as utils

import matplotlib.pyplot as plt
import numpy as np 
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
        self.beta_1 = nn.Parameter(torch.tensor(beta))
        self.alpha_1 = nn.Parameter(torch.tensor(alpha))
        self.threshold = nn.Parameter(torch.tensor(0.8))
        self.scale_fator = nn.Parameter(torch.tensor(1.0))
        self.lif1 = utils.sLIFin(0.1, 30, -50, 0.0002, 0.1, 0.1, -20)#snn.Synaptic(beta=self.beta_1, alpha=self.alpha_1, init_hidden=False,  learn_alpha=True, learn_beta=True, learn_threshold=False, 
                        #         threshold=self.threshold)
        self.tanh = nn.Tanh()
        
    def _forward(self, x):
        syn1, mem1 = self.lif1.init_synaptic()

        #record the cur, syn1 and mem for the chosen unit
        input_rec = []
        spikes_rec = []
        mem_rec = []
        for step in range(num_steps):
                cur1 = x[:, step] 
                spk1, syn1, mem1 = self.lif1(cur1 * self.scale_fator, torch.tensor(0), syn1, mem1)
                mem1 = self.tanh(mem1)
                #record the cur, syn1 and mem for the chosen unit======================================
                input_rec.append(cur1)
                spikes_rec.append(spk1)
                mem_rec.append(mem1)


        return  torch.stack(input_rec, dim=0), torch.stack(spikes_rec, dim=0), torch.stack(mem_rec, dim=0)

    def forward(self, x):
        cur_daq, spikes_daq, mem_daq = self._forward(x)

        #write and read from DAQ
        #the curr and spikes are in the format of [num_steps, batch_size, num_units]
        #we want to write the spikes of the chosen unit to the DAQ
        out_mem = []
        spk_ = []
        for x in range(cur_daq.shape[1]):
            temp_cur = cur_daq[:, x]
            
            resp, spikes = lif_model(temp_cur.detach().cpu().numpy() * 1000, dt=dt)
            _x_proxy = torch.tensor(resp, dtype=torch.float32)
            out_mem.append(torch.tensor(_x_proxy, dtype=torch.float32, device=device))
            spk_.append(torch.tensor(spikes, dtype=torch.float32, device=device))

        out_mem = torch.stack(out_mem, dim=1)
        out_mem = (out_mem - -70)/(-40 - -70)
        spk_ = torch.stack(spk_, dim=1)
        return mem_daq,  out_mem, spikes_daq, spk_
        
    
#daq = ni.init_ni(0.1, 0.1, 1/0.05)

# Define Network
network = proxy_net()


# Define Loss
error = nn.MSELoss()
error2 = utils.EMD(dt=dt, duration=0.25, bin_size=0.002)

# Define Optimizer
learning_rate = 5e-2
momentum = 0.9
optimizer = torch.optim.Adam(network.parameters(), lr=learning_rate)
lr_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=1000, eta_min=1e-9)
loss_hist = []
for _ in range(600):
    # Forward pass
    # Initialize input
    input = torch.randn(batch_size, num_steps, device=device)
    #at time step 100, set the input to 1
    input[:, 200:-200] += 0.25
    output, _output_proxy, spk, spk_proxy = network(input)
    #loss is between the output and the proxy
    loss = error(output, _output_proxy)
    loss += torch.mean(error2(spk, spk_proxy))

    loss_hist.append(loss.detach().cpu().numpy())
    # Backward pass
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    if loss < 5.0:
        lr_scheduler.step()
    print(loss)
    #plot the output and the proxy
    plt.figure(num=1)
    plt.clf()
    plt.plot(output.detach().cpu().numpy()[:,0], label="output")
    plt.plot(_output_proxy.detach().cpu().numpy()[:,0], label="proxy")
    #plot the spikes
    plt.legend()
    plt.pause(0.01)
    plt.figure(num=3)
    plt.clf()
    plt.plot(spk.detach().cpu().numpy()[:,0], label="output")
    plt.plot(spk_proxy.detach().cpu().numpy()[:,0] + 1, label="proxy")
    plt.pause(0.01)
    plt.figure(num=2)
    plt.clf()
    plt.plot(loss_hist)
    plt.pause(0.01)

#clean up the DAQ
#ni.clean_up()
#plot the losst hist
plt.figure()
plt.plot(loss_hist)
plt.show()