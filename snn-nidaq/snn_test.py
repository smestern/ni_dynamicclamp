import torch, torch.nn as nn
import snntorch as snn
import sys
import faulthandler
faulthandler.enable() #to debug seg faults and timeouts
batch_size = 128
data_path='./data/mnist'
sys.path.append('/home/smestern/Dropbox/RTXI/ni_interface')
import ni_generic as ni
pA = 1e-12
device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
from torch.utils.data import DataLoader
from torchvision import datasets, transforms
import numpy as np
import matplotlib.pyplot as plt
import nidaqmx

from snntorch import surrogate
from snntorch import utils

import snntorch.functional as SF
# Network Architecture
num_inputs = 28*28
num_hidden = 100
num_outputs = 10

# Temporal Dynamics
num_steps = int(0.25 * (1/0.0001)) 
beta = 0.08 #membrane potential decay rate in seconds 
alpha = 0.05 #synaptic decay rate in seconds


# Define Network
class NiDAQ_NET(nn.Module):
    def __init__(self, daq):
        super().__init__()

        # Initialize layers
        self.fc1 = nn.Linear(num_inputs, num_hidden)
        self.beta_1 = nn.Parameter(torch.tensor(beta))
        self.alpha_1 = nn.Parameter(torch.tensor(alpha))

        self.lif1 = snn.Synaptic(beta=self.beta_1, alpha=self.alpha_1, init_hidden=False, reset_mechanism='zero', learn_alpha=True, learn_beta=True)
        self.lif2 = snn.Synaptic(beta=beta, alpha=0.8, reset_mechanism='zero')
        
        self.layers = [self.fc1, self.lif1, self.lif2]

        #intiliaze DAQ
        #unit to replace, NOTE: parameterize this later
        self.unit_idx = 0
        self.daq = daq


    def _forward(self, x):
        syn1, mem1 = self.lif1.init_synaptic()
        syn2, mem2 = self.lif2.init_synaptic()

        spk2_rec = []  # Record the output trace of spikes
        mem2_rec = []  # Record the output trace of membrane potential

        #record the cur, syn1 and mem for the chosen unit
        input_rec = []
        spikes_rec = []
        mem_rec = []
        for step in range(num_steps):
                cur1 = self.fc1(x)
                spk1, syn1, mem1 = self.lif1(cur1, syn1, mem1)
                cur2 = spk1 #torch.zeros((batch_size, num_outputs), device=device)
                spk2, syn2, mem2 = self.lif2(cur2, syn2, mem2)

                spk2_rec.append(spk2)
                mem2_rec.append(mem2)

                #record the cur, syn1 and mem for the chosen unit======================================
                input_rec.append(cur1)
                spikes_rec.append(spk1)
                mem_rec.append(mem1)


        return torch.stack(spk2_rec, dim=0), torch.stack(mem2_rec, dim=0), torch.stack(input_rec, dim=0), torch.stack(spikes_rec, dim=0), torch.stack(mem_rec, dim=0)

    def forward(self, x):
        spk_rec, mem_rec, cur_daq, spikes_daq, mem_daq = self._forward(x)

        #write and read from DAQ
        #the curr and spikes are in the format of [num_steps, batch_size, num_units]
        #we want to write the spikes of the chosen unit to the DAQ
        out_mem = []
        for x in range(cur_daq.shape[1]):
            temp_cur = cur_daq[:, x, self.unit_idx]
            temp_spikes = spikes_daq[:, x, self.unit_idx]
            adq = self.write_daq(temp_spikes+temp_cur, num_steps)
            out_mem.append(torch.tensor(adq, dtype=torch.float32, device=device))

        out_mem = torch.stack(out_mem, dim=1)

        #for now compute the error here
        loss_daq = torch.nn.functional.mse_loss(out_mem, mem_daq[:, :, self.unit_idx])

        self.out_mem = out_mem
        return spk_rec, mem_rec, loss_daq

    def write_daq(self, input, samps):
        #write to DAQ
        v_out = ni.loop_clamp(input.detach().cpu().numpy()*pA * 20)*1000
        #min max scale between -100mV and 100mV
        v_out = np.clip(v_out, -100, 100)
        v_out = (v_out - -70)/(100 - -70)
        return v_out

def print_batch_accuracy(data, targets, train=False):
    _, output,_ = net(data.view(batch_size, -1))
    acc = SF.accuracy_rate(output, targets, population_code=True, num_classes=10)

    if train:
        print(f"Train set accuracy for a single minibatch: {acc*100:.2f}%")
    else:
        print(f"Test set accuracy for a single minibatch: {acc*100:.2f}%")

def train_printer():
    print(f"Epoch {epoch}, Iteration {iter_counter}")
    print(f"Train Set Loss: {loss_hist[counter]:.2f}")
    print(f"Test Set Loss: {test_loss_hist[counter]:.2f}")
    print_batch_accuracy(data, targets, train=True)
    print_batch_accuracy(test_data, test_targets, train=False)
    print("\n")


daq = ni.init_ni(0.1, 0.1, 1/0.05)

loss = nn.CrossEntropyLoss()

# Define a transform
transform = transforms.Compose([
            transforms.Resize((28, 28)),
            transforms.Grayscale(),
            transforms.ToTensor(),
            transforms.Normalize((0,), (1,))])

mnist_train = datasets.MNIST(data_path, train=True, download=True, transform=transform)
mnist_test = datasets.MNIST(data_path, train=False, download=True, transform=transform)

# Create DataLoaders
train_loader = DataLoader(mnist_train, batch_size=batch_size, shuffle=True, drop_last=True)
test_loader = DataLoader(mnist_test, batch_size=batch_size, shuffle=True)

beta = 0.9  # neuron decay rate
spike_grad = surrogate.fast_sigmoid()

#  Initialize Network
net = NiDAQ_NET(daq).to(device)


optimizer = torch.optim.Adam(net.parameters(), lr=1e-2, betas=(0.9, 0.999))
loss = SF.loss.ce_count_loss(population_code=True, num_classes=10)

num_epochs = 5
loss_hist = []
test_loss_hist = []
counter = 0

# Outer training loop
for epoch in range(num_epochs):
    iter_counter = 0
    train_batch = iter(train_loader)

    # Minibatch training loop
    for data, targets in train_batch:
        data = data.to(device)
        targets = targets.to(device)

        # forward pass
        net.train()
        spk_rec, mem_rec, loss_daq = net(data.view(batch_size, -1))

        
        #loss_val = loss(spk_rec, targets)/4
        loss_val = loss_daq
        # Gradient calculation + weight update
        optimizer.zero_grad()
        loss_val.backward()
        optimizer.step()

        # Store loss history for future plotting
        loss_hist.append(loss_val.item())

        # Test set
        with torch.no_grad():
            net.eval()
            test_data, test_targets = next(iter(test_loader))
            test_data = test_data.to(device)
            test_targets = test_targets.to(device)

            # Test set forward pass
            test_spk, test_mem, _ = net(test_data.view(batch_size, -1))

            # Test set loss
            test_loss = loss(test_spk, test_targets)
            test_loss_hist.append(test_loss.item())

            # Print train/test loss/accuracy
            if counter % 10 == 0:
                train_printer()
                plot_mem = test_mem[:, 0, :40].detach().cpu().numpy().T
                plt.figure(figsize=(15, 5), num=0)
                plt.clf()
                for neuron in range(plot_mem.shape[0]):
                    if neuron != 0:
                        continue
                    spike_times = np.where(test_spk[:, 0, neuron].detach().cpu().numpy() == 1)[0]
                    plt.plot(plot_mem[neuron], c='k')
                    plt.plot(net.out_mem[:,0], c='r')
                    plt.twinx()

                    if len(spike_times) > 0:
                        plt.scatter(spike_times, plot_mem[neuron, spike_times], s=2, c='r')

                    
                plt.xlabel('Time Steps')
                plt.ylabel('')
                plt.title('Membrane Potential of First 40 Neurons in the First Layer')
                plt.pause(0.1)

            counter += 1
            iter_counter +=1
            #plot the membrane potential of the first 40 neurons in the first layer
            
            
