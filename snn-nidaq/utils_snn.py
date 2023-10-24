import ot
import torch





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
            _coords = torch.arange(x.shape[0], device=x.device)
        else:
            _coords = torch.arange(x.shape[0], device=x.device)
        losses = []
        for i in range(x.shape[1]):
            losses.append(ot.wasserstein_1d(_coords, _coords, x[:, i], y[:, i],  p=2))
        error = torch.mean(torch.stack(losses, dim=0))
        return  error

    def bin(self, x):
        #bin the spikes, in a rolling window
        #x is the spikes
        window_size = int(self.bin_size/self.dt)
        step_size = int(window_size/2)
        x = torch.sum(x.unfold(0,window_size,step_size), dim=0)
        return x