import sys
import os
#import ctypes
from ctypes import CDLL, POINTER, cdll
import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import faulthandler
faulthandler.enable(file=sys.stdout) #to debug seg faults and timeouts
#get the path to the current file
path = os.path.dirname(os.path.abspath(__file__))

libc = cdll.LoadLibrary("/home/smestern/Dropbox/RTXI/ni_interface/interface_c.so")

step_clamp_loop = libc.run_step_loop
step_clamp_loop.restype = None

step_clamp = libc.step_clamp
step_clamp.restype = ctypes.c_double

def loop_clamp(I):
    #create a pointer to the numpy input array
    size = I.shape
    I = np.ctypeslib.as_ctypes(I.astype(np.float64))
    #create a numpy array to store the output
    V_out = np.ctypeslib.as_ctypes(np.zeros(size, dtype=np.float64))
    #create a pointer to the numpy output array
    size = ctypes.c_int(size[0])
    #call the function in the shared library
    step_clamp_loop(I, V_out, size)
    return np.ctypeslib.as_array(V_out)


def init_ni(dt, scalefactor_in=0.1, scalefactor_out=1/0.5):

    #convert the input to the correct type
    dt = ctypes.c_double(dt)
    scalefactor_in = ctypes.c_double(scalefactor_in)
    scalefactor_out = ctypes.c_double(scalefactor_out)


    #call the function in the shared library
    libc.init_ni(dt, scalefactor_in, scalefactor_out)
    return

#init_ni(0.1)
if __name__=="__main__":
    import matplotlib.pyplot as plt
    daq = init_ni(0.1, 0.1, 1/0.05)
    #try writing zeros
    for x in range(10):
        #test = step_clamp(ctypes.c_double(0.1*x), ctypes.c_double(-10*x))
        #print(test)
        pass
    #make a sine wave
    I = np.sin(np.arange(0, 10, 1e-4))*30
    V_out = loop_clamp(I*1e-12)
    plt.plot(V_out)
    print(V_out)
    libc.clean_up()
    plt.show()