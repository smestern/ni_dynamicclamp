from brian2 import *
from brian2.equations.equations import SUBEXPRESSION, PARAMETER, DIFFERENTIAL_EQUATION
import os
import numpy
from brian2.utils.topsort import topsort
from collections import namedtuple

#Create the equation object to bypass the default equations check

