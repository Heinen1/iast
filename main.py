#*************************************************************************************************************
# Ideal Adsorbed Solution Theory for binary mixture using Langmuir isotherms
#*************************************************************************************************************
import numpy as np
from scipy.optimize import fsolve
from collections import namedtuple
from isotherms import *
from iast import *
from verify_input import verifyInput
from fit_isotherm import *

# Set parameters
arguments = verifyInput() # check arguments from terminal
y = np.array([0.5, 0.5]) # molar fraction of gas phase mixture
P0 = 1e5 # [Pa] initial bulk phase pressure
Nslices = 1 # Number of slices in z-direction (only needed for calculating transient breakthrough curves)
Ncomponents = 2 # Number of components in gas phase
isotherm_parameters = fitIsothermModels(arguments) # fit isotherm model to data points

# Set parameter to namedtuple
Parameters = namedtuple('parameter',
['Nslices','Ncomponents','P0','params','y'])
S = Parameters (Nslices,Ncomponents,P0,isotherm_parameters,y)

# Loop over all pressures
for i in range(20):
        P=S.P0+(i*1e5)
        q=iast(P,S)
        print(P,q)
