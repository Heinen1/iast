#*************************************************************************************************************
# Ideal Adsorbed Solution Theory for binary mixture using Langmuir isotherms
#*************************************************************************************************************
import numpy as np
from scipy.optimize import fsolve
from collections import namedtuple
from isotherms import *
from iast import *

# Set parameters
isotherm_type = 1 # 1 = Langmuir
isotherm_parameters = np.array([[6.84703,1.17537e-07],[6.37321,4.85485e-07]]) # Langmuir [qmax1,b1],[qmax2,b2]
y = np.array([0.5, 0.5]) # molar fraction of gas phase mixture
P0 = 1e5 # [Pa] initial bulk phase pressure
Nslices = 1 # Number of slices in z-direction (only needed for calculating transient breakthrough curves)
Ncomponents = 2 # Number of components in gas phase

# Set parameter to namedtuple
Parameters = namedtuple('parameter',
['isotherm_type','Nslices','Ncomponents','P0','params','y'])
S = Parameters (isotherm_type,Nslices,Ncomponents,P0,isotherm_parameters,y)

# Loop over all pressures
for i in range(2):
        P=S.P0+(i*1e5)
        q=iast(P,S)
        print(P,q)
