#*************************************************************************************************************
# Ideal Adsorbed Solution Theory for binary mixture using Langmuir isotherms
#*************************************************************************************************************
import numpy as np
from scipy.optimize import fsolve
from collections import namedtuple
from isotherms import *
from iast import *
from fit_isotherm import *

fileName1 = "Results.dat-ZIF-8-433K-methane" 
fileName2 = "Results.dat-ZIF-8-433K-ethane"
isothermType = 0 
nComponents = 2
y1 = 0.5 
pressureInitial = 1e5 
pressureFinal = 2e6 
nPressureIntervals = 20 

deltaPressure = (pressureFinal-pressureInitial)/nPressureIntervals 
isothermParameters = fitIsothermModels(isothermType,nComponents,fileName1,fileName2) 

Parameters = namedtuple('parameter',
['pressureInitial','isothermParameters','y','nComponents','isothermType'])
S = Parameters (pressureInitial,isothermParameters,np.array([y1,1-y1]),nComponents,isothermType)

for i in range(nPressureIntervals):
        P = S.pressureInitial+(i*deltaPressure)
        q = iast(P,S)
        print(P,q)
