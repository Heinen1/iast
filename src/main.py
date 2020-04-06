#
# Ideal Adsorbed Solution Theory for binary mixture using Langmuir isotherms
#

import numpy as np

from collections import namedtuple
from iast import iast
from fit_isotherm import fitIsotherm

# define parameters: number of components in mixture, gas phase composition,
# initial and final pressure, number of pressure points
n_components = 2
y1 = 0.5
y2 = 1 - y1
pressure_initial = 1e5
pressure_final = 2e6
n_pressure_intervals = 20
pressure = np.linspace(pressure_initial, pressure_final, n_pressure_intervals)

# define RASPA output files from monte carlo simulations that contain
# the pressure-loading datapoints
# fit isotherms to Langmuit model
file_name1 = "Results.dat-ZIF-8-433K-methane"
file_name2 = "Results.dat-ZIF-8-433K-ethane"
isotherm_parameters = fitIsotherm(file_name1, file_name2)

# put all parameters in tuple for easier handling of input parameters
Parameters = namedtuple('parameter', ['pressure_initial',
                                      'isotherm_parameters', 'y',
                                      'n_components',
                                      'n_pressure_intervals'])
S = Parameters(pressure_initial, isotherm_parameters, np.array([y1, y2]),
               n_components, n_pressure_intervals)

# perform actual iast
loading = iast(pressure, S)
print(pressure)
print(loading)
