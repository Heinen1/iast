import numpy as np

from sys import exit
from os import path
from isotherms import langmuirFit
from scipy.optimize import curve_fit

# The first columns contains the pressures in [Pa]
column_pressure = 1
# The 8th columns contains the absolute loading in [mol/kg]
column_loading = 8
# The first 22 lines of the RASPA output file should be skiped
skip_n_lines = 22
# Langmuir isotherm has 2 parameters
number_of_isotherm_parameters = 2


# Fit the RASPA calculations to an isotherm model (Langmuir only for now)
def fitIsotherm(isotherm_type, *files_names):
    coefficients = []

    for file_name in files_names:
        checkIfFileExists(file_name)
        pressure_loading_list = getDataPoints(file_name)
        coefficients.append(doFit(pressure_loading_list))

    return coefficients


# If the file does not exists, kill the program
def checkIfFileExists(file_name):
    file_exists = path.isfile(file_name)
    if not file_exists:
        print("File ", file_name, " does not exist")
        exit()


# Open the RASPA output file that contains the pressure and the loadings of
# the single component adsorption isotherms
def getDataPoints(file_name):
    i = 0
    pressure_and_loading = list()

    with open(file_name, "r") as file_content:
        for line in file_content:
            if i > skip_n_lines:
                pressure_and_loading.append(getPressureAndLoadingOnly(line))
            i += 1

    return pressure_and_loading


# Remove white spaces and split the line into a list
# Get the elements in the list that corresponds to the pressure and the loading
def getPressureAndLoadingOnly(line):
    line_split = line.strip().split()
    return [line_split[column_pressure-1], line_split[column_loading-1]]


# Convert list with pressure and loadings to a numpy array
# Separate pressures from loadings and set initial values
# Fit the isotherm to the datapoints using the initial values from above
def doFit(pressure_loading_list):
    pressure_loading_array = np.array(pressure_loading_list, dtype=float)
    n = len(pressure_loading_list)

    pressure = pressure_loading_array[:, 0]
    loading = pressure_loading_array[:, 1]

    init_values = [loading[n-1], loading[1]/pressure[1]]
    popt, pcov = curve_fit(langmuirFit, pressure, loading, p0=init_values)

    return popt
