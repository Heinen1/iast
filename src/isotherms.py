import numpy as np


# Convert input parameters to numpy array.
def langmuir(p, params):
    np_params = np.array(params)
    np_p = np.array(p)

    langmuir = langmuirFit(np_p, np_params[:, 0], np_params[:, 1])

    return langmuir


# Actual Langmuir model with x= pressure, qstat = saturation capacity and
# b = equilibrium constant.
def langmuirFit(x, qsat, b):
    return qsat*x*b / (1 + b*x)


# For IAST, the Langmuir devided by the pressure is integrated. For numericall
# integrating, this should be specified explicitly.
def langmuirNumerical(x, params):
    qstat = params[0]
    b = params[1]

    return langmuirFit(x, qstat, b) / x
