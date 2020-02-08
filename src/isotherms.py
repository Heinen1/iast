import numpy as np


# Convert input parameters to numpy array.
def langmuir(p, params):
    np_params = np.array(params)
    np_p = np.array(p)

    langmuir = langmuirFit(np_p, np_params[:, 0], np_params[:, 1])

    return langmuir


def langmuirFit(x, qsat, b):
    return qsat*x*b / (1 + b*x)
