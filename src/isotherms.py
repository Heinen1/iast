import numpy as np


def langmuir_fit(x, qsat, b):
    x = np.asarray(x, dtype=float)
    qsat = np.asarray(qsat, dtype=float)
    b = np.asarray(b, dtype=float)
    return (qsat * x * b) / (1.0 + b * x)


def langmuir(p, params):
    params = np.asarray(params, dtype=float)
    if params.ndim != 2 or params.shape[1] != 2:
        raise ValueError("Langmuir parameters must have shape (n_components, 2)")
    return langmuir_fit(p, params[:, 0], params[:, 1])


def langmuirNumerical(x, params):
    qstat = float(params[0])
    b = float(params[1])
    return langmuir_fit(x, qstat, b) / x


# Backwards-compatible aliases retained for the original project API.
langmuirFit = langmuir_fit
