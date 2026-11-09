from __future__ import annotations

import numpy as np
from scipy import integrate
from scipy.optimize import fsolve

try:
    from .isotherms import langmuir, langmuirNumerical
except ImportError:  # pragma: no cover - direct script execution
    from isotherms import langmuir, langmuirNumerical

LOWER_LIMIT = 1e-5


def iast(P, S):
    """Compute the binary IAST mixture loading for each pressure point."""
    pressure = np.asarray(P, dtype=float)
    composition = np.asarray(S["y"], dtype=float)
    n_components = S["n_components"]
    isotherm_parameters = np.asarray(S["isotherm_parameters"], dtype=float)

    if pressure.ndim != 1:
        raise ValueError("Pressure must be a 1D array")
    if composition.shape != (n_components,):
        raise ValueError("Gas composition length must match the number of components")
    if isotherm_parameters.shape != (n_components, 2):
        raise ValueError("Isotherm parameters must have shape (n_components, 2)")

    x = np.zeros((n_components, pressure.size))
    x0 = setInitialValuesForSolver(pressure, S)
    solution = fsolve(
        iastLangmuirNumerical,
        x0,
        args=(pressure, composition, isotherm_parameters),
    )

    x[0, :] = solution
    x[1, :] = 1.0 - x[0, :]

    pi0 = np.outer(composition, pressure) / x
    q = langmuir(pi0.T, isotherm_parameters)
    q_total = 1.0 / np.sum(x.T / q, axis=1)
    return q_total * x


def setInitialValuesForSolver(P, S):
    product = np.outer(P, S["y"])
    q_guess = langmuir(product, S["isotherm_parameters"])
    guess = q_guess[:, 0] / np.sum(q_guess, axis=1)
    if S["y"][0] < LOWER_LIMIT and S["y"][1] < LOWER_LIMIT:
        return np.zeros_like(guess)
    return guess


def iastLangmuirNumerical(x, pressure, composition, isotherm):
    isotherm1 = isotherm[0]
    isotherm2 = isotherm[1]
    f_total = np.zeros_like(pressure, dtype=float)

    pi01 = (composition[0] * pressure) / x
    pi02 = (composition[1] * pressure) / (1.0 - x)

    for idx, (pi0_1, pi0_2) in enumerate(zip(pi01, pi02)):
        f1, _ = integrate.quad(langmuirNumerical, 0.0, pi0_1, args=(isotherm1,))
        f2, _ = integrate.quad(langmuirNumerical, 0.0, pi0_2, args=(isotherm2,))
        f_total[idx] = f1 - f2
    return f_total


def iastLangmuirAnalytical(x, pressure, composition, isotherm):
    isotherm1 = isotherm[0]
    isotherm2 = isotherm[1]
    return (
        isotherm1[0] * np.log(1.0 + isotherm1[1] * composition[0] * pressure / x)
        - isotherm2[0] * np.log(1.0 + isotherm2[1] * composition[1] * pressure / (1.0 - x))
    )


compute_iast = iast
