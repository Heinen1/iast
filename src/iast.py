import numpy as np

from scipy.optimize import fsolve
from scipy import integrate
from isotherms import langmuir
from isotherms import langmuirNumerical

lower_limit = 0.00001


# Perform iast calculation. Using fsolve() to find the root of the equation.
# For a binary mixture this is the the equiality of the integrated isotherms
# with a unknow upper limit.
def iast(P, S):
    q = np.zeros([S.n_components, S.n_pressure_intervals])
    x = np.zeros([S.n_components, S.n_pressure_intervals])
    qtot = np.zeros(S.n_pressure_intervals)

    x0 = setInitialValuesForSolver(P, S)
    data = [P, S.y, S.isotherm_parameters]

    # Integrating can be done numerically or analytically for Langmuir. Use:
    # iastLangmuirNumerical
    # iastLangmuirAnalytical
    f = fsolve(iastLangmuirNumerical, x0, args=data)

    x[0, :] = f
    x[1, :] = 1 - x[0, :]

    pi0 = np.outer(S.y, P) / x

    q = langmuir(pi0.T, S.isotherm_parameters)

    qtot = 1.0/np.sum(x.T/q, axis=1)
    q_mix = qtot*x

    return q_mix


# The fsolver requires initial conditions. Set as initial condition,
# the partial pressure of component 1 divided by the sum of all partial
# pressures. If gas phase fraction is almost zero (lower limit), set the
# initial value at that fraction to zero.
def setInitialValuesForSolver(P, S):
    q_guess = np.zeros([S.n_components, S.n_pressure_intervals])
    quess = np.zeros(S.n_pressure_intervals)

    product = np.outer(P, S.y)
    q_guess = langmuir(product, S.isotherm_parameters)
    quess = q_guess[:, 0]/(np.sum(q_guess, axis=1))

    return np.where(S.y[0] < lower_limit and S.y[1] < lower_limit, 0, quess)


# Integrate the langmuir isotherm / p numerically with upper limit set as
# p_i^0 = y_i * P / x_i, with x_i the unknown value.
# Cannot specify list / array as upper integration limit for integrate.quad!
def iastLangmuirNumerical(x, data):
    P, y, isotherm = data[0], data[1], data[2]

    isotherm1 = isotherm[0]
    isotherm2 = isotherm[1]

    f_total = np.zeros(len(P))

    pi01 = (y[0] * P) / x
    pi02 = (y[1] * P) / (1 - x)

    for i, pi0 in enumerate(zip(pi01, pi02)):
        f1, err1 = integrate.quad(langmuirNumerical, 0, pi0[0],
                                  args=(isotherm1))
        f2, err2 = integrate.quad(langmuirNumerical, 0, pi0[1],
                                  args=(isotherm2))
        f_total[i] = f1 - f2

    return f_total


# Analystical expression for Langmuir isotherm, divided by the pressure,
# that is used in IAST calculation.
def iastLangmuirAnalytical(x, data):
    P, y, isotherm = data[0], data[1], data[2]

    isotherm1 = isotherm[0]
    isotherm2 = isotherm[1]

    f = isotherm1[0]*np.log(1+isotherm1[1]*y[0]*P/x) - isotherm2[0]*np.log(
                    1+isotherm2[1]*((y[1]*P)/(1 - x)))

    return f
