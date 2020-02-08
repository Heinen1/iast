import numpy as np

from scipy.optimize import fsolve
from scipy import integrate
from isotherms import langmuir

lower_limit = 0.00001


# Perform iast calculation. Using fsolve() to find the root of the equation.
# For a binary mixture this is the the equiality of the integrated isotherms
# with a unknow upper limit.
def iast(P, S):
    q_guess = np.zeros(S.n_components)
    q = np.zeros(S.n_components)
    x = np.zeros(S.n_components)

    q_guess = langmuir(P*S.y, S.isothermParameters)
    quess = q_guess[0]/(np.sum(q_guess))

    x0 = np.where(S.y[0] < lower_limit and S.y[1] < lower_limit, 0, quess)
    data = [P, S.y, S.isothermParameters]

    f = fsolve(iastLangmuirAnalytical, x0, args=data)

    x[0] = f
    x[1] = 1 - x[0]

    q = langmuir(P*(S.y/x), S.isothermParameters)

    qtot = 1.0/np.sum(x/q)
    q_mix = qtot*x

    return q_mix


# Analystical expression for Langmuir isotherm, divided by the pressure,
# that is used in IAST calculation.
def iastLangmuirAnalytical(x, data):
    P, y, isotherm = data[0], data[1], data[2]

    isotherm1 = isotherm[0]
    isotherm2 = isotherm[1]

    f = isotherm1[0]*np.log(1+isotherm1[1]*y[0]*P/x) - isotherm2[0]*np.log(
                    1+isotherm2[1]*((y[1]*P)/(1 - x)))
    return f