import numpy as np
from scipy.optimize import fsolve
from isotherms import *

def iast(P,S):
	q_guess = np.zeros(S.Ncomponents)
	q = np.zeros(S.Ncomponents)
	x = np.zeros(S.Ncomponents)

	q_guess = langmuir(P*S.y,S.params)
	quess = q_guess[0]/(np.sum(q_guess))

	x0 = np.where(S.y[0]<0.00001 and S.y[1]<0.00001,0,quess)
	data = [P,S.y,S.params]

	f = fsolve(iast_langmuir_analytical,x0,args=data)

	x[0] = f
	x[1] = 1 - x[0]
	
	q = langmuir(P*(S.y/x),S.params)

	qtot = 1.0/np.sum(x/q) 
	q_mix = qtot * x

	return q_mix

def iast_langmuir_analytical(x,data):
	P,y,isotherm = data[0], data[1], data[2]
	isotherm1 = isotherm[0]
	isotherm2 = isotherm[1]

	f = isotherm1[0]*np.log(1+isotherm1[1]*y[0]*P/x)-isotherm2[0]*np.log(1+isotherm2[1]*((y[1]*P)/(1-x)))
	return f

