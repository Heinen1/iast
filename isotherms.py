import numpy as np

isotherm_models = ['langmuir', 'freundlich'] 

def langmuir(p,params):
	n = len(p)
	langmuir = np.zeros(n)
	for i in range(n):
		langmuir[i] = params[i][0]*params[i][1]*p[i]/(1+(params[i][1]*p[i]))
	return langmuir 

def langmuir_fit(x,qsat,b):
	return qsat*x*b/(1+b*x)
