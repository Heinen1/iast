import numpy as np

def langmuir(p,params):
	n = len(p)
	langmuir = np.zeros(n)
	for i in range(n):
		langmuir[i] = params[i][0]*params[i][1]*p[i]/(1+(params[i][1]*p[i]))
	return langmuir 
