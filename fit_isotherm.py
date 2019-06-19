import sys
import os
import numpy as np
from isotherms import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

col_pressure = 1
col_loading = 8
skip_number_of_lines = 22
number_of_adsorbates = 2
number_of_isotherm_parameters = 2

def fitIsothermModels(arguments):
	print("Isotherm model: ",isotherm_models[int(arguments[1])])
	Ncomponents = len(arguments)-2
	coefficients = np.zeros([Ncomponents,number_of_isotherm_parameters])
	pressure_loading_1 = storeDataPoints(arguments[2])
	pressure_loading_2 = storeDataPoints(arguments[3])
	coefficients[0] = fitIsotherm(pressure_loading_1)
	coefficients[1] = fitIsotherm(pressure_loading_2)
	return coefficients

def storeDataPoints(file_name):
	i=0
	pressure_loading = list()        
	with open (file_name, "r") as myfile:
		for line in myfile:
			if i > skip_number_of_lines:
				line_strip = line.strip()
				line_split = line_strip.split()
				line_segments = [line_split[col_pressure-1], line_split[col_loading-1]]
				pressure_loading.append(line_segments) 
			i = i + 1
	return pressure_loading 

def fitIsotherm(pressure_loading):
	n = len(pressure_loading)
	p = np.zeros(n) 
	q = np.zeros(n) 
	for i in range(n):
		p[i] = pressure_loading[i][0] 
		q[i] = pressure_loading[i][1] 
	# Initial guess values Langmuir isotherm: qsat = loading last point, b = Henry coefficient loading / pressure in low pressure regime
	init_values = [q[n-1], q[1]/p[1]]
	popt, pcov = curve_fit(langmuir_fit, p, q, p0=init_values)
	return popt

