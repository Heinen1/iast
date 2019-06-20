import sys
import os
import numpy as np
from isotherms import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

columnPressure = 1
columnLoading = 8
skipNumberOfLines = 22
number_of_isotherm_parameters = 2

def fitIsothermModels(isothermType,nComponents,fileName1,fileName2):
	coefficients = np.zeros([nComponents,number_of_isotherm_parameters])

	checkIfFileExists(fileName1)
	checkIfFileExists(fileName2)

	pressure_loading_1 = storeDataPoints(fileName1)
	pressure_loading_2 = storeDataPoints(fileName2)

	coefficients[0] = fitIsotherm(pressure_loading_1)
	coefficients[1] = fitIsotherm(pressure_loading_2)
	return coefficients

def checkIfFileExists(fileName):
	file_exists = os.path.isfile(fileName)
	if file_exists == False:
		print("File ",fileName," does not exist")
		exit()

def storeDataPoints(fileName):
	i=0
	pressureLoading = list()        
	with open (fileName, "r") as myfile:
		for line in myfile:
			if i > skipNumberOfLines:
				lineStrip = line.strip()
				lineSplit = lineStrip.split()
				lineSegments = [lineSplit[columnPressure-1], lineSplit[columnLoading-1]]
				pressureLoading.append(lineSegments) 
			i = i + 1
	return pressureLoading 

def fitIsotherm(pressureLoading):
	n = len(pressureLoading)
	pressure = np.zeros(n) 
	loading = np.zeros(n) 
	for i in range(n):
		pressure[i] = pressureLoading[i][0] 
		loading[i] = pressureLoading[i][1] 
	init_values = [loading[n-1], loading[1]/pressure[1]]
	popt, pcov = curve_fit(langmuirFit, pressure, loading, p0=init_values)
	return popt

