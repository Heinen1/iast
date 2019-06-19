import sys
import os
import numpy as np
from isotherms import *

expected_number_of_arguments = 4

def verifyInput():
	checkNumberOfArguments()
	checkIfIsothermModelExists()
	checkIfFileExists(sys.argv[2])
	checkIfFileExists(sys.argv[3])
	return sys.argv 

def checkIfFileExists(file_name):
        file_exists = os.path.isfile(file_name)
        if file_exists == False:
                print("File ",file_name," does not exist")
                exit()
                
def checkNumberOfArguments():
	if len(sys.argv) < expected_number_of_arguments:
		print("Number of expected arguments is ",expected_number_of_arguments-1)
		exit()

def checkIfIsothermModelExists():
        if int(sys.argv[1]) > len(isotherm_models):
                print("This isotherm model does not exist.");
                exit()
