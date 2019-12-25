# Milena Mathew
# OD
# Determines orbital elements from observations of asteroid

import numpy as np
from math import *
from ODFunctions import *

################## USER INPUT #####################

tolerance = 1 * 10**(-12) #float(input("What is the name of the file containing the reference stars? "))
fileName = "1992JE Input.txt" #input("What is the name of the file containing the reference stars? ")
resultsFileName = "Sample OD Results.txt" #input("What is the name of the file you would like to save the results to?")

######################### COMBINATIONS OF OBSERVATIONS ###################
observations = openInput(fileName)
numObservations = len(observations)
obs1 = np.copy(observations[0])
obs3 = np.copy(observations[-1])
middleObservation = np.arange(1, numObservations-1)
prelimOrbitalElements = np.array([])
middleTimes = np.array([])
for num in middleObservation:
    obs2 = np.copy(observations[num])
    sun2 = np.array([obs2[3], obs2[4], obs2[5]])
    print("\nUsing observations 1,", num + 1, "and", str(numObservations) + ":")
    resultsMoG = MethodOfGauss(obs1, obs2, obs3, sun2, tolerance, "y")
    prelimOrbitalElements = np.append(prelimOrbitalElements, resultsMoG[0])
    middleTimes = np.append(middleTimes, resultsMoG[1])

# OUTPUTS RESULTS TO A FILE THAT CAN BE USED FOR DIFFERENTIAL CORRECTION/MONTE CARLO/GENERAL USE #
if len(prelimOrbitalElements)/6 == 1:
    whichData = 0
else:
    whichData = int(input("Which set of results would you like to save? " + str(np.arange((len(prelimOrbitalElements)/6)))))

#Write to file
resultsFile = open(resultsFileName, "w+")
resultsFile.write("# RESULTS FILE OF ORBITAL DETERMINATION \n")
resultsFile.write("# ----------- DATA FROM RELEVANT OBSERVATIONS (FIRST, SECOND, THIRD) ----------- \n")
resultsFile.write("# RA (radians), DEC (radians), TIME (JD), SUN VECTOR (x,y,z [rectangular equatorial]) \n")
resultsFile.write("\n")
for item in observations[0]:
  resultsFile.write("%s\n" % item)
resultsFile.write("\n")
for item in observations[whichData+1]:
  resultsFile.write("%s\n" % item)
resultsFile.write("\n")
for item in observations[-1]:
  resultsFile.write("%s\n" % item)
resultsFile.write("\n")
resultsFile.write("# LIGHT CORRECTED TIME OF MIDDLE OBSERVATION (JD) \n")
resultsFile.write(str(middleTimes[whichData]))
resultsFile.write("\n")
resultsFile.write("# ORBITAL ELEMENTS (a, e, I, Omega, w, M) IN AU AND RADIANS WHERE APPLICABLE \n")
if whichData == 0:
    for item in prelimOrbitalElements:
        resultsFile.write("%s\n" % item)
else:
    for item in prelimOrbitalElements[whichData]:
        resultsFile.write("%s\n" % item)
resultsFile.write("\n")
resultsFile.write("# RA AND DEC UNCERTAINTIES (IN ARCSECONDS) TO USE IN MONTE CARLO PROGRAM \n")
resultsFile.write("# FORMAT: RA1, DEC1, RA2, DEC2, RA3, DEC3 \n")
resultsFile.close()
