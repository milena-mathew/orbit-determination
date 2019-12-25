# Milena Mathew
# Monte Carlo Error Analysis
# Varies initial alphas and decs to determine error (Gaussian distribution)

import numpy as np
from math import *
import matplotlib.pyplot as plt
from ODFunctions import *
from time import time, sleep

########### POSSIBLE USER INPUT #################
tolerance = 1*10**(-12)
numberToSample = 5000
fileName = "1992JE Results.txt"
n_bins = 50 # Number of bins used in histogram display
##################################################

# Reads input file which has appropriate observations, errors etc.
inputText = np.loadtxt(fileName)

obs1 = np.array([inputText[0], inputText[1], inputText[2], inputText[3], inputText[4], inputText[5]]) #grabs observation information
obs2 = np.array([inputText[6], inputText[7], inputText[8], inputText[9], inputText[10], inputText[11]])
sun2 = np.array([obs2[3], obs2[4], obs2[5]])
obs3 = np.array([inputText[12], inputText[13], inputText[14], inputText[15], inputText[16], inputText[17]])
variedOrbitalElements = np.array([])

errors = np.array([[radians(inputText[25]/3600), radians(inputText[26]/3600)], # in arcseconds
                   [radians(inputText[27]/3600), radians(inputText[28]/3600)],
                   [radians(inputText[29]/3600), radians(inputText[30]/3600)]])

deltaMatrix1 = np.array([np.random.normal(obs1[0], errors[0,0], numberToSample), np.random.normal(obs1[1], errors[0,1], numberToSample)]) # creates arrays of variations [alpha, dec]
deltaMatrix2 = np.array([np.random.normal(obs2[0], errors[1,0], numberToSample), np.random.normal(obs2[1], errors[1,1], numberToSample)])
deltaMatrix3 = np.array([np.random.normal(obs3[0], errors[2,0], numberToSample), np.random.normal(obs3[1], errors[2,1], numberToSample)])
count = 0

t0 = time() #Keeps track of how long this takes
for num in range(numberToSample):
    count += 1
    obs1[0] = deltaMatrix1[0, num]
    obs1[1] = deltaMatrix1[1, num]
    obs2[0] = deltaMatrix2[0, num]
    obs2[1] = deltaMatrix2[1, num]
    obs3[0] = deltaMatrix3[0, num]
    obs3[1] = deltaMatrix3[1, num]
    variedOrbitalElements = np.append(variedOrbitalElements, MethodOfGauss(obs1, obs2, obs3, sun2, tolerance, "n", "n", "n"))
    print(count) #prints what iteration it's on to keep track of progress
t1 = time()
print(numberToSample, "iterations took", t1-t0, "seconds")

a = variedOrbitalElements[0::6] # grabs data for analysis 
e = variedOrbitalElements[1::6]
I = variedOrbitalElements[2::6]
Omega = variedOrbitalElements[3::6]
w = variedOrbitalElements[4::6]
M = variedOrbitalElements[5::6]

a = a[~np.isnan(a)]
e = e[~np.isnan(e)]
I = I[~np.isnan(I)]
Omega = Omega[~np.isnan(Omega)]
w = w[~np.isnan(w)]
M = M[~np.isnan(M)]

aError = np.nanstd(a)
eError = np.nanstd(e)
IError = np.nanstd(I)
OmegaError = np.nanstd(Omega)
wError = np.nanstd(w)
MError = np.nanstd(M)

aMedian = np.nanmedian(a)
eMedian = np.nanmedian(e)
IMedian = np.nanmedian(I)
OmegaMedian = np.nanmedian(Omega)
wMedian = np.nanmedian(w)
MMedian = np.nanmedian(M)

# Histograms for data analysis
aplot = plt.figure()
plt.hist(a, bins= n_bins)
plt.xlabel("a")
plt.ylabel("Frequency")

eplot = plt.figure()
plt.hist(e, bins= n_bins)
plt.xlabel("e")
plt.ylabel("Frequency")

Iplot = plt.figure()
plt.hist(I, bins= n_bins)
plt.xlabel("I")
plt.ylabel("Frequency")

Omegaplot = plt.figure()
plt.hist(Omega, bins= n_bins)
plt.xlabel("Omega")
plt.ylabel("Frequency")

wPlot = plt.figure()
plt.hist(w, bins= n_bins)
plt.xlabel("w")
plt.ylabel("Frequency")

MPlot = plt.figure()
plt.hist(M, bins= n_bins)
plt.xlabel("M")
plt.ylabel("Frequency")

# Saves figures to folder for later use
aplot.savefig("aplot.png")
eplot.savefig("eplot.png")
Iplot.savefig("Iplot.png")
Omegaplot.savefig("Omegaplot.png")
wPlot.savefig("wplot.png")
MPlot.savefig("Mplot.png")

print("The median value of a is", aMedian, "AU")
print("The median value of e is", eMedian)
print("The median value of I is", degrees(IMedian), "degrees")
print("The median value of Omega is", degrees(OmegaMedian), "degrees")
print("The median value of w is", degrees(wMedian), "degrees")
print("The median value of M is", degrees(MMedian), "degrees")

print("The error of a is", aError, "AU")
print("The error of e is", eError)
print("The error of I is", degrees(IError), "degrees")
print("The error of Omega is", degrees(OmegaError), "degrees")
print("The error of w is", degrees(wError), "degrees")
print("The error of M is", degrees(MError), "degrees")

plt.show()
