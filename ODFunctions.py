##########################################
##### Functions required for the OD ######
#####    MILENA MATHEW, SSP 2018    ######
##########################################

import numpy as np
from math import *

### CONSTANTS USED IN FUNCTIONS ###
k = 0.01720209895
c = 173.145
eEarth = radians(23.4352)

#########################################

def JD(year, month, day, hour): # Converts civil date and time (UT) to a Julian date
    J0 = 367 * year - floor(1.75 * (year + floor((month + 9) / 12))) + floor(275 * month / 9) + day + 1721013.5
    JD = J0 + hour/24
    return JD

def gaussianTime(t1, t2, t3): # Converts a Julian date to Gaussian time
    tau1 = k * (t1 - t2)
    tau3 = k * (t3 - t2)
    tau = tau3 - tau1
    results = [tau1, tau3, tau]
    return results

def rhohat(alpha, dec): # Generates a unit vector in the direction of the earth-object vector
    rhohat = np.array([cos(alpha)*cos(dec), sin(alpha)*cos(dec), sin(dec)])
    return rhohat

def initialFandG(r2, times): # Calculates initial values for f and g; times are in gaussian days
    u2 = 1/(r2**3)
    f1 = 1 - 0.5*u2*times[0]**2
    f3 = 1 - 0.5*u2*times[1]**2
    g1 = times[0] - (u2*times[0]**3)/6
    g3 = times[1] - (u2*times[1]**3)/6
    fandg = np.array([f1, f3, g1, g3, r2])
    return fandg

#Scalar Equation of Lagrange
def scalarEquationOfLagrange(obs1, obs2, obs3): # Solves the Scalar Equation of Lagrange to determine initial value of r2
    times = gaussianTime(obs1[2], obs2[2], obs3[2])
    # Generates earth-object unit vectors
    rho1 = rhohat(obs1[0], obs1[1])
    rho2 = rhohat(obs2[0], obs2[1])
    rho3 = rhohat(obs3[0], obs3[1])
    # Grabs sun vectors from initial data
    sun1 = np.array([obs1[3], obs1[4], obs1[5]])
    sun2 = np.array([obs2[3], obs2[4], obs2[5]])
    sun3 = np.array([obs3[3], obs3[4], obs3[5]])
    # Calculates constants needed for the Scalar Equation of Lagrange
    A1 = times[1]/times[2]
    A3 = -times[0]/times[2]
    B1 = A1 * ((times[2]**2)-(times[1]**2))/6
    B3 = A3 * ((times[2]**2)-(times[0]**2))/6
    D0 = np.dot(rho1, np.cross(rho2, rho3))
    D21 = np.dot(np.cross(rho1, sun1), rho3)
    D22 = np.dot(np.cross(rho1, sun2), rho3)
    D23 = np.dot(np.cross(rho1, sun3), rho3)

    A = (A1*D21-D22+A3*D23)/(-D0)
    B = (B1*D21+B3*D23)/(-D0)
    E = -2*np.dot(rho2, sun2)
    F = (np.linalg.norm(sun2))**2
    a = -(A**2+A*E+F)
    b = -(2*A*B+B*E)
    c = -B**2
    # Determines the roots and throws out imaginary ones
    r2possible = np.roots([1, 0, a, 0, 0, b, 0, 0, c])
    r2possible = r2possible[np.isreal(r2possible)].real
    reasonableRoots = np.array([])
    # Throws out roots that result initially in a negative value for rho
    for root in r2possible:
        reasonableRho = A + B/(root**3)
        if root > 0 and reasonableRho > 0:
            reasonableRoots = np.append(reasonableRoots, root)
    return reasonableRoots

def findRho(obs1, obs2, obs3, fandg): # Determines rho for each observation given fs and gs
    # Unit vector
    rhoHat1 = rhohat(obs1[0], obs1[1])
    rhoHat2 = rhohat(obs2[0], obs2[1])
    rhoHat3 = rhohat(obs3[0], obs3[1])

    sun1 = np.array([obs1[3], obs1[4], obs1[5]])
    sun2 = np.array([obs2[3], obs2[4], obs2[5]])
    sun3 = np.array([obs3[3], obs3[4], obs3[5]])

    D0 = np.dot(rhoHat1, np.cross(rhoHat2, rhoHat3))
    D11 = np.dot(np.cross(sun1, rhoHat2), rhoHat3)
    D12 = np.dot(np.cross(sun2, rhoHat2), rhoHat3)
    D13 = np.dot(np.cross(sun3, rhoHat2), rhoHat3)
    D21 = np.dot(np.cross(rhoHat1, sun1), rhoHat3)
    D22 = np.dot(np.cross(rhoHat1, sun2), rhoHat3)
    D23 = np.dot(np.cross(rhoHat1, sun3), rhoHat3)
    D31 = np.dot(rhoHat1, np.cross(rhoHat2, sun1))
    D32 = np.dot(rhoHat1, np.cross(rhoHat2, sun2))
    D33 = np.dot(rhoHat1, np.cross(rhoHat2, sun3))

    f1 = fandg[0]
    f3 = fandg[1]
    g1 = fandg[2]
    g3 = fandg[3]

    denom = f1*g3-g1*f3
    d1 = -f3/denom
    d3 = f1/denom
    c1 = g3/denom
    c2 = -1
    c3 = -g1/denom
    # Magnitude of earth-object vector
    pmag1 = (c1*D11+c2*D12+c3*D13)/(c1*D0)
    pmag2 = (c1*D21+c2*D22+c3*D23)/(c2*D0)
    pmag3 = (c1*D31+c2*D32+c3*D33)/(c3*D0)
    # Complete earth-object vector
    rho1 = pmag1 * rhoHat1
    rho2 = pmag2 * rhoHat2
    rho3 = pmag3 * rhoHat3

    rho = np.array([rho1, rho2, rho3, d1, d3])
    return rho

def findR2dot(obs1, obs2, obs3, rho): # Determines r2dot from rho values
    sun1 = np.array([obs1[3], obs1[4], obs1[5]])
    sun2 = np.array([obs2[3], obs2[4], obs2[5]])
    sun3 = np.array([obs3[3], obs3[4], obs3[5]])
    r1 = rho[0] - sun1
    r3 = rho[2] - sun3
    r2dot = rho[3]*r1+rho[4]*r3
    return r2dot

def fAndg(time, r2, r2dot): # Determines new f and g values using the Taylor series
    u = 1/(np.linalg.norm(r2)**3)
    z = np.dot(r2, r2dot)/(np.linalg.norm(r2)**2)
    q = np.dot(r2dot, r2dot)/(np.linalg.norm(r2)**2) - u
    f = 1 - 0.5*u*(time**2)+0.5*u*z*(time**3)+ (1/24)*(3*u*q - 15*u*(z**2)+u**2)*(time**4)
    g = time - (1/6)*u*(time**3)+0.25*u*z*(time**4)
    return f, g

def rotateInv(vector, angle, axis): # Rotates vectors, opposite direction from normal
    x = np.linalg.inv(np.array([[1,0,0],
                  [0, cos(angle), -1*sin(angle)],
                  [0, sin(angle), cos(angle)]]))
    y = np.linalg.inv(np.array([[cos(angle),0,sin(angle)],
                  [0, 1, 0],
                  [-sin(angle), 0, cos(angle)]]))
    z = np.linalg.inv(np.array([[cos(angle),-sin(angle),0],
                  [sin(angle), cos(angle), 0],
                  [0, 0, 1]]))
    if axis=="x":
        result = np.dot(x, vector)
    if axis=="y":
        result = np.dot(y, vector)
    if axis=="z":
        result = np.dot(z, vector)
    return result

def findQuadrant(sine, cosine): # Removes quadrant ambiguities
    if cosine > 0 and sine > 0:
        return asin(sine)
    if cosine < 0 and sine > 0:
        return acos(cosine)
    if cosine < 0 and sine < 0:
        return pi - asin(sine)
    if cosine > 0 and sine < 0:
        return 2*pi + asin(sine)

def babyOD(r, rdot, middleTime, toPrint = "n", reportMiddle = "y"): # Determines orbital elements from r and rdot
    # Rotates to rectangular ecliptic coordinates
    r = rotateInv(r, eEarth, "x")
    rdot = rotateInv(rdot, eEarth, "x")
    a = ((2/np.linalg.norm(r))-(np.dot(rdot, rdot)))**-1
    e = (1-(np.linalg.norm(np.cross(r, rdot))**2/a))**0.5
    h = k * np.cross(r, rdot)
    I = acos(h[2]/np.linalg.norm(h))

    sinomega = h[0]/(np.linalg.norm(h)*sin(I))
    cosomega = -1*h[1]/(np.linalg.norm(h)*sin(I))
    Omega = findQuadrant(sinomega, cosomega)

    sinwtheta = r[2]/(np.linalg.norm(r)*sin(I))
    coswtheta = (cos(Omega)**-1)*((r[0]/np.linalg.norm(r))+cos(I)*sinwtheta*sin(Omega))
    wtheta = findQuadrant(sinwtheta, coswtheta)
    theta = asin((np.dot(r, rdot)*((a*(1-e**2))**0.5))/(e*np.linalg.norm(r)))
    w = wtheta - theta

    if theta > 0:
        E = acos((1/e)*(1-(np.linalg.norm(r)/a)))
    if theta < 0:
        E = (2*pi)-acos((1/e)*(1-(np.linalg.norm(r)/a)))

    M2 = E - e*sin(E)
    mu = k**2
    time = JD(2018, 7, 22, 6) # year, month, day, hour
    n = (mu/(a**3))**0.5
    M = M2 + n*(time - middleTime) # Precesses orbital elements to date defined in time variable
    if toPrint == "y": # Prints output if wanted
        print("a =", a, "AU")
        print("e =", e)
        print("I =", degrees(I), "degrees")
        print("Omega =", degrees(Omega), "degrees")
        print("w =", degrees(w), "degrees")
        print("M =", degrees(M), "degrees")
    orbitalElements = np.array([a, e, I, Omega, w, M])
    if reportMiddle == "n":
        return orbitalElements
    return orbitalElements, middleTime

    ################ ITERATIONS ###################
def iteratingMoG(tolerance, sun2, obs1, obs2, obs3, fandg, toPrint="n", returnR = "n", reportMiddle = "y"): # Iteration portion of the Method of Gauss
    rho = findRho(obs1, obs2, obs3, fandg)
    results = findR2dot(obs1, obs2, obs3, rho)
    r2dot = results[0]
    r2 = rho[1] - sun2
    numIterations = 0 # Initializes variables for use in while loop
    prevr2 = 0
    originalTimes = (obs1[2], obs2[2], obs3[2])
    times = gaussianTime(originalTimes[0], originalTimes[1], originalTimes[2])
    middleTime = 0 # stores light corrected JD time from middle observation
    while abs(np.linalg.norm(r2)-np.linalg.norm(prevr2)) > tolerance:
        numIterations = numIterations + 1 # Keeps track of iterations
        prevr2 = r2
        for x in range(0, 3): # Corrects for light travel time
            times[x] = originalTimes[x] - (np.linalg.norm(rho[x])/c)
        middleTime = times[1]
        times = gaussianTime(times[0], times[1], times[2]) # Converts light corrected time to gaussian days
        f1 = fAndg(times[0], r2, r2dot)[0] # Inefficiently calculates new f and g values
        g1 = fAndg(times[0], r2, r2dot)[1]
        f3 = fAndg(times[1], r2, r2dot)[0]
        g3 = fAndg(times[1], r2, r2dot)[1]
        denom = (f1*g3)-(g1*f3)
        c1 = g3/denom
        c2 = -1
        c3 = -g1/denom
        fandg = [f1, f3, g1, g3]
        rho = findRho(obs1, obs2, obs3, fandg)
        r2 = rho[1] - sun2
        r2dot = findR2dot(obs1, obs2, obs3, rho)
        if numIterations > 1000: # Leaves loop if it seems the calculation likely will not converge
            if toPrint == "y":
                print("The number of iterations has exceeded 1000 and the calculation has not converged. The current value for r2 is", r2, ". Ending iterations.")
                return np.array([float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')])
            else:
                return np.array([float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN'), float('NaN')])
    if toPrint == "y":
        print("The asteroid's position vector for the central observation is", rotateInv(r2, eEarth, "x"), "AU",
          "and the velocity vector is", rotateInv(r2dot, eEarth, "x"), "AU/day")
        print("The range to the asteroid in the central observation is", np.linalg.norm(rho[1]), "AU")
    if returnR == "y": # returns r and rdot instead of orbital elements
        return r2, r2dot
    babyODel = babyOD(r2, r2dot, middleTime, toPrint, reportMiddle)
    return babyODel

def MethodOfGauss(obs1, obs2, obs3, sun2, tolerance, toPrint="n", returnR = "n", reportMiddle = "y"):
    times = gaussianTime(obs1[2], obs2[2], obs3[2])
    reasonableRoots = scalarEquationOfLagrange(obs1, obs2, obs3)
    allOrbitalElements = np.array([])
    for root in reasonableRoots: # Tests all reasonable roots
        if toPrint == "y":
            print("Calculations for reasonable root:", root)
        fandg = initialFandG(root, times)
        allOrbitalElements = np.append(allOrbitalElements, iteratingMoG(tolerance, sun2, obs1, obs2, obs3, fandg, toPrint, returnR, reportMiddle))
    return allOrbitalElements # returns an array for orbital elements of all the orbital elements for all reasonable calculations

def openInput(fileName): # Opens the default input file
    with open(fileName, 'r') as data:
        inputText = data.read()
        inputText = inputText.replace(':',' ')
    inputText = inputText.split()
    numObservations = int(len(inputText)/15)
    observations = np.zeros([numObservations, 6])

    for num in range(numObservations):
        beginCount = (num)*15
        year = float(inputText[beginCount])
        month = float(inputText[beginCount+1])
        day = float(inputText[beginCount+2])
        hour = float(inputText[beginCount+3]) + float(inputText[beginCount+4])/60 + float(inputText[beginCount+5])/3600
        time = JD(year, month, day, hour)
        RA = radians((float(inputText[beginCount+6])+float(inputText[beginCount+7])/60+float(inputText[beginCount+8])/3600)*15)
        dec = radians(float(inputText[beginCount+9])+float(inputText[beginCount+10])/60+float(inputText[beginCount+11])/3600)
        sunx = float(inputText[beginCount+12])
        suny = float(inputText[beginCount+13])
        sunz = float(inputText[beginCount+14])
        observations[num] = np.array([RA, dec, time, sunx, suny, sunz])
    return observations
