"""Astronomy calculations: orbital_elements."""

import math as m
import numpy as np

def getMomentum(r, v):
    return np.cross(r, v)

def getAngle(cosine, sine):
    initialAcos = m.acos(cosine)
    initialAsin = m.asin(sine)
    if (abs(cosine**2 + sine**2 - 1) > 0.1):
        return "error"
    else:
        result = 0
        if (sine == 0 and cosine == 1):
            result = 0
        elif (sine == 1 and cosine == 0):
            result = m.pi*0.5
        elif (sine == 0 and cosine == -1):
            result = m.pi
        elif (sine == -1 and cosine == 0):
            result = 1.5*m.pi
        elif (sine > 0 and cosine > 0):
            result = abs(initialAsin)
        elif (sine > 0 and cosine < 0):
            result = initialAcos
        elif (sine < 0 and cosine < 0):
            result = 1*m.pi + abs(initialAsin)
        elif (sine < 0 and cosine > 0):
            result = 2*m.pi-abs(initialAsin)
    return m.degrees(result)

def getPerihelion(rDotVec, hVec):
    return np.cross(-1*hVec, rDotVec) - np.linalg.norm(rDotVec)

def getA(rVec, rDotVec):
    return 1/((2/np.linalg.norm(rVec))-np.linalg.norm(rDotVec)**2)

def getInclination(hVec):
    return np.degrees(np.arctan(np.sqrt(hVec[0]**2 + hVec[1]**2)/hVec[2]))

def getEccentricity(hVec, a):
    return np.sqrt(1-(np.linalg.norm(hVec)**2/a))

def getLongitudeOfAscendingNode(hVec, iDegrees):
    cosOmega = -hVec[1]/(np.linalg.norm(hVec)*np.sin(np.radians(iDegrees)))
    sinOmega = hVec[0]/(np.linalg.norm(hVec)*np.sin(np.radians(iDegrees)))
    return getAngle(cosOmega, sinOmega)

def getTrueAnomaly(a, e, hVec, rVec, rDotVec):
    rMag = np.linalg.norm(rVec)
    hMag = np.linalg.norm(hVec)
    sinV = ((a*(1 - e**2))/(e*hMag))*((np.dot(rVec, rDotVec))/rMag)
    # print("sinV:", sinV)
    cosV = (1/e)*(((a*(1 - e**2))/rMag) - 1)
    # print("cosV:", cosV)
    return getAngle(cosV, sinV)

def getArgOfPerihelion(tAnomaly, rVec, LAN, i):
    x = rVec[0]
    y = rVec[1]
    z = rVec[2]
    rMag = np.linalg.norm(rVec)
    omega = np.radians(LAN)

    cosU = (x*np.cos(omega) + y*np.sin(omega))/rMag
    sinU = z/(rMag*np.sin(np.radians(i)))
    # print(sinU)

    u = getAngle(cosU, sinU)
    initialArg = np.radians(u)-np.radians(tAnomaly)
    while initialArg < 0:
        initialArg = initialArg + 2*np.pi
    while initialArg > 2*np.pi:
        initialArg = initialArg - 2*np.pi
    return np.degrees(initialArg)

def getMeanAnomaly(rVec, trueAnomaly, a, e):
    r = np.linalg.norm(rVec)
    E0 = np.arccos((1/e)*(1-(r/a)))
    if(np.radians(trueAnomaly) > np.pi):
        return 360 - np.degrees(E0 - e*np.sin(E0))
    return np.degrees(E0 - e*np.sin(E0))

def getTimeOfLastPerihelion(meanAnomaly, a, julianDate):
    k = 0.01720209895
    return -1*((np.radians(meanAnomaly)/np.sqrt((k**2)/(a**3))) - julianDate)
