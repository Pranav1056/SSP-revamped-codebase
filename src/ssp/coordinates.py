"""Astronomy calculations: coordinates."""

import numpy as np
from .kepler import newtonRaphson
from .orbital_elements import getAngle

def julian(isLocalTime, yr, month, day, time):
    term1 = (month+9)//12
    term2 = 7*(term1 + yr)//4
    julianAt0 = 367*yr - term2 + 275*month//9 + day + 1721013.5
    UT = float(time[0:2])+ float(time[3:5])/60 + float(time[6:])/3600
    if isLocalTime:
        UT = UT + 4
    return julianAt0  + UT/24

def convertToDecimalHours(stringToConvert, hasSign, toDegrees):
    coef = 1
    sum = 0
    if(hasSign):
        if (stringToConvert[0] == "-"):
            coef = -1
        sum = coef*(float(stringToConvert[1:3]) + float(stringToConvert[4:6])/60 + float(stringToConvert[7:])/3600)
    else:
        sum = float(stringToConvert[:2]) + float(stringToConvert[3:5])/60 + float(stringToConvert[6:])/3600
    if(toDegrees):
        sum = sum*360/24
    return sum

def convertToDecDegrees(stringToConvert, hasSign, wasHours):
    if (not wasHours):
        coef = 1
        sum = 0
        if(hasSign):
            if (stringToConvert[0] == "-"):
                coef = -1
            sum = coef*(float(stringToConvert[1:3]) + float(stringToConvert[4:6])/60 + float(stringToConvert[7:])/3600)
        else:
            sum = float(stringToConvert[:2]) + float(stringToConvert[3:5])/60 + float(stringToConvert[6:])/3600
        return sum
    else: 
        coef = 1
        sum = 0
        if(hasSign):
            if (stringToConvert[0] == "-"):
                coef = -1
            sum = coef*((float(stringToConvert[1:3])*15) + float(stringToConvert[4:6])*0.25 + float(stringToConvert[7:])*360/86400)
        else:
            sum = (float(stringToConvert[:2])*15) + float(stringToConvert[3:5])*0.25 + float(stringToConvert[6:])*360/86400
        return sum

def getRaAndDec(a, eccentricity, timePerihelion, inclination, LAN, argOfPerihelion, currentTimeJulian, earthSunVec):
    k = 0.01720209895
    n = k*np.sqrt(1/(a**3))
    initialM = n*(currentTimeJulian - timePerihelion)
    E = newtonRaphson(initialM, eccentricity, 0.000000001) #may change tolerance

    inclination = np.radians(inclination)
    LAN = np.radians(LAN)
    argOfPerihelion = np.radians(argOfPerihelion)

    x = a*np.cos(E) - a*eccentricity
    y = a*np.sqrt(1-eccentricity**2)*np.sin(E)
    z = 0

    posVec = np.array([[x], [y], [z]])
    matrix1 = np.array([[np.cos(LAN), -np.sin(LAN), 0], [np.sin(LAN), np.cos(LAN), 0], [0, 0, 1]])
    matrix2 = np.array([[1, 0, 0], [0, np.cos(inclination), -np.sin(inclination)], [0, np.sin(inclination), np.cos(inclination)]])
    matrix3 = np.array([[np.cos(argOfPerihelion), -np.sin(argOfPerihelion), 0], [np.sin(argOfPerihelion), np.cos(argOfPerihelion), 0], [0, 0, 1]])

    newVector = np.matmul(matrix1, np.matmul(matrix2, np.matmul(matrix3, posVec)))
    epsilon = np.radians(23.44)
    rotationMatrix = np.array([[1, 0, 0], [0, np.cos(epsilon), -np.sin(epsilon)], [0, np.sin(epsilon), np.cos(epsilon)]])
    equitorialPosition = np.matmul(rotationMatrix, newVector)
    rangeVec = equitorialPosition + earthSunVec
    # print("rangevec: ", rangeVec)
    pHat = rangeVec/np.linalg.norm(rangeVec)
    # print(pHat)

    Dec = np.arcsin(pHat[2][0])
    RA = getAngle(pHat[0][0]/np.cos(Dec), pHat[1][0]/np.cos(Dec))
    return RA, np.degrees(Dec)
