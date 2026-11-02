import numpy as np
from odlib import *

file = open("D:\Ji_OD\JiInput.txt").read().split()
# print(file)
# *************UNITS IN AU & DAYS******************
realE = float(file[9])
# print(realE)
realIN = float(file[11])
# print(realIN)
realOM = float(file[13])
# print(realOM)
realW = float(file[15])
# print(realW)
realTA = float(file[17])
# print(realTA)
realA = float(file[19])
# print(realA)
realMA = float(file[21])
# print(realMA)
realTp = float(file[23])
# print(realTp)

positionVec = np.array([float(file[1]), float(file[2]), float(file[3])])
velocityVec = np.array([float(file[5]), float(file[6]), float(file[7])])
velocityVec = velocityVec*(365.2563835/(2*np.pi))

hVec = getMomentum(positionVec, velocityVec)
print("************************************************************")
print("momentum:", getMomentum(positionVec, velocityVec))

perihelion = getPerihelion(velocityVec, hVec)
print("************************************************************")
print("perihelion:", perihelion)

semiMajor = getA(positionVec, velocityVec)
print("************************************************************")
print("semi major axis:", semiMajor)
print("actual:", realA)
print("percent error:", 100*abs(realA-semiMajor)/realA)

inclination = getInclination(hVec)
print("************************************************************")
print("inclination:", inclination)
print("actual:", realIN)
print("percent error:", 100*abs(realIN-inclination)/realIN)

eccentricity = getEccentricity(hVec, semiMajor)
print("************************************************************")
print("eccentricity:", eccentricity)
print("actual:", realE)
print("percent error:", 100*abs(realE-eccentricity)/realE)

LAN = getLongitudeOfAscendingNode(hVec, inclination)
print("************************************************************")
print("longitude of ascending node:", LAN)
print("actual:", realOM)
print("percent error:", 100*abs(realOM - LAN)/realOM)

TA = getTrueAnomaly(semiMajor, eccentricity, hVec, positionVec, velocityVec)
print("************************************************************")
print("true anomaly: ", TA)
print("actual:", realTA)
print("percent error:", 100*abs(realTA-TA)/realTA)

argOfPerihelion = getArgOfPerihelion(TA, positionVec, LAN, inclination)
print("************************************************************")
print("argument of perihelion:", argOfPerihelion)
print("actual:", realW)
print("percent error:", 100*abs(realW-argOfPerihelion)/realW)

meanAnomaly = getMeanAnomaly(positionVec, semiMajor, eccentricity)
print("************************************************************")
print("mean anomaly:", meanAnomaly)
print("actual: ", realMA)
print("percent error:", 100*abs(realMA-meanAnomaly)/realMA)

timeOfLastPerihelion = getTimeOfLastPerihelion(meanAnomaly, semiMajor, julian(False, 2018, 7, 13, "00:00:00"))
print("************************************************************")
print("time of last perihelion:", timeOfLastPerihelion)
print("actual: ", real‘p 
prigt,b0ercGnu$gsror:",!±2!bshbd+,Ùp-4imeGFL#stPeri(eliol)/rÂalTv)çzinv∏"*)j+*j*
**:j(*™**j&´**
®****.:8*.****
j*Û/*´*™**"***
*™**;**i