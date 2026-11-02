"""Original three-observation orbit-determination workflow."""
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("input_file")
args = parser.parse_args()

import numpy as np
from ssp.odlib import *

with open(args.input_file, encoding="utf-8") as input_file:
    file = input_file.read().split()
tausInitial = []
k = 0.0172020989484
for i in range(0,3):
    tausInitial.append(julian(False, int(file[0+9*i]), int(file[1+9*i]), int(file[2+9*i]), file[3+9*i]))
timesInitialJD = [tausInitial[0], tausInitial[1], tausInitial[2]]
taus = [k*(tausInitial[0]-tausInitial[1]), k*(tausInitial[2]-tausInitial[1]), k*(tausInitial[2]-tausInitial[0])] #tau1, tau3, tau
# print("taus:",taus)

sunVec1 = [float(file[6]), float(file[7]), float(file[8])]
sunVec2 = [float(file[15]), float(file[16]), float(file[17])]
sunVec3 = [float(file[24]), float(file[25]), float(file[26])]
# print(sunVec1)
# print(sunVec2)
# print(sunVec3)

rho1HatVec = [] #source: pg7 of OD packet
rho1HatVec.append(np.cos(np.radians(convertToDecDegrees(file[4], False, True)))*np.cos((np.radians(convertToDecDegrees(file[5], True, False)))))
rho1HatVec.append(np.sin(np.radians(convertToDecDegrees(file[4], False, True)))*np.cos((np.radians(convertToDecDegrees(file[5], True, False)))))
rho1HatVec.append((np.sin(np.radians(convertToDecDegrees(file[5], True, False)))))
# print(convertToDecDegrees(file[5], True, False))
# print(rho1Vec)

rho2HatVec = [] #source: pg7 of OD packet
rho2HatVec.append(np.cos(np.radians(convertToDecDegrees(file[13], False, True)))*np.cos((np.radians(convertToDecDegrees(file[14], True, False)))))
rho2HatVec.append(np.sin(np.radians(convertToDecDegrees(file[13], False, True)))*np.cos((np.radians(convertToDecDegrees(file[14], True, False)))))
rho2HatVec.append((np.sin(np.radians(convertToDecDegrees(file[14], True, False)))))

# print(rho2Vec)

rho3HatVec = [] #source: pg7 of OD packet
rho3HatVec.append(np.cos(np.radians(convertToDecDegrees(file[22], False, True)))*np.cos((np.radians(convertToDecDegrees(file[23], True, False)))))
rho3HatVec.append(np.sin(np.radians(convertToDecDegrees(file[22], False, True)))*np.cos((np.radians(convertToDecDegrees(file[23], True, False)))))
rho3HatVec.append((np.sin(np.radians(convertToDecDegrees(file[23], True, False)))))
# print(rho2Vec)

Ds = []
Ds.append(np.dot(rho1HatVec, np.cross(rho2HatVec, rho3HatVec))) #add d0
Ds.append(np.dot(np.cross(rho1HatVec, sunVec1), rho3HatVec)) #add d21
Ds.append(np.dot(np.cross(rho1HatVec, sunVec2), rho3HatVec)) #add d22
Ds.append(np.dot(np.cross(rho1HatVec, sunVec3), rho3HatVec)) #add d23

initialRoots, initialRhos = scaleEqLagrange(taus, sunVec2, rho2HatVec, Ds)
# print(initialRoots, initialRhos)

print("Roots:", initialRoots)
print("Rhos:", initialRhos)
userChoice = int(input("enter index for root: "))
initialRoot = initialRoots[userChoice]
initialRho = initialRhos[userChoice]
print("You picked: Rhos:", initialRho, "Roots:", initialRoot)

mu = 1
u = mu/(initialRoot**3)

D0 = Ds[0]
D1s = [np.dot(np.cross(sunVec1, rho2HatVec), rho3HatVec), 
       np.dot(np.cross(sunVec2, rho2HatVec), rho3HatVec), 
       np.dot(np.cross(sunVec3, rho2HatVec), rho3HatVec)]
D2s = [Ds[1], Ds[2], Ds[3]]
D3s = [np.dot(rho1HatVec, np.cross(rho2HatVec, sunVec1)), 
       np.dot(rho1HatVec, np.cross(rho2HatVec, sunVec2)), 
       np.dot(rho1HatVec, np.cross(rho2HatVec, sunVec3))]

f1, g1 = initialFg(initialRoot, taus[0])
f3, g3 = initialFg(initialRoot, taus[1])
print("fs: ", f1, g1, f3, g3)

c1 = g3/(f1*g3-g1*f3)
c2 = -1
c3 = -g1/(f1*g3-g1*f3)

print("cs:", c1, c2, c3)
rhoMagsVec = []
rhoMagsVec.append((c1*D1s[0]+c2*D1s[1]+c3*D1s[2])/(c1*D0))
rhoMagsVec.append((c1*D2s[0]+c2*D2s[1]+c3*D2s[2])/(c2*D0))
rhoMagsVec.append((c1*D3s[0]+c2*D3s[1]+c3*D3s[2])/(c3*D0))

r1Vec = rhoMagsVec[0]*np.array(rho1HatVec) - sunVec1
r2Vec = rhoMagsVec[1]*np.array(rho2HatVec) - sunVec2
r3Vec = rhoMagsVec[2]*np.array(rho3HatVec) - sunVec3

d1 = -f3/(f1*g3 - f3*g1)
d3 = f1/(f1*g3 - f3*g1)

r2dot = d1*r1Vec + d3*r3Vec

rhoGuess = rhoMagsVec[1]
print("rhomags:", rhoMagsVec)
num = 0
while (abs(rhoGuess-initialRho) > 1e-12):
    num = num+1
    initialRho = rhoGuess
    newTimes = []
    for i in range(0, 3):
        newTimes.append(timesInitialJD[i] - (rhoMagsVec[i]/173.14463))
    taus = [k*(newTimes[0]-newTimes[1]), k*(newTimes[2]-newTimes[1]), k*(newTimes[2]-newTimes[0])]

    f1, g1 = fg(taus[0], r2Vec, r2dot, 5)
    f3, g3 = fg(taus[1], r2Vec, r2dot, 5)

    c1 = g3/(f1*g3-g1*f3)
    c2 = -1
    c3 = -g1/(f1*g3-g1*f3)

    rhoMagsVec = []
    rhoMagsVec.append((c1*D1s[0]+c2*D1s[1]+c3*D1s[2])/(c1*D0))
    rhoMagsVec.append((c1*D2s[0]+c2*D2s[1]+c3*D2s[2])/(c2*D0))
    rhoMagsVec.append((c1*D3s[0]+c2*D3s[1]+c3*D3s[2])/(c3*D0))

    r1Vec = rhoMagsVec[0]*np.array(rho1HatVec) - sunVec1
    r2Vec = rhoMagsVec[1]*np.array(rho2HatVec) - sunVec2
    r3Vec = rhoMagsVec[2]*np.array(rho3HatVec) - sunVec3

    d1 = -f3/(f1*g3 - f3*g1)
    d3 = f1/(f1*g3 - f3*g1)

    r2dot = d1*r1Vec + d3*r3Vec
    rhoGuess = rhoMagsVec[1]

    # print("Guess:", rhoGuess)
    # print("Diff:", abs(rhoGuess-initialRho))
    
print(num, "iterations")
print("r2dot =", r2dot)
print("r2 =", r2Vec)

epsilon = -np.radians(23.43829194)
rotationMatrix = np.array([[1, 0, 0], [0, np.cos(epsilon), -np.sin(epsilon)], [0, np.sin(epsilon), np.cos(epsilon)]])

velocityVec = np.matmul(rotationMatrix,r2dot)
positionVec = np.matmul(rotationMatrix, r2Vec)

realE = 0.5317493696688808
# print(realE)
realIN = 3.875284337689891
# print(realIN)
realOM = 238.8360495504369
# print(realOM)
realW = 107.8214075204697
# print(realW)
realA = 2.726573219901208
# print(realA)
realMA = 341.378679
# print(realMA)
realTp = 2457853.36766
# print(realTp)

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

argOfPerihelion = getArgOfPerihelion(TA, positionVec, LAN, inclination)
print("************************************************************")
print("argument of perihelion:", argOfPerihelion)
print("actual:", realW)
print("percent error:", 100*abs(realW-argOfPerihelion)/realW)

meanAnomaly = getMeanAnomaly(positionVec, TA, semiMajor, eccentricity)
print("************************************************************")
print("mean anomaly:", meanAnomaly)
print("actual: ", realMA)
print("percent error:", 100*abs(realMA-meanAnomaly)/realMA)

timeOfLastPerihelion = getTimeOfLastPerihelion(meanAnomaly, semiMajor, julian(False, 2018, 7, 13, "00:00:00"))
print("************************************************************")
print("time of last perihelion:", timeOfLastPerihelion)
print("actual: ", realTp)
print("percent error:", 100*abs(realTp-timeOfLastPerihelion)/realTp)
print("************************************************************")
