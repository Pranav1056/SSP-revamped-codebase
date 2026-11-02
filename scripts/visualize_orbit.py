import argparse
parser = argparse.ArgumentParser()
parser.add_argument("elements_file")
args = parser.parse_args()

import txaio
txaio.use_asyncio()
from vpython import *
import numpy as np
from math import sin, cos

file = open("D:\elements.txt", "r")
elements = file.readlines()

a = float(elements[0])
e = float(elements[1])
M = radians(float(elements[2]))
Oprime = radians(float(elements[3]))
iprime = radians(float(elements[4]))
wprime = radians(float(elements[5]))

def solvekep(M):
    Eguess = M
    Mguess = Eguess - e*sin(Eguess)
    Mtrue = M
    while abs(Mguess - Mtrue) > 1e-004:
        Mguess = Eguess - e*sin(Eguess)
        Eguess = Eguess - (M - (Eguess - e*sin(Eguess))) / (e*cos(Eguess)-1)
    return Eguess

sqrtmu = 0.01720209895
mu = sqrtmu**2
time = 0
dt = .05
period = sqrt(4*pi**2*a**3/mu)
r1ecliptic = vector(0, 0, 0)
Mtrue = 2*pi/period*(time) + M
Etrue = solvekep(Mtrue)

x = a*cos(Etrue) - a*e
y = a*np.sqrt(1-e**2)*np.sin(Etrue)
z = 0
posVec = np.array([[x], [y], [z]])
print(posVec)

matrix1 = np.array([[cos(Oprime), -sin(Oprime), 0], [sin(Oprime), cos(Oprime), 0], [0, 0, 1]])
matrix2 = np.array([[1, 0, 0], [0, cos(iprime), -sin(iprime)], [0, sin(iprime), cos(iprime)]])
matrix3 = np.array([[cos(wprime), -sin(wprime), 0], [sin(wprime), cos(wprime), 0], [0, 0, 1]])

np.matmul(matrix1, matrix2)

newVector = np.matmul(matrix1, np.matmul(matrix2, np.matmul(matrix3, posVec)))
asteroid = sphere(pos=r1ecliptic*150, radius=(15), color=color.white)
asteroid.trail = curve(color=color.white)
sun = sphere(pos=vector(0,0,0), radius=(50), color=color.yellow)
epsilon = np.radians(23.44)

rotationMatrix = np.array([[1, 0, 0], [0, cos(epsilon), -sin(epsilon)], [0, sin(epsilon), cos(epsilon)]])
equitorialV = np.matmul(rotationMatrix, newVector)
finalPos = vector(equitorialV[0][0], equitorialV[1][0], equitorialV[2][0])
print(finalPos)
while (True):
    rate(200)
    time = time + 1
    Mtrue = 2*pi/period*(time) + M
    Etrue = solvekep(Mtrue)
    x = a*cos(Etrue) - a*e
    y = a*np.sqrt(1-e**2)*np.sin(Etrue)
    z = 0
    posVec = np.array([[x], [y], [z]])

    newVector = np.matmul(matrix1, np.matmul(matrix2, np.matmul(matrix3, posVec)))
    equitorialV = np.matmul(rotationMatrix, newVector)
    finalPos = vector(equitorialV[0][0], equitorialV[1][0], equitorialV[2][0])
    asteroid.pos = finalPos*150
    asteroid.trail.append(pos=asteroid.pos)  
