"""Astronomy calculations: kepler."""

import numpy as np
from scipy import optimize
from .orbital_elements import getA, getEccentricity, getMomentum

def newtonRaphson(m, e, tol):
    oldEtest = m
    ETest = oldEtest - (oldEtest - e*np.sin(oldEtest) - m)/(-e*np.cos(oldEtest) + 1)
    iterationNum = 0
    while (abs(ETest - oldEtest) > tol):
        iterationNum = iterationNum + 1
        oldEtest = ETest
        ETest = oldEtest - (oldEtest - e*np.sin(oldEtest) - m)/(-e*np.cos(oldEtest) + 1)
    # print("E: ", ETest)
    # print("Tolerance: ", tol)
    # print("Iterations: ", iterationNum)
    return ETest

def fg(tau, r2, r2dot, order):
    r = np.linalg.norm(r2)
    a = getA(r2, r2dot)
    mu = 1
    n = np.sqrt(mu/(a**3))
    e = getEccentricity(getMomentum(r2, r2dot), a)
    sign = (np.dot(r2, r2dot)/(n*(a**2)))*np.cos(n*tau - (np.dot(r2, r2dot))/(n*(a**2)))
    + (1-r/a)*np.sin(n*tau - (np.dot(r2, r2dot)/(n*(a**2))))
    if (sign > 0):
        x01 = n*tau + 0.85*e - (np.dot(r2, r2dot)/(n*(a**2)))
    else:
        x01 = n*tau - 0.85*e - (np.dot(r2, r2dot)/(n*(a**2)))
    lambdaThing = lambda x: x - (1-(r)/a)*np.sin(x) + (np.dot(r2, r2dot)/(n*(a**2)))*(1-np.cos(x)) - (n*tau)
    E = optimize.newton(lambdaThing, x0 = x01, tol=1E-12)
    if order == 5:
        return 1-(a/(r))*(1-np.cos(E)), tau + (1/n)*(np.sin(E) - E)
    u = mu/(r**3)
    z = (np.dot(r2, r2dot)/(r**2))
    q = (np.dot(r2dot, r2dot)/(r**2)) - u
    tau = tau
    if(order == 3):
        f1 = 1 - 0.5*u*(tau**2) + 0.5*u*z*(tau**3)
        g1 = tau - (1/6)*mu*(tau**3)/r**3 
    if(order == 4):
        f1 = 1 - 0.5*mu*(tau**2) + 0.5*u*z*(tau**3) + (1/24)*(3*u*q - 15*u*(z**2) + u**2)
        g1 = tau - (1/6)*mu*(tau**3)/r**3  + (1/4)*u*z*(tau**4)
    return f1, g1

def initialFg(r, tau):
    mu = 1.0
    u = mu/r**3
    f1 = 1 - 0.5*u*(tau**2)
    g1 = tau - (1/6)*mu*(tau**3)/r**3 
    return f1, g1

def scaleEqLagrange(tauVec, sunVec, rhoHatVec, dsVec):
    mu = 1
    tau = tauVec[2]
    tau1 = tauVec[0]
    tau3 = tauVec[1]
    A1 = tau3/tau
    B1 = (A1/6)*((tau**2) - (tau3**2))
    A3 = -tau1/tau
    B3 = (A3/6)*((tau**2) - (tau1**2))

    A = (A1*dsVec[1] - dsVec[2] + A3*dsVec[3])/(-dsVec[0])
    B = (B1*dsVec[1] + B3*dsVec[3])/(-dsVec[0])

    E = (np.dot(rhoHatVec, sunVec))*-2
    F = (np.linalg.norm(sunVec))**2
    a = -((A**2) + A*E  + F)
    b = -mu*(2*A*B + B*E)
    c = -(mu**2)*(B**2)

    # print([c, 0, 0, b[0], 0, 0, a[0], 0, 1])
    # print(a, b, c)

    answers = np.polynomial.polynomial.polyroots([c, 0, 0, b, 0, 0, a, 0, 1])
    rootFinal = []
    for root in answers:
        if np.imag(root) == 0 and np.real(root) > 0:
            rootFinal.append(np.real(root))
    rhos = []
    for root2 in rootFinal:
        rhos.append(A + (mu*B)/(root2**3))

    return rootFinal, rhos
